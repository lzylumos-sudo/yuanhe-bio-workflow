"""Remote compute host staging (plan §05 Mode B). The remote never reads this machine's paths:
inputs are copied to <remote_root>/incoming/<stage_id>/ with a SHA-256 manifest verified on the remote side;
results come back into a NEW local directory and are verified against the remote manifest before use.

Only ssh/rsync argv arrays are used (no shell strings built from user input). Credentials stay in ~/.ssh.
"""
from __future__ import annotations
import json, re, shlex, subprocess
from pathlib import Path

from .hashing import sha256_file, atomic_write_json

SAFE = re.compile(r"^[A-Za-z0-9_.@:-]+$")


class RemoteError(Exception):
    pass


class Remote:
    def __init__(self, host: str, root: str):
        if not SAFE.match(host):
            raise ValueError(f"unsafe ssh host alias: {host!r}")
        if not root.startswith("/") and not root.startswith("~"):
            raise ValueError("remote root must be absolute or ~-relative")
        self.host, self.root = host, root

    def _ssh(self, script: str, timeout=120) -> str:
        p = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15", self.host, "bash", "-lc", shlex.quote(script)],
                           capture_output=True, text=True, timeout=timeout)
        if p.returncode != 0:
            raise RemoteError(f"ssh rc={p.returncode}: {p.stderr[-400:]}")
        return p.stdout

    def remote_abs(self) -> str:
        return self._ssh(f"cd {shlex.quote(self.root) if not self.root.startswith('~') else self.root} && pwd").strip()

    def stage(self, local_dir: Path, stage_id: str) -> dict:
        """Copy a directory of inputs + manifest; verify hashes ON THE REMOTE. Idempotent for identical content."""
        if not SAFE.match(stage_id):
            raise ValueError("unsafe stage id")
        local_dir = Path(local_dir)
        files = sorted(p for p in local_dir.rglob("*") if p.is_file() and not p.name.startswith("._"))
        manifest = {str(p.relative_to(local_dir)): sha256_file(p) for p in files}
        mf = local_dir / "STAGE_MANIFEST.json"; atomic_write_json(mf, manifest)
        dest = f"{self.remote_abs()}/incoming/{stage_id}"
        self._ssh(f"mkdir -p {shlex.quote(dest)}")
        r = subprocess.run(["rsync", "-rt", "--exclude=._*", "-e", "ssh -o BatchMode=yes", f"{local_dir}/", f"{self.host}:{dest}/"],
                           capture_output=True, text=True, timeout=3600)
        if r.returncode != 0:
            raise RemoteError("rsync failed: " + r.stderr[-400:])
        out = self._ssh(f"cd {shlex.quote(dest)} && sha256sum -c --quiet <(python3 -c \"import json;[print(h+'  '+k) for k,h in json.load(open('STAGE_MANIFEST.json')).items()]\") && echo VERIFIED")
        if "VERIFIED" not in out:
            raise RemoteError("remote hash verification failed")
        return {"stage_id": stage_id, "remote_dir": dest, "n_files": len(manifest), "verified": True}

    def fetch(self, remote_dir: str, local_dest: Path) -> dict:
        """Pull results to a NEW local dir; compare against a remote-generated manifest. Never merges into an existing dir."""
        local_dest = Path(local_dest)
        if local_dest.exists():
            raise FileExistsError(f"refusing to fetch into existing dir: {local_dest}")
        self._ssh(f"cd {shlex.quote(remote_dir)} && find . -type f ! -name RESULT_MANIFEST.sha256 -print0 | sort -z | xargs -0 sha256sum > RESULT_MANIFEST.sha256 && wc -l < RESULT_MANIFEST.sha256")
        local_dest.mkdir(parents=True)
        r = subprocess.run(["rsync", "-rt", "-e", "ssh -o BatchMode=yes", f"{self.host}:{remote_dir}/", f"{local_dest}/"],
                           capture_output=True, text=True, timeout=7200)
        if r.returncode != 0:
            raise RemoteError("rsync fetch failed: " + r.stderr[-400:])
        bad, n = [], 0
        for line in (local_dest / "RESULT_MANIFEST.sha256").read_text().splitlines():
            h, _, rel = line.partition("  ")
            n += 1
            p = local_dest / rel
            if not p.is_file() or sha256_file(p) != h:
                bad.append(rel)
        return {"files": n, "mismatch": bad, "verified": not bad}
