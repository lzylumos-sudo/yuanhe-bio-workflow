"""Engine discovery & device check (plan §05, §22). Never downloads a wrong-arch binary, never fakes success."""
from __future__ import annotations
import platform, shutil, subprocess

SUPPORTED_ARCH = {"x86_64", "amd64"}


def _run(argv, timeout=20):
    try:
        p = subprocess.run(argv, capture_output=True, text=True, timeout=timeout)
        return p.returncode, (p.stdout + p.stderr).strip()
    except Exception as e:   # noqa
        return -1, str(e)


def find_engine(name: str, extra_dirs=()):
    p = shutil.which(name)
    if p:
        return p
    for d in extra_dirs:
        cand = shutil.which(name, path=str(d))
        if cand:
            return cand
    return None


def engine_status(name: str, extra_dirs=()) -> dict:
    arch = platform.machine().lower()
    exe = find_engine(name, extra_dirs)
    if arch not in SUPPORTED_ARCH:
        return {"engine": name, "status": "unsupported_arch", "arch": arch, "path": exe,
                "note": "x86-64 Linux is the first-release target; use a remote x86 compute host (SSH tunnel)."}
    if not exe:
        return {"engine": name, "status": "not_installed", "arch": arch}
    rc, out = _run([exe, "--version"])
    return {"engine": name, "status": "installed_unvalidated" if rc == 0 else "broken", "path": exe,
            "version": out.splitlines()[0] if out else "", "arch": arch}


def gpu_status() -> dict:
    from .. import settings as _s
    exe = _s.load().nvidia_smi      # path comes from config/env/PATH; a WSL location must be configured, not assumed
    if not exe:
        return {"status": "no_nvidia_smi"}
    rc, out = _run([exe, "--query-gpu=name,memory.total,driver_version", "--format=csv,noheader"])
    return {"status": "visible" if rc == 0 else "error", "detail": out, "note": "visible != CNN verified; run a real GNINA job"}
