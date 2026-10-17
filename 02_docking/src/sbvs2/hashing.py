from __future__ import annotations
import hashlib, json, os, tempfile
from pathlib import Path


def sha256_file(path: str | os.PathLike) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_json(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def hash_obj(obj) -> str:
    return sha256_bytes(canonical_json(obj).encode("utf-8"))


def atomic_write_bytes(path: str | os.PathLike, data: bytes) -> None:
    """Temp file + fsync + rename in the same directory (plan §25)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(data)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def atomic_write_json(path, obj) -> None:
    atomic_write_bytes(path, (json.dumps(obj, indent=1, ensure_ascii=False, sort_keys=True) + "\n").encode("utf-8"))


def _is_os_junk(p: Path) -> bool:
    """macOS AppleDouble (._*) and Finder files that exFAT/macOS create behind our back; not project content."""
    return p.name.startswith("._") or p.name in (".DS_Store", "Thumbs.db")


def write_checksums(root: str | os.PathLike, name: str = "checksums.sha256") -> int:
    root = Path(root)
    lines = []
    for p in sorted(root.rglob("*")):
        if p.is_file() and p.name != name and not _is_os_junk(p):
            lines.append(f"{sha256_file(p)}  {p.relative_to(root).as_posix()}")
    atomic_write_bytes(root / name, ("\n".join(lines) + "\n").encode())
    return len(lines)


def verify_checksums(root: str | os.PathLike, name: str = "checksums.sha256") -> list[str]:
    """Return list of problems (empty == all verified)."""
    root = Path(root)
    problems = []
    for line in (root / name).read_text().splitlines():
        digest, rel = line.split("  ", 1)
        p = root / rel
        if not p.is_file():
            problems.append(f"missing: {rel}")
        elif sha256_file(p) != digest:
            problems.append(f"hash mismatch: {rel}")
    return problems
