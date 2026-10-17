"""Unified runtime configuration (plan §05). One place resolves every location; nothing in the source assumes
a particular machine, drive letter, /opt/sbvs, a user's home env, or a WSL GPU path.

Precedence (high -> low): explicit argument > environment variable > ``sbvs.toml`` (SBVS_CONFIG or <home>/sbvs.toml) > default.
All values are recorded in each run manifest so a run can be explained later.
"""
from __future__ import annotations
import os, platform, shutil, sys
from dataclasses import dataclass, field, asdict
from pathlib import Path

try:                                    # py3.11+: stdlib
    import tomllib
except ModuleNotFoundError:             # pragma: no cover
    tomllib = None

ENV_KEYS = {
    "home": "SBVS_HOME", "projects": "SBVS_PROJECT_ROOT", "scratch": "SBVS_SCRATCH", "db": "SBVS_DB",
    "cache": "SBVS_CACHE", "gnina": "SBVS_GNINA", "vina": "SBVS_VINA", "nvidia_smi": "SBVS_NVIDIA_SMI",
    "python": "SBVS_PYTHON", "gnina_lib_dirs": "SBVS_GNINA_LIB_DIRS", "entry": "SBVS_ENTRY",
}
ENTRIES = ("native_ubuntu", "windows_wsl", "mac_dev", "unknown")


@dataclass(frozen=True)
class Settings:
    home: Path
    projects: Path
    scratch: Path            # high-frequency temp: should be on LOCAL disk, not the archive volume
    db: Path                 # SQLite queue/state DB: local disk (exFAT/network FS are unsafe for SQLite locking)
    cache: Path
    gnina: str | None
    vina: str | None
    nvidia_smi: str | None
    python: str
    gnina_lib_dirs: tuple
    entry: str
    sources: dict = field(default_factory=dict, compare=False)   # key -> where the value came from

    def to_record(self) -> dict:
        d = asdict(self)
        return {k: (str(v) if isinstance(v, Path) else list(v) if isinstance(v, tuple) else v) for k, v in d.items()}


def _toml(path: Path) -> dict:
    if tomllib is None or not path.is_file():
        return {}
    with path.open("rb") as f:
        return tomllib.load(f)


def detect_entry() -> str:
    """Best-effort label only. It is NOT proof that a platform was validated (see STATUS maturity matrix)."""
    sysname = platform.system()
    if sysname == "Darwin":
        return "mac_dev"
    if sysname == "Linux":
        rel = platform.release().lower()
        if "microsoft" in rel or "wsl" in rel or os.environ.get("WSL_DISTRO_NAME"):
            return "windows_wsl"
        return "native_ubuntu"
    return "unknown"


def _which(name, extra=()):
    p = shutil.which(name)
    if p:
        return p
    for d in extra:
        c = shutil.which(name, path=str(d))
        if c:
            return c
    return None


def load(home: str | os.PathLike | None = None, cfg_file: str | os.PathLike | None = None, env=None) -> Settings:
    env = os.environ if env is None else env
    h = Path(home or env.get("SBVS_HOME") or Path(__file__).resolve().parents[2]).resolve()
    cf = Path(cfg_file or env.get("SBVS_CONFIG") or h / "sbvs.toml")
    t = _toml(cf).get("paths", {}) | _toml(cf).get("engines", {}) | _toml(cf).get("runtime", {})
    src: dict = {}

    def pick(key, default):
        e = env.get(ENV_KEYS[key])
        if e:
            src[key] = "env:" + ENV_KEYS[key]; return e
        if key in t:
            src[key] = f"toml:{cf.name}"; return t[key]
        src[key] = "default"; return default

    local = Path(env.get("XDG_STATE_HOME") or Path.home() / ".local" / "state") / "sbvs2"
    projects = Path(pick("projects", h / "projects")).expanduser().resolve()
    scratch = Path(pick("scratch", local / "scratch")).expanduser().resolve()
    db = Path(pick("db", local / "sbvs.sqlite3")).expanduser().resolve()
    cache = Path(pick("cache", local / "cache")).expanduser().resolve()
    bin_dirs = [h / "vendor", Path(sys.executable).parent]
    gnina = pick("gnina", None) or _which("gnina", bin_dirs)
    vina = pick("vina", None) or _which("vina", bin_dirs)
    smi = pick("nvidia_smi", None) or _which("nvidia-smi")
    libs = pick("gnina_lib_dirs", "")
    libs = tuple(x for x in (libs.split(os.pathsep) if isinstance(libs, str) else libs) if x)
    entry = pick("entry", None) or detect_entry()
    if entry not in ENTRIES:
        raise ValueError(f"unknown entry {entry!r}; expected one of {ENTRIES}")
    return Settings(h, projects, scratch, db, cache, gnina, vina, smi, str(pick("python", sys.executable)), libs, entry, src)


def is_network_or_exfat_unsafe_for_sqlite(path: Path) -> bool:
    """SQLite needs POSIX locks; refuse exFAT/NTFS-on-mount/SMB for the DB (plan: DB lives on local ext4/APFS)."""
    try:
        import subprocess
        out = subprocess.run(["df", "-T" if platform.system() == "Linux" else "-h", str(path.parent)], capture_output=True, text=True, timeout=5).stdout
        low = out.lower()
        return any(x in low for x in ("exfat", "ntfs", "fuseblk", "smb", "cifs", "nfs", "drvfs", "9p"))
    except Exception:
        return False
