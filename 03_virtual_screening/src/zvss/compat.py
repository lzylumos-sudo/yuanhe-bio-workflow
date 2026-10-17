"""Environment compatibility layer.

Two host facts this project has to live with on the build machine:

1. The working directory contains CJK characters
   (``C:\\Users\\user\\Desktop\\元核生物``).  RDKit's C++ layer opens files
   through narrow ``char*`` paths, so any RDKit call that takes a *filename*
   fails there with ``OSError: Bad input file``.  Solution: every RDKit file
   read goes through :func:`read_mol_blocks` (Python reads bytes, RDKit parses a
   block) and every RDKit write goes through :func:`write_sdf` (RDKit makes the
   text, Python writes it).

2. Some RDKit data files (``BaseFeatures.fdef`` used by PoseBusters/SUCOS, the
   Crippen tables, ...) are loaded by C++ from ``RDConfig.RDDataDir``, i.e. from
   the interpreter's own site-packages path -- which is CJK here.  Solution:
   :func:`bootstrap` mirrors ``rdkit/Data`` into an ASCII scratch directory and
   repoints ``RDConfig.RDDataDir`` at the mirror *before* posebusters imports.

Nothing here changes scientific behaviour; it only relocates byte I/O.
"""

from __future__ import annotations

import os
import shutil
import sys
import tempfile

_BOOTSTRAPPED = False
_STATE: dict[str, str | None] = {"scratch": None, "rdkit_data": None}


def _is_ascii(path: str) -> bool:
    try:
        path.encode("ascii")
        return True
    except UnicodeEncodeError:
        return False


def _writable_ascii_dir(candidate: str) -> bool:
    try:
        os.makedirs(candidate, exist_ok=True)
        probe = os.path.join(candidate, ".zvss_write_probe")
        with open(probe, "w", encoding="utf-8") as handle:
            handle.write("ok")
        os.remove(probe)
        return True
    except OSError:
        return False


def _temp_roots() -> list[str]:
    """Every plausible platform temp root, before writability is decided.

    ``tempfile.gettempdir()`` is *not* trustworthy here: when the real TEMP is
    not writable it silently falls back to the current directory, which is the
    CJK project path.  The environment variables themselves are the ground truth.
    """
    roots: list[str] = []
    for key in ("TMPDIR", "TEMP", "TMP"):
        value = os.environ.get(key)
        if value:
            roots.append(value)
    local = os.environ.get("LOCALAPPDATA")
    if local:
        roots.append(os.path.join(local, "Temp"))
    roots.append(tempfile.gettempdir())
    seen: set[str] = set()
    unique = []
    for root in roots:
        normalised = os.path.normpath(root)
        if normalised not in seen:
            seen.add(normalised)
            unique.append(normalised)
    return unique


def _session_scratch_candidates() -> list[str]:
    """ASCII scratch candidates, including per-session temp folders.

    Sandboxes frequently allow writes only inside one private temp folder, so
    each temp root is probed directly and then through its children.
    """
    candidates: list[str] = []
    for root in _temp_roots():
        candidates.append(os.path.join(root, "zvss-scratch"))
        if not os.path.isdir(root):
            continue
        try:
            children = sorted(os.listdir(root))
        except OSError:
            continue
        for child in children:
            path = os.path.join(root, child)
            if not os.path.isdir(path) or not _is_ascii(path):
                continue
            if child.startswith(("dsh-", "zvss-", "tmp", "Temp", "pip-")):
                candidates.append(os.path.join(path, "zvss-scratch"))
    return candidates


def scratch_root() -> str:
    """Return a writable ASCII directory used for staging and RDKit data.

    Order of preference: explicit ``ZVSS_SCRATCH``; the platform temp roots and
    their per-session sub-folders; a sibling of the temp root.  The last resort
    is a project-local scratch directory, which still works for pure-Python code
    but disables RDKit file APIs -- callers test :func:`ascii_scratch_available`
    and degrade honestly instead of pretending the check ran.
    """
    if _STATE["scratch"]:
        return _STATE["scratch"]  # type: ignore[return-value]
    candidates: list[str] = []
    env = os.environ.get("ZVSS_SCRATCH")
    if env:
        candidates.append(env)
    candidates.extend(_session_scratch_candidates())
    for root in _temp_roots():
        if _is_ascii(root):
            candidates.append(os.path.join(os.path.dirname(root.rstrip("\\/")), "zvss-scratch"))
    for candidate in candidates:
        if candidate and _is_ascii(candidate) and _writable_ascii_dir(candidate):
            _STATE["scratch"] = candidate
            _STATE["scratch_source"] = "environment"
            return candidate
    fallback = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), ".scratch")
    os.makedirs(fallback, exist_ok=True)
    _STATE["scratch"] = fallback
    _STATE["scratch_source"] = "project_fallback_non_ascii"
    return fallback


def ascii_scratch_available() -> bool:
    """True when RDKit's C++ file readers can be pointed at ASCII paths."""
    return _is_ascii(scratch_root())



def configure_temp(scratch: str | None = None) -> str:
    """Point TMP/TEMP at an ASCII staging directory so temp-file APIs behave."""
    root = scratch or scratch_root()
    stage = os.path.join(root, "tmp")
    os.makedirs(stage, exist_ok=True)
    for key in ("TMP", "TEMP", "TMPDIR"):
        os.environ[key] = stage
    tempfile.tempdir = stage
    return stage


def rdkit_data_dir() -> str:
    """Mirror ``rdkit/Data`` into the ASCII scratch and return the mirror path."""
    if _STATE["rdkit_data"]:
        return _STATE["rdkit_data"]  # type: ignore[return-value]
    from rdkit import RDConfig  # imported lazily; rdkit import is slow

    source = RDConfig.RDDataDir
    if _is_ascii(source):
        _STATE["rdkit_data"] = source
        return source
    root = scratch_root()
    target = os.path.join(root, "rdkit-data")
    sentinel = os.path.join(target, "BaseFeatures.fdef")
    if not os.path.exists(sentinel):
        if os.path.isdir(target):
            shutil.rmtree(target, ignore_errors=True)
        shutil.copytree(source, target)
    _STATE["rdkit_data"] = target
    return target


def bootstrap() -> dict[str, str | None]:
    """Idempotent process-wide setup.  Safe to call from every entry point."""
    global _BOOTSTRAPPED
    if _BOOTSTRAPPED:
        return dict(_STATE)
    scratch = scratch_root()
    _STATE["scratch"] = scratch
    configure_temp(scratch)
    try:
        from rdkit import RDConfig

        data = rdkit_data_dir()
        RDConfig.RDDataDir = data
        _STATE["rdkit_data"] = data
    except Exception as exc:  # pragma: no cover - rdkit missing is a hard error elsewhere
        print(f"[zvss.compat] rdkit data relocation skipped: {exc}", file=sys.stderr)
    _BOOTSTRAPPED = True
    return dict(_STATE)


def report() -> dict[str, object]:
    """Diagnostics for ``zvss doctor``."""
    info: dict[str, object] = {
        "python": sys.version.split()[0],
        "executable": sys.executable,
        "cwd": os.getcwd(),
        "cwd_is_ascii": _is_ascii(os.getcwd()),
        "scratch": scratch_root(),
        "scratch_is_ascii": _is_ascii(scratch_root()),
        "rdkit_data_dir": _STATE.get("rdkit_data"),
        "tmpdir": tempfile.gettempdir(),
    }
    return info
