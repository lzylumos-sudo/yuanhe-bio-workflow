"""Versioned input automation; the audited legacy simulation backend is retained."""
from ._platform_version import __version__


def governance():
    """Locate shared governance in a checkout or a packaged membrane export."""
    import sys
    from pathlib import Path
    root=Path(__file__).resolve().parents[3]
    shared=root/'00_workbench'
    if (shared/'yuanhe').is_dir() and str(shared) not in sys.path:sys.path.insert(0,str(shared))
    try:
        from yuanhe import records
        return records
    except ModuleNotFoundError as exc:
        if exc.name!='yuanhe':raise
        return None
