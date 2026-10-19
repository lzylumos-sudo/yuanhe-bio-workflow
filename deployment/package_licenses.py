"""License files shared by independently distributed source archives."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def license_files():
    files = {name: ROOT / name for name in (
        'LICENSE', 'LICENSE_SCOPE.txt', 'THIRD_PARTY_NOTICES.md',
        'THIRD_PARTY_NOTICES.txt',
    )}
    # Keep the paths referenced in the notices valid inside exported archives.
    for folder in ('02_docking/docs/licenses', '03_virtual_screening',
                   '00_workbench/webui/static/vendor'):
        base = ROOT / folder
        for path in base.rglob('*'):
            if not path.is_file() or path.is_symlink():
                continue
            if (folder.endswith('docs/licenses') or
                path.name.startswith(('LICENSE', 'COPYING', 'NOTICE'))):
                files[path.relative_to(ROOT).as_posix()] = path
    return files
