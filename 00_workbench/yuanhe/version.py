"""One platform version; Git and environment identities are separate."""
from pathlib import Path
import os
import platform
import re
import subprocess
import sys
from importlib import metadata

ROOT = next((p for p in Path(__file__).resolve().parents if (p / 'VERSION').is_file() and (p / 'configs/module_registry.yaml').is_file()), Path(__file__).resolve().parents[2])
PRODUCT = "Yuanhe Bio Platform"
SUBTITLE = "AI-Driven Drug Discovery Platform"

def software_version():
    value = (ROOT / 'VERSION').read_text(encoding='utf-8').strip()
    if not re.fullmatch(r'(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)(?:-[0-9A-Za-z.-]+)?', value):
        raise ValueError('Invalid platform SemVer')
    return value

def git_identity(root=ROOT):
    def git(*args):
        try:
            p = subprocess.run(['git', '-C', str(root), *args], capture_output=True, text=True, timeout=10)
            return p.stdout.strip() if p.returncode == 0 else None
        except (OSError, subprocess.TimeoutExpired):
            return None
    commit = git('rev-parse', 'HEAD')
    status = git('status', '--porcelain', '--untracked-files=normal')
    return {'commit': commit or 'unknown', 'branch': git('branch', '--show-current') or 'unknown',
            'dirty': None if status is None else bool(status)}

def identity():
    return {'product': PRODUCT, 'subtitle': SUBTITLE, 'software_version': software_version(),
            'git': git_identity(), 'api_version': '1', 'release_status': 'development'}

def environment():
    packages = {}
    for name in ('rdkit', 'torch', 'openmm', 'parmed', 'openmmforcefields', 'openff-toolkit',
                 'MDAnalysis', 'meeko', 'numpy', 'PyYAML', 'psi4'):
        try: packages[name] = metadata.version(name)
        except metadata.PackageNotFoundError: packages[name] = 'unknown'
    return {'os': platform.platform(), 'python': sys.version, 'executable': sys.executable,
            'cpu': platform.processor() or platform.machine(), 'cpu_count': os.cpu_count(),
            'packages': packages, 'gpu': 'unknown', 'cuda': 'unknown',
            'container_image_digest': os.environ.get('YUANHE_IMAGE_DIGEST') if re.fullmatch(r'(?:[^@]+@)?sha256:[0-9a-f]{64}',os.environ.get('YUANHE_IMAGE_DIGEST','')) else 'unknown'}
