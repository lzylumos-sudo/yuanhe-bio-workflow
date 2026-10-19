"""Verify actual standalone archives, hashes, license references and entrypoints."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import zipfile
import pytest

ROOT = Path(__file__).resolve().parents[2]

@pytest.mark.parametrize('script,module,entry', [
    ('package_component.py', 'docking', ['sbvs.cli']),
    ('package_component.py', 'screening', ['zvss.cli']),
    ('package_md.py', None, ['scientific.md_cli']),
    ('package_evidence.py', None, ['scientific.peptide_contacts']),
])
def test_standalone_package(tmp_path, script, module, entry):
    archive = tmp_path/'package.zip'
    cmd = [sys.executable, str(ROOT/'deployment'/script), '--output', str(archive)]
    if module:
        cmd += ['--module', module]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        manifest = json.loads(z.read('MANIFEST.json'))
        for name, digest in manifest['files'].items():
            assert hashlib.sha256(z.read(name)).hexdigest() == digest
        for name in ['LICENSE', 'LICENSE_SCOPE.txt', 'THIRD_PARTY_NOTICES.txt']:
            assert z.read(name) == (ROOT/name).read_bytes()
        assert '02_docking/docs/licenses/ccsb-scripps_AutoDock-Vina_LICENSE' in z.namelist()
        z.extractall(tmp_path/'unpacked')
    import os
    env = dict(os.environ)
    env['PYTHONPATH'] = str(tmp_path/'unpacked'/('module/src' if module else '.'))
    subprocess.run([sys.executable, '-m', *entry, '--help'],
                   cwd=tmp_path/'unpacked', env=env, check=True, capture_output=True, timeout=30)
    # A second invocation must preserve the existing package.
    original = archive.read_bytes()
    assert subprocess.run(cmd, capture_output=True).returncode != 0
    assert archive.read_bytes() == original
