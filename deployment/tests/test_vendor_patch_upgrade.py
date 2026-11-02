"""Fresh and previously deployed vendor patches converge to the same verified source."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'deployment'))
import apply_vendor_patches as patches


@pytest.mark.parametrize('state', ['original', 'previous', 'current'])
def test_patch_install_and_upgrade_are_idempotent(tmp_path, monkeypatch, state):
    row = json.loads((patches.PATCHES / 'manifest.json').read_text())[0]
    source = Path(os.environ.get('PEPTIDE2MOL_TEST_SOURCE', str(
        ROOT / '01_molecule_generation/vendor/Peptide2Mol' / row['file'])))
    if not source.is_file(): pytest.skip('Pinned sampling fixture not downloaded')
    vendor = tmp_path / '01_molecule_generation/vendor/Peptide2Mol'
    target = vendor / row['file']; target.parent.mkdir(parents=True); target.write_bytes(source.read_bytes())
    assert patches.sha(target) == row['after_sha256']
    if state != 'current':
        name = row['patch'] if state == 'original' else row['upgrades'][0]['patch']
        subprocess.run(['patch', '--batch', '-R', '-p1', '-i', str(patches.PATCHES / name)], cwd=vendor, check=True)
        expected = row['before_sha256'] if state == 'original' else row['upgrades'][0]['before_sha256']
        assert patches.sha(target) == expected
    monkeypatch.setattr(patches, 'ROOT', tmp_path)
    patches.main('Peptide2Mol'); patches.main('Peptide2Mol')
    assert patches.sha(target) == row['after_sha256']
    target.write_text(target.read_text() + '# unknown source\n')
    with pytest.raises(SystemExit, match='revision differs'): patches.main('Peptide2Mol')
