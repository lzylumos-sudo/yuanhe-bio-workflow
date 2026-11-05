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
    rows = [row for row in json.loads((patches.PATCHES / 'manifest.json').read_text())
            if row['vendor'] == 'Peptide2Mol']
    vendor = tmp_path / '01_molecule_generation/vendor/Peptide2Mol'
    targets = []
    for row in rows:
        source = ROOT / '01_molecule_generation/vendor/Peptide2Mol' / row['file']
        if row['file'].endswith('/model.py'):
            source = Path(os.environ.get('PEPTIDE2MOL_TEST_SOURCE', str(source)))
        if not source.is_file(): pytest.skip('Pinned patch fixture not downloaded: ' + row['file'])
        target = vendor / row['file']
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(source.read_bytes())
        assert patches.sha(target) == row['after_sha256']
        if state != 'current':
            upgrade = row.get('upgrades', [])[0] if state == 'previous' and row.get('upgrades') else None
            name = upgrade['patch'] if upgrade else row['patch']
            subprocess.run(['patch', '--batch', '-R', '-p1', '-i', str(patches.PATCHES / name)], cwd=vendor, check=True)
            assert patches.sha(target) == (upgrade['before_sha256'] if upgrade else row['before_sha256'])
        targets.append((target, row))
    monkeypatch.setattr(patches, 'ROOT', tmp_path)
    patches.main('Peptide2Mol'); patches.main('Peptide2Mol')
    for target, row in targets:
        assert patches.sha(target) == row['after_sha256']
        original = target.read_bytes()
        target.write_bytes(original + b'# unknown source\n')
        with pytest.raises(SystemExit, match='revision differs'): patches.main('Peptide2Mol')
        target.write_bytes(original)
