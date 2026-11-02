"""A damaged bundled executable cannot pass the original official hash pin."""
import hashlib
import json
from pathlib import Path
import shutil
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'deployment'))
from verify_vina_bundle import verified_binary


def test_official_binary_identity():
    binary, spec = verified_binary()
    assert len(binary) == 4056088
    assert binary.startswith(b'\x7fELF')
    assert hashlib.sha256(binary).hexdigest() == 'f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644'
    assert spec['version'] == '1.2.7'


@pytest.mark.parametrize('damaged', ['archive', 'license_path'])
def test_damaged_distribution_is_rejected(tmp_path, damaged):
    manifest = ROOT / 'deployment/vina-bundle.json'
    spec = json.loads(manifest.read_text())
    paths = ['deployment/vina-bundle.json', spec['archive'], spec['license_path']]
    for name in paths:
        target = tmp_path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, target)
    (tmp_path / spec[damaged]).write_bytes(b'corrupt data')
    with pytest.raises(ValueError, match='hash mismatch'):
        verified_binary(tmp_path)
