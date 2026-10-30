"""Model reuse never imports an old installation's modified vendor code."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'deployment'))
import restore_models


def git(path, *args):
    return subprocess.check_output(['git', '-C', str(path), *args], stderr=subprocess.PIPE).decode().strip()


@pytest.fixture
def installation(tmp_path):
    root = tmp_path / 'new-platform'
    (root / 'deployment').mkdir(parents=True)
    old = tmp_path / 'old-platform' / '01_molecule_generation'
    pins = {}
    for key, name, field in [('peptide2mol', 'Peptide2Mol', 'checkpoint'),
                             ('reinvent4', 'REINVENT4', 'prior')]:
        upstream = tmp_path / 'upstream' / name
        upstream.mkdir(parents=True)
        (upstream / 'LICENSE').write_text('synthetic fixture license')
        (upstream / 'model.py').write_text('PINNED_SOURCE = True\n')
        (upstream / 'uv.lock').write_text('synthetic locked dependencies\n')
        git(upstream, 'init')
        git(upstream, 'add', '.')
        git(upstream, '-c', 'user.name=Fixture', '-c', 'user.email=fixture@example.invalid',
            'commit', '-m', 'Pinned source')
        commit = git(upstream, 'rev-parse', 'HEAD')
        directory = old / 'vendor' / name
        directory.mkdir(parents=True)
        weights = b'public synthetic model weights ' + name.encode()
        relative = f'vendor/{name}/weights/model.bin'
        pins[key] = {'repo_commit': commit, 'source_url': str(upstream), field: relative,
                     field + '_sha256': hashlib.sha256(weights).hexdigest()}
        (directory / 'weights').mkdir()
        (old / relative).write_bytes(weights)
        (directory / 'LICENSE').write_text('old local license placeholder')
        (directory / 'model.py').write_text('OLD_LOCAL_MODIFICATION = True\n')
        (directory / '.env').write_text('SYNTHETIC_OLD_SETTING=never-copy')
        (directory / 'old-wrapper.py').write_text('old wrapper must not survive')
    (root / 'deployment/model_dependencies.json').write_text(json.dumps(pins))
    return root, old, pins


@pytest.mark.parametrize('local_git_cache', [False, True])
def test_restore_rebuilds_exact_source_and_reuses_only_verified_weights(installation, local_git_cache):
    root, old, pins = installation
    if local_git_cache:
        for key, name in [('peptide2mol', 'Peptide2Mol'), ('reinvent4', 'REINVENT4')]:
            source = old / 'vendor' / name
            git(source, 'init')
            git(source, 'fetch', pins[key]['source_url'], pins[key]['repo_commit'])
            # There is no checkout: the dirty old files must still be ignored.
            pins[key]['source_url'] = str(root / 'unavailable-upstream')
        (root / 'deployment/model_dependencies.json').write_text(json.dumps(pins))
    restore_models.restore(old.parent, root)
    for key, name, field in [('peptide2mol', 'Peptide2Mol', 'checkpoint'),
                             ('reinvent4', 'REINVENT4', 'prior')]:
        vendor = root / '01_molecule_generation/vendor' / name
        assert git(vendor, 'rev-parse', 'HEAD') == pins[key]['repo_commit']
        assert (vendor / 'model.py').read_text() == 'PINNED_SOURCE = True\n'
        assert (vendor / 'uv.lock').read_text() == 'synthetic locked dependencies\n'
        assert (vendor / 'LICENSE').read_text() == 'synthetic fixture license'
        assert not (vendor / '.env').exists()
        assert not (vendor / 'old-wrapper.py').exists()
        assert (root / '01_molecule_generation' / pins[key][field]).read_bytes() == (old / pins[key][field]).read_bytes()
        assert (old / 'vendor' / name / 'model.py').read_text() == 'OLD_LOCAL_MODIFICATION = True\n'
    with pytest.raises(SystemExit, match='refusing overwrite'):
        restore_models.restore(old.parent, root)


def test_all_weights_are_checked_before_any_source_restore(installation):
    root, old, pins = installation
    (old / pins['reinvent4']['prior']).write_bytes(b'wrong weights')
    with pytest.raises(SystemExit, match='mismatched source weights'):
        restore_models.restore(old.parent, root)
    assert not (root / '01_molecule_generation/vendor').exists()
