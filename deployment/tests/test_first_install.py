"""First install needs no old project; failed prerequisites never start installation."""
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'deployment'))
import fetch_generation_models as fetch
import setup_workbench as setup


def test_model_plan_has_pins_and_is_read_only(tmp_path):
    destination = tmp_path / 'new'
    record = fetch.plan(destination)
    pins = json.loads((ROOT / 'deployment/model_dependencies.json').read_text())
    for model, key, field in zip(record['models'], ('peptide2mol', 'reinvent4'), ('checkpoint', 'prior')):
        assert model['sha256'] == pins[key][field + '_sha256']
        assert model['commit'] == pins[key]['repo_commit']
        assert model['url'].startswith('https://')
    assert not destination.exists()
    destination.mkdir()
    with pytest.raises(ValueError, match='NEW'):
        fetch.plan(destination)


def test_verified_download_and_overwrite_refusal(tmp_path, monkeypatch):
    content = b'verified public model fixture'
    entry = dict(name='fixture', url='https://example.invalid/model', bytes=len(content),
                 sha256=hashlib.sha256(content).hexdigest())
    monkeypatch.setattr(fetch.urllib.request, 'urlopen', lambda *a, **k: io.BytesIO(content))
    target = tmp_path / 'model'
    fetch.download(entry, target)
    assert target.read_bytes() == content
    assert not target.with_name('model.part').exists()
    with pytest.raises(ValueError, match='overwrite'):
        fetch.download(entry, target)
    assert target.read_bytes() == content


@pytest.mark.parametrize('content', [b'wrong!', b'too long to accept', b'short'])
def test_invalid_model_never_becomes_a_weight(tmp_path, monkeypatch, content):
    monkeypatch.setattr(fetch.urllib.request, 'urlopen', lambda *a, **k: io.BytesIO(content))
    target = tmp_path / 'model'
    entry = dict(name='fixture', url='https://example.invalid/model', bytes=6, sha256='0' * 64)
    with pytest.raises(ValueError, match='size|SHA-256'):
        fetch.download(entry, target)
    assert not target.exists()


def test_source_fetch_failure_retains_receipt_without_false_success(tmp_path, monkeypatch):
    destination = tmp_path / 'new'
    record = fetch.plan(destination)
    monkeypatch.setattr(fetch.shutil, 'which', lambda name: '/git')
    def fail(argv, **kwargs):
        raise subprocess.CalledProcessError(128, argv)
    monkeypatch.setattr(fetch.subprocess, 'run', fail)
    with pytest.raises(subprocess.CalledProcessError):
        fetch.execute(record)
    receipt = json.loads((destination / 'model-bundle.json').read_text())
    assert receipt['status'] == 'failed' and receipt['verified'] == []
    assert receipt['inference'] == 'not_run'


@pytest.mark.parametrize('change', [{'port': True}, {'port': '8777'}, {'host': ''}, {'unexpected': 1}])
def test_site_config_rejects_ambiguous_or_extra_fields(tmp_path, change):
    value = json.loads((ROOT / 'deployment/site.example.json').read_text())
    value.update(change)
    config = tmp_path / 'site.json'
    config.write_text(json.dumps(value))
    with pytest.raises(ValueError):
        setup.load_config(config)


def test_failed_preflight_collects_errors_and_blocks_execute(tmp_path, monkeypatch, capsys):
    value = json.loads((ROOT / 'deployment/site.example.json').read_text())
    value['prefix'] = str(tmp_path / 'never-created')
    value['model_source'] = str(tmp_path / 'missing-models')
    config = tmp_path / 'site.json'
    config.write_text(json.dumps(value))
    monkeypatch.setattr(setup.platform, 'system', lambda: 'unsupported')
    monkeypatch.setattr(setup.shutil, 'which', lambda name: None)
    def failed_probe(argv, **kwargs):
        return subprocess.CompletedProcess(argv, 1, '', 'missing')
    monkeypatch.setattr(setup.subprocess, 'run', failed_probe)
    monkeypatch.setattr(setup.deploy_workbench, 'execute', lambda *a: pytest.fail('installation must not start'))
    report = setup.preflight(value)
    assert not report['ready']
    failed = {c['check'] for c in report['checks'] if not c['ok']}
    assert {'platform', 'driver', 'gromacs', 'model_bundle', 'installation_plan'} <= failed
    monkeypatch.setattr(sys, 'argv', ['setup', '--config', str(config), '--execute'])
    assert setup.main() == 1
    assert not Path(value['prefix']).exists()


def test_complete_preflight_and_install_dispatch_use_same_plan(tmp_path, monkeypatch):
    # Small hashes here are local fixtures; production manifest is not modified.
    root = tmp_path / 'checkout'
    (root / 'deployment').mkdir(parents=True)
    source = tmp_path / 'bundle'
    pins = {}
    for key, name, field in (('peptide2mol', 'Peptide2Mol', 'checkpoint'), ('reinvent4', 'REINVENT4', 'prior')):
        model = source / 'vendor' / name
        model.mkdir(parents=True)
        (model / 'LICENSE').write_text('fixture license')
        weight = model / 'weight'
        weight.write_bytes(b'fixture')
        pins[key] = {field: weight.relative_to(source).as_posix(),
                     field + '_sha256': hashlib.sha256(b'fixture').hexdigest()}
    (root / 'deployment/model_dependencies.json').write_text(json.dumps(pins))
    value = dict(prefix=str(tmp_path / 'new'), model_source=str(source), gmx='/gmx', host='127.0.0.1', port=8777)
    monkeypatch.setattr(setup.platform, 'system', lambda: 'Linux')
    monkeypatch.setattr(setup.platform, 'machine', lambda: 'x86_64')
    monkeypatch.setattr(setup.shutil, 'which', lambda name: '/tool')
    monkeypatch.setattr(setup.subprocess, 'run', lambda argv, **kw:
                        subprocess.CompletedProcess(argv, 0, 'NVIDIA\nGPU support: CUDA', ''))
    record = {'source_commit': 'a' * 40}
    monkeypatch.setattr(setup.deploy_workbench, 'plan', lambda *a: record)
    report = setup.preflight(value, root)
    assert report['ready'] and report['plan'] is record
    assert report['scientific_acceptance'] == 'not_run'
    config = tmp_path / 'site.json'
    config.write_text(json.dumps(value))
    monkeypatch.setattr(setup, 'preflight', lambda config: report)
    calls = []
    monkeypatch.setattr(setup.deploy_workbench, 'execute', calls.append)
    monkeypatch.setattr(sys, 'argv', ['setup', '--config', str(config)])
    assert setup.main() == 0 and not calls
    monkeypatch.setattr(sys, 'argv', ['setup', '--config', str(config), '--execute'])
    assert setup.main() == 0 and calls == [record]


def test_release_requires_browser_regression():
    source = (ROOT / '.github/workflows/release.yml').read_text()
    assert 'engineering.yml predocking-models.yml evidence-ui.yml' in source
