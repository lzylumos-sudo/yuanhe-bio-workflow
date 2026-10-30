"""Versioned deployment settings fail closed without losing legacy compatibility."""
import json
import hashlib
from pathlib import Path
import subprocess
import sys
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'deployment'))
import setup_workbench as setup

def legacy():
    return dict(prefix='/new', model_source='/models', gmx='/gmx', host='127.0.0.1', port=8778)

def test_legacy_machine_settings_remain_readable(tmp_path):
    path = tmp_path / 'site.json'
    path.write_text(json.dumps(legacy()))
    assert setup.load_config(path) == legacy()

@pytest.mark.parametrize('metadata', [
    {'schema_version':'yuanhe.site/99', 'software_version':'1.0.0'},
    {'schema_version':'yuanhe.site/2'},
    {'software_version':'1.0.0'},
    {'schema_version':'yuanhe.site/2', 'software_version':None},
])
def test_invalid_schema_never_silently_downgrades(tmp_path, metadata):
    path = tmp_path / 'site.json'
    path.write_text(json.dumps({**legacy(), **metadata}))
    with pytest.raises(ValueError): setup.load_config(path)

def test_wrong_release_blocks_install_even_if_machine_plan_passes(tmp_path, monkeypatch):
    (tmp_path / 'VERSION').write_text('0.11.0-dev.4\n')
    monkeypatch.setattr(setup.deploy_workbench, 'plan', lambda *a: {'source_commit':'a'*40})
    monkeypatch.setattr(setup.shutil, 'which', lambda *a: None)
    monkeypatch.setattr(setup.subprocess, 'run', lambda *a, **kw: (_ for _ in ()).throw(OSError('fixture')))
    config = {**legacy(), 'schema_version':'yuanhe.site/2', 'software_version':'0.11.0-dev.3'}
    report = setup.preflight(config, tmp_path)
    assert not report['ready']
    assert not next(c for c in report['checks'] if c['check']=='release_version')['ok']

def test_current_templates_and_browser_match_platform_version():
    version = (ROOT / 'VERSION').read_text().strip()
    for name in ('site.example.json', 'site.campus-upgrade.json'):
        cfg = setup.load_config(ROOT / 'deployment' / name)
        assert cfg['software_version'] == version
        assert cfg['schema_version'] == 'yuanhe.site/2'
    browser = (ROOT / '00_workbench/webui/static/js/platform-version.js').read_text()
    assert json.dumps(version) in browser

def test_service_rejects_mixed_frontend_even_with_correct_version_and_commit():
    import verify_service as service
    version = (ROOT / 'VERSION').read_text().strip()
    endpoints = {'/api/health': {'ok': True, 'version': version},
                 '/api/v1/versions': {'software_version': version, 'git': {'commit': 'a'*40, 'dirty': False}},
                 '/api/studio/diagnostics': {'status': 'passed', 'all_runtime_checks_pass': True}}
    static = ROOT / '00_workbench/webui/static'
    def asset(name): return (static / name.lstrip('/')).read_bytes()
    ok = service.verify('', version, 'a'*40, endpoints.__getitem__, asset)
    assert ok['frontend_verified'] and len(ok['frontend_sha256']) == 7
    assert ok['scientific_acceptance'] == 'not_run'
    def stale(name): return b'old frontend' if name.endswith('/module-jobs.js') else asset(name)
    with pytest.raises(ValueError, match='Frontend differs'):
        service.verify('', version, 'a'*40, endpoints.__getitem__, stale)


@pytest.mark.parametrize('git_status', [
    {'commit': 'a'*40},
    {'commit': 'a'*40, 'dirty': None},
    {'commit': 'a'*40, 'dirty': True},
    {'commit': 'unknown', 'dirty': True},
])
def test_service_never_treats_unverified_git_status_as_clean(git_status):
    import verify_service as service
    version = (ROOT / 'VERSION').read_text().strip()
    endpoints = {'/api/health': {'ok': True, 'version': version},
                 '/api/v1/versions': {'software_version': version, 'git': git_status,
                                     'source': {'verified': True, 'commit': 'a'*40}},
                 '/api/studio/diagnostics': {'status': 'passed'}}
    with pytest.raises(ValueError, match='source commit'):
        service.verify('', version, 'a'*40, endpoints.__getitem__)


@pytest.mark.parametrize('dirty', [None, True])
def test_install_plan_requires_a_verified_clean_git_tree(tmp_path, monkeypatch, dirty):
    import deploy_workbench as deployment
    monkeypatch.setattr(deployment, 'git_identity', lambda root: {'commit': 'a'*40, 'dirty': dirty})
    with pytest.raises(ValueError, match='clean Git checkout'):
        deployment.source_identity(tmp_path)


def test_release_mismatch_alone_blocks_an_otherwise_ready_install(tmp_path, monkeypatch):
    root = tmp_path / 'checkout'
    (root / 'deployment').mkdir(parents=True)
    version = '0.11.0-dev.4'
    (root / 'VERSION').write_text(version)
    bundle = tmp_path / 'models'
    pins = {}
    for key, name, field in (('peptide2mol', 'Peptide2Mol', 'checkpoint'),
                             ('reinvent4', 'REINVENT4', 'prior')):
        directory = bundle / 'vendor' / name
        directory.mkdir(parents=True)
        (directory / 'LICENSE').write_text('synthetic fixture license')
        weight = directory / 'weight'
        weight.write_bytes(b'synthetic weight fixture')
        pins[key] = {field: weight.relative_to(bundle).as_posix(),
                     field + '_sha256': hashlib.sha256(weight.read_bytes()).hexdigest()}
    (root / 'deployment/model_dependencies.json').write_text(json.dumps(pins))
    monkeypatch.setattr(setup.platform, 'system', lambda: 'Linux')
    monkeypatch.setattr(setup.platform, 'machine', lambda: 'x86_64')
    monkeypatch.setattr(setup.shutil, 'which', lambda *a: '/synthetic-tool')
    monkeypatch.setattr(setup.subprocess, 'run', lambda argv, **kw:
                        subprocess.CompletedProcess(argv, 0, 'NVIDIA\nGPU support: CUDA', ''))
    monkeypatch.setattr(setup.deploy_workbench, 'plan', lambda *a: {'source_commit': 'a'*40})
    config = {**legacy(), 'model_source': str(bundle), 'schema_version': 'yuanhe.site/2',
              'software_version': version}
    assert setup.preflight(config, root)['ready'] is True
    report = setup.preflight({**config, 'software_version': '0.11.0-dev.3'}, root)
    assert report['ready'] is False
    assert [c['check'] for c in report['checks'] if not c['ok']] == ['release_version']
