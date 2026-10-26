"""Exercise fresh configuration, preservation, install failure and service identity."""
import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'deployment'))
import configure_workbench as configuration
import deploy_workbench as deployment
import verify_service as service


def assets(root):
    from yuanhe_preflight.config import defaults
    cfg = defaults(root)
    cfg['engine']['toxicity']['python'] = '/installation/admet/bin/python'
    cfg['engine']['synthesis'].update(scscore_path='/installation/scscore.gz', scscore_sha256='1' * 64)
    for key in ('model', 'templates', 'stock'):
        cfg['retro'].update({key + '_path': '/installation/' + key, key + '_version': 'version1',
                             key + '_license': 'test-license'})
    cfg['retro']['python'] = '/installation/retro/bin/python'
    return cfg


def test_fresh_install_uses_current_schema_and_preserves_science(tmp_path):
    from research.config import defaults, resolve
    initial = defaults()
    result = configuration.configured({}, tmp_path, Path('/tools/gmx'), assets(tmp_path), tmp_path)
    assert resolve(result['research']) == result['research']
    assert result['predocking'] == result['research']['predocking']
    assert 'catalogs' not in result['predocking']
    assert result['predocking']['engine']['version'] == '3.0.0'
    for key in ('md_production_ns', 'md_replicas', 'initial', 'refine', 'peptide_attempts'):
        assert result['research'][key] == initial[key]
    assert result['chem_python'] == result['screen_python'] == str(tmp_path / 'envs/runtime/bin/python')


def test_asset_binding_does_not_reset_existing_thresholds_or_search_budget(tmp_path):
    from research.config import defaults
    cfg = {'research': defaults()}
    cfg['research']['md_limit'] = 2
    cfg['research']['wetlab_min'] = cfg['research']['wetlab_max'] = 2
    pc = cfg['research']['predocking']
    pc['retro']['time_limit_seconds'] = 30
    pc['engine']['synthesis']['sa_max'] = 5
    previous = copy.deepcopy(cfg)
    result = configuration.configured(cfg, tmp_path, Path('/tools/gmx'), assets(tmp_path), tmp_path)
    assert cfg == previous
    assert result['research']['md_limit'] == 2
    assert result['predocking']['engine']['synthesis']['sa_max'] == 5
    assert result['predocking']['retro']['time_limit_seconds'] == 30


def setup_config(tmp_path):
    root = tmp_path / 'checkout'
    prefix = tmp_path / 'installation'
    (root / 'deployment').mkdir(parents=True)
    (prefix / 'predocking-assets').mkdir(parents=True)
    (prefix / 'predocking-assets/predocking-config.json').write_text(json.dumps(assets(root)))
    return root, prefix


def test_legacy_runtime_is_rejected_without_writing(tmp_path, monkeypatch):
    root, prefix = setup_config(tmp_path)
    old = (ROOT / '00_workbench/configs/research_dual_blackwell.json').read_text()
    dest = root / 'deployment/runtime.json'
    dest.write_text(json.dumps({'research': json.loads(old)}))
    before = dest.read_bytes()
    monkeypatch.setattr(configuration.subprocess, 'run', lambda *a, **k:
                        subprocess.CompletedProcess(a, 0, 'GPU support: CUDA', ''))
    with pytest.raises(ValueError, match='incompatible'):
        configuration.write_config(prefix, Path('/gmx'), root)
    assert dest.read_bytes() == before
    assert not (prefix / 'evidence').exists()


def test_configure_backs_up_existing_valid_runtime(tmp_path, monkeypatch):
    root, prefix = setup_config(tmp_path)
    dest = root / 'deployment/runtime.json'
    dest.write_text('{"structure":"data/structures/9nu4.pdb"}\n')
    before = dest.read_bytes()
    monkeypatch.setattr(configuration.subprocess, 'run', lambda *a, **k:
                        subprocess.CompletedProcess(a, 0, 'GPU support: CUDA', ''))
    configuration.write_config(prefix, Path('/gmx'), root)
    assert next((prefix / 'evidence').glob('runtime.before-*.json')).read_bytes() == before
    assert json.loads(dest.read_text())['research']['md_production_ns'] == 500


def test_configure_rejects_cpu_gromacs_before_touching_runtime(tmp_path, monkeypatch):
    root, prefix = setup_config(tmp_path)
    monkeypatch.setattr(configuration.subprocess, 'run', lambda *a, **k:
                        subprocess.CompletedProcess(a, 0, 'GPU support: disabled', ''))
    with pytest.raises(ValueError, match='CUDA'):
        configuration.write_config(prefix, Path('/gmx'), root)
    assert not (root / 'deployment/runtime.json').exists()


def test_plan_is_read_only_and_includes_synthesis_before_configuration(tmp_path, monkeypatch):
    root = tmp_path / 'checkout'
    root.mkdir()
    models = tmp_path / 'licensed'
    models.mkdir()
    gmx = tmp_path / 'gmx'
    gmx.write_text('test')
    gmx.chmod(0o700)
    monkeypatch.setattr(deployment, 'source_identity', lambda root: 'a' * 40)
    prefix = tmp_path / 'installation'
    result = deployment.plan(prefix, models, gmx, '127.0.0.1', 8777, root)
    steps = [name for name, argv in result['commands']]
    assert steps.index('synthesis') < steps.index('configure') < steps.index('check')
    assert not prefix.exists()
    prefix.mkdir()
    with pytest.raises(ValueError, match='NEW'):
        deployment.plan(prefix, models, gmx, '127.0.0.1', 8777, root)


def test_install_failure_records_stage_and_never_hands_off_service(tmp_path, monkeypatch):
    monkeypatch.setattr(deployment.platform, 'system', lambda: 'Linux')
    monkeypatch.setattr(deployment.platform, 'machine', lambda: 'x86_64')
    monkeypatch.setattr(deployment.shutil, 'which', lambda name: '/test/' + name)
    calls = []
    def run(argv, **kwargs):
        calls.append(argv)
        return subprocess.CompletedProcess(argv, 0 if '--version' in argv else 7,
                                           'GPU support: CUDA', '')
    monkeypatch.setattr(deployment.subprocess, 'run', run)
    prefix = tmp_path / 'new'
    record = {'prefix': str(prefix), 'gmx': '/gmx', 'commands': [('runtime', ['install']), ('next', ['never'])]}
    with pytest.raises(RuntimeError, match='runtime'):
        deployment.execute(record, tmp_path)
    state = json.loads((prefix / 'evidence/deployment.json').read_text())
    assert state['status'] == 'failed' and state['steps'][0]['returncode'] == 7
    assert ['never'] not in calls
    assert not (prefix / 'evidence/start-workbench.sh').exists()


def test_service_verification_rejects_old_same_version_and_dirty_code():
    values = {'/api/health': {'ok': True, 'version': '0.11.0-dev.3'},
              '/api/v1/versions': {'software_version': '0.11.0-dev.3', 'git': {'commit': 'b'*40, 'dirty': False}},
              '/api/studio/diagnostics': {'status': 'idle'}}
    with pytest.raises(ValueError, match='source commit'):
        service.verify('', '0.11.0-dev.3', 'a'*40, values.__getitem__)
    values['/api/v1/versions']['git'] = {'commit': 'a'*40, 'dirty': True}
    with pytest.raises(ValueError, match='modified'):
        service.verify('', '0.11.0-dev.3', 'a'*40, values.__getitem__)


def test_verified_package_can_identify_service_without_claiming_science():
    values = {'/api/health': {'ok': True, 'version': '0.11.0-dev.3'},
              '/api/v1/versions': {'software_version': '0.11.0-dev.3', 'git': {'commit': 'unknown'},
                                   'source': {'commit': 'a'*40, 'verified': True}},
              '/api/studio/diagnostics': {'status': 'idle'}}
    result = service.verify('', '0.11.0-dev.3', 'a'*40, values.__getitem__)
    assert result['service_verified'] and not result['readiness_passed']
    assert result['scientific_acceptance'] == 'not_run'
    values['/api/v1/versions']['source']['verified'] = False
    with pytest.raises(ValueError):
        service.verify('', '0.11.0-dev.3', 'a'*40, values.__getitem__)


def test_package_source_identity_detects_tampering(tmp_path):
    import hashlib
    from yuanhe.version import source_identity
    version = tmp_path / 'VERSION'
    version.write_text('0.11.0-dev.3\n')
    manifest = {'schema': 'yuanhe.source-package/1', 'git_commit': 'a'*40,
                'software_version': '0.11.0-dev.3', 'files': {'VERSION': hashlib.sha256(version.read_bytes()).hexdigest()}}
    (tmp_path / 'SOURCE_MANIFEST.json').write_text(json.dumps(manifest))
    unknown = {'commit': 'unknown', 'dirty': None}
    assert source_identity(tmp_path, unknown)['verified']
    version.write_text('0.11.0-dev.2\n')
    assert not source_identity(tmp_path, unknown)['verified']
