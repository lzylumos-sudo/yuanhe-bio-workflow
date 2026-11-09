"""The doctor diagnoses the deployed configuration and remains read-only."""
import importlib.util
import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'00_workbench'))
sys.path.insert(0, str(ROOT/'02_docking/src'))
spec = importlib.util.spec_from_file_location('preflight_doctor', ROOT/'deployment/preflight_doctor.py')
doctor = importlib.util.module_from_spec(spec); spec.loader.exec_module(doctor)
from webui import workflow
from yuanhe_preflight.config import defaults


def configured(root):
    cfg = defaults(root)
    for key in ('model_path', 'templates_path', 'stock_path'):
        target = root/'assets'/key; target.parent.mkdir(exist_ok=True)
        target.write_text(key)
        cfg['retro'][key] = str(target.relative_to(root))
    return cfg


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def test_default_uses_deployed_research_config(tmp_path, monkeypatch):
    cfg = configured(tmp_path)
    runtime = tmp_path/'deployment/runtime.json'
    save(runtime, {'research':{'predocking':cfg}, 'predocking':defaults(tmp_path)})
    before = runtime.read_bytes()
    monkeypatch.setattr(doctor, 'ROOT', tmp_path)
    monkeypatch.setattr(workflow, 'ROOT', tmp_path)
    monkeypatch.setattr(sys, 'argv', ['preflight_doctor.py'])
    monkeypatch.setattr(doctor.importlib.util, 'find_spec', lambda name: object())
    # Loading the runtime uses the same research-over-top-level precedence as the UI.
    resolved, source = doctor.load_config()
    assert source == str(runtime)
    assert resolved['retro']['model_path'] == str(tmp_path/'assets/model_path')
    assert doctor.main() == 0
    assert runtime.read_bytes() == before


def test_top_level_runtime_predocking_is_honored(tmp_path, monkeypatch):
    cfg = configured(tmp_path)
    save(tmp_path/'deployment/runtime.json', {'predocking':cfg})
    monkeypatch.setattr(doctor, 'ROOT', tmp_path)
    monkeypatch.setattr(workflow, 'ROOT', tmp_path)
    resolved, _ = doctor.load_config()
    assert resolved['retro']['stock_path'] == str(tmp_path/'assets/stock_path')


@pytest.mark.parametrize('wrapper', ['standalone', 'runtime', 'research', 'resolved'])
def test_explicit_config_forms_and_cwd_independence(tmp_path, monkeypatch, wrapper):
    cfg = configured(tmp_path)
    wrapped = {'standalone':cfg, 'runtime':{'predocking':cfg},
               'research':{'research':{'predocking':cfg}},
               'resolved':{'config':{'research':{'predocking':cfg}}}}[wrapper]
    path = tmp_path/'custom.json'; save(path, wrapped)
    monkeypatch.setattr(doctor, 'ROOT', tmp_path)
    monkeypatch.chdir(tmp_path.parent)
    resolved, source = doctor.load_config(path)
    assert resolved['retro']['templates_path'] == str(tmp_path/'assets/templates_path')
    assert source == str(path)


def test_no_runtime_uses_current_workbench_defaults(tmp_path, monkeypatch):
    monkeypatch.setattr(doctor, 'ROOT', tmp_path)
    monkeypatch.setattr(workflow, 'ROOT', tmp_path)
    cfg, source = doctor.load_config()
    assert source == 'workbench research defaults'
    assert cfg['engine']['version'] == workflow.defaults()['predocking']['engine']['version']


def test_invalid_explicit_config_does_not_fall_back(tmp_path):
    path = tmp_path/'bad.json'; save(path, {})
    with pytest.raises(ValueError, match='config_fields'):
        doctor.load_config(path)


def test_missing_explicit_config_does_not_fall_back(tmp_path):
    with pytest.raises(FileNotFoundError):
        doctor.load_config(tmp_path/'missing.json')


def test_missing_deployed_asset_still_fails(tmp_path, monkeypatch, capsys):
    cfg = configured(tmp_path)
    (tmp_path/'assets/stock_path').unlink()
    save(tmp_path/'deployment/runtime.json', {'predocking':cfg})
    monkeypatch.setattr(doctor, 'ROOT', tmp_path)
    monkeypatch.setattr(workflow, 'ROOT', tmp_path)
    monkeypatch.setattr(sys, 'argv', ['preflight_doctor.py'])
    monkeypatch.setattr(doctor.importlib.util, 'find_spec', lambda name: object())
    assert doctor.main() == 1
    report = json.loads(capsys.readouterr().out)
    assert [item['item'] for item in report['missing']] == ['stock_path']
