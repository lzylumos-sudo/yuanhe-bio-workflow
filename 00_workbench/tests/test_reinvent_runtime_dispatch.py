"""Public/synthetic launch probes; no model weights, sampling, or GPU execution."""
from contextlib import nullcontext
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / '00_workbench'))
sys.path.insert(0, str(ROOT / '01_molecule_generation'))


@pytest.fixture
def launcher(tmp_path, monkeypatch):
    if os.name == 'nt':
        pytest.skip('Linux env shebang contract')
    bindir = tmp_path / 'selected/bin'; bindir.mkdir(parents=True)
    (bindir / 'python').symlink_to(sys.executable)
    executable = bindir / 'reinvent'
    executable.write_text('#!/usr/bin/env python\nimport json,os,sys\n'
                          'print(json.dumps({"python":sys.executable,"home":os.environ.get("PYTHONHOME"),'
                          '"path":os.environ.get("PYTHONPATH"),"gpu":os.environ.get("CUDA_VISIBLE_DEVICES")}))\n')
    executable.chmod(0o755)
    monkeypatch.setenv('PYTHONHOME', '/invalid/old-environment')
    monkeypatch.setenv('PYTHONPATH', '/invalid/old-packages')
    return executable


def probe_launcher(executable, environment):
    # Mirror the standalone adapter's documented overlay boundary and actually run
    # a harmless console script. A wrong PATH or restored PYTHONHOME fails this.
    process = subprocess.run([str(executable)], env={**os.environ, **environment},
                             capture_output=True, text=True, check=True)
    observed = json.loads(process.stdout)
    assert Path(observed['python']).parent == executable.parent
    assert observed['home'] == observed['path'] == ''
    return observed


def test_independent_reinvent_uses_selected_runtime(launcher, tmp_path, monkeypatch):
    pytest.importorskip('rdkit')
    from scientific import module_adapters
    from src.adapters import registry
    observed = []
    def run(**kwargs):
        observed.append(probe_launcher(launcher, kwargs['environment']))
        return {'report': {'completed': 1, 'valid': 0, 'unique': 0}}
    monkeypatch.setattr(registry, 'backend', lambda name: SimpleNamespace(run=run))
    result = module_adapters.execute(tmp_path, {
        'module': 'reinvent4', 'config': {'reinvent_python': str(launcher.parent/'python'),
                                       'reinvent_executable': str(launcher)},
        'parameters': {'budget': 1, 'seed': 17, 'scaffold': 'c1ccccc1[*:1]'}, 'inputs': {},
    }, {}, tmp_path/'out')
    assert len(observed) == 1 and result['report']['valid'] == 0


def test_legacy_cpu_reinvent_uses_selected_runtime(launcher, tmp_path, monkeypatch):
    chem = pytest.importorskip('rdkit.Chem')
    import pipeline_bridge
    from src.adapters import registry
    seed_dir = tmp_path/'peptide2mol'; seed_dir.mkdir()
    with chem.SDWriter(str(seed_dir/'public.sdf')) as writer:
        writer.write(chem.MolFromSmiles('Cc1ccccc1'))
    (seed_dir/'seeds.json').write_text(json.dumps([{'id': 'public', 'raw_sdf': 'peptide2mol/public.sdf'}]))
    out = tmp_path/'out'; out.mkdir()
    class ProbeComplete(Exception):
        pass
    def run(**kwargs):
        probe_launcher(launcher, kwargs['environment'])
        raise ProbeComplete()
    monkeypatch.setattr(registry, 'backend', lambda name: SimpleNamespace(run=run))
    with pytest.raises(ProbeComplete):
        pipeline_bridge.stage_reinvent4(tmp_path, {
            'reinvent_python': str(launcher.parent/'python'), 'reinvent_executable': str(launcher),
            'core_selections': [{'seed_id': 'public', 'key_atom_indices': [2]}], 'budget': 1, 'seed': 17,
        }, out)


@pytest.mark.parametrize('device', ['cpu', 'cuda'])
def test_research_reinvent_keeps_runtime_and_device_policy(launcher, tmp_path, monkeypatch, device):
    from research import generation
    monkeypatch.setattr(generation, 'disk_guard', lambda *args: None)
    monkeypatch.setattr(generation, 'detect', lambda path: {'gpus': [{'allocatable': True, 'uuid': 'synthetic-gpu'}]})
    monkeypatch.setattr(generation, 'gpu_lease', lambda gpu: nullcontext())
    inputs = {'executable': str(launcher)}
    seen = []
    def run(**kwargs):
        seen.append(kwargs)
        return probe_launcher(launcher, kwargs['environment'])
    result = generation.sample_reinvent(run, inputs, tmp_path, {'reinvent_device': device})
    assert result['gpu'] == ('' if device == 'cpu' else 'synthetic-gpu')
    assert seen[0]['device'] == ('cpu' if device == 'cpu' else 'cuda:0')
    assert inputs == {'executable': str(launcher)}  # No full environment enters a persisted batch signature.
