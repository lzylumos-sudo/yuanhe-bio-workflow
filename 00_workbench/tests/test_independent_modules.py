"""Independent module contracts; mocked engines are labelled and never scientific evidence."""
import base64
import copy
import json
from pathlib import Path
import sys
from unittest.mock import Mock
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'00_workbench'))
from webui import module_jobs as jobs, workflow as w
from scientific import prepared_md


@pytest.fixture
def isolated(tmp_path, monkeypatch):
    monkeypatch.setattr(jobs, 'JOBS', tmp_path/'module_jobs')
    monkeypatch.setattr(w, 'RUNS', tmp_path/'formal_runs')
    monkeypatch.setattr(w, '_active', {})
    return tmp_path


def upload(name, content):
    return {'name': name, 'base64': base64.b64encode(content).decode()}


def request(kind):
    spec = jobs.MODULES[kind]
    p = {f['key']: copy.deepcopy(f['default']) for f in spec['fields']}
    for key in ['review_reason', 'orientation_source']:
        if key in p:
            p[key] = 'Public fixture engineering review'
    if kind == 'reinvent4':
        p['scaffold'] = 'c1ccccc1[*:1]'
    if kind == 'interaction':
        p['ligand'] = 'BEN'
    files = {}
    for role, extensions in spec['files'].items():
        if extensions.endswith('?'):
            continue
        extension = extensions.split(',')[0]
        content = b'ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n' if extension == '.pdb' else b'unit-test-fixture'
        files[role] = upload(role+extension, content)
    return {'module': kind, 'parameters': p, 'files': files, 'acknowledged': True}


@pytest.mark.parametrize('kind', list(jobs.MODULES))
def test_every_module_saves_without_any_upstream(isolated, kind):
    result = jobs.create(request(kind))
    assert result['module'] == kind and not result['operations']
    assert result['candidate_qualified'] is False
    assert not w.RUNS.exists()
    assert jobs.verified(result['id'])[1]['schema'] == 'yuanhe.module-job/1'


@pytest.mark.parametrize('mutation,match', [
    (lambda r: r.update(acknowledged=False), 'Acknowledge'),
    (lambda r: r.update(module='shell'), 'Unknown'),
    (lambda r: r['parameters'].update(executable='/bin/sh'), 'Unknown'),
    (lambda r: r.update(files={}), 'required'),
    (lambda r: r['files']['receptor'].update(name='../run.py'), 'format'),
    (lambda r: r['files']['receptor'].update(base64='%%%'), 'encoding'),
    (lambda r: r['parameters'].update(normal=[0,0,0]), 'nonzero'),
    (lambda r: r['parameters'].update(production_ns=float('nan')), 'Finite'),
    (lambda r: r['parameters'].update(seed=True), 'Finite'),
])
def test_invalid_inputs_rejected_before_persisting(isolated, mutation, match):
    value = request('md'); mutation(value)
    with pytest.raises(ValueError, match=match):
        jobs.create(value)
    assert not jobs.JOBS.exists()


def test_input_tampering_and_job_ids(isolated):
    result = jobs.create(request('md'))
    d = jobs.folder(result['id'])
    (d/'inputs/receptor.pdb').write_text('changed')
    with pytest.raises(ValueError, match='input changed'):
        jobs.get(result['id'])
    for ident in ['../other', 'YH-1', 'MOD-'+('a'*32)+'/../other']:
        with pytest.raises(ValueError):
            jobs.folder(ident)


def test_source_change_blocks_start(isolated, monkeypatch):
    result = jobs.create(request('docking_control'))
    monkeypatch.setattr(w, 'source_snapshot', lambda: {})
    with pytest.raises(ValueError, match='Source changed'):
        jobs.launch(result['id'], {})


def test_md_needs_build_and_explicit_review(isolated):
    job = jobs.create(request('md')); d = jobs.folder(job['id'])
    with pytest.raises(ValueError, match='built system'):
        jobs.launch(job['id'], {'action': 'run', 'system_id': 'MD_123456789abc'})
    w.atomic(d/'md_systems/MD_123456789abc/system.json', {'id': 'MD_123456789abc', 'status': 'built'})
    with pytest.raises(ValueError, match='Review'):
        jobs.launch(job['id'], {'action': 'run', 'system_id': 'MD_123456789abc'})
    with pytest.raises(ValueError, match='analysis requires'):
        jobs.launch(job['id'], {'action': 'analyze', 'system_id': 'MD_123456789abc', 'system_reviewed': True, 'review_reason': 'test'})


def test_missing_environment_does_not_leave_running_job(isolated, monkeypatch):
    defaults = w.defaults()
    defaults['chem_python'] = '/missing/module/python'
    monkeypatch.setattr(w, 'defaults', lambda: defaults)
    lease = Mock(); monkeypatch.setattr(w, 'acquire_compute', lambda: lease)
    result = jobs.create(request('docking_control'))
    with pytest.raises(ValueError, match='environment'):
        jobs.launch(result['id'], {})
    assert jobs.get(result['id'])['operations'] == []
    lease.close.assert_called_once()


def test_mixed_reinvent_environment_rejected_before_worker(isolated, monkeypatch):
    from yuanhe import runtime
    cfg = w.defaults(); cfg['chem_python'] = sys.executable
    monkeypatch.setattr(w, 'defaults', lambda: cfg)
    monkeypatch.setattr(runtime, 'reinvent_environment_check', lambda _: {
        'ok': False, 'detail': 'REINVENT launcher uses an old environment'})
    lease = Mock(); monkeypatch.setattr(w, 'acquire_compute', lambda: lease)
    spawn = Mock(); monkeypatch.setattr(w, 'spawn_job', spawn)
    job = jobs.create(request('reinvent4'))
    with pytest.raises(ValueError, match='old environment'):
        jobs.launch(job['id'], {})
    assert jobs.get(job['id'])['operations'] == []
    spawn.assert_not_called()
    lease.close.assert_called_once()


def test_reference_contract_is_not_candidate_evidence(tmp_path):
    from scientific.reference_docking import reference_inputs
    from yuanhe_preflight.gate import validate_evidence, QualificationError
    receptor = ROOT/'02_docking/data/prepared/demo_1drf/receptor_prepared.pdb'
    ligand = ROOT/'02_docking/data/reference/FOL_crystal.sdf'
    proof = reference_inputs(receptor, ligand)
    assert proof['candidate_qualified'] is False
    assert proof['schema'] == 'yuanhe.reference-control/1'
    changed = tmp_path/'ligand.sdf'; changed.write_bytes(ligand.read_bytes()+b'changed')
    with pytest.raises(ValueError, match='pinned public'):
        reference_inputs(receptor, changed)
    with pytest.raises((ValueError, QualificationError)):
        validate_evidence(proof, 'CCO')


def test_docking_stops_if_strict_preflight_fails(tmp_path, monkeypatch):
    from scientific import module_adapters as adapters
    out = tmp_path/'out'; out.mkdir()
    monkeypatch.setattr(adapters, 'molecule', lambda _: (None, 'CCO'))
    monkeypatch.setattr(adapters, 'preflight', lambda *a, **k: {'risk': {'passed': 0}})
    import sbvs2.docking.vina_official as vina
    engine = Mock(side_effect=AssertionError('Must not dock')); monkeypatch.setattr(vina, 'run_vina_job', engine)
    result = adapters.execute(tmp_path, {'module':'docking','config':{},'parameters':{},
        'inputs': {'ligand':{'path':'x.sdf'},'receptor':{'path':'r.pdb'}}}, {}, out)
    assert result['docked'] == 0 and result['status'] == 'blocked_by_strict_preflight'
    engine.assert_not_called()


@pytest.mark.parametrize('damage', ['unfinished', 'step', 'checkpoint', 'lincs', 'time', 'missing_trajectory'])
def test_prepared_md_does_not_accept_partial_output(tmp_path, monkeypatch, damage):
    # Explicit command double exercises completion validation; no real MD is claimed.
    engine = tmp_path/'gmx'; engine.write_bytes(b'command-double')
    source = tmp_path/'input.tpr'; source.write_bytes(b'tpr-double')
    def capture(argv, folder, name, timeout=120):
        text = 'nsteps = 10\ninit-step = 0\ndt = 0.002\ntinit = 0\nintegrator = md\nnstxout-compressed = 5\n'
        if name == 'checkpoint_audit.txt':
            text = 'step = '+('5' if damage=='checkpoint' else '10')
        (folder/name).write_text(text)
        return text
    def execute(argv, cwd, **kwargs):
        log = 'Step Time\n'+('5' if damage=='step' else '10')+' '+('0.01' if damage=='time' else '0.02')+'\n'
        if damage != 'unfinished':
            log += 'Finished mdrun\n'
        if damage == 'lincs':
            log += 'LINCS WARNING\n'
        (cwd/'production.log').write_text(log)
        for ext in ['gro','cpt','xtc','edr']:
            if damage == 'missing_trajectory' and ext=='xtc':
                continue
            (cwd/('production.'+ext)).write_bytes(b'command-double')
        return Mock(returncode=0)
    monkeypatch.setattr(prepared_md, 'capture', capture)
    monkeypatch.setattr(prepared_md.subprocess, 'run', execute)
    with pytest.raises((RuntimeError, FileNotFoundError)):
        prepared_md.run(tmp_path, {'gmx':str(engine)}, {'device':'cpu','threads':1}, source)
    assert not (tmp_path/'prepared_md/completion.json').exists()

def test_thread_budget_and_clean_environment_reach_worker(isolated, monkeypatch):
    cfg = w.defaults(); cfg.update(md_python=sys.executable, chem_python=sys.executable)
    monkeypatch.setattr(w, 'defaults', lambda: cfg)
    monkeypatch.setenv('PYTHONHOME', '/wrong/python')
    monkeypatch.setenv('OMP_NUM_THREADS', '99')
    monkeypatch.setattr(w, 'acquire_compute', lambda: Mock())
    spawn = Mock(); monkeypatch.setattr(w, 'spawn_job', spawn)
    payload = request('md_tpr'); payload['parameters']['threads'] = 1
    job = jobs.create(payload)
    jobs.launch(job['id'], {'action':'run'})
    args = spawn.call_args.args
    env = args[-1]
    assert env['OMP_NUM_THREADS'] == env['YUANHE_MD_THREADS'] == '1'
    assert 'PYTHONHOME' not in env
    args[-2].close()


@pytest.mark.skipif(sys.platform != 'linux', reason='Linux worker process lifecycle')
def test_stop_cleans_detached_engine_before_worker_exit(tmp_path):
    import os, signal, subprocess, time
    script = tmp_path/'worker.py'
    script.write_text("""
import os,sys,subprocess,time
sys.path.insert(0,sys.argv[1])
from scientific.worker_processes import install
install()
p=subprocess.Popen([sys.executable,'-c','import time;time.sleep(90)'],start_new_session=True)
print(p.pid,flush=True)
time.sleep(90)
""")
    p = subprocess.Popen([sys.executable,str(script),str(ROOT/'00_workbench')],stdout=subprocess.PIPE,text=True,start_new_session=True)
    child = int(p.stdout.readline().strip())
    try:
        os.killpg(p.pid, signal.SIGTERM)
        assert p.wait(timeout=15) != 0
        path = Path(f'/proc/{child}/stat')
        assert not path.exists() or path.read_text().split(')')[1].split()[0] == 'Z'
    finally:
        if p.poll() is None:
            os.killpg(p.pid, signal.SIGKILL); p.wait()
        try:
            os.kill(child,signal.SIGKILL)
        except ProcessLookupError:
            pass

@pytest.mark.parametrize('scaffold', ['c1ccccc1[*]', '[*:1]CC[*:2]', 'invalid', 'CCO', '[*:1]'])
def test_reinvent_attachment_is_validated_before_job_creation(isolated, scaffold):
    payload=request('reinvent4');payload['parameters']['scaffold']=scaffold
    with pytest.raises(ValueError,match='labelled attachment'):
        jobs.create(payload)
    assert not jobs.JOBS.exists()
