"""Real Linux processes exercise scheduling; model inference is not simulated acceptance."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from webui import workflow as w, studio, server
from webui.scheduling import acquire, receptor_environment

pytestmark = pytest.mark.skipif(sys.platform != 'linux', reason='Linux/WSL process-group and flock contract')


def eventually(check, timeout=8):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if check():
            return
        time.sleep(.02)
    assert check(), 'Timed out waiting for worker state'


@pytest.mark.parametrize('first,second', [('compute','receptor'), ('receptor','compute')])
@pytest.mark.parametrize('same_run', [True, False])
def test_independent_lanes_and_bounded_capacity(tmp_path, first, second, same_run):
    run = tmp_path/'a'
    with acquire(tmp_path, run, first):
        with acquire(tmp_path, run if same_run else tmp_path/'b', second):
            for lane in (first, second):
                with pytest.raises(RuntimeError, match='capacity'):
                    acquire(tmp_path, tmp_path/'third', lane)
            with pytest.raises(RuntimeError, match='background'):
                acquire(tmp_path)
    with acquire(tmp_path):
        for lane in ('compute','receptor'):
            with pytest.raises(RuntimeError, match='background'):
                acquire(tmp_path, run, lane)
    with acquire(tmp_path, run, 'compute'):
        pass


def test_failed_acquisition_releases_partial_locks(tmp_path):
    with acquire(tmp_path, tmp_path/'a', 'compute'):
        with pytest.raises(RuntimeError):
            acquire(tmp_path, tmp_path/'b', 'compute', stage='docking')
        # The failed contender must not leave its review lock behind.
        with acquire(tmp_path, tmp_path/'b', review='docking'):
            pass
    with acquire(tmp_path):
        pass


def test_all_leases_survive_parent_close_and_reach_grandchild(tmp_path):
    ui = str(Path(w.__file__).resolve().parents[1])
    child = "from scientific.artifacts import lease_fds; import subprocess,sys; fds=lease_fds(); print(len(fds),flush=True); subprocess.run([sys.executable,'-c','import sys; sys.stdin.read()'],pass_fds=fds)"
    lease = acquire(tmp_path, tmp_path/'run', 'compute', stage='peptide2mol')
    proc = subprocess.Popen([sys.executable, '-c', child], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            text=True, pass_fds=lease.fds(), env={**os.environ, 'PYTHONPATH': ui, **lease.environment()})
    try:
        assert int(proc.stdout.readline()) == len(lease.fds())
        lease.close()
        with pytest.raises(RuntimeError):
            acquire(tmp_path, tmp_path/'other', 'compute')
        with pytest.raises(RuntimeError):
            acquire(tmp_path, tmp_path/'run', review='peptide2mol')
        with acquire(tmp_path, tmp_path/'run', 'receptor'):
            pass
    finally:
        proc.communicate('', timeout=8)
        lease.close()
    with acquire(tmp_path):
        pass


WORKER = '''import argparse,json,os,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--run-id');p.add_argument('--stage');p.add_argument('--operation');a=p.parse_args()
d=Path(__file__).parent/'runs'/a.run_id
out=d/a.stage if a.stage else d/'operations'/a.operation
record=out/('stage.json' if a.stage else 'operation.json')
(out/'ready').write_text('ready')
print('worker '+(a.stage or a.operation),flush=True)
while not (out/'release').exists():time.sleep(.02)
if a.operation:
 (out/'receptor.pdb').write_text('public test receptor placeholder')
 import hashlib
 result={'receptor':str((out/'receptor.pdb').relative_to(d)),'sha256':hashlib.sha256((out/'receptor.pdb').read_bytes()).hexdigest(),'pH':7.4}
 (out/'result.json').write_text(json.dumps(result))
r=json.loads(record.read_text())
if r['status']=='running':
 r.update(status='completed',ended=time.time())
 temp=record.with_suffix('.done');temp.write_text(json.dumps(r));temp.replace(record)
'''


@pytest.fixture
def runtime(tmp_path, monkeypatch):
    monkeypatch.setattr(w, 'ROOT', tmp_path)
    monkeypatch.setattr(w, 'UI', tmp_path)
    monkeypatch.setattr(w, 'RUNS', tmp_path/'runs')
    monkeypatch.setattr(w, '_active', {})
    for script in ('pipeline_bridge.py','studio_worker.py'):
        (tmp_path/script).write_text(WORKER)
    for rid in ('a','b'):
        d=w.run_dir(rid)
        cfg=w.defaults();cfg.update(research={},predocking={},chem_python=sys.executable,prep_python=sys.executable)
        w.atomic(d/'config.json',cfg)
        w.atomic(d/'run.json',{'config_sha256':w.digest(d/'config.json')})
        w.atomic(d/'prepare/stage.json',{'status':'completed','config_sha256':w.digest(d/'config.json'),'artifacts':{}})
    yield tmp_path
    for rid in ('a','b'):
        w.cancel(rid)
    eventually(lambda:not w._active)


def launch_pair(rid='a', operation_run=None):
    w.start(rid,'peptide2mol')
    op=studio.launch(operation_run or rid,'receptor',{})['id']
    out=w.run_dir(operation_run or rid)/'operations'/op
    eventually(lambda:(out/'ready').exists() and (w.run_dir(rid)/'peptide2mol/ready').exists())
    return op,out


def test_generation_and_receptor_can_overlap_apply_and_finish_independently(runtime):
    op,out=launch_pair()
    d=w.run_dir('a');stage=d/'peptide2mol'
    frozen=w.digest(stage/'review.json')
    assert w.status('a')['stages'][1]['status']=='running'
    assert studio.operations('a')[0]['status']=='running'
    assert len(w._active)==2
    with pytest.raises(RuntimeError):w.start('b','prepare')
    with pytest.raises(RuntimeError):studio.launch('b','receptor',{})
    with pytest.raises(RuntimeError):w.review('a','peptide2mol',{'peptide_python':'changed'})
    assert w.review('a','docking',{'docking_box':{'center':[0,0,0],'size':[20]*3}})['saved']
    (out/'release').touch()
    eventually(lambda:('a','operation',op) not in w._active)
    assert ('a','stage','peptide2mol') in w._active
    assert studio.apply_receptor('a',{'operation':op,'reviewed':True,'box':{'center':[0,0,0],'size':[20]*3}})['saved']
    assert w.digest(stage/'review.json')==frozen
    assert (stage/'process.log').read_text().strip()=='worker peptide2mol'
    (stage/'release').touch()
    eventually(lambda:not w._active)
    assert w.read(stage/'stage.json')['status']=='completed'
    assert w.read(out/'operation.json')['status']=='completed'


@pytest.mark.parametrize('selector', ['operation','stage','all'])
def test_cancel_is_scoped_and_stop_all_reaps_both_groups(runtime, selector):
    op,out=launch_pair()
    pstage=w._active[('a','stage','peptide2mol')]
    pop=w._active[('a','operation',op)]
    if selector=='operation':
        code,_,body=server.dispatch('POST','/api/studio/run/a/cancel-operation',{}, {'operation':op})
        assert code==200 and json.loads(body)['cancelled']
        eventually(lambda:pop.poll() is not None)
        assert pstage.poll() is None
        assert w.read(w.run_dir('a')/'peptide2mol/stage.json')['status']=='running'
    elif selector=='stage':
        code,_,body=server.dispatch('POST','/api/project/run/a/cancel',{}, {'stage':'peptide2mol'})
        assert code==200 and json.loads(body)['cancelled']
        eventually(lambda:pstage.poll() is not None)
        assert pop.poll() is None
        assert w.read(out/'operation.json')['status']=='running'
    else:
        assert w.cancel('a')['cancelled']
        eventually(lambda:pstage.poll() is not None and pop.poll() is not None)
        assert w.read(out/'operation.json')['status']=='cancelled'


def test_cancel_one_run_leaves_other_run_alive(runtime):
    op,out=launch_pair(operation_run='b')
    proc=w._active[('b','operation',op)]
    w.cancel('a')
    assert proc.poll() is None and studio.operations('b')[0]['status']=='running'


def test_restart_detects_each_worker_and_locks_still_block_duplicates(runtime):
    op,out=launch_pair()
    original=w._active
    w._active={}
    try:
        assert w.status('a')['stages'][1]['status']=='running'
        assert studio.operations('a')[0]['status']=='running'
        with pytest.raises(RuntimeError):studio.launch('a','receptor',{})
        w.cancel('a',operation=op)
        eventually(lambda:original[('a','operation',op)].poll() is not None)
        assert original[('a','stage','peptide2mol')].poll() is None
    finally:
        w._active=original
        # A watcher may have already finished while the registry was replaced.
        for key,proc in list(original.items()):
            if proc.poll() is not None:original.pop(key,None)


def test_stale_operation_is_not_hidden_by_live_generation(runtime):
    w.start('a','peptide2mol')
    w.atomic(w.run_dir('a')/'operations/old/operation.json',{'status':'running','kind':'receptor','pid':2147483647,'process_start':'old'})
    assert studio.operations('a')[0]['status']=='interrupted'


def test_failed_receptor_launch_releases_slot_without_stopping_generation(runtime):
    w.start('a','peptide2mol')
    d=w.run_dir('b');cfg=w.read(d/'config.json');cfg['prep_python']='/missing/python'
    w.atomic(d/'config.json',cfg);w.atomic(d/'run.json',{'config_sha256':w.digest(d/'config.json')})
    w.atomic(d/'prepare/stage.json',{'status':'completed','config_sha256':w.digest(d/'config.json'),'artifacts':{}})
    with pytest.raises(FileNotFoundError):studio.launch('b','receptor',{})
    with w.acquire_compute('b',lane='receptor'):pass
    assert w._active[('a','stage','peptide2mol')].poll() is None


def test_maintenance_and_export_wait_for_both_lanes(runtime):
    op,out=launch_pair()
    with pytest.raises(ValueError):studio.export_run('a',{})
    with pytest.raises(ValueError):studio.diagnostics_start({})
    with pytest.raises(RuntimeError):studio.launch('b','environment',{})


def test_cpu_receptor_has_no_cuda_and_bounded_threads():
    env=receptor_environment()
    assert env['CUDA_VISIBLE_DEVICES']=='' and env['OPENMM_DEFAULT_PLATFORM']=='CPU'
    assert 1<=int(env['OPENMM_CPU_THREADS'])<=2
    assert env['OMP_NUM_THREADS']==env['MKL_NUM_THREADS']==env['OPENMM_CPU_THREADS']


@pytest.mark.parametrize('value', ['', '../other', 'a/b'])
def test_invalid_operation_cancel_never_cancels_whole_run(runtime,value):
    with pytest.raises(ValueError):w.cancel('a',operation=value)


def test_null_stage_cancel_is_rejected(runtime):
    with pytest.raises(ValueError):server.api_project_cancel({}, {'stage':None}, 'a')
