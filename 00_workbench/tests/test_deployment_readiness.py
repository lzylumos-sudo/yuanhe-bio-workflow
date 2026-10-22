"""Deployment regressions; no target inference, docking or production MD."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import pytest
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / '00_workbench'))
from webui import readiness, server, studio, workflow


def test_creation_without_site_packages(tmp_path):
    script = r'''
import sys
from pathlib import Path
root,temporary=map(Path,sys.argv[1:])
sys.path.insert(0,str(root/'00_workbench'))
from webui import workflow as w
from yuanhe_preflight.config import defaults,validate
validate(defaults(root))
w.ROOT=temporary;w.RUNS=temporary/'runs'
rules=temporary/'02_docking/src/yuanhe_preflight/configs/risk_rules.v1.json'
rules.parent.mkdir(parents=True)
rules.write_bytes((root/'02_docking/src/yuanhe_preflight/configs/risk_rules.v1.json').read_bytes())
target=temporary/'03_virtual_screening/configs/targets/ecmurj.yaml'
target.parent.mkdir(parents=True)
target.write_bytes((root/'03_virtual_screening/configs/targets/ecmurj.yaml').read_bytes())
(temporary/'input.pdb').write_text('HEADER TEST\nEND\n')
run=w.create({'structure':'input.pdb'})
assert (w.run_dir(run['id'])/'config.json').is_file()
assert 'rdkit' not in sys.modules
assert not w._active
print('created without RDKit')
'''
    p = subprocess.run([sys.executable, '-S', '-c', script, str(ROOT), str(tmp_path)],capture_output=True,text=True)
    assert p.returncode == 0, p.stderr
    assert 'created without RDKit' in p.stdout


def test_relative_predocking_paths_are_frozen_against_project(tmp_path):
    from yuanhe_preflight.config import absolute_paths, defaults
    cfg=defaults(ROOT);cfg['route_reviews']=[{'search_artifact':'private/route.json'}]
    cfg['retro'].update(python='env/bin/python',model_path='private/model.onnx')
    fixed=absolute_paths(cfg,tmp_path)
    assert fixed['route_reviews'][0]['search_artifact']==str(tmp_path/'private/route.json')
    assert fixed['retro']['model_path']==str(tmp_path/'private/model.onnx')
    assert fixed['retro']['python']==str(tmp_path/'env/bin/python')
    assert cfg['route_reviews'][0]['search_artifact']=='private/route.json'


def test_readiness_reports_all_missing_before_generation(tmp_path,monkeypatch):
    monkeypatch.setattr(readiness,'probe',lambda *a,**k:{'ok':False,'error':'missing executable'})
    seen=[]
    result=readiness.collect(workflow.defaults(),tmp_path,publish=lambda r:seen.append(len(r['checks'])))
    assert result['status']=='blocked' and not result['all_runtime_checks_pass']
    assert {'ui_interpreter','chemistry','membrane_and_analysis','synthesis_path','gromacs_cuda_build','risk_rules','peptide_cuda'}<=result['checks'].keys()
    assert not result['checks']['synthesis_path']['ok']
    assert len(seen)>10


def test_check_uses_screen_interpreter_and_md_tools(monkeypatch):
    calls=[]
    def probe(argv,root,**kwargs):
        calls.append(argv);return {'ok':False,'error':'test inventory'}
    monkeypatch.setattr(readiness,'probe',probe)
    cfg=workflow.defaults();cfg.update(chem_python='/test/chem/python',screen_python='/test/screen/python',md_python='/test/md/python')
    readiness.collect(cfg,ROOT)
    assert any(a[0]=='/test/screen/python' and 'posebusters' in a[-1] for a in calls)
    assert any(a[0]=='/test/md/python' and 'MMPBSA.py' in a[-1] and 'tleap' in a[-1] for a in calls)


def test_md_environment_does_not_depend_on_service_shell(tmp_path,monkeypatch):
    prefix=tmp_path/'md';(prefix/'bin').mkdir(parents=True);(prefix/'dat/leap').mkdir(parents=True)
    monkeypatch.setenv('PATH','/usr/bin');monkeypatch.setenv('PYTHONPATH','/wrong/scientific/environment');monkeypatch.setenv('PYTHONHOME','/wrong/python')
    env=readiness.runtime_env(prefix/'bin/python')
    assert env['PATH'].split(os.pathsep)[0]==str(prefix/'bin')
    assert env['AMBERHOME']==str(prefix)
    assert 'PYTHONPATH' not in env and 'PYTHONHOME' not in env


def test_prepare_refuses_missing_deployment_before_structure_work(tmp_path,monkeypatch):
    import pipeline_bridge
    def blocked(*args):raise RuntimeError('Deployment readiness blocked: chemistry')
    monkeypatch.setattr(readiness,'require_ready',blocked)
    with pytest.raises(RuntimeError,match='readiness blocked'):pipeline_bridge.stage_prepare(tmp_path,{},tmp_path/'prepare')
    assert not (tmp_path/'prepare').exists()


def test_independent_diagnostics_without_run(tmp_path,monkeypatch):
    monkeypatch.setattr(workflow,'ROOT',tmp_path);monkeypatch.setattr(workflow,'_active',{})
    monkeypatch.setattr(studio,'_diagnostics_state',{'status':'idle','checks':{}})
    class Lease:
        closed=False
        def close(self):self.closed=True
    lease=Lease();monkeypatch.setattr(workflow,'acquire_compute',lambda:lease)
    def collect(cfg,root,ui_python,publish):
        publish({'status':'blocked','checks':{'ui_interpreter':{'ok':False,'required':True}},'all_runtime_checks_pass':False})
    monkeypatch.setattr(readiness,'collect',collect)
    code,_,body=server.dispatch('POST','/api/studio/diagnostics',{},{});assert code==202
    deadline=time.monotonic()+3
    while not lease.closed and time.monotonic()<deadline:time.sleep(.01)
    code,_,body=server.dispatch('GET','/api/studio/diagnostics',{},None)
    data=json.loads(body)
    assert code==200 and data['status']=='blocked'
    assert lease.closed and Path(data['report']).is_file()
    assert not (tmp_path/'runs').exists()


def test_diagnostic_cannot_accept_arbitrary_commands():
    with pytest.raises(ValueError,match='Unknown diagnostics setting'):studio.diagnostics_start({'command':'anything'})


def test_checked_launcher_fails_early_without_rdkit(monkeypatch,capsys):
    spec=importlib.util.spec_from_file_location('checked_launcher',ROOT/'00_workbench/run_webui.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    monkeypatch.setenv('YUANHE_UI_PYTHON',sys.executable)
    monkeypatch.setattr(sys,'argv',['run_webui.py','--no-browser'])
    monkeypatch.setattr(readiness,'ui_check',lambda *a:{'ok':False,'error':"No module named rdkit"})
    assert mod.main()==2
    assert 'WORKBENCH_DEPLOYMENT.md' in capsys.readouterr().err
