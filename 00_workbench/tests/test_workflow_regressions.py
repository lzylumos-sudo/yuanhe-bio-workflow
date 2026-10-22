"""Engineering regressions only. No target docking or model inference."""
import importlib,json,os,sys,threading,urllib.request,urllib.error
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from webui import workflow,server,jobs

@pytest.fixture
def sandbox(tmp_path,monkeypatch):
 # Include the actual packaged rules; run creation must still snapshot them.
 rules=Path(__file__).resolve().parents[2]/'02_docking/src/yuanhe_preflight/configs/risk_rules.v1.json'
 target=tmp_path/'02_docking/src/yuanhe_preflight/configs/risk_rules.v1.json'
 target.parent.mkdir(parents=True);target.write_bytes(rules.read_bytes())
 profile=Path(__file__).resolve().parents[2]/'03_virtual_screening/configs/targets/ecmurj.yaml'
 local_profile=tmp_path/'03_virtual_screening/configs/targets/ecmurj.yaml'
 local_profile.parent.mkdir(parents=True);local_profile.write_bytes(profile.read_bytes())
 monkeypatch.setattr(workflow,'ROOT',tmp_path);monkeypatch.setattr(workflow,'RUNS',tmp_path/'runs')
 (tmp_path/'input.pdb').write_text('HEADER    TEST ONLY\nEND\n')
 return tmp_path

def test_path_escape_and_identical_chains(sandbox):
 with pytest.raises(ValueError):workflow.create({'structure':'../outside.pdb'})
 with pytest.raises(ValueError):workflow.create({'structure':'input.pdb','peptide_chain':'A'})
 with pytest.raises(ValueError):workflow.run_dir('../escape')

def test_unknown_and_bad_budget(sandbox):
 with pytest.raises(ValueError):workflow.create({'structure':'input.pdb','arbitrary_command':'sh'})
 with pytest.raises(ValueError):workflow.create({'structure':'input.pdb','budget':True})

def test_create_does_not_launch(sandbox):
 r=workflow.create({'structure':'input.pdb'});d=workflow.run_dir(r['id'])
 assert (d/'input.pdb').read_bytes()==(sandbox/'input.pdb').read_bytes()
 assert all(s['status']=='pending' for s in workflow.status(r['id'])['stages'])
 assert not workflow._active

def test_upstream_hash_tampering(sandbox):
 r=workflow.create({'structure':'input.pdb'});d=workflow.run_dir(r['id']);(d/'prepare').mkdir();p=d/'prepare/x';p.write_text('valid')
 workflow.atomic(d/'prepare/stage.json',{'status':'completed','config_sha256':workflow.digest(d/'config.json'),'artifacts':{'prepare/x':workflow.digest(p)}})
 assert workflow.require(d,'prepare')['status']=='completed'
 p.write_text('changed')
 with pytest.raises(RuntimeError,match='modified'):workflow.require(d,'prepare')

def test_block_out_of_order(sandbox):
 rid=workflow.create({'structure':'input.pdb'})['id']
 with pytest.raises(RuntimeError,match='(?i)upstream'):workflow.start(rid,'docking')
 assert not (workflow.run_dir(rid)/'docking').exists()

def test_interrupted_is_not_success(sandbox):
 rid=workflow.create({'structure':'input.pdb'})['id'];d=workflow.run_dir(rid)
 workflow.atomic(d/'prepare/stage.json',{'status':'running'})
 assert workflow.status(rid)['stages'][0]['status']=='interrupted'

def test_symlink_escape(tmp_path):
 outside=tmp_path/'secret';outside.write_text('do not expose');base=tmp_path/'public';base.mkdir();(base/'link').symlink_to(outside)
 with pytest.raises(server.HttpError):server._safe_join(str(base),'link')

def test_jobs_path_escape():
 for path in ['/tmp/escape','../outside','\\\\host\\outside']:
  with pytest.raises(ValueError):jobs._validate('pipeline',{'out':path})

def test_log_cursor_survives_truncation(monkeypatch):
 monkeypatch.setattr(jobs,'MAX_LINES',2);j={'lines':[]}
 for x in range(5):jobs._append(j,str(x))
 assert [x['i'] for x in j['lines']]==[3,4] and j['next_line']==5

def test_cancel_route_accepts_body(monkeypatch):
 monkeypatch.setattr(jobs,'cancel',lambda x:True)
 status,_,payload=server.dispatch('POST','/api/jobs/abc/cancel',{}, {})
 assert status==200 and json.loads(payload)['cancelled']

def test_http_origin_body_and_host_guards():
 httpd=server.Server(('127.0.0.1',0),server.Handler);threading.Thread(target=httpd.serve_forever,daemon=True).start()
 url=f'http://127.0.0.1:{httpd.server_address[1]}'
 try:
  assert json.load(urllib.request.urlopen(url+'/api/health'))['ok']
  for body,headers,code in [(b'{}',{'Origin':'https://attacker.example'},403),(b'[]',{},400),(b'{}',{'Host':'attacker.example'},403)]:
   req=urllib.request.Request(url+'/api/project/runs',data=body,headers={'Content-Type':'application/json',**headers})
   with pytest.raises(urllib.error.HTTPError) as e:urllib.request.urlopen(req)
   assert e.value.code==code
 finally:httpd.shutdown();httpd.server_close()

def test_review_cannot_change_completed_stage(sandbox):
 rid=workflow.create({'structure':'input.pdb'})['id'];d=workflow.run_dir(rid)
 assert workflow.review(rid,'peptide2mol',{'peptide_python':'/tmp/python'})['saved']
 workflow.atomic(d/'peptide2mol/stage.json',{'status':'completed'})
 with pytest.raises(RuntimeError):workflow.review(rid,'peptide2mol',{'peptide_python':'other'})

def test_review_can_fix_blocked_stage(sandbox):
 rid=workflow.create({'structure':'input.pdb'})['id'];d=workflow.run_dir(rid)
 workflow.atomic(d/'peptide2mol/stage.json',{'status':'blocked'})
 assert workflow.review(rid,'peptide2mol',{'peptide_python':'/tmp/python'})['saved']
 with pytest.raises(ValueError):workflow.review(rid,'peptide2mol',{'command':'sh'})

def test_failed_retry_archives_previous_attempt(sandbox,monkeypatch):
 rid=workflow.create({'structure':'input.pdb'})['id'];d=workflow.run_dir(rid)
 workflow.atomic(d/'prepare/stage.json',{'status':'failed','error':'old failure'})
 (d/'prepare/old.log').write_text('previous attempt')
 class Finished:
  pid=12345
  def wait(self):return 1
 monkeypatch.setattr(workflow.subprocess,'Popen',lambda *a,**k:Finished())
 workflow.start(rid,'prepare')
 assert len(list((d/'attempts').glob('prepare_*/old.log')))==1

def test_pose_viewer_uses_recorded_relative_sdf(sandbox):
 from rdkit import Chem
 from rdkit.Chem import AllChem
 rid=workflow.create({'structure':'input.pdb'})['id'];d=workflow.run_dir(rid);stage=d/'peptide2mol';stage.mkdir()
 mol=Chem.AddHs(Chem.MolFromSmiles('CO'));assert AllChem.EmbedMolecule(mol,randomSeed=7)==0
 with Chem.SDWriter(str(stage/'seed.sdf')) as w:w.write(mol)
 workflow.atomic(stage/'seeds.json',[{'id':'seed_00000','raw_sdf':'peptide2mol/seed.sdf'}])
 code,_,body=server.api_project_pose({'id':['seed_00000']},rid,'peptide2mol')
 assert code==200 and len(json.loads(body)['atoms'])==2
