"""GUI service and membrane contracts, without target inference or MD."""
import importlib,json,sys,threading,hashlib
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from webui import workflow as w,studio,server
from scientific import membrane

@pytest.fixture
def case(tmp_path,monkeypatch):
 monkeypatch.setattr(w,'ROOT',tmp_path);monkeypatch.setattr(w,'RUNS',tmp_path/'runs')
 (tmp_path/'test.pdb').write_text('HEADER TEST\nATOM      1  CA  ALA A   3       1.000   2.000   3.000  1.00 20.00           C\nEND\n')
 rid=w.create({'structure':'test.pdb'})['id'];return tmp_path,rid,w.run_dir(rid)

def complete(d,stage):
 w.atomic(d/stage/'stage.json',{'status':'completed','config_sha256':w.digest(d/'config.json'),'artifacts':{}})

def test_upload_snapshots_and_rejects_nonpdb(case):
 root,rid,d=case
 got=studio.upload({'name':'../../input.pdb','text':(root/'test.pdb').read_text()})
 assert (root/got['structure']).is_file() and got['structure'].startswith('inputs/')
 with pytest.raises(ValueError):studio.upload({'name':'a.py','text':'ATOM  '})
 with pytest.raises(ValueError):studio.upload({'name':'a.pdb','text':'not a structure'})

def test_upload_limit(case):
 with pytest.raises(ValueError):studio.upload({'name':'a.pdb','text':'ATOM  '+'x'*(8*1024*1024)})

def test_finite_box_and_screen_budget(case):
 _,rid,d=case
 with pytest.raises(ValueError):w.review(rid,'docking',{'docking_box':{'center':[0,float('nan'),0],'size':[20]*3}})
 with pytest.raises(ValueError):w.review(rid,'screening',{'screening_budget':'../../other'})
 assert w.review(rid,'screening',{'screening_budget':'full'})['saved']

def test_compute_lock_is_cross_instance(case):
 lease=w.acquire_compute()
 try:
  with pytest.raises(RuntimeError):w.acquire_compute()
 finally:lease.close()
 again=w.acquire_compute();again.close()

def test_core_numbers_and_save(case):
 from rdkit import Chem
 from rdkit.Chem import AllChem
 _,rid,d=case;complete(d,'peptide2mol')
 mol=Chem.AddHs(Chem.MolFromSmiles('CCO'));AllChem.EmbedMolecule(mol,randomSeed=17)
 with Chem.SDWriter(str(d/'peptide2mol/seed.sdf')) as f:f.write(mol)
 w.atomic(d/'peptide2mol/seeds.json',[{'id':'seed_0','smiles':'CCO','raw_sdf':'peptide2mol/seed.sdf'}])
 info=studio.core(rid,'seed_0');assert [a['index'] for a in info['atoms']]==[1,2,3]
 with pytest.raises(ValueError):studio.save_core(rid,{'seed_id':'seed_0','indices':[0],'reason':'test'})
 with pytest.raises(ValueError):studio.save_core(rid,{'seed_id':'seed_0','indices':[1],'reason':''})
 assert studio.save_core(rid,{'seed_id':'seed_0','indices':[1,2],'reason':'Preserve continuous carbon core'})['selections'][0]['key_atom_indices']==[1,2]

def test_relative_receptor_hash_guard(case):
 _,rid,d=case;out=d/'operations/op1';out.mkdir(parents=True);(out/'receptor.pdb').write_text('test')
 w.atomic(out/'operation.json',{'status':'completed'});w.atomic(out/'result.json',{'receptor':'operations/op1/receptor.pdb','sha256':w.digest(out/'receptor.pdb'),'pH':7.4})
 assert studio.apply_receptor(rid,{'operation':'op1','reviewed':True,'box':{'center':[0,0,0],'size':[20,20,20]}})['saved']
 assert not Path(w.read(d/'reviews/docking.json')['docking_receptor']).is_absolute()
 (out/'receptor.pdb').write_text('changed')
 with pytest.raises(ValueError,match='changed'):studio.apply_receptor(rid,{'operation':'op1','reviewed':True,'box':{}})

@pytest.mark.parametrize('normal',[[0,0,1],[0,0,-1],[1,0,0],[1,2,3]])
def test_membrane_rotation_preserves_geometry(normal):
 r=membrane.rotation(normal);n=np.array(normal)/np.linalg.norm(normal)
 assert np.allclose(r@n,[0,0,1]);assert np.allclose(r.T@r,np.eye(3));assert np.linalg.det(r)==pytest.approx(1)

@pytest.mark.parametrize('normal',[[0,0,0],[float('nan'),0,1],[0,1]])
def test_invalid_membrane_rotation(normal):
 with pytest.raises(ValueError):membrane.rotation(normal)

def test_restraints_are_local_indices(tmp_path):
 top=tmp_path/'test.top';top.write_text('[ moleculetype ]\nProtein 3\n[ atoms ]\n1 CT 1 ALA CA 1 0 12.01\n2 H 1 ALA H 1 0 1.008\n[ bonds ]\n1 2 1\n[ moleculetype ]\nLIG 3\n[ atoms ]\n1 C 1 LIG C1 1 0 12.01\n2 H 1 LIG H1 1 0 1.008\n[ system ]\ntest\n[ molecules ]\nProtein 1\nLIG 1\n')
 membrane.add_restraints(top);s=top.read_text()
 assert s.count('[ position_restraints ]')==2
 assert s.count('1 1 POSRES_FC')==2 and '2 1 POSRES_FC' not in s
 assert s.rfind('#endif')<s.index('[ system ]')

def test_equilibration_and_production_protocol():
 p={'temperature_K':310,'production_ns':2,'equilibration_ns':1}
 assert 'gen-vel = yes' in membrane.mdp('nvt',p,17)
 assert 'gen-vel = no' in membrane.mdp('production',p,17)
 assert 'pcoupltype = semiisotropic' in membrane.mdp('production',p,17)
 assert 'POSRES' not in membrane.mdp('production',p,17)
 assert 'POSRES_FC=100' in membrane.mdp('npt_soft',p,17)
 assert 'nsteps = 1000000' in membrane.mdp('production',p,17)

def test_md_refuses_changed_input_before_engine(case,monkeypatch):
 _,rid,d=case;folder=d/'md_systems/MD_test';folder.mkdir(parents=True);(folder/'system.top').write_text('before')
 w.atomic(folder/'system.json',{'status':'built','files':{'system.top':w.digest(folder/'system.top')}})
 (folder/'system.top').write_text('after')
 monkeypatch.setattr(membrane,'execute',lambda *a,**k:pytest.fail('engine must not execute'))
 with pytest.raises(ValueError,match='changed'):membrane.run(d,{}, {'system_id':'MD_test','system_reviewed':True},d)

def test_md_requires_review_and_safe_identifier(case):
 _,rid,d=case
 with pytest.raises(ValueError):membrane.system_dir(d,'../escape')
 folder=d/'md_systems/MD_test';w.atomic(folder/'system.json',{'status':'built'})
 with pytest.raises(ValueError,match='(?i)review'):membrane.run(d,{}, {'system_id':'MD_test'},d)

def test_no_md_or_receptor_work_before_prerequisite(case):
 _,rid,d=case
 with pytest.raises(RuntimeError,match='(?i)upstream'):studio.launch(rid,'md_build',{})
 assert not (d/'operations').exists()

def test_environment_rejects_arbitrary_keys(case):
 with pytest.raises(ValueError):studio.environment_save({'shell_command':'rm'})
 assert studio.environment_save({'md_python':'/opt/example/python'})['saved']

def test_resumed_tpr_tamper_is_rejected(case,monkeypatch):
 _,rid,d=case;folder=d/'md_systems/MD_test';folder.mkdir(parents=True)
 (folder/'system.top').write_text('stable');(folder/'em.tpr').write_text('changed')
 w.atomic(folder/'system.json',{'status':'built','files':{'system.top':w.digest(folder/'system.top')},'parameters':{'temperature_K':310,'production_ns':1,'equilibration_ns':1}})
 w.atomic(folder/'progress.json',{'phases':{'em':{'status':'running','tpr_sha256':'different'}}})
 calls=[];monkeypatch.setattr(membrane,'execute',lambda argv,*a,**k:calls.append(argv))
 with pytest.raises(ValueError,match='TPR'):membrane.run(d,{'gmx':sys.executable,'seed':17},{'system_id':'MD_test','system_reviewed':True},d)
 assert calls==[[sys.executable,'--version']]

def test_export_records_omissions_and_avoids_recursive_archives(case):
 import zipfile
 _,rid,d=case
 (d/'production.xtc').write_bytes(b'trajectory')
 got=studio.export_run(rid,{})
 path=d/'exports'/got['url'].split('/')[-1].split('?')[0]
 with zipfile.ZipFile(path) as z:
  m=json.loads(z.read('EXPORT_MANIFEST.json'))
  assert m['omitted'][0]['path']=='production.xtc'
  assert 'input.pdb' in m['files']
 got=studio.export_run(rid,{'full':True})
 path=d/'exports'/got['url'].split('/')[-1].split('?')[0]
 with zipfile.ZipFile(path) as z:
  assert 'production.xtc' in z.namelist()
  assert not any(n.startswith('exports/') for n in z.namelist())

def test_export_blocked_by_other_ui_compute(case):
 _,rid,d=case;lease=w.acquire_compute()
 try:
  with pytest.raises(RuntimeError):studio.export_run(rid,{})
 finally:lease.close()

def test_success_exit_does_not_mean_completed_simulation(tmp_path):
 (tmp_path/'production.gro').write_text('structure')
 (tmp_path/'production.xtc').write_bytes(b'partial')
 (tmp_path/'production.log').write_text('           Step           Time\n            500          1.000\nFinished mdrun\n')
 with pytest.raises(RuntimeError,match='before the planned'):membrane.check_completion(tmp_path,'production','nsteps = 1000\n')
 (tmp_path/'production.log').write_text('           Step           Time\n           1000          2.000\nFinished mdrun\n')
 membrane.check_completion(tmp_path,'production','nsteps = 1000\n')
 (tmp_path/'em.log').write_text('Steepest Descents did not converge in 5000 steps')
 with pytest.raises(RuntimeError,match='convergence'):membrane.check_completion(tmp_path,'em','nsteps = 5000\n')

def test_cancel_orphan_record_allows_explicit_retry(case):
 _,rid,d=case
 w.atomic(d/'prepare/stage.json',{'status':'running','pid':2147483647,'process_start':'never'})
 assert w.status(rid)['stages'][0]['status']=='interrupted'
 assert w.cancel(rid)['cancelled']
 assert w.read(d/'prepare/stage.json')['status']=='cancelled'

def test_live_record_rejects_pid_reuse(case,monkeypatch):
 _,rid,d=case
 monkeypatch.setattr(w,'process_identity',lambda pid:{'pid':pid,'process_start':'new','argv':['python',str(w.UI/'pipeline_bridge.py'),'--run-id',rid]})
 assert not w.live_record({'pid':1234,'process_start':'old'},rid)

def test_review_update_preserves_bound_receptor_hash(case):
 _,rid,d=case
 w.review(rid,'docking',{'receptor_review':{'sha256':'abc','reviewed':True}})
 w.review(rid,'docking',{'vina':'/bin/vina'})
 assert w.read(d/'reviews/docking.json')['receptor_review']['sha256']=='abc'

def test_md_checkpoint_resume_skips_finished_phases(case,monkeypatch):
 _,rid,d=case;folder=d/'md_systems/MD_test';folder.mkdir(parents=True)
 (folder/'system.top').write_text('stable')
 w.atomic(folder/'system.json',{'id':'MD_test','status':'built','files':{'system.top':w.digest(folder/'system.top')},'parameters':{'temperature_K':310,'production_ns':.01,'equilibration_ns':.1}})
 calls=[];interrupted=[False]
 def engine(argv,cwd,*args,**kwargs):
  calls.append(argv)
  if 'grompp' in argv:(folder/argv[argv.index('-o')+1]).write_text('tpr');return
  if 'mdrun' not in argv:return
  phase=argv[argv.index('-deffnm')+1]
  (folder/(phase+'.cpt')).write_text('checkpoint')
  if phase=='npt' and not interrupted[0]:
   interrupted[0]=True;(folder/(phase+'.log')).write_text('partial');raise RuntimeError('simulated interruption')
  recipe=(folder/(phase+'.mdp')).read_text()
  import re
  n=re.search(r'nsteps = (\d+)',recipe).group(1)
  (folder/(phase+'.log')).write_text('converged to Fmax < 1000' if phase=='em' else f'    Step    Time\n    {n}    10.000\n')
  (folder/(phase+'.gro')).write_text('final')
  if phase=='production':(folder/'production.xtc').write_text('trajectory')
 monkeypatch.setattr(membrane,'execute',engine)
 monkeypatch.setattr(membrane,'checkpoint_probe',lambda *a,**k:{'status':'mocked for orchestration only'})
 c={'gmx':sys.executable,'seed':17};p={'system_id':'MD_test','system_reviewed':True}
 with pytest.raises(RuntimeError,match='simulated'):membrane.run(d,c,p,d)
 assert membrane.run(d,c,p,d)['status']=='completed'
 runs=[a for a in calls if 'mdrun' in a]
 assert sum(a[a.index('-deffnm')+1]=='em' for a in runs)==1
 assert any('-cpi' in a and 'npt.cpt' in a and '-append' in a for a in runs)
 assert all(v['status']=='completed' for v in w.read(folder/'progress.json')['phases'].values())

def test_energy_parser_preserves_units_and_rejects_nan(tmp_path):
 x=tmp_path/'e.xvg';header='@ s0 legend "Temperature"\n@ s1 legend "Pressure"\n@ s2 legend "Potential"\n'
 x.write_text(header+'1000 310 1 -200\n')
 assert membrane.energy_series(x)==[{'time_ns':1,'Temperature':310,'Pressure':1,'Potential':-200}]
 x.write_text(header+'1000 nan 1 -200\n')
 with pytest.raises(ValueError,match='nonfinite'):membrane.energy_series(x)

def test_structure_reader_accepts_overflow_residue_ids(tmp_path):
 p=tmp_path/'large.pdb'
 p.write_text('HETATM    1  O   HOH WA000       1.000   2.000   3.000  1.00  0.00           O\n')
 assert studio.structure(p)[0]['resid']=='A000'
