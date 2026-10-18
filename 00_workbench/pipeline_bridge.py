#!/usr/bin/env python3
"""Real stage adapters. No surrogate molecules, scores, or success markers."""
from __future__ import annotations
import argparse,csv,json,math,os,shutil,subprocess,sys,time
from pathlib import Path
from webui.workflow import ROOT,UI,STAGES,atomic,digest,read,require,run_dir,lease_fds,gnina_environment
GEN=ROOT/'01_molecule_generation'
SCREEN=ROOT/'03_virtual_screening'
sys.path.insert(0,str(GEN))

class Blocked(RuntimeError):pass

def executable(value,label):
 path=shutil.which(value) if value else None
 if not path and value and Path(value).is_file():path=str(Path(value).resolve())
 if not path or not os.access(path,os.X_OK):raise Blocked(label+' Executable is missing or cannot run')
 return path

def command(argv,cwd,timeout,env=None):
 print(json.dumps({'argv':list(map(str,argv)),'cwd':str(cwd)},ensure_ascii=False),flush=True)
 p=subprocess.run(list(map(str,argv)),cwd=cwd,env=env,timeout=timeout,check=False,pass_fds=lease_fds())
 if p.returncode:raise RuntimeError(f'Subprocess failed with exit code {p.returncode}')

def csv_write(path,rows,fields):
 with Path(path).open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)

def stage_prepare(d,cfg,out):
 from src.structure import prepare_structure
 if digest(d/'input.pdb')!=cfg['input_sha256']:raise RuntimeError('Input structure hash mismatch')
 result=prepare_structure(d/'input.pdb',out,structure_id=Path(cfg['structure']).stem,
 receptor_chain=cfg['receptor_chain'],peptide_chain=cfg['peptide_chain'],receptor_start=cfg['receptor_range'][0],
 receptor_end=cfg['receptor_range'][1],crop_distance_A=float(cfg['pocket_radius']),source_url='local project input')
 atomic(out/'preparation.json',result['report'])
 return {'complexes':1,'receptor_atoms':result['report']['receptor']['range_atom_count'],'pocket_atoms':result['report']['receptor']['prepared_atom_count']}

def stage_peptide2mol(d,cfg,out):
 from src.adapters import peptide2mol
 from concurrent.futures import ThreadPoolExecutor
 mc=read(GEN/'configs/project.json')['models']['peptide2mol']
 if cfg.get('generation_recovery_run'):
  # Explicit UI-created continuation of this run's interrupted real inference.
  # The original failed stage, raw stdout/stderr and molecules remain immutable.
  parent=run_dir(cfg['generation_recovery_run']);oldcfg=read(parent/'config.json',{})
  if oldcfg.get('input_sha256')!=cfg['input_sha256']:raise ValueError('Recovery source PDB hash differs from current input')
  require(parent,'prepare')
  if digest(parent/'prepare/receptor_prepared.pdb')!=digest(d/'prepare/receptor_prepared.pdb'):raise ValueError('Recovery source pocket differs')
  previous=read(parent/'peptide2mol/stage.json',{})
  if previous.get('status') not in {'failed','cancelled','interrupted'}:raise ValueError('Recovery requires a stopped failed or interrupted sampling stage')
  import re
  rows=[];seen=set();audit=[];rejections=[];started=0;completed=0
  for model in sorted((parent/'peptide2mol').glob('model_*')):
   manifest=read(model/'manifest.json')
   if not manifest:raise ValueError('Original inference manifest missing; recovery refused')
   if manifest.get('input_sha256')!=digest(d/'prepare/receptor_prepared.pdb') or manifest.get('checkpoint',{}).get('actual_sha256')!=mc['checkpoint_sha256']:raise ValueError('Model or pocket provenance verification failed')
   inference=manifest['commands'][-1];stdout=inference.get('stdout','');stderr=inference.get('stderr','')
   nstart=len(re.findall(r'(?<!\d)0/1000',stderr));nend=stdout.count('Success:')+stdout.count('Incomplete molecule:')+stdout.count('Reconstruction error encountered.')
   if not nstart or nstart<nend:raise ValueError('Sampling attempt count cannot be verified')
   started+=nstart;completed+=nend
   for line in stdout.splitlines():
    if line.startswith(('Incomplete molecule:','Reconstruction error encountered.')):rejections.append({'source_batch':model.name,'reason':line})
   source=model/'model_output_SDF';archive=out/'raw_sdf'/model.name
   collected,reasons=peptide2mol._collect_sdfs(source,archive,16)
   for i,x in enumerate(collected):
    if x['valid'] and x['smiles'] not in seen:
     seen.add(x['smiles']);rows.append({'id':f'recovered_{model.name}_{i:05d}','smiles':x['smiles'],'raw_sdf':str(Path(x['raw_sdf']).relative_to(d))})
    else:rejections.append({'source_batch':model.name,'reason':x.get('failure_reason') or 'duplicate'})
   audit.append({'parent_run':parent.name,'batch':model.name,'manifest_sha256':digest(model/'manifest.json'),'original_status':manifest.get('status'),'started_attempts':nstart,'completed_attempts':nend,'checkpoint_sha256':mc['checkpoint_sha256'],'raw_sdf_sha256':{p.name:digest(p) for p in archive.glob('*.sdf')}})
  if started>16:raise ValueError('The 16-attempt budget has been exceeded; do not sample further')
  atomic(out/'recovery_audit.json',{'completion_mode':'validated_partial_output_from_interrupted_inference','parent_run':parent.name,'original_stage_status':previous.get('status'),'started_attempts':started,'completed_attempts':completed,'interrupted_attempts':started-completed,'additional_sampling':0,'batches':audit})
  atomic(out/'rejected.json',rejections);atomic(out/'seeds.json',rows);csv_write(out/'molecules.csv',rows,['id','smiles','raw_sdf'])
  if not rows:raise RuntimeError('No valid molecule in completed partial outputs')
  return {'requested':oldcfg['budget'],'sampled':completed,'attempted':started,'interrupted':started-completed,'remaining':len(rows),'recovered':len(rows)}
 py=executable(cfg['peptide_python'],'Peptide2Mol Python')
 workers=min(int(os.environ.get('YUANHE_GENERATION_WORKERS','4')),cfg['budget'])
 workers=max(1,min(workers,4));quotas=[cfg['budget']//workers+(i<cfg['budget']%workers) for i in range(workers)]
 atomic(out/'resource_plan.json',{'workers':workers,'sample_quotas':quotas,'seeds':[cfg['seed']+i for i in range(workers)],'total_attempts':sum(quotas),'diffusion_steps':'model default, unchanged'})
 def sample(i):
  return peptide2mol.run(pocket_pdb=d/'prepare/receptor_prepared.pdb',checkpoint=GEN/mc['checkpoint'],checkpoint_sha256=mc['checkpoint_sha256'],output_dir=out/f'model_{i:02d}',budget=quotas[i],seed=cfg['seed']+i,executable=py,repo=GEN/'vendor/Peptide2Mol')
 with ThreadPoolExecutor(max_workers=workers) as pool:results=list(pool.map(sample,range(workers)))
 rows=[];seen=set();rejected=[];sampled=0
 for i,result in enumerate(results):
  candidates=read(out/f'model_{i:02d}/candidates.json')['candidates'];sampled+=result['report']['completed']
  for reason,count in result['report'].get('failure_reasons',{}).items():
   rejected.append({'source_batch':i,'reason':reason,'count':count})
  for n,x in enumerate(candidates):
   cid=f'seed_{i:02d}_{n:05d}'
   if x['valid'] and x['smiles'] not in seen:
    seen.add(x['smiles']);rows.append({'id':cid,'smiles':x['smiles'],'raw_sdf':str(Path(x['raw_sdf']).relative_to(d))})
   else:rejected.append({'id':cid,'reason':x.get('failure_reason') or 'duplicate_across_batches','raw_sdf':x['raw_sdf']})
 atomic(out/'rejected.json',rejected);atomic(out/'seeds.json',rows);csv_write(out/'molecules.csv',rows,['id','smiles','raw_sdf'])
 if not rows:raise RuntimeError('No valid unique generated molecules')
 return {'requested':cfg['budget'],'sampled':sampled,'rejected':sampled-len(rows),'remaining':len(rows)}

def stage_reinvent4(d,cfg,out):
 from src.adapters import reinvent4
 from src.cores import derive_design_core
 from rdkit import Chem
 from rdkit.Chem import AllChem
 selections=cfg.get('core_selections') or []
 if not selections:raise Blocked('Reviewed core selection required: in core_selections enter seed_id and key_atom_indices(from 1 onward). First inspect Peptide2Mol seed structures.')
 seeds={r['id']:r for r in read(d/'peptide2mol/seeds.json',[])};mc=read(GEN/'configs/project.json')['models']['reinvent4']
 exe=executable(cfg['reinvent_executable'],'REINVENT4');rows=[];failures=[];seen=set();core_records=[];sampled=0
 for n,selection in enumerate(selections):
  seed=seeds.get(selection.get('seed_id'))
  if not seed:raise ValueError('Core references an unknown seed')
  mol=Chem.SDMolSupplier(str(d/seed['raw_sdf']),removeHs=True)[0]
  if mol is None:raise ValueError('Seed SDF could not be parsed')
  keys=selection.get('key_atom_indices',[])
  if not keys or any(not isinstance(x,int) or x<1 or x>mol.GetNumAtoms() for x in keys):raise ValueError('Invalid core atom indices')
  for atom in mol.GetAtoms():atom.SetAtomMapNum(atom.GetIdx()+1)
  core=derive_design_core(mol,keys,seed_id=seed['id'],max_attachment_points=1)
  coremol=Chem.MolFromSmiles(core.design_core)
  for atom in coremol.GetAtoms():
   if atom.GetAtomicNum()!=0:atom.SetAtomMapNum(0)
  scaffold=Chem.MolToSmiles(coremol,isomericSmiles=True)
  core_records.append({'seed_id':seed['id'],'core':scaffold,'key_atom_indices':keys,'mapped_core':core.design_core})
  result=reinvent4.run(scaffold_smiles=scaffold,checkpoint=GEN/mc['prior'],checkpoint_sha256=mc['prior_sha256'],
  output_dir=out/f'core_{n:04d}',budget=cfg['budget'],seed=cfg['seed']+n,executable=exe,repo=GEN/'vendor/REINVENT4',device='cpu')
  sampled+=result['report']['completed']
  for c in read(out/f'core_{n:04d}/candidates.json')['candidates']:
   smiles=c['canonical_smiles'];base=Chem.MolFromSmiles(smiles)
   for atom in base.GetAtoms():atom.SetAtomMapNum(0)
   smiles=Chem.MolToSmiles(base,isomericSmiles=True)
   if smiles in seen:continue
   seen.add(smiles);cid='CMP'+__import__('hashlib').sha256(smiles.encode()).hexdigest()[:20]
   mol=Chem.AddHs(base);params=AllChem.ETKDGv3();params.randomSeed=cfg['seed']
   if AllChem.EmbedMolecule(mol,params)!=0:
    failures.append({'id':cid,'reason':'conformer_generation_failed'});continue
   mol.SetProp('_Name',cid);p=out/'ligands'/f'{cid}.sdf';p.parent.mkdir(exist_ok=True)
   with Chem.SDWriter(str(p)) as w:w.write(mol)
   rows.append({'id':cid,'smiles':smiles,'seed_id':seed['id'],'core':scaffold,'sdf':str(p.relative_to(d)),
   'route_status':'not_evaluated','protonation':'as_input_not_pH_enumerated','provenance_manifest':str((out/f'core_{n:04d}/manifest.json').relative_to(d))})
 atomic(out/'cores.json',core_records);atomic(out/'failures.json',failures);atomic(out/'candidates.json',rows)
 csv_write(out/'molecules.csv',rows,['id','smiles','seed_id','core','route_status','protonation','sdf','provenance_manifest'])
 if not rows:raise RuntimeError('No valid unique molecules available for docking')
 return {'sampled':sampled,'unique':len(seen),'failed':len(failures),'remaining':len(rows)}

def stage_docking(d,cfg,out):
 sys.path.insert(0,str(ROOT/'02_docking/src'))
 from sbvs2.docking.run_gnina import run_gnina_job
 from sbvs2.docking.protocol import Protocol
 from sbvs2.docking.box import Box
 exe=executable(cfg['gnina'],'GNINA')
 if cfg['gnina_device']=='gpu':
  smi=shutil.which('nvidia-smi') or '/usr/lib/wsl/lib/nvidia-smi'
  probe=subprocess.run([smi,'-L'],capture_output=True,text=True,timeout=20)
  if probe.returncode or not probe.stdout.strip():raise Blocked('GPU docking requested but no NVIDIA GPU detected')
 # Explicit preparation is a scientific input, not guessed from raw coordinates.
 rec=Path(cfg.get('docking_receptor') or '')
 if not rec.is_absolute():rec=d/rec
 if not cfg.get('docking_receptor') or not rec.is_file():raise Blocked('Configure a protonated,/residue-reviewed full receptor docking_receptor; The local model pocket is not used directly for docking')
 review=cfg.get('receptor_review',{})
 if review and (not review.get('reviewed') or digest(rec)!=review.get('sha256')):raise ValueError('Receptor differs from the reviewed file; prepare and review it again')
 expected=read(d/'prepare/preparation.json')['receptor']['range_atom_count']
 atoms=[l for l in rec.read_text().splitlines() if l.startswith(('ATOM  ','HETATM')) and l[76:78].strip().upper()!='H']
 if any(l[21:22]==cfg['peptide_chain'] for l in atoms):raise ValueError('Docking receptor still contains the peptide chain')
 if len(atoms)<expected*.8:raise ValueError('Too few receptor heavy atoms; possible cropped pocket or incorrect file')
 b=cfg.get('docking_box')
 if not isinstance(b,dict) or len(b.get('center',[]))!=3 or len(b.get('size',[]))!=3:raise Blocked('Explicitly reviewed values required for docking_box.center and size, units Å')
 box=Box(tuple(map(float,b['center'])),tuple(map(float,b['size'])),'manual')
 if box.validate():raise ValueError('Invalid or oversized docking box: '+','.join(box.validate()))
 receptor=out/'receptor.pdb';shutil.copyfile(rec,receptor);atomic(out/'box.json',b)
 proto=Protocol.preset('standard',seeds=(cfg['seed'],));rows=list(read(d/'reinvent4/candidates.json',[]));results=[];successful=[]
 from rdkit import Chem
 seen={Chem.MolToSmiles(Chem.MolFromSmiles(r['smiles']),isomericSmiles=True) for r in rows}
 for seed in read(d/'peptide2mol/seeds.json',[]):
  canon=Chem.MolToSmiles(Chem.MolFromSmiles(seed['smiles']),isomericSmiles=True)
  if canon not in seen:
   seen.add(canon);rows.append({'id':seed['id'],'smiles':canon,'sdf':seed['raw_sdf'],'source_stage':'peptide2mol'})
 if len(rows)>16:raise ValueError('This bounded workflow accepts at most 16 unique docking inputs')
 atomic(out/'input_molecules.json',rows)
 for row in rows:
  result=run_gnina_job(exe,receptor,d/row['sdf'],box,proto,cfg['seed'],out/'jobs'/row['id'],device=cfg['gnina_device'],cpu=int(os.environ.get('YUANHE_DOCKING_THREADS','24')),timeout=cfg['timeout_seconds'],nvidia_smi=shutil.which('nvidia-smi'),env_extra=gnina_environment())
  entry={**row,'status':result['status'],'error':result.get('error')}
  if result['status']=='ok':
   best=max(result['poses'],key=lambda x:x['cnn_pose_score']);entry.update(best);successful.append(row['id'])
  results.append(entry);atomic(out/'progress.json',{'input':len(rows),'processed':len(results),'successful':len(successful),'failed':len(results)-len(successful)})
 atomic(out/'results.json',results);csv_write(out/'molecules.csv',results,['id','smiles','status','empirical_score_kcal_mol','cnn_pose_score','cnn_affinity_pK','error'])
 if not successful:raise RuntimeError('All molecules failed docking; inspect each job_record.json')
 return {'input':len(rows),'remaining':len(successful),'failed':len(rows)-len(successful),'poses':sum(len(read(out/'jobs'/x/'job_record.json')['poses']) for x in successful)}

def stage_screening(d,cfg,out):
 py=executable(cfg.get('screen_python') or cfg['chem_python'],'Virtual screening Python')
 # Separate stage-local database prevents cross-run counts and candidate contamination.
 workspace=out/'workspace';workspace.mkdir();shutil.copytree(SCREEN/'configs',workspace/'configs')
 host={'paths':{'structures_dir':str(ROOT)},'tools':{'gnina':[executable(cfg['gnina'],'GNINA')], 'gnina_settings':{'cnn_model':'crossdock_default2018_1_3','device':cfg['gnina_device'],'exhaustiveness':16,'num_modes':9,'cpu':4,'cnn_scoring':'rescore','environment':gnina_environment()}},
 'executors':{'md_openmm_linux':{'mode':'disabled','note':'GROMACS used separately; membrane system not auto-prepared'}}}
 if cfg['target_profile']=='generic':
  from scientific.target_profile import bind_generic
  host['paths']['structures_dir']=str(bind_generic(workspace,d,cfg))
 atomic(workspace/'host.json',host)
 handoff=out/'docking_input';handoff.mkdir()
 shutil.copyfile(d/'docking/receptor.pdb',handoff/'receptor.pdb')
 for record in read(d/'docking/results.json',[]):
  if record['status']!='ok':continue
  target=handoff/record['id'];target.mkdir()
  for name in ('poses.sdf','gnina.log'):
   source=d/'docking/jobs'/record['id']/name
   if source.is_file():shutil.copyfile(source,target/name)
 env=os.environ.copy();env['PYTHONPATH']=str(SCREEN/'src')
 import tempfile
 env['ZVSS_SCRATCH']=tempfile.mkdtemp(prefix='yuanhe_zvss_')
 command([py,'-m','zvss.cli','--root',workspace,'run','--input',handoff,'--run-id','screen',
 '--target',workspace/'configs/targets'/f"{cfg['target_profile']}.yaml",'--budget',workspace/'configs/budgets'/f"{cfg.get('screening_budget','quick')}.yaml",
 '--host',workspace/'host.json','--hints',json.dumps({'receptor_file':str(d/'docking/receptor.pdb'),'receptor_sha256':digest(d/'docking/receptor.pdb'),'box':cfg['docking_box'],'seed':cfg['seed']})],workspace,cfg['timeout_seconds'],env)
 state=read(workspace/'runs/screen/run_state.json',{})
 if state.get('errors'):raise RuntimeError('Screening module reported an error')
 files=list((workspace/'runs/screen').rglob('candidates.csv'))
 if not files:raise RuntimeError('Screening did not produce candidates.csv')
 with files[0].open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
 from rdkit import Chem
 selected=out/'selected_poses';selected.mkdir()
 for row in rows:
  pid=row.get('pose_id','')
  if not pid or '/' in pid or '\\' in pid:continue
  matches=list(workspace.rglob(pid+'.sdf'))
  if not matches:continue
  mol=next((m for m in Chem.SDMolSupplier(str(matches[0]),removeHs=False) if m is not None),None)
  if mol is None:continue
  target=selected/(pid+'.sdf');shutil.copyfile(matches[0],target)
  row['sdf']=str(target.relative_to(d));row['smiles']=Chem.MolToSmiles(Chem.RemoveHs(mol),isomericSmiles=True)
 atomic(out/'candidates.json',rows)
 csv_write(out/'molecules.csv',rows,list(dict.fromkeys(k for row in rows for k in row)))
 if rows:
  from scientific.peptide_contacts import calculate
  evidence_dir=out/'peptide_evidence';evidence_dir.mkdir()
  calculate(d,cfg,{},evidence_dir)
 return {'input':read(d/'docking/stage.json')['counts']['remaining'],'remaining':len({x.get('compound_id') for x in rows if x.get('compound_id')}),'formal':len({x.get('compound_id') for x in rows if x.get('qualification')=='formal'}),'preliminary':len({x.get('compound_id') for x in rows if x.get('qualification')!='formal'})}

def stage_gromacs(d,cfg,out):
 if cfg.get('md_system_ids'):
  py=executable(cfg.get('md_python'),'membrane simulation Python')
  command([py,str(UI/'md_stage.py'),'--run-id',d.name,'--output',out],out,cfg['timeout_seconds'])
  systems=read(out/'systems.json',[])
  return {'simulations_completed':len(systems),'remaining':len({s['compound_id'] for s in systems}),'scientifically_validated':0}
 exe=executable(cfg['gmx'],'GROMACS')
 if not cfg.get('gromacs_system_reviewed'):raise Blocked('Membrane system, force field and equilibration review required; a OpenMM solution workflow is not a substitute for GROMACS membrane simulation')
 tpr=Path(cfg.get('gromacs_tpr') or '')
 if not cfg.get('gromacs_tpr') or not tpr.is_file():raise Blocked('A constructed, equilibrated and reviewed system is required for production.tpr; Build and select a system in the membrane wizard, or import a reviewed TPR')
 shutil.copyfile(tpr,out/'production.tpr')
 command([exe,'--version'],out,60)
 command([exe,'mdrun','-s','production.tpr','-deffnm','production'],out,cfg['timeout_seconds'])
 if not (out/'production.gro').is_file() or not (out/'production.log').is_file():raise RuntimeError('GROMACS Final structure or log missing')
 atomic(out/'analysis_status.json',{'status':'not_evaluated','note':'Trajectory analysis requires receptor-specific/and ligand groups; external TPR mode only runs simulation; use the membrane wizard for integrated analysis'})
 return {'simulations_completed':1,'scientifically_validated':0}

def main():
 p=argparse.ArgumentParser();p.add_argument('--run-id',required=True);p.add_argument('--stage',choices=[x[0] for x in STAGES],required=True);a=p.parse_args()
 d=run_dir(a.run_id);out=d/a.stage;cfg=read(d/'config.json');rec=read(out/'stage.json',{})
 try:
  if rec.get('status')!='running':raise RuntimeError('Stage must be started by the workflow manager')
  if digest(d/'config.json')!=read(d/'run.json')['config_sha256']:raise RuntimeError('Configuration changed')
  order=[x[0] for x in STAGES];i=order.index(a.stage)
  if i:require(d,order[i-1])
  for prior in order[:i+1]:
   cfg.update(read(d/prior/'review.json',{}))
  counts=globals()['stage_'+a.stage](d,cfg,out)
  artifacts={str(x.relative_to(d)):digest(x) for x in out.rglob('*') if x.is_file() and x.name not in {'stage.json','process.log'} and not x.name.endswith(('.sqlite3-wal','.sqlite3-shm'))}
  current=read(out/'stage.json',{})
  if current.get('status')=='cancelled':return 1
  atomic(out/'stage.json',{**rec,'status':'completed','counts':counts,'artifacts':artifacts,'ended':time.time(),'scientific_validation':'not_evaluated'})
  return 0
 except Exception as e:
  import traceback;traceback.print_exc()
  if read(out/'stage.json',{}).get('status')!='cancelled':atomic(out/'stage.json',{**rec,'status':'blocked' if isinstance(e,Blocked) else 'failed','error':str(e),'ended':time.time()})
  return 1
if __name__=='__main__':raise SystemExit(main())
