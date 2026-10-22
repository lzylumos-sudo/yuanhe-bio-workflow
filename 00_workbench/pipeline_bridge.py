#!/usr/bin/env python3
"""Real stage adapters. No surrogate molecules, scores, or success markers."""
from __future__ import annotations
import argparse,csv,json,os,shutil,subprocess,sys,time
from pathlib import Path
from webui.workflow import ROOT,UI,STAGES,atomic,digest,read,require,run_dir,lease_fds
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
 from webui.readiness import require_ready
 require_ready(cfg,ROOT,d/'deployment_readiness.json')
 from src.structure import prepare_structure
 if digest(d/'input.pdb')!=cfg['input_sha256']:raise RuntimeError('Input structure hash mismatch')
 result=prepare_structure(d/'input.pdb',out,structure_id=Path(cfg['structure']).stem,
 receptor_chain=cfg['receptor_chain'],peptide_chain=cfg['peptide_chain'],receptor_start=cfg['receptor_range'][0],
 receptor_end=cfg['receptor_range'][1],crop_distance_A=float(cfg['pocket_radius']),source_url='local project input')
 atomic(out/'preparation.json',result['report'])
 return {'complexes':1,'receptor_atoms':result['report']['receptor']['range_atom_count'],'pocket_atoms':result['report']['receptor']['prepared_atom_count']}

def stage_peptide2mol(d,cfg,out):
 from src.adapters.registry import backend
 peptide2mol=backend("peptide2mol")
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
 workers=max(1,min(workers,4)) if not cfg.get('research') else min(cfg['research']['generation_workers'],cfg['budget'])
 if cfg.get('research'):
  # Serial model allocation until measured Peptide2Mol memory requirements are available.
  # Each worker receives a single leased GPU through its subprocess environment.
  from research.resources import detect
  workers=min(workers,max(1,len([g for g in detect(d)['gpus'] if g['allocatable']])))
 quotas=[cfg['budget']//workers+(i<cfg['budget']%workers) for i in range(workers)]
 atomic(out/'resource_plan.json',{'workers':workers,'sample_quotas':quotas,'seeds':[cfg['seed']+i for i in range(workers)],'total_attempts':sum(quotas),'diffusion_steps':'model default, unchanged'})
 def sample(i):
  kwargs=dict(pocket_pdb=d/'prepare/receptor_prepared.pdb',checkpoint=GEN/mc['checkpoint'],checkpoint_sha256=mc['checkpoint_sha256'],output_dir=out/f'model_{i:02d}',budget=quotas[i],seed=cfg['seed']+i,executable=py,repo=GEN/'vendor/Peptide2Mol')
  if cfg.get('research'):
   from research.generation import guarded_sample
   return guarded_sample(peptide2mol.run,kwargs,d,cfg['research'],'peptide2mol')
  return peptide2mol.run(**kwargs)
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
 from src.adapters.registry import backend
 reinvent4=backend("reinvent4")
 from src.cores import derive_design_core
 from rdkit import Chem
 from rdkit.Chem import AllChem
 selections=cfg.get('core_selections') or []
 rc=cfg.get('research')
 if rc and len(selections)>rc['core_limit']:raise ValueError('Reviewed core limit exceeded')
 if not selections:raise Blocked('Reviewed core selection required: in core_selections enter seed_id and key_atom_indices(from 1 onward). First inspect Peptide2Mol seed structures.')
 seeds={r['id']:r for r in read(d/'peptide2mol/seeds.json',[])};mc=read(GEN/'configs/project.json')['models']['reinvent4']
 exe=executable(cfg['reinvent_executable'],'REINVENT4');rows=[];failures=[];seen=set();core_records=[];sampled=0
 if rc:(out/'lineage.jsonl').write_text('')
 quotas=[rc['reinvent_attempts']//len(selections)+(i<rc['reinvent_attempts']%len(selections)) for i in range(len(selections))] if rc else [cfg['budget']]*len(selections)
 atomic(out/'budget_plan.json',{'quotas':quotas,'total_attempts':sum(quotas),'scope':'whole_task' if rc else 'historical_per_core'})
 for n,selection in enumerate(selections):
  if not quotas[n]:continue
  if rc and len(rows)>=rc['generated_unique_limit']:break
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
  kwargs=dict(scaffold_smiles=scaffold,checkpoint=GEN/mc['prior'],checkpoint_sha256=mc['prior_sha256'],output_dir=out/f'core_{n:04d}',budget=quotas[n],seed=cfg['seed']+n,executable=exe,repo=GEN/'vendor/REINVENT4')
  if rc:
   from research.generation import guarded_sample
   result=guarded_sample(reinvent4.run,kwargs,d,rc,'reinvent4')
  else:result=reinvent4.run(**kwargs,device='cpu')
  core_records[-1].update(attempts=result['report']['completed'],valid=result['report'].get('valid'),unique=result['report'].get('unique'))
  sampled+=result['report']['completed']
  for c in read(out/f'core_{n:04d}/candidates.json')['candidates']:
   smiles=c['canonical_smiles'];base=Chem.MolFromSmiles(smiles)
   for atom in base.GetAtoms():atom.SetAtomMapNum(0)
   smiles=Chem.MolToSmiles(base,isomericSmiles=True)
   if rc:
    from research.qualification import identity
    _,smiles,_=identity(smiles);base=Chem.MolFromSmiles(smiles)
   if rc:
    with (out/'lineage.jsonl').open('a') as lineage:lineage.write(json.dumps({'compound_id':identity(smiles)[0],'core_index':n,'seed_id':seed['id'],'sample_record':c,'manifest':str((out/f'core_{n:04d}/manifest.json').relative_to(d))})+'\n')
   if smiles in seen:continue
   seen.add(smiles);cid=identity(smiles)[0] if rc else 'CMP'+__import__('hashlib').sha256(smiles.encode()).hexdigest()[:20]
   mol=Chem.AddHs(base);params=AllChem.ETKDGv3();params.randomSeed=cfg['seed']
   if AllChem.EmbedMolecule(mol,params)!=0:
    failures.append({'id':cid,'reason':'conformer_generation_failed'});continue
   mol.SetProp('_Name',cid);p=out/'ligands'/f'{cid}.sdf';p.parent.mkdir(exist_ok=True)
   with Chem.SDWriter(str(p)) as w:w.write(mol)
   if rc and len(rows)>=rc['generated_unique_limit']:continue
   rows.append({'id':cid,'smiles':smiles,'seed_id':seed['id'],'core':scaffold,'sdf':str(p.relative_to(d)),
   'route_status':'not_evaluated','protonation':'as_input_not_pH_enumerated','provenance_manifest':str((out/f'core_{n:04d}/manifest.json').relative_to(d))})
 atomic(out/'cores.json',core_records);atomic(out/'failures.json',failures);atomic(out/'candidates.json',rows)
 csv_write(out/'molecules.csv',rows,['id','smiles','seed_id','core','route_status','protonation','sdf','provenance_manifest'])
 if not rows:raise RuntimeError('No valid unique molecules available for docking')
 return {'sampled':sampled,'unique':len(seen),'failed':len(failures),'remaining':len(rows)}

def stage_risk(d,cfg,out):
 if cfg.get('research'):
  from research.qualification import prepare_risk
  return prepare_risk(ROOT,d,cfg,out)
 from yuanhe_preflight.pipeline import risk_stage
 if cfg.get('imported_candidates'):
  from yuanhe_preflight.cli import import_candidates
  rows=import_candidates(Path(cfg['imported_candidates']),out/'imported')
 else:
  rows=list(read(d/'reinvent4/candidates.json',[]))
  rows += [{**r,'sdf':r['raw_sdf'],'source_stage':'peptide2mol'} for r in read(d/'peptide2mol/seeds.json',[])]
 return risk_stage(rows,cfg['predocking'],out,base=d,batch_id=d.name)

def stage_availability(d,cfg,out):
 from yuanhe_preflight.pipeline import availability_stage
 require(d,'risk')
 return availability_stage(read(d/'risk/records.json',[]),cfg['predocking'],out)

def stage_docking(d,cfg,out):
 if cfg.get('research'):
  from research.docking import run
  return run(ROOT,d,cfg,out,executable)
 sys.path.insert(0,str(ROOT/'02_docking/src'))
 from sbvs2.docking.vina_official import run_vina_job, validate_resources
 from sbvs2.docking.protocol import Protocol
 from sbvs2.docking.box import Box
 require(d,'risk');require(d,'availability')
 from yuanhe_preflight.gate import require_ligand
 rows=read(d/'availability/qualified.json',[])
 for row in rows:require_ligand(row['source_sdf'],row['qualification_evidence'])
 if not rows:
  atomic(out/'results.json',[]);csv_write(out/'molecules.csv',[],['id','smiles','status'])
  return {'input':0,'remaining':0,'failed':0,'poses':0,'shortfall':read(d/'risk/summary.json',{}).get('input',0)}
 exe=executable(cfg['vina'],'AutoDock Vina')
 validate_resources(int(cfg['docking_cpu']), int(cfg['docking_workers']), int(cfg['docking_total_threads']))
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
 if b.get('units','angstrom') not in {'angstrom','A','Å'}:raise ValueError('Docking box units must be angstrom')
 box=Box(tuple(map(float,b['center'])),tuple(map(float,b['size'])),'manual')
 if box.validate():raise ValueError('Invalid or oversized docking box: '+','.join(box.validate()))
 receptor=out/'receptor.pdb';shutil.copyfile(rec,receptor);atomic(out/'box.json',b)
 proto=Protocol.preset('standard',seeds=(cfg['seed'],),exhaustiveness=cfg['exhaustiveness'],num_modes=cfg['num_modes'],energy_range=cfg['energy_range']);results=[];successful=[]
 atomic(out/'input_molecules.json',rows)
 from concurrent.futures import ThreadPoolExecutor, as_completed
 seeds=cfg.get('docking_seeds') or [cfg['seed']]
 if not seeds or len(set(seeds))!=len(seeds) or any(type(s) is not int or not 0<s<2**31 for s in seeds):raise ValueError('Distinct positive integer docking seeds required')
 tasks=[(row,int(seed)) for row in rows for seed in seeds]
 def dock(task):
  row,seed=task;job_id=row['id']+'_seed'+str(seed)
  result=run_vina_job(exe,receptor,d/row['sdf'],box,proto,seed,out/'jobs'/job_id,cpu=int(cfg['docking_cpu']),timeout=cfg['timeout_seconds'],qualification_evidence=row['qualification_evidence'],meta={'predocking':{k:row.get(k) for k in ('compound_id','state_id','source_core','generation_batch','labels','qualification_sha256','qualification_evidence')}})
  entry={**row,'job_id':job_id,'seed':seed,'status':result['status'],'engine':'vina','engine_version':'1.2.7','error':result.get('error')}
  if result['status']=='ok':entry.update(min(result['poses'],key=lambda x:x['vina_affinity_kcal_mol']))
  return entry
 with ThreadPoolExecutor(max_workers=int(cfg['docking_workers'])) as pool:
  for future in as_completed([pool.submit(dock,task) for task in tasks]):
   entry=future.result();results.append(entry)
   if entry['status']=='ok':successful.append(entry['job_id'])
   atomic(out/'progress.json',{'input':len(tasks),'processed':len(results),'successful':len(successful),'failed':len(results)-len(successful)})
 results.sort(key=lambda r:(r['id'],r['seed']))
 atomic(out/'results.json',results);csv_write(out/'molecules.csv',results,['id','compound_id','state_id','smiles','status','vina_affinity_kcal_mol','engine','labels','source_core','generation_batch','qualification_evidence','qualification_sha256','error'])
 if not successful:raise RuntimeError('All molecules failed docking; inspect each job_record.json')
 return {'input':len(rows),'shortfall':max(0,read(d/'risk/summary.json',{}).get('input',len(rows))-len(rows)),'remaining':len({r['id'] for r in results if r['status']=='ok'}),'failed':len(tasks)-len(successful),'poses':sum(len(read(out/'jobs'/x/'job_record.json')['poses']) for x in successful)}

def stage_screening(d,cfg,out):
 if cfg.get('research') and not cfg.get('_screen_receptor'):
  from research.screening import run
  return run(stage_screening,d,cfg,out)
 docking_records=read(d/'docking/results.json')
 if docking_records == []:
  atomic(out/'candidates.json',[])
  csv_write(out/'molecules.csv',[],['compound_id','state_id','pose_id','qualification'])
  return {'input':0,'remaining':0,'formal':0}
 screen_receptor=Path(cfg['_screen_receptor']['path']) if cfg.get('_screen_receptor') else d/'docking/receptor.pdb'
 py=executable(cfg.get('screen_python') or cfg['chem_python'],'Virtual screening Python')
 # Separate stage-local database prevents cross-run counts and candidate contamination.
 workspace=out/'workspace';workspace.mkdir();shutil.copytree(SCREEN/'configs',workspace/'configs')
 if cfg.get('research'):
  import yaml
  rc=cfg['research'];budget=yaml.safe_load((workspace/'configs/budgets/full.yaml').read_text())
  budget.update(budget_id='research',redock_max_parents=0,deep_max_parents=rc['review_limit'],md_max_parents=rc['md_limit'],md_primary_slots=rc['md_limit'],md_exploratory_slots=0,not_run=['redock'],seeds=[cfg['seed']+i for i in range(rc['refine']['seeds'])])
  budget['qualification'].update(redock_seeds_required=rc['refine']['seeds'],display_primary_slots=rc['md_limit'])
  jobs=[]
  for record in read(d/'docking/results.json',[]):
   if record.get('receptor_sha256')!=cfg['_screen_receptor']['sha256']:continue
   poses=d/'docking/jobs'/record['job_id']/'poses.sdf'
   jobs.append({'smiles':record['smiles'],'seed':record['seed'],'status':record['status'],'poses':str(poses),'poses_sha256':digest(poses) if poses.is_file() else None})
  evidence=workspace/'refinement_seeds.json';atomic(evidence,{'jobs':jobs,'receptor_sha256':cfg['_screen_receptor']['sha256'],'cluster_rmsd_A':budget['cluster_rmsd_A']})
  budget['qualification']['imported_seed_evidence']={'path':str(evidence),'sha256':digest(evidence)}
  (workspace/'configs/budgets/research.yaml').write_text(yaml.safe_dump(budget))
  cfg={**cfg,'screening_budget':'research'}
 host={'paths':{'structures_dir':str(ROOT)},'tools':{'vina':[executable(cfg['vina'],'AutoDock Vina')], 'vina_settings':{'engine_version':'1.2.7','exhaustiveness':cfg['exhaustiveness'],'num_modes':cfg['num_modes'],'cpu':cfg['docking_cpu'],'energy_range':cfg['energy_range']}},
 'executors':{'md_openmm_linux':{'mode':'disabled','note':'GROMACS used separately; membrane system not auto-prepared'}}}
 if cfg['target_profile']=='generic':
  from scientific.target_profile import bind_generic
  host['paths']['structures_dir']=str(bind_generic(workspace,d,cfg))
 atomic(workspace/'host.json',host)
 handoff=out/'docking_input';handoff.mkdir()
 shutil.copyfile(screen_receptor,handoff/'receptor.pdb')
 for record in read(d/'docking/results.json',[]):
  if cfg.get('_screen_receptor') and record.get('receptor_sha256')!=cfg['_screen_receptor']['sha256']:continue
  if record['status']!='ok':continue
  target=handoff/(record.get('job_id') or record['id']);target.mkdir()
  for name in ('poses.sdf','poses.sdf.qualification.json','job_record.json'):
   source=d/'docking/jobs'/(record.get('job_id') or record['id'])/name
   if source.is_file():shutil.copyfile(source,target/name)
 env=os.environ.copy();env['PYTHONPATH']=str(SCREEN/'src')
 import tempfile
 env['ZVSS_SCRATCH']=tempfile.mkdtemp(prefix='yuanhe_zvss_')
 command([py,'-m','zvss.cli','--root',workspace,'run','--input',handoff,'--run-id','screen',
 '--target',workspace/'configs/targets'/f"{cfg['target_profile']}.yaml",'--budget',workspace/'configs/budgets'/f"{cfg.get('screening_budget','quick')}.yaml",
 '--host',workspace/'host.json','--hints',json.dumps({'engine':'vina','engine_version':'1.2.7','receptor_file':str(screen_receptor),'receptor_sha256':digest(screen_receptor),'box':cfg['docking_box'],'seed':cfg['seed']})],workspace,cfg['timeout_seconds'],env)
 state=read(workspace/'runs/screen/run_state.json',{})
 if state.get('errors'):raise RuntimeError('Screening module reported an error')
 files=list((workspace/'runs/screen').rglob('candidates.csv'))
 if not files:raise RuntimeError('Screening did not produce candidates.csv')
 with files[0].open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
 from rdkit import Chem
 evidence_index=read(d/'availability/identity_evidence.json',{})
 for row in rows:
  evidence=evidence_index.get(row.get('compound_id'))
  if evidence:row.update({k:v for k,v in evidence.items() if k not in ('compound_id','state_id')})
 atomic(out/'predocking_handoff.json',{'scope':'acquisition evidence only; MD parameters and final wet-lab review unchanged','compounds':evidence_index})
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
 if cfg.get('research'):
  from research.qualification import identity
  proofs={x['compound_id']:x for x in map(json.loads,(d/'docking/qualification/qualification.jsonl').read_text().splitlines())}
  kept=[];seen_parents=set()
  for row in rows:
   if not row.get('smiles'):continue
   cid=identity(row['smiles'])[0]
   if cid not in proofs or not proofs[cid]['screening_qualified']:continue
   if cid not in seen_parents and len(seen_parents)>=cfg['research']['review_limit']:continue
   seen_parents.add(cid);row['research_compound_id']=cid;row['research_evidence']=proofs[cid];row['research_receptor']=cfg['_screen_receptor'];kept.append(row)
  rows=kept
 atomic(out/'candidates.json',rows)
 csv_write(out/'molecules.csv',rows,list(dict.fromkeys(k for row in rows for k in row)))
 if rows:
  from scientific.peptide_contacts import calculate
  evidence_dir=out/'peptide_evidence';evidence_dir.mkdir()
  calculate(d,cfg,{},evidence_dir)
 return {'input':read(d/'docking/stage.json')['counts']['remaining'],'remaining':len({x.get('compound_id') for x in rows if x.get('compound_id')}),'formal':len({x.get('compound_id') for x in rows if x.get('qualification')=='formal'}),'preliminary':len({x.get('compound_id') for x in rows if x.get('qualification')!='formal'})}

def stage_gromacs(d,cfg,out):
 if cfg.get('research') and not cfg.get('md_system_ids'):raise Blocked('Research MD requires reviewed membrane systems; external TPR bypass disabled')
 if cfg.get('md_system_ids'):
  py=executable(cfg.get('md_python'),'membrane simulation Python')
  timeout=cfg['timeout_seconds'] if not cfg.get('research') else cfg['research']['md_task_timeout_seconds']*max(1,len(cfg['md_system_ids']))
  from webui.readiness import runtime_env
  command([py,str(UI/'md_stage.py'),'--run-id',d.name,'--output',out],out,timeout,env=runtime_env(py))
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
 from webui.workflow import stages_for
 try:
  if rec.get('status')!='running':raise RuntimeError('Stage must be started by the workflow manager')
  if digest(d/'config.json')!=read(d/'run.json')['config_sha256']:raise RuntimeError('Configuration changed')
  order=[x[0] for x in stages_for(cfg)];i=order.index(a.stage)
  if i and not (a.stage=='risk' and (read(d/'reviews/risk.json',{}).get('imported_candidates') or cfg.get('imported_candidates'))):require(d,order[i-1])
  if (d/'run_manifest.json').is_file():
   from yuanhe.configuration import fingerprint
   from yuanhe.records import atomic as record_atomic
   from yuanhe.version import environment
   resolved=read(out/'resolved_config.json')
   if not resolved or resolved['sha256']!=fingerprint(resolved['config']):raise ValueError('Resolved stage configuration integrity failure')
   if read(out/'resolved_config.yaml')!=resolved:raise ValueError('Resolved configuration formats differ')
   from yuanhe.configuration import stage_config
   expected=stage_config(cfg,[(prior,read(d/prior/'review.json',{})) for prior in order[:i+1]])
   if expected['sha256']!=resolved['sha256']:raise ValueError('Dispatch configuration does not match frozen run and reviews')
   cfg=resolved['config']
   record_atomic(out/'execution_environment.json',environment())
   from yuanhe.contracts import validate_stage
   manifest=read(d/'run_manifest.json')
   contract={'run_id':d.name,'stage_id':a.stage,'config_hash':fingerprint(cfg),'input_artifacts':manifest['input_artifacts']+manifest['output_artifacts']}
   record_atomic(out/'stage_input.json',validate_stage(contract,'stage-input'))
  else:
   for prior in order[:i+1]:
    cfg.update(read(d/prior/'review.json',{}))
  counts=globals()['stage_'+a.stage](d,cfg,out)
  artifacts={str(x.relative_to(d)):digest(x) for x in out.rglob('*') if x.is_file() and x.name not in {'stage.json','process.log'} and not x.name.endswith(('.sqlite3-wal','.sqlite3-shm'))}
  if (d/'run_manifest.json').is_file():
   from yuanhe.records import artifact
   output={'run_id':d.name,'stage_id':a.stage,'config_hash':fingerprint(cfg),'output_artifacts':[artifact(d/name,d,a.stage,content_hash=sha) for name,sha in artifacts.items()],'scientific_review':'not_reviewed'}
   record_atomic(out/'stage_output.json',validate_stage(output,'stage-output'))
   artifacts[str((out/'stage_output.json').relative_to(d))]=digest(out/'stage_output.json')
  current=read(out/'stage.json',{})
  if current.get('status')=='cancelled':return 1
  atomic(out/'stage.json',{**rec,'status':'completed','counts':counts,'artifacts':artifacts,'ended':time.time(),'scientific_validation':'not_evaluated'})
  from yuanhe.records import sync,event
  sync(d,order)
  if (d/'run_manifest.json').is_file():event(d,'stage_completed','WORKBENCH',a.stage,'Artifacts recorded; scientific review remains separate')
  return 0
 except Exception as e:
  import traceback;traceback.print_exc()
  if read(out/'stage.json',{}).get('status')!='cancelled':atomic(out/'stage.json',{**rec,'status':'blocked' if isinstance(e,Blocked) else 'failed','error':str(e),'ended':time.time()})
  from yuanhe.records import sync
  sync(d,[x[0] for x in stages_for(cfg)])
  return 1
if __name__=='__main__':raise SystemExit(main())
