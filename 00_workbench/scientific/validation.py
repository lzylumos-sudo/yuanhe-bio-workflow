"""Bounded validation operations, with scientific claims separated from execution."""
from __future__ import annotations
import json,shutil,sys,subprocess,uuid,os,math
from pathlib import Path
from . import artifacts as a

MD_QUALITY_PROTOCOL={'version':'single-500ns-quality/2','scope':'provisional trajectory quality and pose retention; not binding, activity or equilibrium proof','min_replicas':1,'min_production_ns':500,'temperature_mean_tolerance_K':10,'max_final_backbone_rmsd_A':3,'max_final_ligand_rmsd_A':3,'minimum_protein_contact_fraction':0.8,'max_ligand_half_mean_difference_A':1,'max_between_replica_ligand_mean_std_A':1}

def replicas(d,c,p,out):
 from .membrane import system_dir
 src=system_dir(d,p.get('system_id'));meta=a.read(src/'system.json',{})
 if meta.get('status') not in {'built','completed'}:raise ValueError('A successfully built system is required')
 if not p.get('reviewed'):raise ValueError('Review common starting coordinates and independent velocity seeds')
 seeds=p.get('seeds',[17])
 if not isinstance(seeds,list) or not 1<=len(seeds)<=5 or len(set(seeds))!=len(seeds) or any(type(s)!=int or not 1<=s<=2147483647 for s in seeds):raise ValueError('Use 1–5 distinct positive integer seeds')
 for name,sha in meta['files'].items():
  if a.digest(a.within(src,name))!=sha:raise ValueError('Built system checksum mismatch')
 rows=[]
 for seed in seeds:
  ident='MD_'+uuid.uuid4().hex[:12];dest=system_dir(d,ident);dest.mkdir(parents=True)
  for name in meta['files']:
   target=dest/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(a.within(src,name),target)
  copy={**meta,'id':ident,'status':'built','parameters':{**meta['parameters'],'replica_seed':seed},'acceptance_protocol':MD_QUALITY_PROTOCOL,'replica':{'parent_system':meta['id'],'parent_system_sha256':a.digest(src/'system.json'),'velocity_seed':seed,'initial_coordinates':'shared; independent velocities, not independent conformations'},'source_revision':a.read(d/'active_source_revision.json')}
  a.atomic(dest/'system.json',copy);rows.append({'system_id':ident,'seed':seed})
 return {'replicas':rows,'next':'Review and run all selected replicas through stage six','scope':'Independent initial velocities from one shared prepared structure'}

def summarize_replicas(d,c,p,out):
 from .membrane import system_dir
 import numpy as np
 ids=p.get('system_ids',[])
 if not ids or len(ids)!=len(set(ids)):raise ValueError('Select at least one completed trajectory')
 rows=[];parents=set();seeds=set();protocols=[]
 completed={r['system_id']:r for r in a.read(d/'gromacs/systems.json',[])}
 for ident in ids:
  folder=system_dir(d,ident);meta=a.read(folder/'system.json',{});rep=meta.get('replica',{})
  if meta.get('status')!='completed' or not rep:raise ValueError('Only completed recorded replicas can be compared')
  recorded=completed.get(ident,{}).get('files',{})
  for relative,sha in recorded.items():
   if a.digest(a.within(d,relative))!=sha:raise ValueError('Recorded simulation provenance mismatch')
  for name in ('analysis/summary.json','analysis/series.json','system.json'):
   path=folder/name;relative=str(path.relative_to(d))
   if not recorded.get(relative) or a.digest(path)!=recorded[relative]:raise ValueError('Recorded analysis provenance mismatch')
  parents.add(rep['parent_system_sha256']);seeds.add(rep['velocity_seed']);protocols.append(meta.get('acceptance_protocol'))
  summary=a.read(folder/'analysis/summary.json',{});series=a.read(folder/'analysis/series.json',[])
  if len(series)<4:raise ValueError('Analyze every replica before comparison')
  for name,sha in summary['source_hashes'].items():
   if a.digest(a.within(folder,name))!=sha:raise ValueError('Analysis provenance mismatch')
  metrics={}
  for key in ('backbone_rmsd_A','ligand_rmsd_A','minimum_contact_A'):
   v=np.array([x[key] for x in series]);n=len(v)//2
   if not np.isfinite(v).all():raise ValueError('Nonfinite trajectory metrics')
   metrics[key]={'mean':float(v.mean()),'first_half_mean':float(v[:n].mean()),'second_half_mean':float(v[n:].mean()),'half_difference':float(v[n:].mean()-v[:n].mean()),'last':float(v[-1])}
  rows.append({'system_id':ident,'compound_id':meta.get('compound_id'),'seed':rep['velocity_seed'],'production_ns':summary.get('production_end_ns',summary['production_span_ns']),'frames':len(series),'metrics':metrics,'last100ns':summary.get('last100ns',{}),'quality':summary['quality_observations'],'temperature_mean_K':summary['thermodynamics']['Temperature']['mean'],'target_temperature_K':meta['parameters']['temperature_K'],'protein_contact_fraction':sum(x['minimum_contact_A']<=4.5 for x in series)/len(series),'analysis_sha256':a.digest(folder/'analysis/summary.json'),'series_sha256':a.digest(folder/'analysis/series.json')})
 if len(parents)!=1 or len(seeds)!=len(ids):raise ValueError('Replicas must share inputs and have distinct seeds')
 result={'replicas':rows,'independent_velocity_repeats':len(rows),'convergence':'not_established','reason':'Block differences and between-replica dispersion are descriptive; short shared-start trajectories do not establish equilibrium or binding affinity','md_validated_molecules':0,'between_replica_mean_std':{k:(float(np.std([r['metrics'][k]['mean'] for r in rows],ddof=1)) if len(rows)>1 else None) for k in rows[0]['metrics']}}
 if any(protocol!=MD_QUALITY_PROTOCOL for protocol in protocols):raise ValueError('Replicate acceptance protocol was absent or differs from this validator')
 criteria=[]
 def criterion(name,passed,observed):criteria.append({'criterion':name,'passed':bool(passed),'observed':observed})
 protocol=protocols[0]
 criterion('independent_velocity_replicas',len(rows)>=protocol['min_replicas'],len(rows))
 for row in rows:
  tag=row['system_id'];q=row['quality']
  criterion(tag+':formal_window',row['last100ns'].get('status')=='computed_requires_review',row['last100ns'])
  criterion(tag+':duration',row['production_ns']>=protocol['min_production_ns']-.001,row['production_ns'])
  criterion(tag+':finite_no_constraint_warning',q['finite_trajectory'] and not q['constraint_warnings'],q)
  criterion(tag+':temperature',abs(row['temperature_mean_K']-row['target_temperature_K'])<=protocol['temperature_mean_tolerance_K'],row['temperature_mean_K'])
  criterion(tag+':backbone_change',row['metrics']['backbone_rmsd_A']['last']<=protocol['max_final_backbone_rmsd_A'],row['metrics']['backbone_rmsd_A']['last'])
  criterion(tag+':ligand_pose_retention',row['metrics']['ligand_rmsd_A']['last']<=protocol['max_final_ligand_rmsd_A'],row['metrics']['ligand_rmsd_A']['last'])
  criterion(tag+':protein_contact_fraction',row['protein_contact_fraction']>=protocol['minimum_protein_contact_fraction'],row['protein_contact_fraction'])
  criterion(tag+':ligand_half_difference',abs(row['metrics']['ligand_rmsd_A']['half_difference'])<=protocol['max_ligand_half_mean_difference_A'],row['metrics']['ligand_rmsd_A']['half_difference'])
 dispersion=result['between_replica_mean_std']['ligand_rmsd_A']
 if dispersion is not None:criterion('between_replica_ligand_dispersion',dispersion<=protocol['max_between_replica_ligand_mean_std_A'],dispersion)
 else:result['between_repeat_status']='not_applicable_single_trajectory'
 result['quality_protocol']=protocol;result['quality_criteria']=criteria
 result['gromacs_records_sha256']=a.digest(d/'gromacs/systems.json')
 result['provisional_md_quality_accepted']=all(x['passed'] for x in criteria)
 result['short_md_quality_passed_compounds']=sorted({r['compound_id'] for r in rows if r['compound_id']}) if result['provisional_md_quality_accepted'] else []
 result['scientific_claim']='Descriptive 500 ns trajectory quality gate; stability requires review and may require extension. Convergence, binding strength and activity remain unestablished.'
 a.atomic(out/'replicate_comparison.json',result);return result

def docking_control(d,c,p,out,root):
 """Export docking package and run its API in a separate interpreter, no UI imports."""
 source=root/'02_docking'
 if not source.exists():source=root/'01_GNINA_SBVS'
 for name in ('src/sbvs2','data/prepared/demo_1drf','data/reference/FOL_crystal.sdf'):
  src=source/name;target=out/'standalone'/name;target.parent.mkdir(parents=True,exist_ok=True)
  if src.is_dir():shutil.copytree(src,target,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
  else:shutil.copy2(src,target)
 a.atomic(out/'standalone/config.json',{'vina':c['vina'],'seeds':[17,41,97],'scope':'1DRF folate crystallographic pose recovery; no activity claim','rmsd_threshold_A':2.0})
 script=out/'standalone/run_control.py'
 script.write_text(CONTROL_SCRIPT)
 env={**os.environ,'PYTHONPATH':str(out/'standalone/src')}
 # Vina CPU execution has no CNN library environment.
 with (out/'standalone/process.log').open('w') as log:
  r=subprocess.run([p.get('python') or c['chem_python'],str(script)],cwd=script.parent,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=7200,pass_fds=a.lease_fds())
 if r.returncode:raise RuntimeError('Independent docking control failed; inspect standalone/process.log')
 result=a.read(out/'standalone/result.json');result['artifact_sha256']={str(x.relative_to(out)):a.digest(x) for x in (out/'standalone').rglob('*') if x.is_file() and '__pycache__' not in x.parts}
 return result

CONTROL_SCRIPT = r'''import json,os,sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from sbvs2.docking.vina_official import run_vina_job
from sbvs2.docking.protocol import Protocol
from sbvs2.docking.box import Box
base=Path(__file__).resolve().parent;c=json.loads((base/'config.json').read_text())
b=json.loads((base/'data/prepared/demo_1drf/box.json').read_text());box=Box(tuple(b['center']),tuple(b['size']),'reference_ligand',approved_reference_pose=True)
refpath=base/'data/reference/FOL_crystal.sdf';ref=Chem.RemoveHs(Chem.SDMolSupplier(str(refpath),removeHs=False)[0])
if ref is None:raise ValueError('Invalid crystallographic reference')
proto=Protocol.preset('standard',seeds=tuple(c['seeds']))
def one(seed):
 job=base/'jobs'/str(seed)
 rec=run_vina_job(c['vina'],base/'data/prepared/demo_1drf/receptor_prepared.pdb',refpath,box,proto,seed,job,cpu=4,timeout=3600)
 if rec['status']!='ok':raise RuntimeError(rec['error'])
 poses=[Chem.RemoveHs(m) for m in Chem.SDMolSupplier(str(job/'poses.sdf'),removeHs=False)]
 rms=[float(rdMolAlign.CalcRMS(m,ref)) for m in poses]
 return {'seed':seed,'rmsd_A':rms,'top1_rmsd_A':rms[0],'best_rmsd_A':min(rms),'top1_recovered':rms[0]<=c['rmsd_threshold_A'],'any_recovered':min(rms)<=c['rmsd_threshold_A'],'cpu_threads':rec['cpu']}
with ThreadPoolExecutor(max_workers=1) as pool:rows=list(pool.map(one,c['seeds']))
result={'scope':c['scope'],'module_import':__import__('sbvs2').__file__,'rmsd':'symmetry-aware heavy atoms in fixed receptor frame; no ligand fitting','threshold_A':c['rmsd_threshold_A'],'runs':rows,'top1_recovery_fraction':sum(r['top1_recovered'] for r in rows)/len(rows),'any_pose_recovery_fraction':sum(r['any_recovered'] for r in rows)/len(rows),'clean_environment':sys.prefix!=sys.base_prefix and 'include-system-site-packages = false' in (Path(sys.prefix)/'pyvenv.cfg').read_text().lower(),'python_executable':sys.executable,'environment_scope':'interpreter isolation is measured; dependencies are recorded separately in the installation lock'}
(base/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False))
'''

def screening_reassessment(d,c,p,out):
 """Recompute pose comparison from immutable GNINA results; never rerun engines."""
 from rdkit import Chem
 from rdkit.Chem import rdMolAlign
 from webui import workflow
 workflow.require(d,'screening')
 base=d/'screening/workspace/runs/screen/stages'
 original=a.read(base/'selection/summary.json');redock=a.read(base/'redock/summary.json')
 candidates=a.read(d/'screening/candidates.json');lookup={x['compound_id']:x for x in candidates}
 result={'scope':'Corrected fixed-frame RMSD reassessment of archived GNINA poses; original stage and previous MD handoffs remain unchanged','source_selection_sha256':a.digest(base/'selection/summary.json'),'source_redock_sha256':a.digest(base/'redock/summary.json'),'source_candidate_sha256':a.digest(d/'screening/candidates.json'),'qualification':{},'formal_qualified':0,'md_validated_molecules':0}
 for compound,q in original['qualification'].items():
  row=lookup[compound];refpath=a.within(d,row['sdf']);ref=Chem.RemoveHs(Chem.SDMolSupplier(str(refpath),removeHs=False)[0]);seeds=[]
  for record in redock['compounds'][compound]['per_seed']:
   posefile=Path(record['file']).resolve()
   if not posefile.is_relative_to(d.resolve()) or not posefile.is_file():raise ValueError('Redocking output outside run')
   rms=[]
   for mol in Chem.SDMolSupplier(str(posefile),removeHs=False):
    if mol is None:raise ValueError('Invalid archived redocked pose')
    probe=Chem.RemoveHs(mol)
    if Chem.MolToSmiles(probe)!=Chem.MolToSmiles(ref):raise ValueError('Redocked molecule identity differs')
    rms.append(float(rdMolAlign.CalcRMS(probe,ref)))
   if not rms:raise ValueError('Empty redocking output')
   seeds.append({'seed':record['seed'],'minimum_rmsd_A':min(rms),'supports_selected_pose':min(rms)<=redock['compounds'][compound]['cluster_rmsd_A'],'source_sha256':a.digest(posefile)})
  support=sum(x['supports_selected_pose'] for x in seeds)
  requirements=[]
  for rule in q['requirements']:
   if rule['requirement']=='redock_seeds':
    required=rule['observed']['required_supporting'];planned=rule['observed']['seeds_planned'];passed=len(seeds)==planned and support>=required
    rule={**rule,'status':'pass' if passed else 'fail','satisfied':passed,'observed':{**rule['observed'],'supporting_seeds_same_cluster':support,'compared_pose_id':row['pose_id']},'detail':'Corrected symmetry-aware fixed-frame comparison of archived poses','source':'versioned reassessment','reason_code':None if passed else 'insufficient_pose_reproduction'}
   if rule['requirement']=='strain':
    observed=dict(rule['observed']);limit=observed.pop('limit',25);observed.update(warning_limit=limit,failure_limit=2*limit);rule={**rule,'observed':observed}
   requirements.append(rule)
  passed=all(x['status']=='pass' for x in requirements);result['formal_qualified']+=int(passed)
  result['qualification'][compound]={**q,'requirements':requirements,'formal_qualified':passed,'qualification':'formal' if passed else 'preliminary','per_seed':seeds,'pose_id':row['pose_id'],'reference_sha256':a.digest(refpath),'failed_requirements':[x['requirement'] for x in requirements if x['status']=='fail'],'blocked_by_missing_evidence':[x['requirement'] for x in requirements if x['status'] not in ('pass','fail')]}
 a.atomic(out/'screening_reassessment.json',result)
 return result

def peptide_evidence(d,c,p,out):
 from webui import workflow
 from .peptide_contacts import calculate
 workflow.require(d,'screening');workflow.require(d,'prepare')
 return calculate(d,c,p,out)


def run(d,c,p,out,root):
 mode=p.get('mode')
 if mode=='docking_control':return docking_control(d,c,p,out,root)
 if mode=='screening_reassessment':return screening_reassessment(d,c,p,out)
 if mode=='peptide_evidence':return peptide_evidence(d,c,p,out)
 if mode=='replicas':return replicas(d,c,p,out)
 if mode=='replica_summary':return summarize_replicas(d,c,p,out)
 raise ValueError('Unknown validation mode')


def effective_candidates(d):
 """Apply only a completed, provenance-matching reassessment; never edit history."""
 rows=a.read(d/'screening/candidates.json',[])
 for path in sorted((d/'operations').glob('*/screening_reassessment.json'),reverse=True):
  if a.read(path.parent/'operation.json',{}).get('status')!='completed':continue
  proof=a.read(path)
  if a.read(path.parent/'result.json')!=proof:raise ValueError('Reassessment result mismatch')
  files={'source_candidate_sha256':d/'screening/candidates.json','source_selection_sha256':d/'screening/workspace/runs/screen/stages/selection/summary.json','source_redock_sha256':d/'screening/workspace/runs/screen/stages/redock/summary.json'}
  if any(a.digest(f)!=proof[k] for k,f in files.items()):raise ValueError('Reassessment inputs changed')
  revised=[]
  for row in rows:
   q=proof['qualification'].get(row['compound_id'])
   if not q or row['pose_id']!=q['pose_id'] or a.digest(a.within(d,row['sdf']))!=q['reference_sha256']:raise ValueError('Reassessment pose mismatch')
   passed=all(r['status']=='pass' for r in q['requirements'])
   if passed!=q['formal_qualified']:raise ValueError('Reassessment verdict mismatch')
   revised.append({**row,'qualification':'formal' if passed else 'preliminary','qualification_provenance':{'artifact':str(path.relative_to(d)),'sha256':a.digest(path),'original':row['qualification']}})
  return sorted(revised,key=lambda r:r['qualification']!='formal')
 return rows
