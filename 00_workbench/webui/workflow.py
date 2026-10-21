"""Persistent orchestration. Scientific subprocesses are isolated from the UI."""
from __future__ import annotations
import hashlib,json,os,re,shutil,subprocess,sys,threading,time,uuid
from pathlib import Path
UI=Path(__file__).resolve().parents[1]
ROOT=UI.parent
sys.path.insert(0,str(ROOT/'02_docking/src'))
RUNS=UI/'workflow_runs'
STAGES=[('prepare','Complex preparation','Separate protein and peptide; retain full receptor and local pocket'),('peptide2mol','Peptide2Mol','Real model sampling and chemistry checks'),('reinvent4','REINVENT4','Expand candidate molecules after core review'),('risk','风险子结构过滤','结构标准化与去重；规则匹配、排除和人工复核'),('availability','可购买／可合成性筛选','精确采购证据优先；本地逆合成完整路线校验'),('docking','AutoDock Vina batch docking','CPU · affinity / kcal mol-1 · lower is better'),('screening','Virtual screening','Quality control, binding-mode and diversity screening'),('gromacs','GROMACS','Membrane minimization, equilibration, production and trajectory analysis')]
_lock=threading.RLock(); _active={}

def digest(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def atomic(path,value):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 tmp=path.with_name(path.name+'.'+uuid.uuid4().hex+'.tmp')
 tmp.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8');tmp.replace(path)

def run_dir(rid):
 if not re.fullmatch(r'[a-zA-Z0-9][a-zA-Z0-9_-]{0,79}',rid):raise ValueError('Invalid run ID')
 return RUNS/rid

def read(path,default=None):
 return json.loads(Path(path).read_text(encoding='utf-8')) if Path(path).is_file() else default

def defaults():
 cfg = {'structure':'data/structures/9nu4.pdb','receptor_chain':'A','peptide_chain':'B','receptor_range':[3,500],
 'budget':20,'generation_recovery_run':'','seed':17,'pocket_radius':6.0,'peptide_python':'','reinvent_executable':'','chem_python':sys.executable,
 'prep_python':'/opt/yuanhe/envs/prep/bin/python','md_python':'/opt/yuanhe/envs/membrane/bin/python','screening_budget':'quick','receptor_review':{},'md_system_ids':[], 'screen_python':'','vina':'vina','docking_engine':'vina','docking_seeds':[],'docking_cpu':4,'docking_workers':1,'docking_total_threads':4,'exhaustiveness':32,'num_modes':9,'energy_range':3.0,'gmx':'','target_profile':'ecmurj',
 'core_selections':[], 'docking_receptor':'','docking_box':None,'gromacs_tpr':'',
 'gromacs_system_reviewed':False,'timeout_seconds':86400}
 from yuanhe_preflight.config import defaults as preflight_defaults
 cfg['predocking']=preflight_defaults(ROOT)
 cfg['imported_candidates']=''
 cfg.update({k:v for k,v in read(ROOT/'deployment/runtime.json',{}).items() if k in cfg})
 return cfg

def source_snapshot():
 roots=['00_workbench/scientific','00_workbench/webui','01_molecule_generation/src','02_docking/src','03_virtual_screening/src','05_molecular_dynamics/modules']
 files=[ROOT/'00_workbench'/n for n in ['pipeline_bridge.py','studio_worker.py','md_stage.py']]
 for root in roots:files.extend((ROOT/root).rglob('*.py'))
 files.extend((ROOT/'03_virtual_screening'/'configs').rglob('*'))
 files.extend((ROOT/'02_docking/src/yuanhe_preflight/configs').rglob('*.json'))
 return {str(p.relative_to(ROOT)):digest(p) for p in sorted(set(files)) if p.is_file()}

def _screen_revision_allowed(name):
 return name.startswith(('03_virtual_screening/src/','03_virtual_screening/configs/')) or name in {'00_workbench/webui/server.py', '00_workbench/webui/workflow.py'}

def _md_revision_allowed(name):
 return _screen_revision_allowed(name) or name.startswith("05_molecular_dynamics/modules/") or name in {str(UI.relative_to(ROOT))+'/'+x for x in ('scientific/membrane.py','scientific/artifacts.py','scientific/md_cli.py','scientific/validation.py','md_stage.py','studio_worker.py','webui/studio.py')}

def continue_md(rid):
 """Audit MD-only changes without recreating or relabeling upstream results."""
 d=run_dir(rid)
 with _lock:
  if _active:raise RuntimeError('A calculation is still active')
  lease=acquire_compute()
  try:
   upstream={k:require(d,k) for k in ('prepare','peptide2mol','reinvent4','docking','screening')}
   initial=read(d/'source_snapshot.json');current=source_snapshot()
   previous=read(d/'active_source_revision.json')
   if previous:
    prev=d/'source_revisions'/previous['file']
    if prev.parent!=d/'source_revisions' or digest(prev)!=previous['sha256']:raise ValueError('Previous revision evidence changed')
   changed={k for k in set(initial)|set(current) if initial.get(k)!=current.get(k)}
   if any(not _md_revision_allowed(k) for k in changed):raise ValueError('Upstream generation or docking code changed; cannot reuse')
   stamp=str(time.time_ns());proof={'scope':'md','reason':'UI-authorized MD continuation; original screening qualification retained',
    'initial_source_sha256':digest(d/'source_snapshot.json'),'source':current,'previous_revision':previous,
    'changes':{k:{'before':initial.get(k),'after':current.get(k)} for k in sorted(changed)},
    'upstream_stage_records':{k:digest(d/k/'stage.json') for k in upstream},'recorded_at':time.time()}
   record=d/'source_revisions'/(stamp+'.json');atomic(record,proof)
   atomic(d/'active_source_revision.json',{'file':record.name,'sha256':digest(record)})
   return {'status':'verified','upstream_stages':list(upstream),'source_revision':record.name,'qualification':'unchanged'}
  finally:lease.close()

def verify_source(d):
 frozen=read(d/'source_snapshot.json')
 revision=read(d/'active_source_revision.json')
 if revision:
  path=d/'source_revisions'/revision['file']
  if path.parent!=d/'source_revisions' or digest(path)!=revision['sha256']:raise ValueError('Invalid source revision evidence')
  proof=read(path)
  if proof['initial_source_sha256']!=digest(d/'source_snapshot.json'):raise ValueError('Original source record changed')
  changed={k for k in set(frozen)|set(proof['source']) if frozen.get(k)!=proof['source'].get(k)}
  if any(not (_md_revision_allowed(k) if proof.get('scope')=='md' else _screen_revision_allowed(k)) for k in changed):raise ValueError('Revision changes an upstream calculation')
  frozen=proof['source']
 if frozen is not None and frozen!=source_snapshot():raise ValueError('Code differs from the recorded run version; preserve outputs and explicitly revise the affected stage')

def rescreen(rid):
 """Archive screening and record a narrowly scoped revision; never rerun models."""
 d=run_dir(rid)
 with _lock:
  if _active:raise RuntimeError('A calculation is still active')
  lease=acquire_compute()
  try:
   if (d/'gromacs').exists() or list((d/'md_systems').glob('*/system.json')):
    raise ValueError('Downstream MD exists; create a separate run instead')
   previous=read(d/'screening/stage.json',{})
   if previous and previous.get('status') not in {'completed','failed','cancelled','blocked'}:
    raise ValueError('Screening must be stopped before revision')
   upstream={key:require(d,key) for key in ('prepare','peptide2mol','reinvent4','docking')}
   initial=read(d/'source_snapshot.json');current=source_snapshot()
   changed={k for k in set(initial)|set(current) if initial.get(k)!=current.get(k)}
   if any(not _screen_revision_allowed(k) for k in changed):raise ValueError('Upstream computation code changed; results cannot be reused')
   stamp=str(time.time_ns());archive=d/'attempts'/('screening_revision_'+stamp)
   proof={'reason':'UI-requested screening rerun after scoped repair; upstream results hash-verified and reused',
          'initial_source_sha256':digest(d/'source_snapshot.json'),'source':current,
          'changes':{k:{'before':initial.get(k),'after':current.get(k)} for k in sorted(changed)},
          'upstream_stage_records':{k:digest(d/k/'stage.json') for k in upstream},
          'previous_screening_stage_sha256':digest(d/'screening/stage.json') if (d/'screening/stage.json').is_file() else None,
          'archive':str(archive.relative_to(d)),'recorded_at':time.time()}
   archive.parent.mkdir(exist_ok=True)
   if (d/'screening').exists():(d/'screening').rename(archive)
   record=d/'source_revisions'/(stamp+'.json');atomic(record,proof)
   atomic(d/'active_source_revision.json',{'file':record.name,'sha256':digest(record)})
  finally:lease.close()
  return start(rid,'screening')


def create(config):
 allowed=set(defaults());unknown=set(config)-allowed
 if unknown:raise ValueError('Unknown configuration field: '+','.join(sorted(unknown)))
 cfg={**defaults(),**config};validate_docking_settings(cfg)
 from yuanhe_preflight.config import validate
 validate(cfg['predocking'])
 source=(ROOT/cfg['structure']).resolve()
 if not source.is_relative_to(ROOT.resolve()) or not source.is_file() or source.suffix.lower()!='.pdb':raise ValueError('Structure must be a project-local PDB')
 if cfg['receptor_chain']==cfg['peptide_chain']:raise ValueError('Receptor and peptide chains must differ')
 for k in ('receptor_chain','peptide_chain'):
  if not isinstance(cfg[k],str) or len(cfg[k])!=1:raise ValueError('PDB Chain identifiers must contain one character')
 for k,low,high in [('budget',1,10000),('seed',0,2147483647),('timeout_seconds',60,604800)]:
  if isinstance(cfg[k],bool) or not isinstance(cfg[k],int) or not low<=cfg[k]<=high:raise ValueError(k+' Out of range')
 rr=cfg['receptor_range']
 if not isinstance(rr,list) or len(rr)!=2 or not all(isinstance(x,int) for x in rr) or rr[0]>rr[1]:raise ValueError('Invalid receptor residue range')
 if cfg['docking_engine']!='vina':raise ValueError('New tasks require official Vina; historical GNINA results remain read-only')
 if cfg['target_profile'] not in ('ecmurj','generic'):raise ValueError('Invalid target configuration')
 rid=time.strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:8];d=run_dir(rid);d.mkdir(parents=True)
 shutil.copyfile(source,d/'input.pdb');cfg['input_sha256']=digest(d/'input.pdb')
 # External executables are configured locally, never inferred from old machine paths.
 atomic(d/'source_snapshot.json',source_snapshot())
 atomic(d/'config.json',cfg);atomic(d/'run.json',{'id':rid,'created':time.time(),'config_sha256':digest(d/'config.json'),'status':'created'})
 return {'id':rid}

def require(d,stage):
 rec=read(d/stage/'stage.json',{})
 if rec.get('status')!='completed':raise RuntimeError('Upstream stage incomplete: '+stage)
 if rec.get('config_sha256')!=digest(d/'config.json'):raise RuntimeError('Configuration changed; create a new run')
 for p,sha in rec.get('artifacts',{}).items():
  path=(d/p).resolve()
  if not path.is_relative_to(d.resolve()) or not path.is_file() or digest(path)!=sha:raise RuntimeError('Upstream artifacts missing or modified: '+p)
 return rec

def process_identity(pid):
 try:
  # Linux /proc starttime prevents a recycled PID from becoming a cancel target.
  text=Path(f'/proc/{pid}/stat').read_text();start=text[text.rfind(')')+2:].split()[19]
  cmd=Path(f'/proc/{pid}/cmdline').read_bytes().split(b'\0')
  return {'pid':pid,'process_start':start,'argv':[x.decode('utf-8','replace') for x in cmd if x]}
 except (OSError,IndexError):return {'pid':pid}

def live_record(rec,rid):
 pid=rec.get('pid')
 if not isinstance(pid,int) or pid<=1 or not rec.get('process_start'):return False
 actual=process_identity(pid)
 if actual.get('process_start')!=rec['process_start']:return False
 argv=actual.get('argv',[])
 if not any(str(UI/x) in argv for x in ('pipeline_bridge.py','studio_worker.py')):return False
 try:return argv[argv.index('--run-id')+1]==rid and os.getpgid(pid)==pid
 except (ValueError,IndexError,ProcessLookupError):return False


def status(rid):
 d=run_dir(rid);cfg=read(d/'config.json')
 if cfg is None:raise FileNotFoundError('Run not found')
 stages=[]
 for key,label,desc in STAGES:
  rec=read(d/key/'stage.json',{})
  if rec.get('status')=='running' and rid not in _active and not live_record(rec,rid):
   rec={**rec,'status':'interrupted','error':'Original process not found after restart. Stop the task to archive its interrupted state, then retry; membrane simulations retain separate checkpoints.'}
  if rec.get('status')=='running' and key in ('risk','availability'):
   rec={**rec,'counts':read(d/key/'progress.json',{})}
  if rec.get('status')=='running' and key=='docking':
   prog=read(d/key/'progress.json',{})
   if prog:rec={**rec,'counts':{'input':prog['input'],'processed':prog['processed'],'remaining':prog['successful'],'failed':prog['failed']}}
  if rec.get('status')=='running' and key=='screening':
   prog=read(d/key/'workspace/runs/screen/live_progress.json',{})
   if prog:rec={**rec,'counts':prog.get('counts',{}),'substage':prog.get('current_stage')}
  if key=='screening' and rec.get('status')=='completed':
   from scientific.validation import effective_candidates
   rows=effective_candidates(d)
   if any(r.get('qualification_provenance') for r in rows):
    rec={**rec,'counts':{**rec.get('counts',{}),'formal':sum(r['qualification']=='formal' for r in rows),'preliminary':sum(r['qualification']!='formal' for r in rows)},'qualification_reassessment':rows[0]['qualification_provenance']}
  if key=='gromacs' and rec.get('status')=='completed':
   for path in sorted((d/'operations').glob('*/replicate_comparison.json'),reverse=True):
    if read(path.parent/'operation.json',{}).get('status')!='completed':continue
    proof=read(path)
    if read(path.parent/'result.json')!=proof or proof.get('gromacs_records_sha256')!=digest(d/'gromacs/systems.json'):raise ValueError('MD quality evidence changed')
    for row in proof['replicas']:
     folder=d/'md_systems'/row['system_id']/'analysis'
     if digest(folder/'summary.json')!=row['analysis_sha256'] or digest(folder/'series.json')!=row['series_sha256']:raise ValueError('MD analysis evidence changed')
    rec={**rec,'counts':{**rec.get('counts',{}),'md_quality_passed':len(proof['short_md_quality_passed_compounds'])},'md_quality_evidence':{'artifact':str(path.relative_to(d)),'sha256':digest(path)}}
    break
  stages.append({'key':key,'label':label,'description':desc,**rec,'status':rec.get('status','pending')})
 return {'id':rid,'config':cfg,'stages':stages,'reviews':{k:read(d/'reviews'/f'{k}.json',{}) for k,_,_ in STAGES}}

def inventory():
 paths={'Peptide2Mol':'01_molecule_generation/vendor/Peptide2Mol/ckpts/PMT_major.ckpt',
 'REINVENT4':'01_molecule_generation/vendor/REINVENT4/priors/libinvent.prior',
 'Vina adapter':'02_docking/src/sbvs2/docking/vina_official.py',
 'Virtual screening':'03_virtual_screening/src/zvss/stages.py',
 'Membrane MD':str(UI.relative_to(ROOT)/'scientific/membrane.py')}
 modules=[{'name':n,'present':(ROOT/p).is_file(),'path':p} for n,p in paths.items()]
 return {'modules':modules,'structures':[str(p.relative_to(ROOT)) for p in sorted(list((ROOT/'data/structures').glob('*.pdb'))+list((ROOT/'inputs').glob('*.pdb')))],
 'runs':[p.name for p in sorted(RUNS.glob('*'),reverse=True) if (p/'run.json').is_file()],
 'defaults':defaults(),'platform':sys.platform,'target_platform':'Windows + WSL2 + NVIDIA',
 'notice':'Code and model availability do not establish real computational acceptance. Historical scaffolds are not counted as model outputs.'}

def acquire_compute():
 import tempfile,fcntl
 key=hashlib.sha256(str(ROOT.resolve()).encode()).hexdigest()[:24]
 fh=open(Path(tempfile.gettempdir())/('yuanhe-compute-'+key+'.lock'),'a+')
 try:fcntl.flock(fh,fcntl.LOCK_EX|fcntl.LOCK_NB)
 except BlockingIOError:
  fh.close();raise RuntimeError('A background process is active, possibly from another interface. Stop it first')
 return fh


def start(rid,stage):
 if stage not in dict((k,l) for k,l,_ in STAGES):raise ValueError('Unknown stage')
 d=run_dir(rid)
 with _lock:
  if _active:raise RuntimeError('A computation is active; to avoid GPU/output conflicts, jobs execute serially')
  verify_source(d)
  cfg=read(d/'config.json');meta=read(d/'run.json')
  if cfg is None or meta is None:raise ValueError('Run not found')
  if digest(d/'config.json')!=meta['config_sha256']:raise ValueError('Configuration changed; copy settings into a new run')
  target=d/stage
  if target.exists():
   previous_record=read(target/'stage.json',{})
   if previous_record.get('status') not in {'failed','blocked','cancelled'}:
    raise RuntimeError('This stage has completed or interrupted records and cannot be overwritten')
  previous=[k for k,_,_ in STAGES];idx=previous.index(stage)
  if idx and not (stage=='risk' and (read(d/'reviews/risk.json',{}).get('imported_candidates') or cfg.get('imported_candidates'))):require(d,previous[idx-1])
  lease=acquire_compute()
  if target.exists():
   archive=d/'attempts'/f'{stage}_{time.time_ns()}';archive.parent.mkdir(exist_ok=True);target.rename(archive)
  target.mkdir();atomic(target/'review.json',read(d/'reviews'/f'{stage}.json',{}));atomic(target/'stage.json',{'status':'running','started':time.time(),'counts':{},'source_revision':read(d/'active_source_revision.json'),'config_sha256':digest(d/'config.json')})
  py=cfg.get('chem_python') or sys.executable
  logfile=(target/'process.log').open('w',encoding='utf-8')
  try:proc=subprocess.Popen([py,str(UI/'pipeline_bridge.py'),'--run-id',rid,'--stage',stage],stdout=logfile,stderr=subprocess.STDOUT,cwd=UI,start_new_session=True,pass_fds=(lease.fileno(),),env={**os.environ,'YUANHE_COMPUTE_FD':str(lease.fileno())})
  except Exception as e:
   logfile.close();lease.close();atomic(target/'stage.json',{'status':'failed','error':str(e)});raise
  _active[rid]=proc
  atomic(target/'stage.json',{**read(target/'stage.json'),**process_identity(proc.pid)})
  def wait():
   code=proc.wait();logfile.close();lease.close()
   with _lock:
    _active.pop(rid,None)
    rec=read(target/'stage.json',{})
    if rec.get('status')=='running':atomic(target/'stage.json',{**rec,'status':'failed','error':f'Process exited {code}, without complete results','ended':time.time()})
  threading.Thread(target=wait,daemon=True).start()
 return {'id':rid,'stage':stage,'pid':proc.pid}

def cancel(rid):
 import signal
 with _lock:
  d=run_dir(rid);paths=[d/key/'stage.json' for key,_,_ in STAGES]+list((d/'operations').glob('*/operation.json'))
  proc=_active.get(rid);records=[(p,read(p,{})) for p in paths]
  pid=proc.pid if proc else next((r['pid'] for _,r in records if r.get('status')=='running' and live_record(r,rid)),None)
  if pid:
   try:os.killpg(pid,signal.SIGTERM)
   except ProcessLookupError:pass
   identity=process_identity(pid)
   def kill_later():
    if proc:
     try:proc.wait(timeout=30);return
     except subprocess.TimeoutExpired:pass
    else:
     # Preserve checkpoint flush time after SIGTERM; identity is rechecked.
     time.sleep(30)
     if not live_record(identity,rid):return
    try:os.killpg(pid,signal.SIGKILL)
    except ProcessLookupError:pass
   threading.Thread(target=kill_later,daemon=True).start()
  else:
   lease=acquire_compute();lease.close()
  changed=False
  for path,record in records:
   if record.get('status')=='running':atomic(path,{**record,'status':'cancelled','ended':time.time()});changed=True
 return {'cancelled':bool(pid) or changed}

def log(rid,stage):
 if stage not in [x[0] for x in STAGES]:raise ValueError('Unknown stage')
 p=run_dir(rid)/stage/'process.log'
 if not p.exists():return {'text':''}
 with p.open('rb') as f:
  f.seek(max(0,p.stat().st_size-64000));return {'text':f.read().decode('utf-8','replace')}

REVIEW_FIELDS={
 'prepare':[],
 'peptide2mol':['peptide_python'],
 'reinvent4':['reinvent_executable','core_selections'],
 'risk':['predocking','imported_candidates'],
 'availability':['predocking'],
 'docking':['vina','docking_seeds','docking_engine','docking_cpu','docking_workers','docking_total_threads','exhaustiveness','num_modes','energy_range','docking_receptor','docking_box','receptor_review'],
 'screening':['screen_python','screening_budget'],
 'gromacs':['gmx','gromacs_tpr','gromacs_system_reviewed','md_system_ids'],
}
def validate_docking_settings(cfg):
 import math
 if cfg.get('docking_engine')!='vina':raise ValueError('New tasks require official AutoDock Vina')
 for key in ('docking_cpu','docking_workers','docking_total_threads','exhaustiveness','num_modes'):
  value=cfg.get(key)
  if type(value)!=int or value<1:raise ValueError(key+' must be a positive integer')
 if cfg['num_modes']>100:raise ValueError('num_modes must not exceed 100')
 if cfg['docking_cpu']*cfg['docking_workers']>cfg['docking_total_threads']:raise ValueError('Concurrent Vina tasks exceed the CPU thread budget')
 value=cfg.get('energy_range')
 if type(value) not in (int,float) or not math.isfinite(value) or value<0:raise ValueError('energy_range must be finite and nonnegative in kcal/mol')
 seeds=cfg.get('docking_seeds',[])
 if not isinstance(seeds,list) or any(type(v)!=int or not 0<=v<=2147483647 for v in seeds) or len(seeds)!=len(set(seeds)):raise ValueError('Docking seeds must be distinct integers from 0 to 2147483647')

def review(rid,stage,payload):
 d=run_dir(rid)
 if 'predocking' in payload:
  from yuanhe_preflight.config import validate
  validate(payload['predocking'])
 if 'screening_budget' in payload and payload['screening_budget'] not in {'quick','redock_trial','full'}:raise ValueError('Invalid screening budget')
 if 'docking_box' in payload and payload['docking_box'] is not None:
  import math
  b=payload['docking_box']
  if not isinstance(b,dict) or any(not isinstance(b.get(k),list) or len(b[k])!=3 or any(type(v) not in (int,float) or not math.isfinite(v) for v in b[k]) for k in ('center','size')):raise ValueError('Box values must be finite')
  if any(not 4<=x<=50 for x in b['size']):raise ValueError('Docking box side lengths must lie within 4–50 Å')
 if stage not in REVIEW_FIELDS or set(payload)-set(REVIEW_FIELDS[stage]):raise ValueError('Invalid stage parameters')
 if not (d/'run.json').is_file():raise ValueError('Run not found')
 if stage=='docking':validate_docking_settings({**read(d/'config.json'),**read(d/'reviews/docking.json',{}),**payload})
 with _lock:
  if rid in _active:raise RuntimeError('Process has not exited')
  if (d/stage).exists() and read(d/stage/'stage.json',{}).get('status') not in {'failed','blocked','cancelled'}:raise RuntimeError('Stage started or completed; settings are immutable')
  atomic(d/'reviews'/f'{stage}.json',{**read(d/'reviews'/f'{stage}.json',{}),**payload})
 return {'saved':True}


def revise_preflight(rid,stage):
 """Explicit UI revision of evidence only. Preserve generation and old attempts."""
 if stage not in ('risk','availability'):raise ValueError('Invalid preflight stage')
 d=run_dir(rid)
 with _lock:
  if _active:raise ValueError('Wait for active work to finish')
  if (d/'docking').exists() or (d/'screening').exists() or (d/'md_systems').exists():raise ValueError('Downstream artifacts exist; create a separate evidence-review run')
  lease=acquire_compute()
  try:
   for key in (('risk','availability') if stage=='risk' else ('availability',)):
    target=d/key
    if target.exists():
     record=read(target/'stage.json',{})
     if record.get('status')=='running':raise ValueError('Stage still running')
     archive=d/'attempts'/(key+'_evidence_revision_'+str(time.time_ns()));archive.parent.mkdir(exist_ok=True);target.rename(archive)
   return {'status':'pending','stage':stage,'generation_unchanged':True}
  finally:lease.close()


def lease_fds():
 try:
  fd=int(os.environ.get('YUANHE_COMPUTE_FD','-1'));os.fstat(fd);return (fd,)
 except (ValueError,OSError):return ()


def gnina_environment():
 # Legacy GNINA libraries are scoped to GNINA, never injected into PyTorch/OpenMM.
 paths=sorted(Path('/opt/sbvs/env/lib/python3.11/site-packages/nvidia').glob('*/lib'))
 paths += [p for p in [Path('/opt/sbvs/env/lib'),Path('/usr/lib/wsl/lib')] if p.is_dir()]
 old=os.environ.get('LD_LIBRARY_PATH','')
 return {'LD_LIBRARY_PATH':os.pathsep.join([str(p) for p in paths]+([old] if old else []))} if paths else {}
