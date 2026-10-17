"""Persistent orchestration. Scientific subprocesses are isolated from the UI."""
from __future__ import annotations
import hashlib,json,os,re,shutil,subprocess,sys,threading,time,uuid
from pathlib import Path
UI=Path(__file__).resolve().parents[1]
ROOT=UI.parent
RUNS=UI/'workflow_runs'
STAGES=[('prepare','复合物准备','蛋白与多肽分离，保留完整受体和局部口袋'),('peptide2mol','Peptide2Mol','真实模型采样与化学检查'),('reinvent4','REINVENT4','审核核心后扩展候选小分子'),('docking','GNINA 批量对接','完整受体 · 独立分子 · 原始 CNN 分数'),('screening','虚拟筛选','质检、结合模式与多样性筛选'),('gromacs','GROMACS','已准备体系的动力学执行与轨迹文件')]
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
 if not re.fullmatch(r'[a-zA-Z0-9][a-zA-Z0-9_-]{0,79}',rid):raise ValueError('无效运行 ID')
 return RUNS/rid

def read(path,default=None):
 return json.loads(Path(path).read_text(encoding='utf-8')) if Path(path).is_file() else default

def defaults():
 cfg = {'structure':'data/structures/9nu4.pdb','receptor_chain':'A','peptide_chain':'B','receptor_range':[5,500],
 'budget':20,'seed':17,'pocket_radius':6.0,'peptide_python':'','reinvent_executable':'','chem_python':sys.executable,
 'screen_python':'','gnina':'','gnina_device':'gpu','gmx':'','target_profile':'ecmurj',
 'core_selections':[], 'docking_receptor':'','docking_box':None,'gromacs_tpr':'',
 'gromacs_system_reviewed':False,'timeout_seconds':86400}
 cfg.update({k:v for k,v in read(ROOT/'deployment/runtime.json',{}).items() if k in cfg})
 return cfg

def create(config):
 allowed=set(defaults());unknown=set(config)-allowed
 if unknown:raise ValueError('未知配置字段: '+','.join(sorted(unknown)))
 cfg={**defaults(),**config};source=(ROOT/cfg['structure']).resolve()
 if not source.is_relative_to(ROOT.resolve()) or not source.is_file() or source.suffix.lower()!='.pdb':raise ValueError('结构必须是项目文件夹内的 PDB')
 if cfg['receptor_chain']==cfg['peptide_chain']:raise ValueError('受体链和多肽链不能相同')
 for k in ('receptor_chain','peptide_chain'):
  if not isinstance(cfg[k],str) or len(cfg[k])!=1:raise ValueError('PDB 链编号必须为一个字符')
 for k,low,high in [('budget',1,10000),('seed',0,2147483647),('timeout_seconds',60,604800)]:
  if isinstance(cfg[k],bool) or not isinstance(cfg[k],int) or not low<=cfg[k]<=high:raise ValueError(k+' 超出范围')
 rr=cfg['receptor_range']
 if not isinstance(rr,list) or len(rr)!=2 or not all(isinstance(x,int) for x in rr) or rr[0]>rr[1]:raise ValueError('受体残基范围无效')
 if cfg['gnina_device'] not in ('cpu','gpu'):raise ValueError('GNINA 设备无效')
 if cfg['target_profile'] not in ('ecmurj','generic'):raise ValueError('靶点配置无效')
 rid=time.strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:8];d=run_dir(rid);d.mkdir(parents=True)
 shutil.copyfile(source,d/'input.pdb');cfg['input_sha256']=digest(d/'input.pdb')
 # External executables are configured locally, never inferred from old machine paths.
 atomic(d/'config.json',cfg);atomic(d/'run.json',{'id':rid,'created':time.time(),'config_sha256':digest(d/'config.json'),'status':'created'})
 return {'id':rid}

def require(d,stage):
 rec=read(d/stage/'stage.json',{})
 if rec.get('status')!='completed':raise RuntimeError('上游阶段未完成: '+stage)
 if rec.get('config_sha256')!=digest(d/'config.json'):raise RuntimeError('配置已改变，请新建运行')
 for p,sha in rec.get('artifacts',{}).items():
  path=(d/p).resolve()
  if not path.is_relative_to(d.resolve()) or not path.is_file() or digest(path)!=sha:raise RuntimeError('上游产物缺失或被修改: '+p)
 return rec

def status(rid):
 d=run_dir(rid);cfg=read(d/'config.json')
 if cfg is None:raise FileNotFoundError('运行不存在')
 stages=[]
 for key,label,desc in STAGES:
  rec=read(d/key/'stage.json',{})
  if rec.get('status')=='running' and rid not in _active:
   rec={**rec,'status':'interrupted','error':'服务已重启；旧进程状态未知。请确认旧进程停止后新建运行，不覆盖原产物。'}
  stages.append({'key':key,'label':label,'description':desc,**rec,'status':rec.get('status','pending')})
 return {'id':rid,'config':cfg,'stages':stages,'reviews':{k:read(d/'reviews'/f'{k}.json',{}) for k,_,_ in STAGES}}

def inventory():
 paths={'Peptide2Mol':'01_molecule_generation/vendor/Peptide2Mol/ckpts/PMT_major.ckpt',
 'REINVENT4':'01_molecule_generation/vendor/REINVENT4/priors/libinvent.prior',
 'GNINA adapter':'02_docking/src/sbvs2/docking/run_gnina.py',
 'Virtual screening':'03_virtual_screening/src/zvss/stages.py',
 'OpenMM':'05_molecular_dynamics/modules/md/openmm_engine.py'}
 modules=[{'name':n,'present':(ROOT/p).is_file(),'path':p} for n,p in paths.items()]
 return {'modules':modules,'structures':[p.relative_to(ROOT).as_posix() for p in sorted((ROOT/'data/structures').glob('*.pdb'))],
 'runs':[p.name for p in sorted(RUNS.glob('*'),reverse=True) if (p/'run.json').is_file()],
 'defaults':defaults(),'platform':sys.platform,'target_platform':'Windows + WSL2 + NVIDIA',
 'notice':'代码与模型存在不等于已通过真实计算验收；本界面不把历史骨架计作模型结果。'}

def start(rid,stage):
 if stage not in dict((k,l) for k,l,_ in STAGES):raise ValueError('未知阶段')
 d=run_dir(rid)
 with _lock:
  if _active:raise RuntimeError('当前已有计算任务；为避免 GPU/输出冲突，执行队列串行运行')
  cfg=read(d/'config.json');meta=read(d/'run.json')
  if cfg is None or meta is None:raise ValueError('运行不存在')
  if digest(d/'config.json')!=meta['config_sha256']:raise ValueError('配置已改变；请复制参数新建运行')
  target=d/stage
  if target.exists():
   previous_record=read(target/'stage.json',{})
   if previous_record.get('status') not in {'failed','blocked','cancelled'}:
    raise RuntimeError('该阶段已有完成或中断记录，不能覆盖')
  previous=[k for k,_,_ in STAGES];idx=previous.index(stage)
  if idx:require(d,previous[idx-1])
  if target.exists():
   archive=d/'attempts'/f'{stage}_{time.time_ns()}';archive.parent.mkdir(exist_ok=True);target.rename(archive)
  target.mkdir();atomic(target/'review.json',read(d/'reviews'/f'{stage}.json',{}));atomic(target/'stage.json',{'status':'running','started':time.time(),'counts':{},'config_sha256':digest(d/'config.json')})
  py=cfg.get('chem_python') or sys.executable
  logfile=(target/'process.log').open('w',encoding='utf-8')
  try:proc=subprocess.Popen([py,str(UI/'pipeline_bridge.py'),'--run-id',rid,'--stage',stage],stdout=logfile,stderr=subprocess.STDOUT,cwd=UI,start_new_session=(os.name!='nt'))
  except Exception as e:
   logfile.close();atomic(target/'stage.json',{'status':'failed','error':str(e)});raise
  _active[rid]=proc
  def wait():
   code=proc.wait();logfile.close()
   with _lock:
    _active.pop(rid,None)
    rec=read(target/'stage.json',{})
    if rec.get('status')=='running':atomic(target/'stage.json',{**rec,'status':'failed','error':f'进程退出 {code}，无完整结果','ended':time.time()})
  threading.Thread(target=wait,daemon=True).start()
 return {'id':rid,'stage':stage,'pid':proc.pid}

def cancel(rid):
 import signal
 with _lock:
  proc=_active.get(rid)
  if not proc:return {'cancelled':False}
  if os.name=='nt':proc.terminate()
  else:os.killpg(proc.pid,signal.SIGTERM)
  def kill_later():
   try:proc.wait(timeout=5)
   except subprocess.TimeoutExpired:
    if os.name=='nt':proc.kill()
    else:
     try:os.killpg(proc.pid,signal.SIGKILL)
     except ProcessLookupError:pass
  threading.Thread(target=kill_later,daemon=True).start()
  d=run_dir(rid)
  for key,_,_ in STAGES:
   p=d/key/'stage.json';rec=read(p,{})
   if rec.get('status')=='running':atomic(p,{**rec,'status':'cancelled','ended':time.time()})
 return {'cancelled':True}

def log(rid,stage):
 if stage not in [x[0] for x in STAGES]:raise ValueError('未知阶段')
 p=run_dir(rid)/stage/'process.log'
 if not p.exists():return {'text':''}
 with p.open('rb') as f:
  f.seek(max(0,p.stat().st_size-64000));return {'text':f.read().decode('utf-8','replace')}

REVIEW_FIELDS={
 'prepare':[],
 'peptide2mol':['peptide_python'],
 'reinvent4':['reinvent_executable','core_selections'],
 'docking':['gnina','gnina_device','docking_receptor','docking_box'],
 'screening':['screen_python'],
 'gromacs':['gmx','gromacs_tpr','gromacs_system_reviewed'],
}
def review(rid,stage,payload):
 d=run_dir(rid)
 if stage not in REVIEW_FIELDS or set(payload)-set(REVIEW_FIELDS[stage]):raise ValueError('无效阶段参数')
 if not (d/'run.json').is_file():raise ValueError('运行不存在')
 with _lock:
  if rid in _active:raise RuntimeError('进程尚未结束')
  if (d/stage).exists() and read(d/stage/'stage.json',{}).get('status') not in {'failed','blocked','cancelled'}:raise RuntimeError('阶段已开始或完成，参数不可修改')
  atomic(d/'reviews'/f'{stage}.json',payload)
 return {'saved':True}
