#!/usr/bin/env python3
"""Real stage adapters. No surrogate molecules, scores, or success markers."""
from __future__ import annotations
import argparse,csv,json,math,os,shutil,subprocess,sys,time
from pathlib import Path
from webui.workflow import ROOT,UI,STAGES,atomic,digest,read,require,run_dir
GEN=ROOT/'01_molecule_generation'
SCREEN=ROOT/'03_virtual_screening'
sys.path.insert(0,str(GEN))

class Blocked(RuntimeError):pass

def executable(value,label):
 path=shutil.which(value) if value else None
 if not path and value and Path(value).is_file():path=str(Path(value).resolve())
 if not path or not os.access(path,os.X_OK):raise Blocked(label+' 可执行程序未配置或不可运行')
 return path

def command(argv,cwd,timeout,env=None):
 print(json.dumps({'argv':list(map(str,argv)),'cwd':str(cwd)},ensure_ascii=False),flush=True)
 p=subprocess.run(list(map(str,argv)),cwd=cwd,env=env,timeout=timeout,check=False)
 if p.returncode:raise RuntimeError(f'子进程失败，退出码 {p.returncode}')

def csv_write(path,rows,fields):
 with Path(path).open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)

def stage_prepare(d,cfg,out):
 from src.structure import prepare_structure
 if digest(d/'input.pdb')!=cfg['input_sha256']:raise RuntimeError('输入结构哈希不匹配')
 result=prepare_structure(d/'input.pdb',out,structure_id=Path(cfg['structure']).stem,
 receptor_chain=cfg['receptor_chain'],peptide_chain=cfg['peptide_chain'],receptor_start=cfg['receptor_range'][0],
 receptor_end=cfg['receptor_range'][1],crop_distance_A=float(cfg['pocket_radius']),source_url='local project input')
 atomic(out/'preparation.json',result['report'])
 return {'complexes':1,'receptor_atoms':result['report']['receptor']['range_atom_count'],'pocket_atoms':result['report']['receptor']['prepared_atom_count']}

def stage_peptide2mol(d,cfg,out):
 from src.adapters import peptide2mol
 mc=read(GEN/'configs/project.json')['models']['peptide2mol']
 py=executable(cfg['peptide_python'],'Peptide2Mol Python')
 result=peptide2mol.run(pocket_pdb=d/'prepare/receptor_prepared.pdb',checkpoint=GEN/mc['checkpoint'],checkpoint_sha256=mc['checkpoint_sha256'],
 output_dir=out/'model',budget=cfg['budget'],seed=cfg['seed'],executable=py,repo=GEN/'vendor/Peptide2Mol')
 candidates=read(out/'model/candidates.json')['candidates'];rows=[]
 for n,x in enumerate(candidates):
  if x['valid']:rows.append({'id':f'seed_{n:05d}','smiles':x['smiles'],'raw_sdf':str(Path(x['raw_sdf']).relative_to(d))})
 atomic(out/'seeds.json',rows);csv_write(out/'molecules.csv',rows,['id','smiles','raw_sdf'])
 return {'requested':cfg['budget'],'sampled':result['report']['completed'],'remaining':len(rows)}

def stage_reinvent4(d,cfg,out):
 from src.adapters import reinvent4
 from src.cores import derive_design_core
 from rdkit import Chem
 from rdkit.Chem import AllChem
 selections=cfg.get('core_selections') or []
 if not selections:raise Blocked('缺少已审核的核心选择：core_selections 中填写 seed_id 和 key_atom_indices（从 1 开始）。先检查 Peptide2Mol 种子的结构。')
 seeds={r['id']:r for r in read(d/'peptide2mol/seeds.json',[])};mc=read(GEN/'configs/project.json')['models']['reinvent4']
 exe=executable(cfg['reinvent_executable'],'REINVENT4');rows=[];failures=[];seen=set();core_records=[];sampled=0
 for n,selection in enumerate(selections):
  seed=seeds.get(selection.get('seed_id'))
  if not seed:raise ValueError('核心引用了不存在的种子')
  mol=Chem.SDMolSupplier(str(d/seed['raw_sdf']),removeHs=True)[0]
  if mol is None:raise ValueError('种子 SDF 无法解析')
  keys=selection.get('key_atom_indices',[])
  if not keys or any(not isinstance(x,int) or x<1 or x>mol.GetNumAtoms() for x in keys):raise ValueError('核心原子编号无效')
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
 if not rows:raise RuntimeError('未生成可进入对接的有效独立分子')
 return {'sampled':sampled,'unique':len(seen),'failed':len(failures),'remaining':len(rows)}

def stage_docking(d,cfg,out):
 sys.path.insert(0,str(ROOT/'02_docking/src'))
 from sbvs2.docking.run_gnina import run_gnina_job
 from sbvs2.docking.protocol import Protocol
 from sbvs2.docking.box import Box
 exe=executable(cfg['gnina'],'GNINA')
 if cfg['gnina_device']=='gpu':
  smi=shutil.which('nvidia-smi') or '/usr/lib/wsl/lib/nvidia-smi'
  probe=subprocess.run([smi,'--query-gpu=name','--format=csv,noheader'],capture_output=True,text=True,timeout=20)
  if probe.returncode or not probe.stdout.strip():raise Blocked('请求 GPU 对接，但 WSL 未检测到 NVIDIA GPU')
 # Explicit preparation is a scientific input, not guessed from raw coordinates.
 rec=Path(cfg.get('docking_receptor') or '')
 if not cfg.get('docking_receptor') or not rec.is_file():raise Blocked('请配置经过质子化/残基审核的完整受体 docking_receptor；不会把局部模型口袋直接用于对接')
 expected=read(d/'prepare/preparation.json')['receptor']['range_atom_count']
 atoms=[l for l in rec.read_text().splitlines() if l.startswith(('ATOM  ','HETATM')) and l[76:78].strip().upper()!='H']
 if any(l[21:22]==cfg['peptide_chain'] for l in atoms):raise ValueError('对接受体仍包含多肽链')
 if len(atoms)<expected*.8:raise ValueError('对接受体重原子数显著不足，疑似局部口袋或错误文件')
 b=cfg.get('docking_box')
 if not isinstance(b,dict) or len(b.get('center',[]))!=3 or len(b.get('size',[]))!=3:raise Blocked('需要明确审核后的 docking_box.center 与 size，单位 Å')
 box=Box(tuple(map(float,b['center'])),tuple(map(float,b['size'])),'manual')
 if box.validate():raise ValueError('对接盒无效或过大: '+','.join(box.validate()))
 receptor=out/'receptor.pdb';shutil.copyfile(rec,receptor);atomic(out/'box.json',b)
 proto=Protocol.preset('standard',seeds=(cfg['seed'],));rows=read(d/'reinvent4/candidates.json',[]);results=[];successful=[]
 for row in rows:
  result=run_gnina_job(exe,receptor,d/row['sdf'],box,proto,cfg['seed'],out/'jobs'/row['id'],device=cfg['gnina_device'],cpu=4,timeout=cfg['timeout_seconds'],nvidia_smi=shutil.which('nvidia-smi'))
  entry={**row,'status':result['status'],'error':result.get('error')}
  if result['status']=='ok':
   best=max(result['poses'],key=lambda x:x['cnn_pose_score']);entry.update(best);successful.append(row['id'])
  results.append(entry);atomic(out/'progress.json',{'input':len(rows),'processed':len(results),'successful':len(successful),'failed':len(results)-len(successful)})
 atomic(out/'results.json',results);csv_write(out/'molecules.csv',results,['id','smiles','status','empirical_score_kcal_mol','cnn_pose_score','cnn_affinity_pK','error'])
 if not successful:raise RuntimeError('所有分子对接失败，详见各 job_record.json')
 return {'input':len(rows),'remaining':len(successful),'failed':len(rows)-len(successful),'poses':sum(len(read(out/'jobs'/x/'job_record.json')['poses']) for x in successful)}

def stage_screening(d,cfg,out):
 py=executable(cfg.get('screen_python') or cfg['chem_python'],'虚拟筛选 Python')
 # Separate stage-local database prevents cross-run counts and candidate contamination.
 workspace=out/'workspace';workspace.mkdir();shutil.copytree(SCREEN/'configs',workspace/'configs')
 host={'paths':{'structures_dir':str(ROOT)},'tools':{'gnina':[executable(cfg['gnina'],'GNINA')]},
 'executors':{'md_openmm_linux':{'mode':'disabled','note':'GROMACS used separately; membrane system not auto-prepared'}}}
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
 '--target',workspace/'configs/targets'/f"{cfg['target_profile']}.yaml",'--budget',workspace/'configs/budgets/quick.yaml',
 '--host',workspace/'host.json','--hints',json.dumps({'receptor_file':str(d/'docking/receptor.pdb'),'box':cfg['docking_box'],'seed':cfg['seed']})],workspace,cfg['timeout_seconds'],env)
 state=read(workspace/'runs/screen/run_state.json',{})
 if state.get('errors'):raise RuntimeError('筛选模块报告错误')
 files=list((workspace/'runs/screen').rglob('candidates.csv'))
 if not files:raise RuntimeError('筛选未产出 candidates.csv')
 with files[0].open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
 atomic(out/'candidates.json',rows);shutil.copyfile(files[0],out/'molecules.csv')
 return {'input':read(d/'docking/stage.json')['counts']['remaining'],'remaining':len({x.get('compound_id') for x in rows if x.get('compound_id')})}

def stage_gromacs(d,cfg,out):
 exe=executable(cfg['gmx'],'GROMACS')
 if not cfg.get('gromacs_system_reviewed'):raise Blocked('缺少膜体系、力场与平衡审核；不会把 OpenMM 水溶液流程冒充 GROMACS 膜模拟')
 tpr=Path(cfg.get('gromacs_tpr') or '')
 if not cfg.get('gromacs_tpr') or not tpr.is_file():raise Blocked('需要已构建、已平衡并审核的 production.tpr；当前项目没有完整膜体系自动构建模块')
 shutil.copyfile(tpr,out/'production.tpr')
 command([exe,'--version'],out,60)
 command([exe,'mdrun','-s','production.tpr','-deffnm','production'],out,cfg['timeout_seconds'])
 if not (out/'production.gro').is_file() or not (out/'production.log').is_file():raise RuntimeError('GROMACS 缺少最终结构或日志')
 atomic(out/'analysis_status.json',{'status':'not_evaluated','note':'轨迹需针对所选受体/配体定义分析组；尚未实现自动膜体系准备与科学终点分析'})
 return {'simulations_completed':1,'scientifically_validated':0}

def main():
 p=argparse.ArgumentParser();p.add_argument('--run-id',required=True);p.add_argument('--stage',choices=[x[0] for x in STAGES],required=True);a=p.parse_args()
 d=run_dir(a.run_id);out=d/a.stage;cfg=read(d/'config.json');rec=read(out/'stage.json',{})
 try:
  if rec.get('status')!='running':raise RuntimeError('Stage must be started by the workflow manager')
  if digest(d/'config.json')!=read(d/'run.json')['config_sha256']:raise RuntimeError('配置已变化')
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
