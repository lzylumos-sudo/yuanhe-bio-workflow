import os,time,uuid,json,shutil,hashlib,traceback
from pathlib import Path
from rdkit import Chem
import yaml
from .paths import ROOT,rootpath,sha,write_json,read_json
from .config import validate
from .preparation.ligands import inputs,prepare,conformers,states as enumerate_states
from .engines import vina,gnina
from .engines.common import Cancelled
from .ranking import aggregate,csvwrite
from .qc import pose_qc
def fingerprint(cfg):
 from importlib.metadata import version
 data={'config':cfg,'receptor':sha(rootpath(cfg['receptor'])),'box':sha(rootpath(cfg['box'])),'engine_binary':sha(Path('/opt/sbvs/vendor')/cfg['engine']),'preparation_packages':{p:version(p) for p in ('rdkit','meeko','numpy')}}
 return hashlib.sha256(json.dumps(data,sort_keys=True).encode()).hexdigest()
def create(source,cfg,out=None):
 cfg=dict(cfg)
 if cfg['engine']=='gnina' and cfg.get('device')=='auto':
  caps=read_json(ROOT/'reports/gnina_capabilities.json') if (ROOT/'reports/gnina_capabilities.json').exists() else {}
  if caps.get('gpu_passed'):cfg.update(device='gpu',device_resolution='automatic: verified GPU')
  elif caps.get('cpu_passed'):cfg.update(device='cpu',device_resolution='automatic fallback to verified GNINA CPU CNN; engine/scoring unchanged')
  else:raise ValueError('GNINA CNN device must pass a smoke test on this host')
 validate(cfg);rid=time.strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:10]
 out=rootpath(out) if out else ROOT/'runs'/rid
 if out.exists():raise ValueError('Run destination exists')
 out.mkdir(parents=True);source=rootpath(source);target=out/('ligands'+source.suffix.lower());shutil.copyfile(source,target)
 cfg=dict(cfg);original_cfg=dict(cfg);snapshot=out/'protocol_inputs';snapshot.mkdir()
 for key,name in [('receptor','receptor.pdbqt'),('receptor_pdb','receptor_prepared.pdb'),('box','box.json')]:
  dst=snapshot/name;shutil.copyfile(rootpath(cfg[key]),dst);cfg[key]=str(dst.relative_to(ROOT))
 if cfg.get('preparation_record'):
  dst=snapshot/'preparation.json';shutil.copyfile(rootpath(cfg['preparation_record']),dst);cfg['preparation_record']=str(dst.relative_to(ROOT))
 write_json(snapshot/'original_config.json',original_cfg)
 m={'run_id':rid,'created':time.time(),'status':'queued','input':target.name,'input_sha256':sha(target),'config':cfg,'protocol_fingerprint':fingerprint(cfg),'completed_entries':{},'failures':[],'counts':{'successful_entries':0,'failed_entries':0},'environment':{'cpu_count':os.cpu_count(),'platform':__import__('platform').platform(),'device':cfg.get('device'),'vendor_lock':read_json(ROOT/'env/vendor.lock.json')},'relative_path':str(out.relative_to(ROOT))}
 write_json(out/'run_manifest.json',m);(out/'config_expanded.yaml').write_text(yaml.safe_dump(cfg,allow_unicode=True));return out
def run(out,resume=False,cancel=lambda:False,update=lambda m:None):
 from filelock import FileLock,Timeout
 lock=FileLock('/opt/sbvs/state/compute.lock')
 while True:
  if cancel():raise Cancelled()
  try:lock.acquire(timeout=.5);break
  except Timeout:continue
 try:return _run(out,resume,cancel,update)
 finally:lock.release()
def _run(out,resume=False,cancel=lambda:False,update=lambda m:None):
 out=rootpath(out);m=read_json(out/'run_manifest.json');cfg=m['config'];validate(cfg)
 if fingerprint(cfg)!=m['protocol_fingerprint']:raise ValueError('Protocol/environment changed; create new run')
 if sha(out/m['input'])!=m['input_sha256']:raise ValueError('Original inputs changed')
 box=read_json(rootpath(cfg['box']));rows=[];failures=[];prepared=[]
 m['status']='preparing';m['started']=time.time();write_json(out/'run_manifest.json',m);update(m)
 recxyz=[]
 for l in rootpath(cfg['receptor_pdb']).read_text().splitlines():
  if l.startswith(('ATOM  ','HETATM')) and l[76:78].strip()!='H':recxyz.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
 import numpy as np
 recxyz=np.asarray(recxyz)
 try:
  for item in inputs(out/m['input']):
   pid=f'P{item["index"]:06d}'
   try:microstates=enumerate_states(item,cfg)
   except Exception as e:microstates=[dict(item,preparation_error=str(e))]
   for si,item in enumerate(microstates,1):
    sid=pid+f'_S{si:02d}'
    try:variants=conformers(item,cfg.get('conformers_generated',1),cfg.get('conformers_retained',1),cfg['seeds'][0]) if cfg.get('conformers_generated',1)>1 else [item]
    except Exception as e:variants=[dict(item,preparation_error=str(e))]
    for ci,item in enumerate(variants,1):
     for seed in cfg['seeds']:
      key=f'{sid}_C{ci:02d}_seed{seed}';dest=out/'entries'/key;dest.mkdir(parents=True,exist_ok=True)
      if cancel():raise Cancelled()
      prior=m['completed_entries'].get(key)
      if resume and prior and all((out/p).is_file() and sha(out/p)==h for p,h in prior['hashes'].items()):
       rows.extend(read_json(dest/'rows.json'));prepared.append(read_json(dest/'preparation.json'));continue
      m['completed_entries'].pop(key,None)
      try:
       info=prepare(item,dest,seed);prepared.append(dict(info,parent_id=pid,state_id=sid,conformer_id=sid+f'_C{ci:02d}',prep_status='success'))
       m['status']='running';write_json(out/'run_manifest.json',m);update(m)
       result=(vina if cfg['engine']=='vina' else gnina).dock(cfg,box,dest,seed,cancel);entryrows=[]
       for i,r in enumerate(result,1):
        mol=r.pop('mol');raw=r.pop('raw_properties',None);poseid=key+f'_pose{i:03d}';posepath=dest/f'pose_{i:03d}.sdf';w=Chem.SDWriter(str(posepath));w.write(mol);w.close()
        flags=info['qc_flags']+pose_qc(mol,info['identity'],box,recxyz)
        if 'chemical_identity_or_stereo_changed' in flags:raise ValueError('CHEMICAL_IDENTITY_CHANGED_ON_EXPORT')
        row=dict(r,run_id=m['run_id'],parent_id=pid,original_name=item['name'],state_id=sid,conformer_id=sid+f'_C{ci:02d}',pose_id=poseid,engine=cfg['engine'],engine_version='1.2.7' if cfg['engine']=='vina' else '1.3.3',scoring_function=cfg.get('empirical_scoring','vina'),seed=seed,cnn_vs=None,prep_status='success',dock_status='success',qc_flags=';'.join(flags),receptor_hash=sha(rootpath(cfg['receptor'])),box_hash=sha(rootpath(cfg['box'])),pose_file=str(posepath.relative_to(out)))
        if cfg.get("operation") in ("score_only","refinement"):
         from .qc import fixed_frame_rmsd
         inputmol=next(iter(Chem.SDMolSupplier(str(dest/"input.sdf"),removeHs=False)));row["input_pose_displacement_rmsd_A"]=fixed_frame_rmsd(mol,inputmol)
        entryrows.append(row)
        if raw:write_json(dest/f'pose_{i:03d}_properties.json',raw)
       write_json(dest/'rows.json',entryrows)
       hashes={str(p.relative_to(out)):sha(p) for p in dest.iterdir() if p.is_file()}
       m['completed_entries'][key]={'hashes':hashes,'completed':time.time()};rows.extend(entryrows)
      except Cancelled:raise
      except Exception as e:
       failures.append({'parent_id':pid,'state_id':sid,'seed':seed,'original_name':item['name'],'error_type':type(e).__name__,'error':str(e)})
       (dest/'failure.txt').write_text(traceback.format_exc())
      m['counts']={'successful_entries':len(m['completed_entries']),'failed_entries':len(failures)};m['failures']=failures
      write_json(out/'run_manifest.json',m);update(m)
  m['status']='partial_failed' if rows and failures else 'completed' if rows else 'failed'
 except Cancelled:m['status']='cancelled'
 except BaseException:
  m['status']='interrupted';write_json(out/'run_manifest.json',m);raise
 states,ranking=aggregate(rows,cfg['engine']);csvwrite(out/'poses.csv',rows);csvwrite(out/'states.csv',states);csvwrite(out/'ranking.csv',ranking);csvwrite(out/'failures.csv',failures,['parent_id','state_id','seed','original_name','error_type','error']);csvwrite(out/'ligand_manifest.csv',prepared)
 m.update(finished=time.time(),elapsed_seconds=time.time()-m['started'],failures=failures,counts={'successful_entries':len(m['completed_entries']),'failed_entries':len(failures),'successful_parents':len(ranking)},artifacts={p.name:sha(p) for p in out.glob('*.csv')})
 write_json(out/'run_manifest.json',m);update(m)
 (out/'REPLAY.txt').write_text('Replay with a NEW run: ./run.sh screen --config '+str((out/'config_expanded.yaml').relative_to(ROOT))+' --ligands '+str((out/m['input']).relative_to(ROOT))+'\nResume only verified completed entries: ./run.sh resume --run '+str(out.relative_to(ROOT))+'\n')
 return m
