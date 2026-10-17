import math
from pathlib import Path
from rdkit import Chem
from .common import execute
from ..paths import rootpath
def value(m,aliases,required=False):
 for k in aliases:
  if m.HasProp(k):
   x=float(m.GetProp(k))
   if not math.isfinite(x):raise ValueError('Nonfinite GNINA property '+k)
   return x
 if required:raise ValueError('GNINA required property absent '+repr(aliases))
 return None
def dock(cfg,box,folder,seed,cancel):
 folder=Path(folder);out=folder/'poses.sdf';device=cfg['device']
 if device=='auto':raise ValueError('Resolve GNINA device using verified capabilities before submitting')
 argv=['/opt/sbvs/vendor/gnina','-r',str(rootpath(cfg['receptor_pdb'])),'-l',str(folder/'input.sdf'),'-o',str(out),'--scoring',cfg.get('empirical_scoring','vina'),'--cnn_scoring','rescore','--cnn',cfg['cnn_model'],'--cpu',str(cfg['cpu_per_job']),'--seed',str(seed)]
 for i,axis in enumerate('xyz'):argv+=['--center_'+axis,str(box['center'][i]),'--size_'+axis,str(box['size'][i])]
 op=cfg.get('operation','redock')
 if op=='score_only':argv+=['--score_only']
 elif op=='refinement':argv+=['--minimize']
 elif op=='redock':argv+=['--exhaustiveness',str(cfg['exhaustiveness']),'--num_modes',str(cfg['num_modes'])]
 else:raise ValueError('Unknown GNINA operation')
 argv+=['--no_gpu'] if device=='cpu' else ['--device',str(cfg.get('gpu_device',0))]
 execute(argv,folder,cfg['per_ligand_timeout_seconds'],cancel)
 result=[]
 for m in Chem.SDMolSupplier(str(out),removeHs=False):
  if m is None:raise ValueError('Invalid GNINA SDF record')
  result.append({'mol':m,'vina_affinity_kcal_mol':None,'empirical_score_kcal_mol':value(m,['minimizedAffinity','affinity'],True),'cnn_pose_score':value(m,['CNNscore','CNN_score'],True),'cnn_affinity_pK':value(m,['CNNaffinity','CNN_affinity'],True),'model_id':cfg['cnn_model'],'raw_properties':m.GetPropsAsDict()})
 if not result:raise ValueError('Empty GNINA output')
 return result
