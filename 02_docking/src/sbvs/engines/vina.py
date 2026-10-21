import re,math
from pathlib import Path
from meeko import PDBQTMolecule,RDKitMolCreate
from rdkit import Chem
from .common import execute
from ..paths import rootpath
def dock(cfg,box,folder,seed,cancel):
 folder=Path(folder);out=folder/'poses.pdbqt'
 from sbvs2.docking.vina_official import vina_version
 executable=cfg.get('executable','vina');vina_version(executable)
 argv=[executable,'--receptor',str(rootpath(cfg['receptor'])),'--ligand',str(folder/'input.pdbqt'),'--scoring','vina']
 for i,axis in enumerate('xyz'):argv+=['--center_'+axis,str(box['center'][i]),'--size_'+axis,str(box['size'][i])]
 for k,v in [('exhaustiveness',cfg['exhaustiveness']),('num_modes',cfg['num_modes']),('energy_range',cfg.get('energy_range',3.0)),('cpu',cfg['cpu_per_job']),('seed',seed),('out',out)]:argv+=['--'+k,str(v)]
 execute(argv,folder,cfg['per_ligand_timeout_seconds'],cancel)
 text=out.read_text();scores=[float(x) for x in re.findall(r'REMARK VINA RESULT:\s+([-+\d.eE]+)',text)]
 if not scores or not all(math.isfinite(x) for x in scores):raise RuntimeError('No finite Vina scores')
 from sbvs2.docking.vina_official import export_sdf
 source=Chem.SDMolSupplier(str(folder/'input.sdf'),removeHs=False)[0]
 export_sdf(out,source,folder/'poses.sdf',[{'vina_affinity':score,'pose_id':i} for i,score in enumerate(scores,1)])
 mols=list(Chem.SDMolSupplier(str(folder/'poses.sdf'),removeHs=False))
 if len(mols)!=len(scores) or any(m is None for m in mols):raise RuntimeError('Incomplete Vina SDF export')
 return [{'mol':m,'vina_affinity_kcal_mol':s,'empirical_score_kcal_mol':None,'cnn_pose_score':None,'cnn_affinity_pK':None,'model_id':None} for m,s in zip(mols,scores)]
