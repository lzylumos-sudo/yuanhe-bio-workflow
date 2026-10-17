import csv,json,math
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation,PDBQTWriterLegacy
from ..paths import write_json
ALLOWED={1,6,7,8,9,15,16,17,35,53}
def inputs(path):
 path=Path(path)
 if path.suffix.lower()=='.sdf':
  for i,m in enumerate(Chem.SDMolSupplier(str(path),removeHs=False),1):
   yield {'index':i,'name':m.GetProp('_Name') if m and m.HasProp('_Name') else f'molecule_{i}','mol':m,'source':'sdf','smiles':Chem.MolToSmiles(Chem.RemoveHs(m)) if m else ''}
 elif path.suffix.lower()=='.csv':
  with path.open(encoding='utf-8-sig',newline='') as f:
   for i,row in enumerate(csv.DictReader(f),1):
    smi=row.get('smiles','').strip();yield {'index':i,'name':row.get('name') or f'molecule_{i}','smiles':smi,'mol':Chem.MolFromSmiles(smi),'source':'smiles'}
 else:
  with path.open(encoding='utf-8-sig') as f:
   for i,line in enumerate(f,1):
    if not line.strip() or line.lstrip().startswith('#'):continue
    a=line.strip().split(maxsplit=1);yield {'index':i,'name':a[1] if len(a)>1 else f'molecule_{i}','smiles':a[0],'mol':Chem.MolFromSmiles(a[0]),'source':'smiles'}
def identity(m):return Chem.MolToSmiles(Chem.RemoveHs(m),isomericSmiles=True)
def states(item,cfg):
 if cfg.get('ligand_states','preserve')=='preserve' or item['mol'] is None:return [item]
 from molscrub import Scrub
 from rdkit.Chem.MolStandardize import rdMolStandardize
 mol=Chem.Mol(item['mol']);fragment=cfg.get('fragment_policy','reject')
 if len(Chem.GetMolFrags(mol))>1:
  if fragment!='largest':raise ValueError('MULTIFRAGMENT_REQUIRES_EXPLICIT_LARGEST_POLICY')
  mol=rdMolStandardize.LargestFragmentChooser().choose(mol)
 scrub=Scrub(ph_low=cfg.get('ph_low',6.4),ph_high=cfg.get('ph_high',8.4),skip_gen3d=True)
 generated={identity(m):m for m in scrub(mol)}
 if not generated:raise ValueError('STATE_ENUMERATION_EMPTY')
 limit=int(cfg.get('max_states',8));selected=sorted(generated)[:limit];result=[]
 for smi in selected:
  m=Chem.Mol(generated[smi]);m.RemoveAllConformers()
  result.append(dict(item,mol=m,source='enumerated',state_enumeration={'tool':'molscrub 0.2.2','ph_low':cfg.get('ph_low',6.4),'ph_high':cfg.get('ph_high',8.4),'input_identity':identity(item['mol']),'state_identity':smi,'generated_states':len(generated),'retained_states':len(selected),'selection':'canonical isomeric SMILES stable sort capped at max_states; not a population ranking','fragment_policy':fragment,'tautomer_policy':'Molscrub rules','stereo_policy':'preserve specified; unspecified not assigned as known'}))
 return result
def conformers(item,count=8,retain=3,seed=20261012):
 if not item['mol']:return [item]
 if not 1<=retain<=count<=8:raise ValueError('Conformer limits require 1<=retain<=count<=8')
 m=Chem.AddHs(Chem.Mol(item['mol']));m.RemoveAllConformers();p=AllChem.ETKDGv3();p.randomSeed=seed;p.pruneRmsThresh=.5
 ids=list(AllChem.EmbedMultipleConfs(m,numConfs=count,params=p))
 if not ids:raise ValueError('EMBED_FAILED')
 if AllChem.MMFFHasAllMoleculeParams(m):scores=AllChem.MMFFOptimizeMoleculeConfs(m,numThreads=1,maxIters=1000);method='MMFF94'
 elif AllChem.UFFHasAllMoleculeParams(m):scores=AllChem.UFFOptimizeMoleculeConfs(m,numThreads=1,maxIters=1000);method='UFF_fallback'
 else:raise ValueError('NO_FORCEFIELD_PARAMETERS')
 ordered=sorted(zip(ids,scores),key=lambda x:(x[1][1],x[0]))[:retain];result=[]
 for cid,(status,energy) in ordered:
  copy=Chem.Mol(m);copy.RemoveAllConformers();copy.AddConformer(Chem.Conformer(m.GetConformer(cid)),assignId=True)
  result.append(dict(item,mol=copy,source='sdf',conformer_generation={'requested':count,'generated':len(ids),'retained':len(ordered),'method':'ETKDGv3_'+method,'seed':seed,'original_conformer_id':cid,'energy_kcal_mol':energy,'optimization_status':status,'dedup':'ETKDG prune RMS 0.5 A','selection':'lowest minimized force-field energy'}))
 return result
def prepare(item,out,seed=20261012):
 if item.get('preparation_error'):raise ValueError(item['preparation_error'])
 out=Path(out);out.mkdir(parents=True,exist_ok=True);m=item['mol']
 if m is None or m.GetNumAtoms()==0:raise ValueError('INVALID_SMILES_OR_MOLECULE')
 Chem.SanitizeMol(m)
 if any(a.GetAtomicNum() not in ALLOWED for a in m.GetAtoms()):raise ValueError('UNSUPPORTED_ELEMENT')
 if len(Chem.GetMolFrags(m))!=1:raise ValueError('MULTIFRAGMENT_REQUIRES_REVIEW')
 before=identity(m);unassigned=Chem.FindMolChiralCenters(m,includeUnassigned=True)
 flags=['unspecified_stereochemistry'] if any(s=='?' for _,s in unassigned) else []
 keep=item['source']=='sdf' and m.GetNumConformers()>0 and m.GetConformer().Is3D()
 m=Chem.AddHs(m,addCoords=keep)
 method='preserved_input_3d'
 if not keep:
  p=AllChem.ETKDGv3();p.randomSeed=int(seed)
  if AllChem.EmbedMolecule(m,p)!=0:raise ValueError('EMBED_FAILED')
  if AllChem.MMFFHasAllMoleculeParams(m):status=AllChem.MMFFOptimizeMolecule(m,maxIters=1000);method='ETKDGv3_MMFF94'
  elif AllChem.UFFHasAllMoleculeParams(m):status=AllChem.UFFOptimizeMolecule(m,maxIters=1000);method='ETKDGv3_UFF_fallback'
  else:raise ValueError('NO_FORCEFIELD_PARAMETERS')
  if status:flags.append('optimization_not_converged')
 xyz=m.GetConformer().GetPositions()
 if not np.isfinite(xyz).all() or (m.GetNumAtoms()>1 and np.max(np.ptp(xyz,axis=0))<.1):raise ValueError('INVALID_3D')
 if identity(m)!=before:raise ValueError('CHEMICAL_IDENTITY_CHANGED')
 m.SetProp('_Name',item['name']);w=Chem.SDWriter(str(out/'input.sdf'));w.write(m);w.close()
 setups=MoleculePreparation().prepare(m)
 if len(setups)!=1:raise ValueError('Unexpected Meeko setup count')
 text,ok,error=PDBQTWriterLegacy.write_string(setups[0])
 if not ok:raise ValueError('MEEKO_PREP_FAILED '+str(error))
 (out/'input.pdbqt').write_text(text)
 info={'name':item['name'],'original_smiles':item['smiles'],'identity':before,'formal_charge':Chem.GetFormalCharge(m),'method':method,'seed':seed,'qc_flags':flags,'state_policy':'preserve','source_type':item['source']}
 if 'conformer_generation' in item:info['conformer_generation']=item['conformer_generation']
 if 'state_enumeration' in item:info['state_enumeration']=item['state_enumeration'];info['state_policy']='enumerate'
 write_json(out/'preparation.json',info);return info
