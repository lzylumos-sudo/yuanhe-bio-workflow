import uuid
from rdkit import Chem
from .paths import ROOT,rootpath,read_json,write_json
from .export import read_rows
from .config import load
from .runner import create
def vina_stage(run,operation='redock',top_k=10,device='cpu'):
 if operation != 'redock':raise ValueError('Unsupported operation')
 if not 1<=top_k<=10000:raise ValueError('Invalid TopK')
 old=rootpath(run);rows=read_rows(old/'ranking.csv')[:top_k]
 if not rows:raise ValueError('No successful parent compounds')
 selection=ROOT/'data/selections'/uuid.uuid4().hex;selection.mkdir(parents=True);sdf=selection/'selected.sdf';w=Chem.SDWriter(str(sdf))
 for row in rows:
  mol=next(iter(Chem.SDMolSupplier(str(old/row['pose_file']),removeHs=False)));mol.SetProp('_Name',row['original_name'])
  if operation=='redock':mol.RemoveAllConformers()
  w.write(mol)
 w.close();selection_info={'source_run':str(old.relative_to(ROOT)),'operation':operation,'basis':'original parent ranking TopK','top_k':top_k,'poses':rows}
 write_json(selection/'selection.json',selection_info)
 cfg=load('research_vina');original=read_json(old/'run_manifest.json')['config']
 for key in ('receptor','receptor_pdb','box'):cfg[key]=original[key]
 cfg.update(device=device,operation=operation)
 out=create(sdf,cfg);write_json(out/'stage_selection.json',selection_info);return out
