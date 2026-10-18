"""Bind generic screening to the uploaded receptor, never an unrelated demo PDB."""
from pathlib import Path
import shutil

def bind_generic(workspace,run_dir,config):
 import yaml
 source=Path(run_dir)/'input.pdb'
 if not source.is_file():raise ValueError('Uploaded input snapshot is missing')
 folder=Path(workspace)/'input_structure';folder.mkdir()
 shutil.copy2(source,folder/'USER_INPUT.pdb')
 profile=Path(workspace)/'configs/targets/generic.yaml'
 data=yaml.safe_load(profile.read_text())
 data['structures']={'primary':{'pdb_id':'USER_INPUT','chain':config['receptor_chain'],'role':'Exact uploaded input; no external accession inferred'}}
 data['construct']={'native_range':config['receptor_range'],'label_native_offset':0,'note':'Identity map in uploaded author numbering; not a UniProt mapping'}
 data['structures_dir']=str(folder)
 data['native_length']=0;data['native_sequence']='';data['uniprot']=''
 profile.write_text(yaml.safe_dump(data,sort_keys=False))
 return folder
