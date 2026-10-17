import numpy as np
from rdkit import Chem
from ..paths import sha,write_json
from ..config import validate_box
def make_box(reference,padding=5,out=None):
 mols=[m for m in Chem.SDMolSupplier(str(reference),removeHs=False) if m]
 if len(mols)!=1:raise ValueError('Reference must contain exactly one valid instance')
 m=Chem.RemoveHs(mols[0]);x=m.GetConformer().GetPositions();center=x.mean(axis=0);size=2*np.max(abs(x-center),axis=0)+2*padding
 b={'center':center.tolist(),'size':size.tolist(),'padding_per_side':padding,'method':'geometric_center_max_deviation','reference_sha256':sha(reference),'units':'angstrom'}
 validate_box(b)
 if not np.all(abs(x-center)<=size/2):raise ValueError('Reference outside box')
 if out:
  from pathlib import Path
  write_json(out,b)
  (Path(out).parent/'vina_conf.txt').write_text('\n'.join(f'{key}_{axis} = {b[key][i]:.8f}' for key in ('center','size') for i,axis in enumerate('xyz'))+'\n')
 return b
