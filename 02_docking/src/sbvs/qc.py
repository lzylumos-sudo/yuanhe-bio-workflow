import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from .preparation.ligands import identity
def pose_qc(m,reference_identity,box,receptor_xyz=None):
 flags=[];heavy=Chem.RemoveHs(m);xyz=heavy.GetConformer().GetPositions()
 if identity(m)!=reference_identity:flags.append('chemical_identity_or_stereo_changed')
 if not np.isfinite(xyz).all():raise ValueError('Nonfinite pose coordinates')
 if np.any(abs(xyz-np.asarray(box['center']))>np.asarray(box['size'])/2+2):flags.append('outside_box_tolerance_2A')
 if receptor_xyz is not None:
  from scipy.spatial import cKDTree
  distances=cKDTree(receptor_xyz).query(xyz)[0]
  if distances.min()<1.2:flags.append('severe_clash_below_1.2A')
 return flags
def fixed_frame_rmsd(probe,reference):
 p=Chem.RemoveHs(probe);r=Chem.RemoveHs(reference)
 if identity(p)!=identity(r):raise ValueError('RMSD chemical identities differ')
 # CalcRMS enumerates symmetry mappings without aligning coordinates.
 return float(rdMolAlign.CalcRMS(p,r,maxMatches=100000))
