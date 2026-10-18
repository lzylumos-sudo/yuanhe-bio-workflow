from rdkit import Chem
from rdkit.Chem import AllChem
from zvss.redock.gnina import _best_rmsd

def molecule():
 m=Chem.AddHs(Chem.MolFromSmiles('CC(=O)Oc1ccccc1'))
 assert AllChem.EmbedMolecule(m,randomSeed=17)==0
 return m

def test_fixed_frame_translation_is_not_fitted_away():
 m=molecule();p=Chem.Mol(m)
 for i in range(p.GetNumAtoms()):
  q=p.GetConformer().GetAtomPosition(i);q.x+=10;p.GetConformer().SetAtomPosition(i,q)
 assert abs(_best_rmsd(p,m)-10)<1e-5
 assert _best_rmsd(m,m)<1e-5

def test_atom_reordering_preserves_symmetry_correct_rmsd():
 m=molecule();p=Chem.RenumberAtoms(m,list(reversed(range(m.GetNumAtoms()))))
 assert _best_rmsd(p,m)<1e-5

def test_different_molecules_do_not_receive_a_comparison():
 assert _best_rmsd(molecule(),Chem.MolFromSmiles('CC')) is None
