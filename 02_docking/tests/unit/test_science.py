import pytest,numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sbvs.config import validate_box
from sbvs.preparation.ligands import prepare,inputs,identity
from sbvs.qc import fixed_frame_rmsd
from sbvs.paths import rootpath,ROOT
def molecule(s):
 m=Chem.AddHs(Chem.MolFromSmiles(s));p=AllChem.ETKDGv3();p.randomSeed=17;assert AllChem.EmbedMolecule(m,p)==0;return m
def test_rmsd_no_alignment():
 m=molecule('CCO');n=Chem.Mol(m);c=n.GetConformer()
 for i in range(n.GetNumAtoms()):v=c.GetAtomPosition(i);c.SetAtomPosition(i,(v.x+5,v.y,v.z))
 assert fixed_frame_rmsd(n,m)==pytest.approx(5,abs=1e-5)
def test_rmsd_symmetry():
 m=molecule('CC(=O)O');n=Chem.RenumberAtoms(m,list(reversed(range(m.GetNumAtoms()))));assert fixed_frame_rmsd(n,m)==pytest.approx(0,abs=1e-5)
@pytest.mark.parametrize('size',[[0,1,2],[-1,2,3],[float('nan'),2,3],[float('inf'),2,3]])
def test_bad_box(size):
 with pytest.raises(ValueError):validate_box({'center':[0,0,0],'size':size})
def test_no_path_escape():
 with pytest.raises(ValueError):rootpath('../outside')
def test_names_and_failures(tmp_path):
 p=tmp_path/'input.smi';p.write_text('CCO Name Space\nCCC Name Space\nINVALID bad\n[Fe] metal\n',encoding='utf-8');items=list(inputs(p));assert len(items)==4 and items[0]['name']=='Name Space'
 assert prepare(items[0],tmp_path/'a')['formal_charge']==0
 with pytest.raises(ValueError,match='INVALID'):prepare(items[2],tmp_path/'bad')
 with pytest.raises(ValueError,match='UNSUPPORTED'):prepare(items[3],tmp_path/'metal')
def test_chirality_charge_preserved(tmp_path):
 for i,s in enumerate(['C[C@H](O)C(=O)O','C[NH3+]']):
  m=Chem.MolFromSmiles(s);r=prepare({'mol':m,'name':'test','smiles':s,'source':'smiles'},tmp_path/str(i));out=next(iter(Chem.SDMolSupplier(str(tmp_path/str(i)/'input.sdf'),removeHs=False)));assert identity(out)==identity(m)
def test_missing_receptor():
 from sbvs.config import load,validate
 import yaml
 c=yaml.safe_load((ROOT/'configs/class_quick.yaml').read_text());c['cpu_per_job']=1;c['receptor']='data/nonexistent_receptor.pdbqt'
 with pytest.raises(ValueError,match='Missing receptor'):validate(c)
def test_external_pdbqt_provenance(tmp_path):
 from sbvs.preparation.receptor import plan
 source=tmp_path/'external.pdbqt';source.write_text('REMARK external receptor fixture\n')
 p=plan(source);assert p['external'] and p['protonation']=='external/unknown' and not p['confirmed']
