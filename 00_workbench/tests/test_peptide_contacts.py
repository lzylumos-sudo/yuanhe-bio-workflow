import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scientific.peptide_contacts import compare

def pdb(path,positions):
 path.write_text(''.join(f'ATOM  {i:5d}  CA  ALA A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n' for i,(x,y,z) in enumerate(positions,1)))

def test_standalone_contacts_and_geometric_controls(tmp_path):
 Chem=pytest.importorskip('rdkit.Chem')
 receptor=tmp_path/'r.pdb';pep=tmp_path/'p.pdb';pdb(receptor,[(0,0,0),(10,0,0),(20,0,0)]);pdb(pep,[(1,0,0)])
 mol=Chem.MolFromSmiles('C');conf=Chem.Conformer(1);conf.SetAtomPosition(0,(1.,0.,0.));mol.AddConformer(conf)
 sdf=tmp_path/'m.sdf';writer=Chem.SDWriter(str(sdf));writer.write(mol);writer.close()
 rows=[dict(compound_id='fixture',pose_id='pose',sdf=str(sdf))];selection={'primary':[{'compound_id':'fixture','raw':{'affinity':-1}}]}
 result=compare(receptor,receptor,pep,rows,selection,tmp_path/'result')
 assert result['controls']['self']['jaccard']==1
 assert result['controls']['translated_100A_beyond_receptor']['jaccard']==0
 assert result['candidates'][0]['recall']==1


def test_coordinate_frame_mismatch_is_refused(tmp_path):
 pytest.importorskip('rdkit')
 r=tmp_path/'r.pdb';original=tmp_path/'o.pdb';p=tmp_path/'p.pdb'
 pdb(original,[(0,0,0),(10,0,0),(20,0,0)]);pdb(r,[(1,0,0),(11,0,0),(21,0,0)]);pdb(p,[(1,0,0)])
 with pytest.raises(ValueError,match='coordinate frames'):compare(r,original,p,[],{},tmp_path/'result')
