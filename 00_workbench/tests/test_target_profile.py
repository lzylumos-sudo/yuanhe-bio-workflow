import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scientific.target_profile import bind_generic

def test_generic_profile_uses_uploaded_chain_and_file(tmp_path):
 yaml=pytest.importorskip('yaml')
 d=tmp_path/'run';d.mkdir();(d/'input.pdb').write_text('original uploaded data')
 w=tmp_path/'workspace';p=w/'configs/targets/generic.yaml';p.parent.mkdir(parents=True);p.write_text('target_id: generic\nstructures: {primary: {pdb_id: 1DRF, chain: A}}\n')
 folder=bind_generic(w,d,{'receptor_chain':'C','receptor_range':[5,200]})
 data=yaml.safe_load(p.read_text())
 assert data['structures']['primary']['pdb_id']=='USER_INPUT'
 assert data['structures']['primary']['chain']=='C'
 assert (folder/'USER_INPUT.pdb').read_text()=='original uploaded data'
 assert data['uniprot']==''
