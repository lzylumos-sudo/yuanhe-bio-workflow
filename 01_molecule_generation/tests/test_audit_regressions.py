import json,importlib.util
from pathlib import Path
import pytest
from src.adapters import peptide2mol,reinvent4
from src.structure import prepare_structure

def test_scaffold_labels_removed_by_official_sampler():
 assert reinvent4.scaffold_identity('c1ccc([*:1])cc1')==reinvent4.scaffold_identity('*c1ccccc1')
 assert reinvent4.scaffold_identity('*c1ccccc1')!=reinvent4.scaffold_identity('*C1CCCCC1')

def test_disconnected_candidates_rejected():
 assert reinvent4._core_retained('c1ccc([*:1])cc1','Cc1ccccc1.C')[1]=='disconnected_molecule'

def test_hydra_path_quoted():
 path='/tmp/中文 项目/a=b/model.ckpt'
 assert json.loads(peptide2mol.hydra_string(path))==path

def test_existing_model_outputs_cannot_be_reused(tmp_path):
 out=tmp_path/'out';out.mkdir();(out/'stale.sdf').write_text('stale')
 with pytest.raises(FileExistsError):peptide2mol.run(pocket_pdb=tmp_path/'x',checkpoint=tmp_path/'x',checkpoint_sha256='0'*64,output_dir=out,budget=1,seed=1)

def test_full_receptor_separate_from_model_pocket(tmp_path):
 root=Path(__file__).resolve().parents[2]
 result=prepare_structure(root/'data/structures/9nu4.pdb',tmp_path,crop_distance_A=6.)
 full=result['receptor_full'].read_text();pocket=result['receptor_prepared'].read_text()
 assert full.count('ATOM')>pocket.count('ATOM')
 assert all(l[21:22]=='A' for l in full.splitlines() if l.startswith('ATOM'))
 assert 'receptor_full.pdb' in result['checksums'].read_text()
