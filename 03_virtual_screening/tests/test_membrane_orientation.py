from pathlib import Path
from types import SimpleNamespace
import json, shutil
import numpy as np
import pytest
from zvss.membrane.orientation import transfer, atom_map
from zvss.select.qualify import _membrane

EVIDENCE = Path(__file__).resolve().parents[1]/'configs/membrane/9nu4'

def test_real_ppm_transfer_preserves_placement_under_rigid_transform(tmp_path):
    baseline = transfer(EVIDENCE, EVIDENCE/'murj.pdb')
    rot = np.array([[0.,-1,0],[1,0,0],[0,0,1]])
    shift = np.array([13.,-7,4])
    lines=[]
    for line in (EVIDENCE/'murj.pdb').read_text().splitlines():
        if line.startswith(('ATOM  ','HETATM')):
            xyz = np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]) @ rot + shift
            line = line[:30]+''.join(f'{v:8.3f}' for v in xyz)+line[54:]
        lines.append(line)
    target=tmp_path/'rotated.pdb';target.write_text('\n'.join(lines)+'\n')
    moved=transfer(EVIDENCE,target)
    assert moved['coverage']==1
    assert moved['rigid_fit_rmsd_A'] < .001
    assert np.allclose(moved['normal'],np.asarray(baseline['normal'])@rot,atol=1e-6)
    probe=np.array([160,154,135])
    before=(probe-baseline['center'])@np.asarray(baseline['normal'])
    after=(probe@rot+shift-moved['center'])@np.asarray(moved['normal'])
    assert after==pytest.approx(before,abs=.002)


def test_tampered_evidence_rejected(tmp_path):
    shutil.copytree(EVIDENCE,tmp_path/'ppm')
    with (tmp_path/'ppm/murjout.pdb').open('a') as f:f.write('REMARK edited\n')
    with pytest.raises(ValueError,match='checksum'):
        transfer(tmp_path/'ppm',EVIDENCE/'murj.pdb')


def test_different_construct_is_not_given_a_membrane_pass(tmp_path):
    rows=(EVIDENCE/'murj.pdb').read_text().splitlines()
    target=tmp_path/'truncated.pdb';target.write_text('\n'.join(rows[:100]))
    with pytest.raises(ValueError,match='construct'):
        transfer(EVIDENCE,target)


@pytest.mark.parametrize('verified,depth,status',[(True,2.4,'pass'),(False,2.4,'not_evaluated'),(True,None,'not_evaluated'),(True,float('nan'),'not_evaluated')])
def test_orientation_requires_frame_verified_measurement(tmp_path,verified,depth,status):
    (tmp_path/'summary.json').write_text(json.dumps({'status':'sourced','orientation_verified':verified,
        'md_capability':{'status':'separate_execution_required'}}))
    rows=[] if depth is None else [{'metric':'insertion_depth','value':depth}]
    ctx=SimpleNamespace(stage_dir=lambda _:str(tmp_path),db=SimpleNamespace(run_id='r',query=lambda *args:rows))
    assert _membrane(ctx,{'pose_id':'pose'})['status']==status
