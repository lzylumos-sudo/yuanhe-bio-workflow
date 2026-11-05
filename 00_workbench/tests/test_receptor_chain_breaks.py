"""Public receptor gaps must never become long covalent peptide bonds."""
import pytest
from studio_worker import preserve_receptor_breaks


def atom(serial, name, residue, x):
    return f'ATOM  {serial:5d} {name:^4s} ALA A{residue:4d}    {x:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00          {name[0]:>2s}\n'


def test_numbered_gap_is_disconnected_without_moving_atoms(tmp_path):
    from openmm.app import PDBFile
    source=tmp_path/'source.pdb';target=tmp_path/'segmented.pdb'
    atoms=[atom(1,'N',215,0),atom(2,'CA',215,1.4),atom(3,'C',215,2.8),
           atom(4,'N',377,12.6),atom(5,'CA',377,14),atom(6,'C',377,15.4)]
    source.write_text(''.join(atoms)+'END\n')
    report=preserve_receptor_breaks(source,target)
    assert report[0]['missing_numbered_positions']==161
    assert report[0]['C_N_distance_A']==pytest.approx(9.8)
    assert [s+'\n' for s in target.read_text().splitlines() if s.startswith('ATOM')]==atoms
    pdb=PDBFile(str(target))
    assert len(list(pdb.topology.chains()))==2
    assert all(a.residue.chain is b.residue.chain for a,b in pdb.topology.bonds())
    assert [r.id for r in pdb.topology.residues()]==['215','377']


def test_numbering_gap_with_real_peptide_geometry_stays_connected(tmp_path):
    source=tmp_path/'source.pdb';target=tmp_path/'segmented.pdb'
    source.write_text(atom(1,'C',215,0)+atom(2,'N',377,1.33)+'END\n')
    assert preserve_receptor_breaks(source,target)==[]
    assert target.read_bytes()==source.read_bytes()


def test_long_bond_without_numbering_gap_requires_coordinate_review(tmp_path):
    source=tmp_path/'source.pdb'
    source.write_text(atom(1,'C',215,0)+atom(2,'N',216,9.8)+'END\n')
    with pytest.raises(ValueError,match='Invalid receptor peptide geometry'):
        preserve_receptor_breaks(source,tmp_path/'segmented.pdb')


def test_existing_ter_is_preserved(tmp_path):
    source=tmp_path/'source.pdb';target=tmp_path/'segmented.pdb'
    source.write_text(atom(1,'C',215,0)+'TER\n'+atom(2,'N',377,9.8)+'END\n')
    assert preserve_receptor_breaks(source,target)==[]
    assert target.read_bytes()==source.read_bytes()
