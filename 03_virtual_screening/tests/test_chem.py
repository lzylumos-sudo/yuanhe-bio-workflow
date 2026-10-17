"""Chemistry layer: identity, recovery, properties, scaffolds, geometry."""

from __future__ import annotations

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from zvss.chem import rdio
from zvss.chem.props import alert_burden, compute_properties
from zvss.chem.recover import recover_from_smiles


def test_canonical_smiles_is_canonical() -> None:
    a = rdio.mol_from_smiles("OCC")
    b = rdio.mol_from_smiles("CCO")
    assert a is not None and b is not None
    assert rdio.canonical_isomeric_smiles(a) == rdio.canonical_isomeric_smiles(b)


def test_properties_flag_that_logd_is_unavailable() -> None:
    mol = rdio.mol_from_smiles("c1ccccc1O")
    assert mol is not None
    panel = compute_properties(mol)
    assert panel.values["heavy_atoms"] == 7
    assert panel.values["molecular_weight"] > 90
    assert panel.methods["logd_available"] == "false"
    assert "NOT logD" in panel.methods["clogp_crippen"]
    assert isinstance(alert_burden(panel), float)


def test_murcko_scaffold_of_phenol_is_benzene() -> None:
    mol = rdio.mol_from_smiles("c1ccccc1O")
    assert mol is not None
    assert rdio.murcko_scaffold(mol) == "c1ccccc1"


def test_recover_from_smiles_restores_identity_from_pdb_block() -> None:
    template = "c1ccc(cc1)C(=O)O"
    mol = rdio.mol_from_smiles(template)
    assert mol is not None
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=7)
    block = Chem.MolToPDBBlock(mol)
    repaired = recover_from_smiles(template, block)
    assert repaired.mol is not None
    assert rdio.canonical_isomeric_smiles(Chem.RemoveHs(repaired.mol)) == rdio.canonical_isomeric_smiles(
        rdio.mol_from_smiles(template)
    )
    assert repaired.atom_map
    assert repaired.status in {"pass", "flag"}


def test_rmsd_of_identical_geometry_is_zero() -> None:
    mol = rdio.mol_from_smiles("CCO")
    assert mol is not None
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=11)
    coords = rdio.coordinates(mol)
    elements = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    assert rdio.symmetry_corrected_rmsd_numpy(coords, coords, elements) == pytest.approx(0.0, abs=1e-6)


def test_rmsd_grows_with_a_real_offset() -> None:
    mol = rdio.mol_from_smiles("CCO")
    assert mol is not None
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=11)
    coords = rdio.coordinates(mol)
    shifted = [(x + 0.5, y, z) for x, y, z in coords]
    elements = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    value = rdio.symmetry_corrected_rmsd_numpy(coords, shifted, elements)
    assert value is not None and value > 0.1
