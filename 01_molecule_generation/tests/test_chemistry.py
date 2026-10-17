"""Chemistry identity, SDF fidelity and export validation (RDKit only)."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.chemistry import (ChemistryError, STANDARDIZATION_VERSION, calculate_properties, compound_id,
                           parse_smiles, pose_id, provenance_id, read_sdf, standardize_smiles,
                           state_id, write_sdf)
from src.export import LIBRARY_FILENAME, REQUIRED_FIELDS, export_library, validate_export

CHECKSUM = "7dc82c4f895efea83764ce25c2fae5787f9f25b7935590e5119bde7199540129"


def embed(smiles: str, seed: int = 7) -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    return mol


@pytest.mark.parametrize("bad", ["", "C1CC", "CC.O", "*CC", "[*:1]CC", "C(C)(C)(C)(C)C", "CC name"])
def test_invalid_input_rejected(bad):
    with pytest.raises(ChemistryError):
        parse_smiles(bad)


def test_ids_are_deterministic_and_atom_order_independent():
    assert compound_id("OCC") == compound_id("C(O)C")
    assert compound_id("CCO") != compound_id("CCN")
    assert state_id("CC(=O)[O-]") != state_id("CC(=O)O")  # protonation state is preserved
    assert standardize_smiles("C[C@H](N)C(=O)O") != standardize_smiles("C[C@@H](N)C(=O)O")


@pytest.mark.parametrize("smiles", [
    "c1ccccc1C#CC", "CC(=O)NC1CC1", "CC(C)C(=O)[O-]", "C[C@H](N)C(=O)[O-]", "C[C@@H](N)C(=O)[O-]",
    "Cc1ncc[nH]1", "O=C1CCCN1c1ccncc1"])
def test_sdf_roundtrip_preserves_graph_charge_stereo_and_coords(tmp_path, smiles):
    mol = embed(smiles)
    path = write_sdf(mol, tmp_path / "m.sdf", {"pose_id": "x"})
    back = read_sdf(path)
    assert len(back) == 1
    assert standardize_smiles(back[0]) == standardize_smiles(smiles)
    assert Chem.GetFormalCharge(back[0]) == Chem.GetFormalCharge(Chem.MolFromSmiles(smiles))
    a = mol.GetConformer().GetPositions()
    b = back[0].GetConformer().GetPositions()
    assert abs(a - b).max() < 1e-3  # SDF stores 4 decimals


def test_write_requires_coordinates(tmp_path):
    with pytest.raises(ChemistryError):
        write_sdf(Chem.MolFromSmiles("CCO"), tmp_path / "m.sdf")


def test_legacy_invalid_sdf_rejected(tmp_path):
    legacy = tmp_path / "legacy.sdf"
    # Old pepwedge format: two header lines, non-integer counts, no bond block.
    legacy.write_text("mol\npepwedge\n 62.  0.  0  0  0  0  0  0  0  0999 V2000\n"
                      "    0.0000    0.0000    0.0000 C   0  0\nM  END\n$$$$\n")
    with pytest.raises(ChemistryError):
        read_sdf(legacy)
    truncated = tmp_path / "trunc.sdf"
    truncated.write_text("")
    with pytest.raises(ChemistryError):
        read_sdf(truncated)


def test_properties_match_rdkit():
    props = calculate_properties("CC(=O)[O-]")
    assert props["formal_charge"] == -1 and props["HBD"] == 0


def test_pose_id_changes_with_coordinates_not_atom_order():
    base = embed("CCO")
    other = embed("CCO", seed=11)
    state = state_id("CCO")
    assert pose_id(base, state) != pose_id(other, state)
    assert pose_id(base, state) == pose_id(Chem.Mol(base), state)


def _library(tmp_path, **overrides):
    smiles = standardize_smiles("CC(C)Cc1ccc(cc1)C(C)C(=O)[O-]")
    cid, sid = compound_id(smiles), state_id(smiles)
    mol = embed(smiles)
    prov = {"model": "REINVENT4", "checkpoint_sha256": CHECKSUM, "run_id": "r1", "seed": "41"}
    prid = provenance_id(prov)
    prov["provenance_id"] = prid
    pid = pose_id(mol, sid, "")
    rel = "sdf/a.sdf"
    write_sdf(mol, tmp_path / rel, {"compound_id": cid, "state_id": sid, "pose_id": pid, "provenance_id": prid})
    row = {"compound_id": cid, "canonical_smiles": smiles, "standardization_version": STANDARDIZATION_VERSION,
           "state_id": sid, "pose_id": pid, "provenance_id": prid, "source": "REINVENT4",
           "route_status": "supported", "geometry_status": "pass", "qc_status": "pass", "sdf_path": rel,
           **{k: v for k, v in calculate_properties(smiles).items()}}
    row.update(overrides)
    return row, prov


def _export(tmp_path, rows, prov):
    result = export_library(rows, tmp_path, {"provenance": [prov]})
    return result


def test_valid_library_passes(tmp_path):
    row, prov = _library(tmp_path)
    result = _export(tmp_path, [row], prov)
    assert result["validation"]["valid"], result["validation"]["errors"]
    assert (tmp_path / LIBRARY_FILENAME).is_file()
    assert validate_export(tmp_path)["row_count"] == 1


def test_duplicate_pose_fails(tmp_path):
    row, prov = _library(tmp_path)
    assert not _export(tmp_path, [row, dict(row)], prov)["validation"]["valid"]


@pytest.mark.parametrize("field", ["canonical_smiles", "route_status", "sdf_path", "MW"])
def test_missing_required_field_fails(tmp_path, field):
    row, prov = _library(tmp_path)
    row[field] = ""
    assert not _export(tmp_path, [row], prov)["validation"]["valid"]


def test_missing_sdf_fails(tmp_path):
    row, prov = _library(tmp_path)
    (tmp_path / row["sdf_path"]).unlink()
    assert not _export(tmp_path, [row], prov)["validation"]["valid"]


def test_placeholder_source_and_unrouted_rejected(tmp_path):
    row, prov = _library(tmp_path, source="demo")
    assert not _export(tmp_path, [row], prov)["validation"]["valid"]
    row, prov = _library(tmp_path, route_status="unresolved")
    assert not _export(tmp_path, [row], prov)["validation"]["valid"]


def test_missing_provenance_fails(tmp_path):
    row, prov = _library(tmp_path)
    prov["checkpoint_sha256"] = ""
    assert not _export(tmp_path, [row], prov)["validation"]["valid"]


def test_required_field_set_matches_plan():
    for field in ("canonical_smiles", "compound_id", "provenance_id", "sdf_path", "route_status"):
        assert field in REQUIRED_FIELDS
