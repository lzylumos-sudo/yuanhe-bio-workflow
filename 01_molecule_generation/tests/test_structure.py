from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "v2"))

from src.hypotheses import derive_hypotheses, generate_hypotheses
from src.structure import (
    Atom,
    crop_complete_residues,
    entity_for_chain,
    parse_pdb,
    prepare_structure,
    terminal_carboxylate_status,
)


def atom_line(serial, name, resname, chain, resseq, x, y, z, *, icode="", altloc="", occupancy=1.0, element=None):
    element = element or next(c for c in name if c.isalpha())
    return (
        f"ATOM  {serial:5d} {name:>4s}{altloc:1s}{resname:>3s} {chain:1s}{resseq:4d}{icode:1s}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{20.0:6.2f}          {element:>2s}  \n"
    )


def pdb_header():
    return (
        "COMPND    MOL_ID: 1;\n"
        "COMPND   2 MOLECULE: LARGE PEPTIDE ENTITY;\n"
        "COMPND   3 CHAIN: B;\n"
        "COMPND   4 MOL_ID: 2;\n"
        "COMPND   5 MOLECULE: SMALL RECEPTOR ENTITY;\n"
        "COMPND   6 CHAIN: A;\n"
    )


def test_parse_preserves_model_chain_residue_and_insertion_code_and_selects_altloc(tmp_path):
    pdb = tmp_path / "models.pdb"
    text = pdb_header()
    text += "MODEL        1\n"
    text += atom_line(1, "CA", "GLY", "A", 10, 1, 0, 0, icode="A", altloc="A", occupancy=0.30)
    text += atom_line(2, "CA", "GLY", "A", 10, 2, 0, 0, icode="A", altloc="B", occupancy=0.70)
    text += atom_line(3, "N", "GLY", "A", 10, 0, 0, 0, icode="A")
    text += "ENDMDL\nMODEL        2\n"
    text += atom_line(4, "CA", "GLY", "A", 10, 9, 0, 0, icode="A")
    text += "ENDMDL\nEND\n"
    pdb.write_text(text)

    first = parse_pdb(pdb, model=1)
    assert first.available_models == (1, 2)
    ca = next(a for a in first.atoms if a.name == "CA")
    assert (ca.model, ca.chain, ca.resseq, ca.icode, ca.altloc, ca.occupancy, ca.xyz[0]) == (1, "A", 10, "A", "B", 0.70, 2.0)
    assert first.discarded_altloc_atoms == 1
    assert parse_pdb(pdb, model=2).atoms[0].xyz[0] == 9.0
    with pytest.raises(ValueError, match="model 3 not present"):
        parse_pdb(pdb, model=3)


def test_entity_assignment_uses_compnd_not_chain_size(tmp_path):
    pdb = tmp_path / "entities.pdb"
    # Chain B is deliberately larger, but COMPND identity remains authoritative.
    text = pdb_header() + atom_line(1, "CA", "ALA", "A", 1, 0, 0, 0)
    for serial in range(2, 10):
        text += atom_line(serial, "CA", "ALA", "B", serial, serial, 0, 0)
    pdb.write_text(text + "END\n")
    structure = parse_pdb(pdb)
    assert entity_for_chain(structure.entities, "A").molecule == "SMALL RECEPTOR ENTITY"
    assert entity_for_chain(structure.entities, "B").molecule == "LARGE PEPTIDE ENTITY"


def test_crop_keeps_complete_receptor_residues():
    receptor = [
        Atom("ATOM", 1, "CA", "", "ALA", "A", 5, "", (0, 0, 0), 1, 20, "C", ""),
        Atom("ATOM", 2, "CB", "", "ALA", "A", 5, "", (20, 0, 0), 1, 20, "C", ""),
        Atom("ATOM", 3, "CA", "", "ALA", "A", 6, "", (30, 0, 0), 1, 20, "C", ""),
    ]
    peptide = [Atom("ATOM", 4, "CA", "", "LEU", "B", 14, "", (1, 0, 0), 1, 20, "C", "")]
    cropped = crop_complete_residues(receptor, peptide, distance_A=6.0)
    assert [(a.resseq, a.name) for a in cropped] == [(5, "CA"), (5, "CB")]


def test_terminal_carboxylate_requires_oxt_and_explicit_status():
    oxygen = Atom("ATOM", 1, "O", "", "ASP", "B", 18, "", (0, 0, 0), 1, 20, "O", "")
    oxt = Atom("ATOM", 2, "OXT", "", "ASP", "B", 18, "", (1, 0, 0), 1, 20, "O", "")
    assert not terminal_carboxylate_status([oxygen], True)["terminal_carboxylate"]
    assert not terminal_carboxylate_status([oxygen, oxt], False)["terminal_carboxylate"]
    assert terminal_carboxylate_status([oxygen, oxt], True)["terminal_carboxylate"]


def test_prepare_actual_9nu4_and_generate_hypotheses(tmp_path):
    source = ROOT / "data/structures/9nu4.pdb"
    assert source.exists(), "actual 9NU4 fixture is required"
    prepared = prepare_structure(
        source,
        tmp_path,
        selected_peptide_residues=[14, 18],
        crop_distance_A=6.0,
        peptide_explicit_c_terminus=False,
    )
    expected = {
        "receptor_full.pdb", "receptor_prepared.pdb", "peptide_reference.pdb", "structures_manifest.csv",
        "residue_map.csv", "preparation_report.json", "checksums.sha256",
    }
    assert expected == {path.name for key, path in prepared.items() if key != "report"}
    assert not (tmp_path / "peptide_reference.sdf").exists()

    receptor = parse_pdb(tmp_path / "receptor_prepared.pdb")
    peptide = parse_pdb(tmp_path / "peptide_reference.pdb")
    assert receptor.atoms and peptide.atoms
    assert {a.chain for a in receptor.atoms} == {"A"}
    assert {a.chain for a in peptide.atoms} == {"B"}
    assert all(5 <= a.resseq <= 500 for a in receptor.atoms)
    assert prepared["report"]["entity_assignment_basis"] == "PDB_COMPND_CHAIN"
    assert prepared["report"]["receptor"]["entity"]["molecule"] == "LIPID II FLIPPASE MURJ"
    assert prepared["report"]["peptide"]["entity"]["molecule"] == "LYSIS PROTEIN"
    assert not prepared["report"]["peptide"]["c_terminus"]["terminal_carboxylate"]

    with (tmp_path / "residue_map.csv").open(newline="") as handle:
        mapping = list(csv.DictReader(handle))
    assert mapping and {row["role"] for row in mapping} == {"receptor", "peptide"}
    assert {"auth_chain", "auth_resseq", "auth_icode", "canonical_index"} <= set(mapping[0])
    manifest = next(csv.DictReader((tmp_path / "structures_manifest.csv").open(newline="")))
    assert manifest["receptor_chain"] == "A" and manifest["peptide_chain"] == "B"
    assert manifest["receptor_range"] == "5-500"
    checksum_lines = (tmp_path / "checksums.sha256").read_text().strip().splitlines()
    assert len(checksum_lines) == 6 and all(len(line.split()[0]) == 64 for line in checksum_lines)

    output = tmp_path / "hypotheses.json"
    result = generate_hypotheses(source, output)
    loaded = json.loads(output.read_text())
    assert loaded == result
    hypotheses = {h["hypothesis_id"]: h for h in result["hypotheses"]}
    assert {"H1", "H2"} <= set(hypotheses)
    assert [f["type"] for f in hypotheses["H1"]["features"]] == ["hydrophobic"]
    assert len(hypotheses["H2"]["features"]) >= 2
    for hypothesis in hypotheses.values():
        assert hypothesis["source_structure_ids"] == ["9NU4"]
        assert hypothesis["frame_id"] and hypothesis["threshold_version"]
        assert hypothesis["residue_mapping_version"]
        assert hypothesis["exclusion_spheres"]
        assert hypothesis["pass_rule"]
        for feature in hypothesis["features"]:
            assert len(feature["center_xyz"]) == 3
            assert feature["evidence_locator"]["source_atoms"]
            assert feature["radius_A"] > 0
            assert "direction" in feature and "confidence" in feature and "required" in feature
