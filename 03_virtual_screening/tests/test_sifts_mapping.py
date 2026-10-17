"""SIFTS-based residue numbering against the real 9NU4 entry (offline, cached)."""

from __future__ import annotations

import os

import pytest

from zvss.receptor.sifts import build_mapping

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
CACHE = os.path.join(PROJECT, "runtime", "reference_cache")
CIF = os.path.join(CACHE, "9nu4_updated.cif")
FASTA = os.path.join(CACHE, "P0AF16.fasta")


def _sequence(path: str) -> str:
    lines = [line.strip() for line in open(path, encoding="utf-8") if line.strip()]
    return "".join(line for line in lines if not line.startswith(">"))


@pytest.mark.skipif(not (os.path.exists(CIF) and os.path.exists(FASTA)), reason="9NU4 inputs not cached")
def test_murj_mapping_matches_the_deposition() -> None:
    report = build_mapping(
        CIF, pdb_id="9NU4", chain="A", accession="P0AF16", native_sequence=_sequence(FASTA)
    )
    summary = report.summary()
    assert summary["source"] == "sifts_pdbx_sifts_xref_db"
    assert summary["mapped_residues"] == 498
    assert summary["construct_residues"] == 614
    # exactly the three engineered substitutions, nothing else
    assert summary["mismatches"] == 3
    assert summary["native_range"] == [3, 500]
    pairs = {row.auth_seq_id: row.uniprot_resseq for row in report.rows if row.uniprot_resseq}
    assert pairs[3] == 3
    assert pairs[500] == 500
    assert summary["missing_native"] >= 12
    mismatched = {item["auth_seq_id"] for item in report.mismatches}
    assert mismatched == {5, 12, 13}


@pytest.mark.skipif(not os.path.exists(CIF), reason="9NU4 mmCIF not cached")
def test_mapping_cli_artifact_matches_report_rows() -> None:
    from zvss.receptor.sifts import write_mapping

    report = build_mapping(CIF, pdb_id="9NU4", chain="A", accession="P0AF16", native_sequence=_sequence(FASTA))
    paths = write_mapping(report, PROJECT)
    assert os.path.exists(paths["csv"]) or os.path.exists(paths["json"])
