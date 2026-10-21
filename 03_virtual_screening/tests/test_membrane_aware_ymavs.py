"""Focused Y-MAVS checks: geometry, pocket, accessibility, evidence ordering."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from zvss.membrane.accessibility import MembraneAccessibilityAnalyzer
from zvss.membrane.compatibility import MembraneCompatibilityAnalyzer
from zvss.membrane.ensemble import MembraneEnsembleAnalyzer
from zvss.membrane.pocket import MembranePocketAnalyzer
from zvss.membrane.regions import MembraneGeometry, MembraneRegionAnalyzer
from zvss.scoring.membrane_aware import MembraneAwareScorer, write_ranking_files
from zvss.select import qualify
from zvss.util import PdbAtom


def atom(serial, resseq, resname, coord, element="C", chain="A", name="CA"):
    return PdbAtom(
        "ATOM  ", serial, name, "", resname, chain, resseq, "",
        float(coord[0]), float(coord[1]), float(coord[2]), element=element,
    )


def geometry() -> MembraneGeometry:
    return MembraneGeometry.from_plane(
        {
            "center": [0.0, 0.0, 0.0],
            "normal": [0.0, 0.0, 1.0],
            "thickness_A": 30.0,
            "source": "synthetic PPM3",
            "method": "explicit test plane",
            "coordinate_frame": "docking_receptor",
            "receptor_sha256": "abc",
        },
        orientation_status="verified_computational",
    )


# --------------------------------------------------------------------------- #
# A. membrane geometry
# --------------------------------------------------------------------------- #
def test_signed_depth_regions_are_geometric_and_do_not_claim_inside_outside(tmp_path):
    geo = geometry()
    analyzer = MembraneRegionAnalyzer(geo, interface_buffer_A=5.0)
    points = np.array([[0, 0, 3.4], [0, 0, 16.0], [0, 0, 25.0], [0, 0, -16.0]])
    result = analyzer.analyze_points(points)
    assert result["region_counts"]["membrane_core"] == 1
    assert result["region_counts"]["membrane_interface"] == 2
    assert result["region_counts"]["outside_membrane_core"] == 1
    assert result["orientation_status"] == "verified_computational"
    descriptor = geo.point_descriptor([0, 0, 3.4], 5.0)
    assert descriptor["signed_depth_A"] == pytest.approx(3.4)
    assert descriptor["region"] == "membrane_core"
    assert "positive" in descriptor["side"]
    path = analyzer.write_geometry(str(tmp_path), atom_rows=[
        {"region": "membrane_core"}, {"region": "membrane_interface"}
    ])
    payload = json.loads(Path(path).read_text())
    assert payload["geometry"]["half_core_thickness_A"] == pytest.approx(15.0)
    assert payload["topology_assignment"] == "not_attempted"
    assert payload["status"] == "evaluated"


def test_residue_and_ligand_outputs_keep_depth_and_region(tmp_path):
    receptor = [
        atom(1, 10, "ALA", (0, 0, 0)),
        atom(2, 10, "ALA", (0, 0, 2)),
        atom(3, 11, "SER", (0, 0, 18)),
    ]
    analyzer = MembraneRegionAnalyzer(geometry(), interface_buffer_A=5.0)
    residues = analyzer.analyze_residues(receptor)
    assert len(residues) == 2
    assert residues[0]["region"] == "membrane_core"
    ligand = analyzer.analyze_ligand(
        np.array([[0, 0, 1.0], [0, 0, 18.0]]),
        pose_id="POS1", compound_id="CMP1", elements=["C", "O"],
    )
    assert ligand["ligand_core_fraction"] == pytest.approx(0.5)
    assert ligand["ligand_interface_fraction"] == pytest.approx(0.5)
    path = analyzer.write_ligand_csv(str(tmp_path), ligand["atom_rows"])
    with open(path, encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2
    assert rows[0]["pose_id"] == "POS1"


# --------------------------------------------------------------------------- #
# B. pocket and accessibility
# --------------------------------------------------------------------------- #
def test_pocket_uses_contacts_and_does_not_fake_volume():
    receptor = [
        atom(1, 10, "ALA", (3.0, 0, 3.0)),
        atom(2, 10, "ALA", (3.5, 0, 3.5)),
        atom(3, 12, "SER", (-3.0, 0, 3.0), element="O"),
        atom(4, 30, "LYS", (0, 0, 25.0), element="N", name="NZ"),
    ]
    analyzer = MembranePocketAnalyzer(
        receptor,
        geometry(),
        MembraneRegionAnalyzer(geometry(), interface_buffer_A=5.0),
        contact_cutoff_A=5.0,
        annotated_segments={"TM2": [("A", 10, "")], "TM7": [("A", 30, "")]},
    )
    pocket = analyzer.analyze(np.array([[0, 0, 3.0], [0.5, 0, 3.2]]), pose_id="POS1", compound_id="CMP1")
    assert pocket["pocket_contact_count"] == 2
    assert pocket["pocket_tm_segments"] == ["TM2"]
    assert pocket["pocket_volume_A3"] is None
    assert pocket["pocket_volume_status"] == "not_evaluated"
    assert pocket["pocket_region"] in {"membrane_embedded", "membrane_interface", "mixed"}


def test_accessibility_returns_real_geometric_path_and_uncertainty_boundary():
    receptor = [atom(1, 1, "ALA", (3, 0, 0)), atom(2, 2, "ALA", (-3, 0, 0))]
    analyzer = MembraneAccessibilityAnalyzer(
        geometry(), grid_spacing_A=1.0, probe_radius_A=1.4, margin_A=4.0, max_grid_points=100_000
    )
    prepared = analyzer.prepare(receptor)
    assert prepared["prepared"] is True
    result = analyzer.analyze(np.array([[0, 0, 2.5], [0, 0, 3.0]]), pose_id="POS1", compound_id="CMP1")
    assert result["evidence_level"] == "geometric_proxy"
    assert result["aqueous_accessibility"] == "accessible"
    assert result["membrane_lateral_accessibility"] == "accessible"
    assert result["paths"]["aqueous"]["length_A"] > 0
    assert result["probe_radius_A"] == 1.4
    assert "not a free-energy" in " ".join(result["limitations"]).lower() or "not a free-energy" in " ".join(result["limitations"]).lower()


def test_compatibility_reports_descriptors_without_scoring_them():
    receptor = [atom(1, 10, "ALA", (3, 0, 0)), atom(2, 12, "SER", (-3, 0, 0), element="O")]
    analyzer = MembraneCompatibilityAnalyzer(
        geometry(), MembraneRegionAnalyzer(geometry(), interface_buffer_A=5.0), receptor
    )
    result = analyzer.analyze(
        np.array([[0, 0, 2.5], [0, 0, 3.0]]),
        pose_id="POS1", compound_id="CMP1", ligand_elements=["C", "N"],
    )
    assert result["steric_compatibility_status"] in {"compatible", "minor_overlap", "clash"}
    assert result["lipid_clash_status"] == "not_evaluated"
    assert "not_binding_affinity" in " ".join(result["warnings"])


# --------------------------------------------------------------------------- #
# C. evidence ordering and output identity
# --------------------------------------------------------------------------- #
class FakeDB:
    def __init__(self, metrics, pose_rows=None):
        self.metrics = metrics
        self.pose_rows_data = pose_rows or {}

    def metric_map(self, scope_type, scope_id, stage):
        return self.metrics.get((scope_type, scope_id, stage), {})

    def query(self, sql, params=()):
        return []

    def one(self, sql, params):
        return self.pose_rows_data.get(params[0])


class FakeStages:
    def __init__(self, payload):
        self.payload = payload

    def get(self, *keys, default=None):
        node = self.payload
        for key in keys:
            if not isinstance(node, dict) or key not in node:
                return default
            node = node[key]
        return node


def fake_ctx(tmp_path, metrics, *, stage_status="sourced", verified=True, pose_rows=None):
    membrane_dir = tmp_path / "membrane"
    membrane_dir.mkdir()
    (membrane_dir / "summary.json").write_text(json.dumps({
        "status": stage_status,
        "orientation_verified": verified,
        "geometry": {
            "status": "evaluated" if verified else "proxy_only",
            "source": "synthetic-test",
            "coordinate_frame": "docking_receptor",
        },
        "input_hashes": {"receptor_sha256": "abc"},
    }))
    (membrane_dir / "ensemble_analysis.json").write_text(json.dumps({"enabled": False, "ensemble_status": "disabled"}))
    target = SimpleNamespace(
        target_id="synthetic",
        screening_mode="membrane_aware",
        membrane_screening={
            "enabled": True,
            "accessibility": {"enabled": True},
            "scoring": {"strategy": "evidence_first"},
        },
        primary=SimpleNamespace(pdb_id="TEST", chain="A"),
    )
    stages = FakeStages({
        "membrane_screening": {
            "enabled": True,
            "accessibility": {"enabled": True},
            "scoring": {"strategy": "evidence_first"},
        }
    })
    ctx = SimpleNamespace(
        run_id="unit",
        stage_dir=lambda *parts: str(tmp_path.joinpath(*parts)),
        reports_dir=lambda: str((tmp_path / "reports")),
        config=SimpleNamespace(target=target, stages=stages, digest=lambda: "digest"),
        db=FakeDB(metrics, pose_rows),
        log=lambda *args, **kwargs: None,
    )
    (tmp_path / "reports").mkdir()
    return ctx


def metric_rows(**values):
    return {key: {"value": value, "value_text": None} for key, value in values.items()}


def metric_text(**values):
    return {key: {"value": None, "value_text": value} for key, value in values.items()}


def test_membrane_aware_ranking_preserves_baseline_and_resolves_ties():
    import tempfile

    with tempfile.TemporaryDirectory() as directory:
        tmp = Path(directory)
        metrics = {
            ("pose", "POS_GOOD", "membrane"): {
                **metric_text(
                    orientation_status="verified_computational",
                    membrane_geometry_status="evaluated",
                    pocket_region="membrane_embedded",
                    aqueous_accessibility="accessible",
                    membrane_lateral_accessibility="uncertain",
                    geometric_accessibility_status="accessible",
                    steric_compatibility_status="compatible",
                ),
                **metric_rows(
                    ligand_center_depth_A=3.0,
                    ligand_core_fraction=1.0,
                    pocket_contact_count=8,
                    geometric_path_length_A=18.0,
                ),
            }
        }
        ctx = fake_ctx(tmp, metrics)
        entries = [
            {"compound_id": "CMP_A", "pose_id": "POS_MISSING", "state_id": "S1", "tier": "A",
             "evidence_score": 0.5, "raw": {"affinity": -8.0, "cnn_pose": 0.6, "cnn_affinity": 7.0},
             "components_missing": []},
            {"compound_id": "CMP_B", "pose_id": "POS_GOOD", "state_id": "S2", "tier": "A",
             "evidence_score": 0.5, "raw": {"affinity": -8.0, "cnn_pose": 0.6, "cnn_affinity": 7.0},
             "components_missing": []},
        ]
        scorer = MembraneAwareScorer(ctx, settings={
            "enabled": True,
            "accessibility": {"enabled": True},
            "scoring": {"strategy": "evidence_first", "allow_pose_reselection": False},
        })
        result = scorer.run(entries)
        ordered = result["entries"]
        assert [entry["compound_id"] for entry in ordered] == ["CMP_B", "CMP_A"]
        by_id = {entry["compound_id"]: entry for entry in ordered}
        assert by_id["CMP_B"]["baseline_rank"] == 2
        assert by_id["CMP_B"]["membrane_aware_rank"] == 1
        assert by_id["CMP_B"]["rank_shift"] == 1
        assert by_id["CMP_A"]["membrane_evidence_status"] == "incomplete"
        assert by_id["CMP_B"]["raw"]["affinity"] == -8.0
        assert by_id["CMP_B"]["membrane_evidence_status"] == "evaluated"

        selection = {
            "qualification": {
                "CMP_A": {"formal_qualified": False, "qualification": "preliminary"},
                "CMP_B": {"formal_qualified": False, "qualification": "preliminary"},
            }
        }
        outputs = write_ranking_files(ctx, {"ranking": ordered, "membrane_aware": {k: v for k, v in result.items() if k != "entries"}}, selection)
        with open(outputs["membrane_aware_ranking"], encoding="utf-8-sig") as handle:
            rows = list(csv.DictReader(handle))
        assert [row["pose_id"] for row in rows] == ["POS_GOOD", "POS_MISSING"]
        assert rows[0]["formal_qualified"] == "False"
        assert rows[0]["baseline_rank"] == "2"
        assert rows[0]["membrane_aware_rank"] == "1"
        assert Path(outputs["membrane_screening_report"]).is_file()


def test_missing_or_proxy_orientation_never_satisfies_membrane_geometry(tmp_path):
    membrane_dir = tmp_path / "membrane"
    membrane_dir.mkdir()
    (membrane_dir / "summary.json").write_text(json.dumps({
        "status": "proxy_only",
        "orientation_verified": False,
        "geometry": {"status": "proxy_only"},
    }))
    ctx = SimpleNamespace(
        stage_dir=lambda *parts: str(tmp_path.joinpath(*parts)),
        db=FakeDB({}),
        config=SimpleNamespace(
            get=lambda key, default=None: {
                "qualification.formal_requirements": ["pose_qa", "membrane_geometry"],
            }.get(key, default),
            target=SimpleNamespace(
                membrane_screening={"enabled": True, "qualification": {"formal_requirements": []}}
            ),
        ),
    )
    requirement = qualify._membrane_geometry(ctx, {"pose_id": "POS1"})
    assert requirement["status"] != "pass"
    assert requirement["satisfied"] is False
    rules = qualify.qualification_rules(ctx)
    assert "membrane_geometry" in rules


def test_ensemble_never_averages_missing_states():
    payload = {
        "enabled": True,
        "ensemble_status": "partial",
        "compounds": {
            "CMP1": {
                "number_of_evaluated_states": 1,
                "number_of_missing_states": 1,
                "state_coverage": 0.5,
                "best_supported_state": "S1",
                "ensemble_status": "partial",
                "state_specific_evidence": {
                    "S1": {"evidence_status": "evaluated", "record": {"gnina_affinity": -9.0, "pose_qc_status": "pass"}},
                    "S2": {"evidence_status": "missing", "record": None},
                },
            }
        },
    }
    # _integrate is the pure state-integrator; no filesystem or score synthesis.
    analyzer = object.__new__(MembraneEnsembleAnalyzer)
    states = [
        {"structure_id": "S1", "conformational_state": "outward", "evidence": {"status": "evaluated", "compounds": {"CMP1": {"gnina_affinity": -9.0}}}},
        {"structure_id": "S2", "conformational_state": "inward", "evidence": {"status": "missing", "compounds": {}}},
    ]
    integrated = analyzer._integrate(states)
    assert integrated["CMP1"]["number_of_evaluated_states"] == 1
    assert integrated["CMP1"]["number_of_missing_states"] == 1
    assert integrated["CMP1"]["state_coverage"] == 0.5
    assert integrated["CMP1"]["best_supported_state"] == "S1"
    assert integrated["CMP1"]["ensemble_status"] == "partial"
    assert integrated["CMP1"]["state_specific_evidence"]["S2"]["record"] is None


def test_ensemble_calpha_alignment_checks_residue_identity_and_rmsd(tmp_path):
    from zvss.util import write_pdb

    primary = tmp_path / "primary.pdb"
    state = tmp_path / "state.pdb"
    atoms = []
    coords = []
    for index in range(12):
        xyz = np.array([10.0 + index * 1.2, 20.0 + (index % 3) * 2.0, 30.0 + (index % 4) * 1.5])
        coords.append(xyz)
        atoms.append(atom(index + 1, 100 + index, "ALA", xyz))
    write_pdb(str(primary), atoms)
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    shift = np.array([4.0, -9.0, 2.5])
    moved = []
    for index, xyz in enumerate(coords):
        point = xyz @ rotation.T + shift
        moved.append(atom(index + 1, 100 + index, "ALA", point))
    write_pdb(str(state), moved)
    analyzer = object.__new__(MembraneEnsembleAnalyzer)
    analyzer.primary_structure_path = str(primary)
    analyzer.rules = {"alignment_min_ca": 12, "alignment_max_rmsd_A": 2.5}
    analyzer.ctx = SimpleNamespace(
        ref_dir=lambda: str(tmp_path),
        config=SimpleNamespace(target=SimpleNamespace(primary=SimpleNamespace(pdb_id="TEST", chain="A"))),
    )
    result = analyzer._map_to_primary(
        str(state), chain="A", primary_pdb="TEST", primary_chain="A", allow_mutations=False
    )
    assert result["status"] == "valid"
    assert result["matched_ca"] == 12
    assert result["rmsd_A"] < 0.01
    assert result["state_to_primary_rotation"]
