"""Membrane orientation plus the Y-MAVS membrane-aware screening sub-analysis.

The orientation-transfer step is unchanged in spirit: a PPM3 placement is used
only after receptor identity, coverage and rigid-fit checks.  When a target opts
into ``screening_mode: membrane_aware``, this stage additionally computes:

* membrane regions for receptor atoms/residues and every docked ligand;
* a ligand-defined membrane pocket description;
* a static geometric accessibility grid;
* membrane-compatibility descriptors and a van der Waals overlap proxy;
* optional per-state ensemble evidence.

All outputs remain in the exact docking-receptor frame.  No sub-analysis is run
when that frame cannot be verified (except explicitly labelled exploration proxies
whose status is never promoted to a formal pass).
"""

from __future__ import annotations

import os
import shutil
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np

from .. import data
from ..constants import AssessmentStatus, Reason
from ..ids import sha256_file
from ..util import atomic_write_json, now_iso, parse_pdb_text, read_text
from .accessibility import MembraneAccessibilityAnalyzer
from .compatibility import MembraneCompatibilityAnalyzer
from .ensemble import MembraneEnsembleAnalyzer
from .pocket import MembranePocketAnalyzer
from .regions import (
    ALGORITHM_VERSION as REGION_ALGORITHM_VERSION,
    MembraneGeometry,
    MembraneRegionAnalyzer,
    REGION_UNKNOWN,
)
from .settings import (
    accessibility_rules,
    compatibility_rules,
    digest as settings_digest,
    geometry_rules,
    is_enabled,
    pocket_rules,
    settings_for,
)


def run(ctx) -> dict[str, Any]:
    profile = ctx.config.target
    membrane = profile.membrane
    settings = settings_for(ctx)
    enabled = is_enabled(ctx)
    summary: dict[str, Any] = {
        "schema": "ymavs-membrane-stage/1",
        "algorithm_version": REGION_ALGORITHM_VERSION,
        "generated_at": now_iso(),
        "screening_mode": profile.screening_mode,
        "membrane_screening_enabled": enabled,
        "membrane_screening_settings": settings,
        "membrane_config": asdict(membrane),
        "usable": False,
        "status": AssessmentStatus.NOT_EVALUATED.value,
        "orientation_verified": False,
        "md_capability": {
            "status": "separate_execution_required",
            "note": "Orientation and membrane geometry do not run MD; the standalone/workbench MD adapter validates its own environment and inputs.",
        },
        "limitations": [
            "membrane orientation is a computational placement unless explicitly declared experimental",
            "membrane depth and accessibility are geometric descriptors, not binding affinities or permeability predictions",
        ],
    }
    settings_hash = settings_digest(settings)
    summary["settings_digest"] = settings_hash
    evidence = membrane.raw.get("evidence_dir")
    plane: dict[str, Any] | None = None
    receptor_atoms, receptor_path = data.receptor_atoms(ctx)
    receptor_sha = sha256_file(receptor_path) if receptor_path and os.path.exists(receptor_path) else ""
    if evidence:
        folder = (Path(profile.source_path).parent / evidence).resolve()
        path = receptor_path
        try:
            from .orientation import transfer

            if not path:
                raise ValueError("docking receptor path is unavailable for orientation transfer")
            plane = transfer(folder, path)
            archive = Path(ctx.stage_dir("membrane")) / "orientation_evidence"
            shutil.copytree(folder, archive, dirs_exist_ok=True)
            summary.update(status="sourced", usable=True, orientation_verified=True)
            summary["orientation_status"] = "verified_computational"
        except (ValueError, KeyError, TypeError, OSError) as exc:
            summary.update(
                status="needs_input",
                reason=str(exc),
                reason_code=str(Reason.MEMBRANE_ORIENTATION_UNSOURCED),
                orientation_status="needs_input",
            )
    else:
        # PCA is only a geometric proxy.  When a docking receptor is available,
        # derive it from that receptor so the plane and the poses share one frame.
        proxy_atoms = receptor_atoms
        proxy_frame = "docking_receptor_pca_proxy"
        if not proxy_atoms:
            proxy_atoms, _record = data.reference_structure(ctx)
            proxy_frame = "reference_structure_pca_proxy"
        proxy = _plane(proxy_atoms)
        if proxy is not None:
            thickness = geometry_rules(settings)["proxy_thickness_A"]
            plane = {
                "center": proxy["center"].tolist(),
                "normal": proxy["normal"].tolist(),
                "thickness_A": thickness,
                "source": "receptor_structure_pca" if receptor_atoms else "reference_structure_pca",
                "atoms": proxy["atoms"],
                "coordinate_frame": proxy_frame,
                "note": (
                    "Geometric proxy only; not sourced membrane placement and not eligible for formal "
                    f"membrane qualification. A configured proxy_thickness_A={thickness:g} A is used solely "
                    "to partition space for exploration."
                ),
            }
            summary["status"] = "proxy_only"
            summary["orientation_status"] = "proxy_only"
    if plane is not None:
        summary["plane"] = plane
    if plane is not None:
        summary["input_hashes"] = {
            "receptor_sha256": receptor_sha,
            "orientation_evidence": plane.get("evidence_hashes") or {},
        }

    geometry = _build_geometry(summary, plane, receptor_sha)
    summary["geometry"] = geometry.to_dict(buffer_A=geometry_rules(settings)["interface_buffer_A"])
    summary["geometry"]["status"] = summary["geometry"].get("status")
    summary["geometry"]["input_hashes"] = summary.get("input_hashes") or {}
    summary["coordinate_frame"] = geometry.coordinate_frame
    summary["orientation_source"] = geometry.source

    region_analyzer = MembraneRegionAnalyzer(
        geometry, interface_buffer_A=geometry_rules(settings)["interface_buffer_A"]
    )
    residue_rows: list[dict[str, Any]] = []
    atom_rows: list[dict[str, Any]] = []
    if receptor_atoms:
        atom_rows = region_analyzer.analyze_atoms(receptor_atoms)
        residue_rows = region_analyzer.analyze_residues(receptor_atoms)
    output_dir = ctx.stage_dir("membrane")
    output_paths: dict[str, str] = {}
    output_paths["membrane_geometry"] = region_analyzer.write_geometry(
        output_dir,
        receptor_path=receptor_path,
        input_hashes=summary.get("input_hashes") or {},
        atom_rows=atom_rows,
    )
    output_paths["residue_membrane_regions"] = region_analyzer.write_residue_csv(output_dir, residue_rows)
    summary["receptor_region_counts"] = _region_counts(atom_rows)
    summary["residues_analyzed"] = len(residue_rows)

    annotated_segments = _annotated_segments(ctx)
    summary["annotated_segments"] = {
        name: {"configured_residue_count": len(residues)} for name, residues in annotated_segments.items()
    }
    pocket_analyzer = MembranePocketAnalyzer(
        receptor_atoms,
        geometry,
        region_analyzer,
        contact_cutoff_A=pocket_rules(settings)["contact_cutoff_A"],
        embedded_core_fraction=pocket_rules(settings)["embedded_core_fraction"],
        interface_fraction=pocket_rules(settings)["interface_fraction"],
        aqueous_outside_fraction=pocket_rules(settings)["aqueous_outside_fraction"],
        mixed_min_fraction=pocket_rules(settings)["mixed_min_fraction"],
        annotated_segments=annotated_segments,
    )
    compatibility_settings = compatibility_rules(settings)
    lipid_xyz, lipid_source, lipid_frame, lipid_hash = _load_explicit_lipids(
        ctx, compatibility_settings, geometry
    )
    compatibility_analyzer = MembraneCompatibilityAnalyzer(
        geometry,
        region_analyzer,
        receptor_atoms,
        contact_cutoff_A=compatibility_settings["contact_cutoff_A"],
        explicit_lipid_xyz=lipid_xyz,
        explicit_lipid_source=lipid_source,
        explicit_lipid_frame=lipid_frame,
        lipid_clash_cutoff_A=compatibility_settings["lipid_clash_cutoff_A"],
    )
    lipid_status = "not_evaluated"
    if lipid_xyz is not None:
        lipid_status = "evaluated" if lipid_frame and lipid_frame == geometry.coordinate_frame else "needs_input"
    summary["explicit_lipids"] = {
        "status": lipid_status,
        "source": lipid_source or None,
        "coordinate_frame": lipid_frame or None,
        "sha256": lipid_hash or None,
        "note": "absence of explicit lipid coordinates is not filled with a membrane plane",
    }
    access_analyzer = None
    accessibility_settings = accessibility_rules(settings)
    accessibility_enabled = enabled and accessibility_settings["enabled"]
    if accessibility_enabled and geometry.usable and receptor_atoms:
        access_analyzer = MembraneAccessibilityAnalyzer(
            geometry,
            grid_spacing_A=accessibility_settings["grid_spacing_A"],
            probe_radius_A=accessibility_settings["probe_radius_A"],
            margin_A=accessibility_settings["margin_A"],
            seed_radius_A=accessibility_settings["seed_radius_A"],
            max_grid_points=accessibility_settings["max_grid_points"],
            interface_buffer_A=geometry_rules(settings)["interface_buffer_A"],
        )
        summary["accessibility_grid"] = access_analyzer.prepare(receptor_atoms)
    else:
        summary["accessibility_grid"] = {
            "prepared": False,
            "status": "disabled" if not accessibility_enabled else "not_evaluated",
            "reason": "accessibility disabled by configuration" if not accessibility_enabled else "membrane geometry is not usable",
        }

    poses = data.pose_rows(ctx)
    per_pose: dict[str, dict[str, Any]] = {}
    ligand_csv_rows: list[dict[str, Any]] = []
    pocket_analyses: list[dict[str, Any]] = []
    accessibility_analyses: list[dict[str, Any]] = []
    compatibility_analyses: list[dict[str, Any]] = []
    for pose in poses:
        pose_result = _analyze_pose(
            ctx,
            pose,
            summary=summary,
            enabled=enabled,
            geometry=geometry,
            region_analyzer=region_analyzer,
            pocket_analyzer=pocket_analyzer,
            access_analyzer=access_analyzer,
            compatibility_analyzer=compatibility_analyzer,
            receptor_sha=receptor_sha,
            settings=settings,
        )
        per_pose[pose["pose_id"]] = pose_result["compact"]
        ligand_csv_rows.extend(pose_result.get("ligand_rows") or [])
        if pose_result.get("pocket"):
            pocket_analyses.append(pose_result["pocket"])
        if pose_result.get("accessibility"):
            accessibility_analyses.append(pose_result["accessibility"])
        if pose_result.get("compatibility"):
            compatibility_analyses.append(pose_result["compatibility"])
    summary["per_pose"] = per_pose
    if enabled:
        per_pose_path = os.path.join(output_dir, "pose_membrane_evidence.json")
        atomic_write_json(
            per_pose_path,
            {
                "schema": "ymavs-pose-membrane-evidence/1",
                "algorithm_version": REGION_ALGORITHM_VERSION,
                "generated_at": now_iso(),
                "note": "Every pose was evaluated separately; ranking and MD still use the single baseline representative pose unless an explicit, evidence-bound reselection implementation is enabled.",
                "coordinate_frame": geometry.coordinate_frame,
                "orientation_status": geometry.orientation_status,
                "poses": per_pose,
            },
        )
        output_paths["pose_membrane_evidence"] = per_pose_path
        output_paths["ligand_membrane_regions"] = region_analyzer.write_ligand_csv(output_dir, ligand_csv_rows)
        output_paths.update(pocket_analyzer.write_outputs(output_dir, pocket_analyses))
        if access_analyzer is not None:
            output_paths.update(access_analyzer.write_outputs(output_dir, accessibility_analyses))
        else:
            access_path = os.path.join(output_dir, "accessibility_analysis.json")
            paths_path = os.path.join(output_dir, "accessibility_paths.json")
            reason = (
                "accessibility disabled by configuration"
                if not accessibility_enabled
                else "membrane geometry or receptor atoms are unavailable"
            )
            atomic_write_json(
                access_path,
                {
                    "schema": "ymavs-accessibility-collection/1",
                    "algorithm_version": "ymavs-accessibility/1.0.0",
                    "generated_at": now_iso(),
                    "status": "not_evaluated",
                    "reason": reason,
                    "method": "grid_connectivity",
                    "evidence_level": "geometric_proxy",
                    "grid": summary.get("accessibility_grid", {}),
                    "count": 0,
                    "analyses": [],
                },
            )
            atomic_write_json(
                paths_path,
                {
                    "schema": "ymavs-accessibility-paths/1",
                    "algorithm_version": "ymavs-accessibility/1.0.0",
                    "status": "not_evaluated",
                    "reason": reason,
                    "paths": [],
                },
            )
            output_paths["accessibility_analysis"] = access_path
            output_paths["accessibility_paths"] = paths_path
        output_paths.update(compatibility_analyzer.write_outputs(output_dir, compatibility_analyses))
        ensemble = MembraneEnsembleAnalyzer(
            ctx,
            settings=settings,
            primary_geometry=geometry,
            primary_receptor_sha256=receptor_sha,
            primary_structure_path=receptor_path,
        )
        ensemble_payload = ensemble.run()
        output_paths["ensemble_analysis"] = ensemble.write_outputs(output_dir, ensemble_payload)
        output_paths["ensemble_state_coverage"] = os.path.join(output_dir, "ensemble_state_coverage.csv")
        summary["ensemble"] = {
            "status": ensemble_payload.get("ensemble_status"),
            "enabled": ensemble_payload.get("enabled"),
            "state_count": ensemble_payload.get("state_count"),
            "states_with_evidence": ensemble_payload.get("states_with_evidence"),
            "path": _relpath(output_paths["ensemble_analysis"], ctx.run_dir),
            "limitations": ensemble_payload.get("limitations") or [],
        }
        summary["pocket_analysis"] = {
            "status": "evaluated" if any(item.get("status") == "evaluated" for item in pocket_analyses) else "needs_input",
            "count": len(pocket_analyses),
            "volume_status": "not_evaluated",
        }
        summary["accessibility_analysis"] = {
            "status": "evaluated" if any(item.get("status") == "evaluated" for item in accessibility_analyses) else "not_evaluated",
            "count": len(accessibility_analyses),
            "grid": summary.get("accessibility_grid", {}),
            "evidence_level": "geometric_proxy",
        }
        summary["membrane_compatibility"] = {
            "status": "evaluated" if any(item.get("status") == "evaluated" for item in compatibility_analyses) else "needs_input",
            "count": len(compatibility_analyses),
            "explicit_lipid_status": summary["explicit_lipids"]["status"],
        }
        summary["membrane_evidence_status"] = _aggregate_evidence_status(summary, enabled)
    else:
        summary["membrane_evidence_status"] = "disabled"
    summary["counts"] = {
        "poses": len(poses),
        "poses_with_region_analysis": sum(1 for item in per_pose.values() if item.get("region")),
        "pocket_analyses": len(pocket_analyses),
        "accessibility_analyses": len(accessibility_analyses),
        "compatibility_analyses": len(compatibility_analyses),
    }
    summary["outputs"] = {name: _relpath(path, ctx.run_dir) for name, path in output_paths.items()}
    summary["algorithm_versions"] = {
        "membrane_regions": REGION_ALGORITHM_VERSION,
        "pocket": "ymavs-pocket/1.0.0",
        "accessibility": "ymavs-accessibility/1.0.0",
        "compatibility": "ymavs-compatibility/1.0.0",
        "ensemble": "ymavs-ensemble/1.0.0",
    }
    summary["config_hash"] = ctx.config.digest()
    summary["settings_digest"] = settings_hash
    summary["finished_at"] = now_iso()
    return summary


def _build_geometry(summary: dict[str, Any], plane: dict[str, Any] | None, receptor_sha: str) -> MembraneGeometry:
    if plane is None:
        geometry = MembraneGeometry.unavailable(reason="no sourced or proxy membrane orientation is available")
        summary["geometry_status"] = "needs_input"
        return geometry
    status = "verified_computational" if summary.get("orientation_verified") else "proxy_only"
    try:
        geometry = MembraneGeometry.from_plane(
            plane,
            orientation_status=status,
            uncertainty=str(summary.get("membrane_config", {}).get("uncertainty") or ""),
        )
    except (ValueError, TypeError) as exc:
        geometry = MembraneGeometry.unavailable(reason=f"invalid membrane plane: {exc}")
        summary["geometry_status"] = "needs_input"
        summary["usable"] = False
        return geometry
    expected = str(plane.get("receptor_sha256") or "")
    if status == "verified_computational":
        if expected and receptor_sha and expected != receptor_sha:
            geometry = MembraneGeometry.unavailable(
                reason="membrane orientation receptor hash does not match the docking receptor"
            )
            summary["orientation_verified"] = False
            summary["usable"] = False
            summary["status"] = "needs_input"
            summary["reason_code"] = str(Reason.MEMBRANE_ORIENTATION_UNSOURCED)
            summary["geometry_status"] = "needs_input"
            return geometry
        if not expected or not receptor_sha:
            geometry = MembraneGeometry.unavailable(
                reason="verified membrane orientation lacks a receptor hash binding"
            )
            summary["orientation_verified"] = False
            summary["usable"] = False
            summary["status"] = "needs_input"
            summary["geometry_status"] = "needs_input"
            return geometry
    summary["geometry_status"] = "evaluated" if geometry.verified else "proxy_only"
    return geometry


def _analyze_pose(
    ctx,
    pose: dict[str, Any],
    *,
    summary: dict[str, Any],
    enabled: bool,
    geometry: MembraneGeometry,
    region_analyzer: MembraneRegionAnalyzer,
    pocket_analyzer: MembranePocketAnalyzer,
    access_analyzer: MembraneAccessibilityAnalyzer | None,
    compatibility_analyzer: MembraneCompatibilityAnalyzer,
    receptor_sha: str,
    settings: dict[str, Any],
) -> dict[str, Any]:
    pose_id = pose["pose_id"]
    compound_id = pose.get("compound_id") or ""
    compact: dict[str, Any] = {"pose_id": pose_id, "compound_id": compound_id}
    mol = data.pose_mol(ctx, pose)
    if mol is None or mol.GetNumConformers() == 0:
        ctx.db.record_assessment("pose", pose_id, "membrane", "insertion_depth", None, AssessmentStatus.NEEDS_INPUT, Reason.MISSING_COORDINATES)
        ctx.db.record_assessment("pose", pose_id, "membrane", "ligand_center_depth_A", None, AssessmentStatus.NEEDS_INPUT, Reason.MISSING_COORDINATES)
        compact["status"] = "needs_input"
        compact["reason"] = "pose has no coordinates"
        return {"compact": compact}
    xyz = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    frame_ok = True
    pose_receptor_sha = str(pose.get("receptor_sha256") or "")
    if receptor_sha and pose_receptor_sha:
        frame_ok = pose_receptor_sha == receptor_sha
    verified = bool(summary.get("orientation_verified")) and frame_ok and geometry.verified
    proxy = geometry.orientation_status == "proxy_only" and frame_ok
    status = AssessmentStatus.PASS if verified else (AssessmentStatus.FLAG if proxy else AssessmentStatus.NEEDS_INPUT)
    reason = None if verified else (Reason.MEMBRANE_ORIENTATION_UNSOURCED if not proxy else None)
    if not frame_ok:
        status = AssessmentStatus.NEEDS_INPUT
        reason = Reason.RECEPTOR_UNBOUND
        compact["frame_status"] = "pose_receptor_sha_mismatch"
        for metric in (
            "insertion_depth", "ligand_center_depth_A", "ligand_core_fraction",
            "pocket_region", "aqueous_accessibility", "membrane_lateral_accessibility",
            "steric_compatibility_status", "membrane_evidence_status",
        ):
            ctx.db.record_assessment(
                "pose", pose_id, "membrane", metric, None, AssessmentStatus.NEEDS_INPUT, reason,
                value_text="not_evaluated" if metric not in {"insertion_depth", "ligand_center_depth_A", "ligand_core_fraction"} else None,
            )
        compact["status"] = "needs_input"
        compact["reason"] = "pose and membrane plane are not in the same verified receptor frame"
        return {"compact": compact}
    depth = None
    if geometry.usable:
        depth = float(np.mean(geometry.signed_depth(xyz)))
        ctx.db.record_assessment(
            "pose", pose_id, "membrane", "insertion_depth", round(depth, 4), status, reason,
            unit="angstrom", note="mean signed heavy-atom depth along the supplied membrane normal",
        )
        orientation_status = geometry.orientation_status
        ctx.db.record_assessment("pose", pose_id, "membrane", "orientation_status", None, status, reason, value_text=orientation_status)
        ctx.db.record_assessment("pose", pose_id, "membrane", "membrane_geometry_status", None, status, reason, value_text=summary.get("geometry_status"))
        ctx.db.record_assessment("pose", pose_id, "membrane", "orientation_source", None, status, reason, value_text=geometry.source or "unknown")
        ctx.db.record_assessment("pose", pose_id, "membrane", "coordinate_frame", None, status, reason, value_text=geometry.coordinate_frame or "unknown")
        compact.update({"region": "evaluated", "signed_depth_A": round(depth, 4)})
    if not enabled:
        compact["status"] = summary.get("status")
        return {"compact": compact}
    ligand = region_analyzer.analyze_ligand(xyz, pose_id=pose_id, compound_id=compound_id, elements=elements)
    ligand_rows = ligand.pop("atom_rows", [])
    ctx.db.record_assessment("pose", pose_id, "membrane", "ligand_center_depth_A", ligand.get("ligand_center_depth_A"), status, reason, unit="angstrom")
    ctx.db.record_assessment("pose", pose_id, "membrane", "ligand_core_fraction", ligand.get("ligand_core_fraction"), status, reason)
    ctx.db.record_assessment("pose", pose_id, "membrane", "ligand_interface_fraction", ligand.get("ligand_interface_fraction"), status, reason)
    ctx.db.record_assessment("pose", pose_id, "membrane", "ligand_outside_core_fraction", ligand.get("ligand_outside_core_fraction"), status, reason)
    ctx.db.record_assessment("pose", pose_id, "membrane", "ligand_region", None, status, reason, value_text=str(ligand.get("ligand_region") or REGION_UNKNOWN))
    pocket = pocket_analyzer.analyze(
        xyz,
        pose_id=pose_id,
        compound_id=compound_id,
        box_center=(pose.get("box") or {}).get("center"),
    )
    for key in ("pocket_region", "pocket_contact_count", "pocket_depth_A", "pocket_contact_centroid_depth_A", "pocket_core_fraction", "pocket_interface_fraction", "pocket_outside_core_fraction", "pocket_extent_A", "pocket_radius_of_gyration_A"):
        value = pocket.get(key)
        if isinstance(value, (int, float)):
            ctx.db.record_assessment("pose", pose_id, "membrane", key, value, status, reason)
        elif value is not None:
            ctx.db.record_assessment("pose", pose_id, "membrane", key, None, status, reason, value_text=str(value))
    ctx.db.record_assessment(
        "pose", pose_id, "membrane", "pocket_contact_residues", None, status, reason,
        value_text=";".join(pocket.get("pocket_contact_residues") or []),
    )
    ctx.db.record_assessment(
        "pose", pose_id, "membrane", "pocket_tm_segments", None, status, reason,
        value_text=";".join(pocket.get("pocket_tm_segments") or []),
    )
    ctx.db.record_assessment(
        "pose", pose_id, "membrane", "pocket_volume_status", None, status, reason,
        value_text=str(pocket.get("pocket_volume_status") or "not_evaluated"),
    )
    pocket_localization_status = (
        "evaluated" if pocket.get("pocket_region") not in (None, REGION_UNKNOWN) and (pocket.get("pocket_contact_count") or 0) > 0 else "not_evaluated"
    )
    ctx.db.record_assessment(
        "pose", pose_id, "membrane", "pocket_localization_status", None, status if pocket_localization_status == "evaluated" else AssessmentStatus.NOT_EVALUATED, reason,
        value_text=pocket_localization_status,
    )
    access = None
    if access_analyzer is not None:
        access = access_analyzer.analyze(xyz, pose_id=pose_id, compound_id=compound_id, pocket_center=pocket.get("pocket_center_A"))
        access_status = _access_assessment_status(access, status)
        ctx.db.record_assessment(
            "pose", pose_id, "membrane", "aqueous_accessibility", None, access_status, reason,
            value_text=str(access.get("aqueous_accessibility")),
        )
        ctx.db.record_assessment(
            "pose", pose_id, "membrane", "membrane_lateral_accessibility", None, access_status, reason,
            value_text=str(access.get("membrane_lateral_accessibility")),
        )
        ctx.db.record_assessment(
            "pose", pose_id, "membrane", "geometric_path_length_A", access.get("geometric_path_length_A"),
            access_status, reason, unit="angstrom",
        )
        ctx.db.record_assessment(
            "pose", pose_id, "membrane", "geometric_accessibility_status", None, access_status, reason,
            value_text=("accessible" if (access.get("aqueous_accessibility") == "accessible" or access.get("membrane_lateral_accessibility") == "accessible") else ("uncertain" if access.get("status") == "uncertain" else "not_evaluated")),
        )
        compact["accessibility"] = access.get("status")
    else:
        access = {
            "status": "not_evaluated",
            "aqueous_accessibility": "not_evaluated",
            "membrane_lateral_accessibility": "not_evaluated",
            "reason": "accessibility disabled or geometry unavailable",
        }
        compact["accessibility"] = "not_evaluated"
    existing = ctx.db.metric_map("pose", pose_id, "mechanism")
    mechanism_metrics = {
        key: (row.get("value") if row.get("value") is not None else row.get("value_text"))
        for key, row in existing.items()
    }
    compatibility = compatibility_analyzer.analyze(
        xyz,
        pose_id=pose_id,
        compound_id=compound_id,
        mol=mol,
        ligand_elements=elements,
        existing_metrics=mechanism_metrics,
    )
    ctx.db.record_assessment("pose", pose_id, "membrane", "steric_compatibility_status", None, status, reason, value_text=str(compatibility.get("steric_compatibility_status")))
    ctx.db.record_assessment("pose", pose_id, "membrane", "vdw_overlap_pairs", compatibility.get("vdw_overlap_pairs"), status, reason)
    ctx.db.record_assessment("pose", pose_id, "membrane", "vdw_severe_overlap_pairs", compatibility.get("vdw_severe_overlap_pairs"), status, reason)
    ctx.db.record_assessment("pose", pose_id, "membrane", "ligand_receptor_min_distance_A", compatibility.get("ligand_receptor_min_distance_A"), status, reason, unit="angstrom")
    ctx.db.record_assessment("pose", pose_id, "membrane", "lipid_clash_status", None, status, reason, value_text=str(compatibility.get("lipid_clash_status")))
    ctx.db.record_assessment("pose", pose_id, "membrane", "lipid_clash_count", compatibility.get("lipid_clash_count"), status, reason)
    ctx.db.record_assessment("pose", pose_id, "membrane", "lipid_minimum_distance_A", compatibility.get("lipid_minimum_distance_A"), status, reason, unit="angstrom")
    ctx.db.record_assessment("pose", pose_id, "membrane", "molecular_descriptors_status", None, status, reason, value_text=str((compatibility.get("molecular_descriptors") or {}).get("status")))
    evidence_status = _pose_evidence_status(summary, geometry, ligand, pocket, access, compatibility)
    ctx.db.record_assessment(
        "pose", pose_id, "membrane", "membrane_evidence_status", None,
        AssessmentStatus.PASS if evidence_status == "evaluated" else (AssessmentStatus.FLAG if evidence_status in {"proxy_only", "incomplete"} else AssessmentStatus.NOT_EVALUATED),
        reason, value_text=evidence_status,
    )
    compact.update(
        {
            "status": "evaluated" if evidence_status == "evaluated" else evidence_status,
            "region": ligand.get("ligand_region"),
            "signed_depth_A": ligand.get("ligand_center_depth_A"),
            "pocket_region": pocket.get("pocket_region"),
            "pocket_contact_count": pocket.get("pocket_contact_count"),
            "pocket_tm_segments": pocket.get("pocket_tm_segments") or [],
            "accessibility": access.get("status"),
            "steric_compatibility": compatibility.get("steric_compatibility_status"),
            "membrane_evidence_status": evidence_status,
        }
    )
    return {
        "compact": compact,
        "ligand_rows": ligand_rows,
        "pocket": pocket,
        "accessibility": access,
        "compatibility": compatibility,
    }


def _pose_evidence_status(
    summary: dict[str, Any],
    geometry: MembraneGeometry,
    ligand: dict[str, Any],
    pocket: dict[str, Any],
    access: dict[str, Any],
    compatibility: dict[str, Any],
) -> str:
    if not geometry.usable:
        return "incomplete"
    if ligand.get("ligand_center_depth_A") is None:
        return "incomplete"
    if pocket.get("pocket_region") in (None, REGION_UNKNOWN) or not pocket.get("pocket_contact_count"):
        return "incomplete"
    if compatibility.get("status") != "evaluated":
        return "incomplete"
    access_required = str((summary.get("accessibility_grid") or {}).get("status") or "") != "disabled"
    if access_required and access.get("status") not in {"evaluated", "uncertain"}:
        return "incomplete"
    if geometry.orientation_status == "proxy_only":
        return "proxy_only"
    if not summary.get("orientation_verified"):
        return "incomplete"
    return "evaluated"


def _aggregate_evidence_status(summary: dict[str, Any], enabled: bool) -> str:
    if not enabled:
        return "disabled"
    poses = summary.get("per_pose") or {}
    if not poses:
        return "incomplete"
    values = {str(item.get("membrane_evidence_status") or "incomplete") for item in poses.values()}
    if "incomplete" in values:
        return "incomplete"
    if summary.get("orientation_verified") and values == {"evaluated"}:
        return "evaluated"
    if summary.get("status") == "proxy_only" or "proxy_only" in values:
        return "proxy_only"
    return "incomplete"


def _access_assessment_status(access: dict[str, Any], base: AssessmentStatus) -> AssessmentStatus:
    if access.get("status") == "not_evaluated":
        return AssessmentStatus.NOT_EVALUATED
    if base == AssessmentStatus.PASS:
        return AssessmentStatus.PASS
    if base == AssessmentStatus.FLAG:
        return AssessmentStatus.FLAG
    return AssessmentStatus.NEEDS_INPUT


def _plane(atoms) -> dict[str, Any] | None:
    if len(atoms) < 12:
        return None
    xyz = np.array([[a.x, a.y, a.z] for a in atoms], dtype=float)
    center = xyz.mean(axis=0)
    _, _, vt = np.linalg.svd(xyz - center, full_matrices=False)
    normal = vt[2]
    norm = np.linalg.norm(normal)
    if not np.isfinite(norm) or norm == 0:
        return None
    return {
        "center": center,
        "normal": normal / norm,
        "atoms": len(atoms),
        "axis_explained_variance": float(_explained(xyz, vt)),
    }


def _explained(xyz: np.ndarray, vt: np.ndarray) -> float:
    centered = xyz - xyz.mean(axis=0)
    total = float((centered ** 2).sum())
    if total == 0:
        return 0.0
    first = float(((centered @ vt[0]) ** 2).sum())
    return first / total


def _annotated_segments(ctx) -> dict[str, list[tuple[str, int, str]]]:
    """Map target-configured annotation residues to author numbering.

    These are annotation sets (for MurJ, TM2/TM7 residues listed by the target
    profile), not automatically detected helices.  Missing mapping means the
    segment is omitted rather than guessed.
    """
    profile = ctx.config.target
    region = getattr(profile, "region", None)
    if region is None:
        return {}
    try:
        from ..mech.interactions import native_to_auth

        segments: dict[str, list[tuple[str, int, str]]] = {}
        for side in region.sides:
            resolved = native_to_auth(ctx, side.native_residues)
            keys = []
            chain = str(getattr(side, "chain", "") or profile.primary.chain or "")
            for native, info in resolved.items():
                auth = info.get("auth_seq_id")
                if auth is not None:
                    keys.append((chain, int(auth), ""))
            if keys:
                segments[side.name] = keys
        return segments
    except Exception:  # noqa: BLE001 - annotation absence must not fabricate segments
        return {}


def _load_explicit_lipids(
    ctx,
    rules: dict[str, Any],
    geometry: MembraneGeometry,
) -> tuple[np.ndarray | None, str, str, str]:
    path_value = rules.get("lipid_coordinates")
    if not path_value:
        return None, "", "", ""
    path = Path(str(path_value))
    if not path.is_absolute():
        path = Path(ctx.config.target.source_path).parent / path
    if not path.exists():
        return None, str(path), str(rules.get("lipid_coordinate_frame") or ""), ""
    frame = str(rules.get("lipid_coordinate_frame") or "")
    try:
        structure = parse_pdb_text(read_text(str(path)), keep_hydrogens=False)
        coordinates = np.asarray([[atom.x, atom.y, atom.z] for atom in structure.atoms if not atom.is_hydrogen()], dtype=float)
    except (OSError, ValueError):
        return None, str(path), frame, sha256_file(str(path))
    return coordinates.reshape((-1, 3)), str(path), frame, sha256_file(str(path))


def _region_counts(atom_rows: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in atom_rows:
        key = str(row.get("region") or REGION_UNKNOWN)
        counts[key] = counts.get(key, 0) + 1
    return counts


def _relpath(path: str, start: str) -> str:
    try:
        return os.path.relpath(path, start).replace(os.sep, "/")
    except ValueError:
        return path.replace(os.sep, "/")
