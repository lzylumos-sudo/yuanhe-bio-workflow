"""Sourced membrane placement with receptor-frame verification; PCA stays a proxy."""

from __future__ import annotations

from typing import Any

import numpy as np

from ..constants import AssessmentStatus, Reason


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
    return {"center": center, "normal": normal / norm, "atoms": len(atoms),
            "axis_explained_variance": float(_explained(xyz, vt))}


def _explained(xyz: np.ndarray, vt: np.ndarray) -> float:
    centered = xyz - xyz.mean(axis=0)
    total = float((centered ** 2).sum())
    if total == 0:
        return 0.0
    first = float(((centered @ vt[0]) ** 2).sum())
    return first / total


def run(ctx) -> dict[str, Any]:
    from pathlib import Path
    import shutil
    from dataclasses import asdict
    from .. import data
    from .orientation import transfer

    profile = ctx.config.target
    membrane = profile.membrane
    summary = {"membrane_config": asdict(membrane), "usable": False,
               "status": "not_evaluated", "md_capability": {
                   "status": "separate_execution_required",
                   "note": "Orientation assessment does not run MD; the standalone/workbench MD adapter validates its own environment and inputs."}}
    evidence = membrane.raw.get("evidence_dir")
    plane = None
    if evidence:
        folder = (Path(profile.source_path).parent / evidence).resolve()
        receptor = data.receptor_record(ctx)
        path = receptor.get("archived") or receptor.get("path") or receptor.get("output")
        try:
            plane = transfer(folder, path)
            archive = Path(ctx.stage_dir("membrane")) / "orientation_evidence"
            shutil.copytree(folder, archive, dirs_exist_ok=True)
            summary.update(status="sourced", usable=True, orientation_verified=True)
        except (ValueError, KeyError, TypeError, OSError) as exc:
            summary.update(status="needs_input", reason=str(exc),
                           reason_code=str(Reason.MEMBRANE_ORIENTATION_UNSOURCED))
    else:
        ref_atoms, record = data.reference_structure(ctx)
        plane = _plane(ref_atoms)
        if plane is not None:
            plane = {"center": plane["center"].tolist(), "normal": plane["normal"].tolist(),
                     "source": "reference_structure_pca", "atoms": plane["atoms"],
                     "note": "Geometric proxy only; not sourced membrane placement"}
            summary["status"] = "proxy_only"
    if plane is not None:
        summary["plane"] = plane
    for pose in data.pose_rows(ctx):
        mol = data.pose_mol(ctx, pose)
        valid = mol is not None and mol.GetNumConformers() > 0 and plane is not None
        if valid and summary["status"] == "sourced":
            valid = pose.get("receptor_sha256") == plane["receptor_sha256"]
        if not valid:
            ctx.db.record_assessment("pose", pose["pose_id"], "membrane", "insertion_depth", None,
                                     AssessmentStatus.NEEDS_INPUT, Reason.MEMBRANE_ORIENTATION_UNSOURCED)
            continue
        xyz = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
        depth = float(np.mean((xyz-np.asarray(plane["center"])) @ np.asarray(plane["normal"])))
        if not np.isfinite(depth):
            raise ValueError("Non-finite membrane depth")
        ctx.db.record_assessment("pose", pose["pose_id"], "membrane", "insertion_depth", round(depth,4), AssessmentStatus.PASS)
        ctx.db.record_assessment("pose", pose["pose_id"], "membrane", "orientation_source", plane["source"], AssessmentStatus.PASS)
    return summary
