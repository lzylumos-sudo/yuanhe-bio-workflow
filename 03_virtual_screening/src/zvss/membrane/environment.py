"""Membrane context: orientation proxy, insertion depth, and MD capability.

No membrane model is built here.  The module computes a transparent geometric
proxy -- the principal axis of the reference transmembrane bundle -- and reports
where a pose sits relative to that plane.  Membrane MD itself is reported as
``blocked_missing_membrane_capability`` until a membrane builder and a
membrane-capable MD engine are actually available.
"""

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
    from .. import data

    profile = ctx.config.target
    membrane = profile.membrane
    summary: dict[str, Any] = {
        "membrane_config": membrane.to_dict() if hasattr(membrane, "to_dict") else {},
        "usable": bool(getattr(membrane, "usable", False)),
        "status": "not_evaluated",
        "md_capability": {
            "status": "blocked_missing_membrane_capability",
            "reason_code": str(Reason.MD_MEMBRANE_CAPABILITY_MISSING),
            "note": (
                "no membrane builder and no membrane-capable MD runner are available in this environment; "
                "membrane MD is reported as blocked, not as passed"
            ),
        },
    }
    ref_atoms, record = data.reference_structure(ctx)
    tm_residues = []
    region = getattr(profile, "region", None)
    for side in list(getattr(region, "sides", []) or []):
        tm_residues.extend(getattr(side, "native_residues", []) or [])
    mapping = record.get("native_to_auth") or {}
    auth = {int(k) for k in mapping.values()} if mapping else set()
    selected = [a for a in ref_atoms if a.resseq in auth and a.resname != "HOH"] if auth else ref_atoms
    plane = _plane(selected)
    if plane is None:
        summary["reason"] = "reference structure missing or too few atoms"
        summary["reason_code"] = str(Reason.REFERENCE_POSE_MISSING)
        for pose in data.pose_rows(ctx):
            ctx.db.record_assessment("pose", pose["pose_id"], "membrane", "insertion_depth", None,
                                     AssessmentStatus.NEEDS_INPUT, Reason.REFERENCE_POSE_MISSING)
        return summary

    summary["plane"] = {
        "center": [round(float(v), 3) for v in plane["center"]],
        "normal": [round(float(v), 4) for v in plane["normal"]],
        "atoms": plane["atoms"],
        "source": "reference_structure_tm_bundle_pca",
        "note": "geometric proxy for the membrane normal; not a built membrane model",
    }
    summary["status"] = "proxy_only"

    for pose in data.pose_rows(ctx):
        mol = data.pose_mol(ctx, pose)
        if mol is None:
            ctx.db.record_assessment("pose", pose["pose_id"], "membrane", "insertion_depth", None,
                                     AssessmentStatus.NOT_EVALUATED, Reason.MISSING_COORDINATES)
            continue
        xyz = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
        depth = float(np.mean((xyz - plane["center"]) @ plane["normal"]))
        ctx.db.record_assessment("pose", pose["pose_id"], "membrane", "insertion_depth", round(depth, 4), AssessmentStatus.PASS)
        ctx.db.record_assessment("pose", pose["pose_id"], "membrane", "orientation_proxy", "tm_bundle_pca",
                                 AssessmentStatus.PASS)
    return summary
