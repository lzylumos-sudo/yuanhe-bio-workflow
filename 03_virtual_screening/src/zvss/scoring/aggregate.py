"""Score integration: one evidence-tiered, fully traceable ranking.

Rules:
* a compound contributes at most one pose to the ranking (its representative);
* every component is normalised across the current candidate set, so the ranking
  is reproducible for a fixed input set;
* components that could not be computed are listed as missing and renormalise
  the weights instead of silently counting as zero.
"""

from __future__ import annotations

import math
from typing import Any

from ..constants import AssessmentStatus, Reason
from ..util import median_or_none, safe_float

DEFAULT_WEIGHTS = {
    "affinity": 0.34,
    "cnn_pose": 0.18,
    "cnn_affinity": 0.18,
    "mechanism": 0.18,
    "pose_quality": 0.12,
}


def _metrics(ctx, scope_type: str, scope_id: str, stage: str) -> dict[str, Any]:
    payload = ctx.db.metric_map(scope_type, scope_id, stage)
    return {key: row.get("value") for key, row in payload.items()}


def representative_pose(poses: list[dict[str, Any]]) -> dict[str, Any]:
    """Deterministic representative: best affinity, then best CNN, then pose_id."""
    def key(pose):
        scores = pose.get("scores") or {}
        affinity = safe_float(scores.get("minimizedAffinity"))
        cnn = safe_float(scores.get("CNNscore"))
        return (
            affinity if affinity is not None else math.inf,
            -(cnn if cnn is not None else -math.inf),
            pose["pose_id"],
        )

    return sorted(poses, key=key)[0]


def run(ctx) -> dict[str, Any]:
    from .. import data

    poses = data.pose_rows(ctx)
    weights = dict(DEFAULT_WEIGHTS)
    weights.update(ctx.config.get("scoring.weights", {}) or {})
    cutoff_affinity = safe_float(ctx.config.get("scoring.affinity_floor")) or 0.0

    by_compound: dict[str, list[dict[str, Any]]] = {}
    for pose in poses:
        by_compound.setdefault(pose["compound_id"], []).append(pose)

    entries: list[dict[str, Any]] = []
    for compound_id in sorted(by_compound):
        group = by_compound[compound_id]
        pose = representative_pose(group)
        scores = pose.get("scores") or {}
        qa = _metrics(ctx, "pose", pose["pose_id"], "pose_qa")
        mech = _metrics(ctx, "pose", pose["pose_id"], "mechanism")
        motion = _metrics(ctx, "pose", pose["pose_id"], "motion")
        props = _metrics(ctx, "compound", compound_id, "properties")
        strain = _metrics(ctx, "compound", compound_id, "strain")

        tier, tier_reasons = _tier(pose, qa, mech, motion)
        entry = {
            "compound_id": compound_id,
            "state_id": pose["state_id"],
            "pose_id": pose["pose_id"],
            "poses_available": len(group),
            "tier": tier,
            "tier_reasons": tier_reasons,
            "raw": {
                "affinity": safe_float(scores.get("minimizedAffinity")),
                "cnn_pose": safe_float(scores.get("CNNscore")),
                "cnn_affinity": safe_float(scores.get("CNNaffinity")),
                "mechanism": safe_float(mech.get("region_contacts_total")),
                "pose_quality": 1.0 - min(1.0, (safe_float(qa.get("min_receptor_distance")) or 0.0) / 4.5),
                "strain": safe_float(strain.get("strain_kcal_mol")),
            },
            "flags": {
                "excessive_strain": bool(strain.get("excessive_strain")),
                "property_alerts": props.get("alert_burden"),
                "single_pose": motion.get("pose_cluster_size") == 1,
            },
        }
        entries.append(entry)

    _normalise(entries, cutoff_affinity)
    available_weights = [w for w in weights if any(e["components"].get(w) is not None for e in entries)]
    total_weight = sum(weights[w] for w in available_weights) or 1.0
    for entry in entries:
        score = 0.0
        used = []
        missing = []
        for component, weight in weights.items():
            value = entry["components"].get(component)
            if value is None:
                missing.append(component)
                continue
            score += (weight / total_weight) * value
            used.append(component)
        entry["evidence_score"] = round(score, 6)
        entry["components_used"] = used
        entry["components_missing"] = missing
        ctx.db.record_assessment("compound", entry["compound_id"], "scoring", "evidence_score", entry["evidence_score"], AssessmentStatus.PASS)
        ctx.db.record_assessment("compound", entry["compound_id"], "scoring", "evidence_tier", entry["tier"], AssessmentStatus.PASS)
        for component, value in entry["components"].items():
            if value is not None:
                ctx.db.record_assessment("compound", entry["compound_id"], "scoring", f"component_{component}", round(value, 6), AssessmentStatus.PASS)

    entries.sort(key=lambda e: (e["tier"], -e["evidence_score"], e["compound_id"]))
    for index, entry in enumerate(entries, start=1):
        entry["rank"] = index
    return {
        "weights_used": {w: weights[w] for w in available_weights},
        "weights_skipped": [w for w in weights if w not in available_weights],
        "compounds": len(entries),
        "tiers": {tier: sum(1 for e in entries if e["tier"] == tier) for tier in ("A", "B", "C", "D")},
        "ranking": entries,
    }


def _tier(pose, qa, mech, motion) -> tuple[str, list[str]]:
    reasons: list[str] = []
    qa_status = str(pose.get("assess_status") or "")
    if qa_status == AssessmentStatus.FAIL.value:
        return "D", [str(Reason.POSE_GEOMETRY_INVALID)]
    contact_count = safe_float(mech.get("contact_count")) or 0.0
    region = safe_float(mech.get("region_contacts_total"))
    has_scores = any(v is not None for v in (pose.get("scores") or {}).values())
    if not has_scores:
        reasons.append(str(Reason.SCORE_MISSING))
    if qa_status == AssessmentStatus.FLAG.value:
        reasons.append("pose_qa_flagged")
    if contact_count <= 0:
        reasons.append(str(Reason.NO_CONTACTS))
        return "D", reasons
    if region is None:
        reasons.append(str(Reason.RESIDUE_MAPPING_MISSING))
    if qa_status == AssessmentStatus.PASS.value and has_scores and (region or 0) > 0:
        return "A", reasons
    if qa_status == AssessmentStatus.PASS.value and has_scores:
        return "B", reasons
    return "C", reasons


def _normalise(entries: list[dict[str, Any]], affinity_floor: float) -> None:
    keys = ("affinity", "cnn_pose", "cnn_affinity", "mechanism", "pose_quality")
    for entry in entries:
        entry["components"] = {}
    for key in keys:
        values = [
            -e["raw"][key] if key == "affinity" and e["raw"][key] is not None else e["raw"][key]
            for e in entries
        ]
        clean = [v for v in values if v is not None]
        if key == "affinity":
            clean = [v for v in clean if v >= -affinity_floor] if affinity_floor else clean
        if len(clean) < 2:
            for entry, value in zip(entries, values):
                if value is not None:
                    entry["components"][key] = 0.5
            continue
        ordered = sorted(clean)
        for entry, value in zip(entries, values):
            if value is None:
                entry["components"][key] = None
                continue
            # percentile rank; ties resolved by position in the sorted list
            rank = _percentile(value, ordered)
            entry["components"][key] = round(rank, 6)
    for entry in entries:
        strain = entry["raw"].get("strain")
        if strain is not None:
            penalty = max(0.0, min(1.0, (strain - 25.0) / 75.0))
            entry["components"]["strain_penalty"] = round(1.0 - penalty, 6)
        else:
            entry["components"]["strain_penalty"] = None


def _percentile(value: float, ordered: list[float]) -> float:
    n = len(ordered)
    below = sum(1 for v in ordered if v < value)
    equal = sum(1 for v in ordered if v == value)
    return (below + 0.5 * equal) / n
