"""Formal candidate qualification, separate from the display shortlist.

A ranking exists to order everything that was measured.  A *qualification* is a
different claim: it says the required evidence for a candidate was actually
produced.  Mixing the two is how a top-of-list compound becomes a "primary
candidate" without the checks that the specification requires.

Every requirement returns one of:

``pass``           measured, and the measurement satisfies the rule
``fail``           measured, and the measurement violates the rule
``not_evaluated``  the check never ran (tool missing, input missing, stage skipped)
``needs_input``    the check needs a configuration value that is absent
``unsupported``    the capability does not exist in this build

Only ``pass`` counts toward formal qualification.  In particular the absence of
membrane-orientation data or of re-dock support is never rounded up to a pass.
"""

from __future__ import annotations

import os
from typing import Any

from ..constants import AssessmentStatus, Reason
from ..util import read_json

PASS = AssessmentStatus.PASS.value
FAIL = AssessmentStatus.FAIL.value
NOT_EVALUATED = AssessmentStatus.NOT_EVALUATED.value
NEEDS_INPUT = AssessmentStatus.NEEDS_INPUT.value
UNSUPPORTED = AssessmentStatus.UNSUPPORTED.value

#: Requirement ids understood by the selector, in report order.
REQUIREMENT_ORDER = ("pose_qa", "target_region", "redock_seeds", "strain", "membrane")

DEFAULT_FORMAL = {
    "pose_qa": {},
    "target_region": {},
    "redock_seeds": {"seeds_required": 3, "supporting_seeds_required": 2},
    "strain": {},
    "membrane": {},
}


def _metric(ctx, scope_type: str, scope_id: str, stage: str) -> dict[str, Any]:
    rows = ctx.db.query(
        "SELECT metric, value, value_text, status, reason_code FROM assessment "
        "WHERE scope_type=? AND scope_id=? AND stage=? AND run_id=?",
        (scope_type, scope_id, stage, getattr(ctx.db, "run_id", "") or ""),
    )
    out: dict[str, Any] = {}
    for row in rows:
        if row.get("value") is not None:
            out[row["metric"]] = row["value"]
        elif row.get("value_text") is not None:
            out[row["metric"]] = row["value_text"]
    return out


def _requirement(
    name: str,
    *,
    status: str,
    observed: Any = None,
    detail: str = "",
    source: str = "",
    reason_code: str | None = None,
) -> dict[str, Any]:
    return {
        "requirement": name,
        "status": status,
        "satisfied": status == PASS,
        "observed": observed,
        "detail": detail,
        "source": source,
        "reason_code": reason_code,
    }


def _pose_qa(ctx, entry) -> dict[str, Any]:
    metrics = _metric(ctx, "pose", entry["pose_id"], "pose_qa")
    run_verdict = metrics.get("qa_status")
    fallback = ctx.db.one("SELECT assess_status, reason_code FROM pose WHERE pose_id=?", (entry["pose_id"],)) or {}
    status = str(run_verdict or fallback.get("assess_status") or NOT_EVALUATED)
    source = "pose_qa assessment for this run" if run_verdict else "pose.assess_status (global fallback)"
    mapping = {AssessmentStatus.PASS.value: PASS, AssessmentStatus.FLAG.value: PASS,
               AssessmentStatus.FAIL.value: FAIL}
    mapped = mapping.get(status, NOT_EVALUATED)
    detail = {
        PASS: "pose passed geometry and chemistry checks"
        if status == PASS else "pose was flagged; retained with its flag recorded",
        FAIL: "pose failed a blocking check",
        NOT_EVALUATED: "pose-level QA never produced a verdict",
    }[mapped]
    return _requirement(
        "pose_qa",
        status=mapped,
        observed=status,
        detail=detail,
        source="pose.assess_status",
        reason_code=pose.get("reason_code") if mapped != PASS else None,
    )


def _target_region(ctx, entry) -> dict[str, Any]:
    profile = ctx.config.target
    region = getattr(profile, "region", None)
    mech = _metric(ctx, "pose", entry["pose_id"], "mechanism")
    contacts = mech.get("region_contacts_total")
    if region is None or not getattr(region, "usable", False):
        return _requirement(
            "target_region",
            status=NOT_EVALUATED,
            observed=None,
            detail="target profile declares no usable local region; no site-specific evidence exists",
            source="configs/targets",
            reason_code=str(Reason.REGION_NOT_CONFIGURED),
        )
    if contacts is None:
        return _requirement(
            "target_region",
            status=NOT_EVALUATED,
            observed=None,
            detail="region is configured but the pose has no region-contact measurement",
            source="mechanism stage",
            reason_code=str(Reason.RESIDUE_MAPPING_MISSING),
        )
    if float(contacts) <= 0:
        return _requirement(
            "target_region",
            status=FAIL,
            observed=float(contacts),
            detail="pose makes no contact with the configured region",
            source="mechanism stage",
            reason_code=str(Reason.NO_CONTACTS),
        )
    return _requirement(
        "target_region",
        status=PASS,
        observed=float(contacts),
        detail="pose contacts residues inside the configured region",
        source="mechanism stage",
    )


def _redock_seeds(ctx, entry, rules: dict[str, Any]) -> dict[str, Any]:
    planned = int(rules.get("seeds_required", 3))
    supporting_required = int(rules.get("supporting_seeds_required", 2))
    stage = read_json(os.path.join(ctx.stage_dir("redock"), "summary.json"), {}) or {}
    per_compound = (stage.get("compounds") or {}).get(entry["compound_id"]) if isinstance(stage.get("compounds"), dict) else None
    metrics = _metric(ctx, "compound", entry["compound_id"], "redock")
    attempted = int(metrics.get("seeds_attempted") or 0)
    succeeded = int(metrics.get("seeds_succeeded") or 0)
    supporting = int(metrics.get("supporting_seeds_same_cluster") or 0)
    status = str(stage.get("status") or "not_run")
    observed = {
        "seeds_planned": planned,
        "seeds_attempted": attempted,
        "seeds_succeeded": succeeded,
        "supporting_seeds_same_cluster": supporting,
        "required_supporting": supporting_required,
        "stage_status": status,
    }
    if status.startswith("blocked"):
        return _requirement(
            "redock_seeds",
            status=NOT_EVALUATED,
            observed=observed,
            detail=f"re-docking did not run: {stage.get('reason_code') or stage.get('status')}",
            source="redock stage",
            reason_code=str(stage.get("reason_code") or Reason.REDOCK_TOOL_UNAVAILABLE),
        )
    if per_compound is None and succeeded == 0:
        return _requirement(
            "redock_seeds",
            status=NOT_EVALUATED,
            observed=observed,
            detail="this compound was never re-docked, so seed support is unknown",
            source="redock stage",
            reason_code=str(Reason.REDOCK_SEED_INCOMPLETE),
        )
    if succeeded < planned:
        return _requirement(
            "redock_seeds",
            status=FAIL,
            observed=observed,
            detail=f"only {succeeded} of {planned} planned seeds completed",
            source="redock stage",
            reason_code=str(Reason.REDOCK_SEED_INCOMPLETE),
        )
    if supporting < supporting_required:
        return _requirement(
            "redock_seeds",
            status=FAIL,
            observed=observed,
            detail=f"{supporting} seeds support one pose cluster, {supporting_required} required",
            source="redock stage",
            reason_code=str(Reason.REDOCK_SEED_INCOMPLETE),
        )
    return _requirement(
        "redock_seeds",
        status=PASS,
        observed=observed,
        detail=f"{supporting} of {planned} seeds reproduce the same pose cluster",
        source="redock stage",
    )


def _strain(ctx, entry) -> dict[str, Any]:
    metrics = _metric(ctx, "compound", entry["compound_id"], "strain")
    value = metrics.get("strain_kcal_mol")
    force_field = metrics.get("force_field")
    if value is None:
        reason = "strain was not measured; either no force field parameterises the molecule or the geometry is degenerate"
        return _requirement(
            "strain",
            status=NOT_EVALUATED,
            observed={"strain_kcal_mol": None, "force_field": force_field},
            detail=reason,
            source="strain stage",
            reason_code=str(Reason.STRAIN_NOT_CONVERGED),
        )
    flag_limit = float(getattr(ctx.config.budget, "strain_max_kcal", 25.0) or 25.0)
    if float(metrics.get("excessive_strain") or 0.0) >= 1.0:
        return _requirement(
            "strain",
            status=FAIL,
            observed={"strain_kcal_mol": float(value), "force_field": force_field, "limit": flag_limit},
            detail="strain is above the failure threshold; the pose is not a usable starting structure",
            source="strain stage",
            reason_code=str(Reason.EXCESSIVE_STRAIN),
        )
    return _requirement(
        "strain",
        status=PASS,
        observed={"strain_kcal_mol": float(value), "force_field": force_field, "limit": flag_limit},
        detail="strain is measurable and below the failure threshold",
        source="strain stage",
    )


def _membrane(ctx, entry) -> dict[str, Any]:
    stage = read_json(os.path.join(ctx.stage_dir("membrane"), "summary.json"), {}) or {}
    capability = stage.get("md_capability") or {}
    capability_status = str(capability.get("status") or "")
    stage_status = str(stage.get("status") or "")
    observed = {"membrane_status": stage_status or None, "md_capability": capability_status or None}
    if not stage:
        return _requirement(
            "membrane",
            status=NOT_EVALUATED,
            observed=observed,
            detail="membrane stage produced no summary in this run",
            source="membrane stage",
            reason_code=str(Reason.MEMBRANE_ORIENTATION_UNSOURCED),
        )
    sourced = stage_status == "sourced" and not capability_status.startswith("blocked")
    if sourced:
        return _requirement(
            "membrane",
            status=PASS,
            observed=observed,
            detail="membrane orientation comes from a configured source",
            source="membrane stage",
        )
    return _requirement(
        "membrane",
        status=NOT_EVALUATED,
        observed=observed,
        detail=(
            "no sourced membrane orientation exists; only a geometric membrane proxy was computed"
            if stage_status == "proxy_only"
            else f"membrane environment is not sourced (stage status: {stage_status or 'unknown'})"
        ),
        source="membrane stage",
        reason_code=str(Reason.MEMBRANE_ORIENTATION_UNSOURCED),
    )


def evaluate(ctx, entry: dict[str, Any], rules: dict[str, Any] | None = None) -> dict[str, Any]:
    rules = rules or DEFAULT_FORMAL
    results = [
        _pose_qa(ctx, entry),
        _target_region(ctx, entry),
        _redock_seeds(ctx, entry, rules.get("redock_seeds") or {}),
        _strain(ctx, entry),
        _membrane(ctx, entry),
    ]
    required = [r for r in results if r["requirement"] in rules]
    unmet = [r for r in required if not r["satisfied"]]
    formal = not unmet
    return {
        "compound_id": entry["compound_id"],
        "pose_id": entry["pose_id"],
        "tier": entry["tier"],
        "formal_qualified": formal,
        "qualification": "formal" if formal else "preliminary",
        "requirements": results,
        "unmet_requirements": [r["requirement"] for r in unmet],
        "blocked_by_missing_evidence": [r["requirement"] for r in unmet if r["status"] in {NOT_EVALUATED, NEEDS_INPUT, UNSUPPORTED}],
        "failed_requirements": [r["requirement"] for r in unmet if r["status"] == FAIL],
        "statement": (
            "all required evidence was produced and satisfied its rule"
            if formal
            else "display-only: required evidence is missing or violated, so this is not a formal candidate"
        ),
    }


def qualification_rules(ctx) -> dict[str, Any]:
    """Rules from the budget file, with the documented defaults as fallback."""
    raw = ctx.config.get("qualification.formal_requirements") or list(REQUIREMENT_ORDER)
    rules: dict[str, Any] = {}
    for name in raw:
        if name == "redock_seeds":
            rules[name] = {
                "seeds_required": int(ctx.config.get("qualification.redock_seeds_required", 3) or 3),
                "supporting_seeds_required": int(
                    ctx.config.get("qualification.supporting_seeds_required", 2) or 2
                ),
            }
        elif name in DEFAULT_FORMAL:
            rules[name] = dict(DEFAULT_FORMAL[name])
    return rules


def display_slots(ctx) -> int:
    return int(ctx.config.get("qualification.display_primary_slots", 2) or 2)
