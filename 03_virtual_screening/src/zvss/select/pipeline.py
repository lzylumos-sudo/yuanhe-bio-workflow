"""Candidate selection: evidence tiers, budget, Pareto front, diversity.

Determinism is enforced here: every sort has an explicit total order ending in
``compound_id``, diversity selection is a greedy max-min over a fixed candidate
list, and the resulting selection carries its own ``selection_version`` so a
report can be reproduced byte-for-byte.
"""

from __future__ import annotations

import math
from typing import Any

from ..chem import rdio
from ..constants import AssessmentStatus, SelectionStatus, Reason
from ..ids import selection_version, sha256_text
from ..util import safe_float
import json
from . import qualify


def run(ctx, ranking: dict[str, Any]) -> dict[str, Any]:
    budget = ctx.config.budget
    entries = ranking["ranking"]
    strategy = ctx.config.get("selection.strategy", "evidence_diversity")
    if strategy not in {"evidence_diversity", "pareto_diversity"}:
        raise ValueError("Unknown selection.strategy: " + str(strategy))
    limits = {
        "max_primary": int(getattr(budget, "md_primary_slots", 8)),
        "max_exploratory": int(getattr(budget, "md_exploratory_slots", 4)),
        "max_poses_per_parent": int(getattr(budget, "poses_per_parent_max", 2)),
        "similarity_ceiling": float(getattr(budget, "similarity_ceiling", 0.7)),
    }
    eligible: list[dict[str, Any]] = []
    exclusions: list[dict[str, Any]] = []
    for entry in entries:
        reason = _exclusion(entry, limits, ctx)
        if reason is None:
            eligible.append(entry)
        else:
            exclusions.append({"compound_id": entry["compound_id"], "pose_id": entry["pose_id"],
                               "tier": entry["tier"], "reason_code": str(reason[0]), "detail": reason[1]})
            ctx.db.set_pose_status(entry["pose_id"], select_status=SelectionStatus.EXCLUDED, reason_code=reason[0])
            ctx.db.record_assessment("compound", entry["compound_id"], "selection", "excluded", 1.0,
                                     status=None, reason_code=reason[0], note=reason[1])

    by_tier = sorted(eligible, key=lambda e: (e["tier"], -e["evidence_score"], e["compound_id"]))

    # Qualification is evaluated on every eligible compound, not only on the
    # shortlist, so the report can show why a well-scoring compound is not a
    # formal candidate.
    rules = qualify.qualification_rules(ctx)
    verdicts = {entry["compound_id"]: qualify.evaluate(ctx, entry, rules) for entry in by_tier}
    formal_pool = [entry for entry in by_tier if verdicts[entry["compound_id"]]["formal_qualified"]]

    primary_pool = formal_pool if strategy == "pareto_diversity" else formal_pool[: max(limits["max_primary"] * 4, 4)]
    primary = _select_policy(primary_pool, limits["max_primary"], diversity=0.6, strategy=strategy)
    formal_ids = {e["compound_id"] for e in primary}
    display_candidates = [
        entry for entry in by_tier if entry["compound_id"] not in formal_ids
    ]
    preliminary = _select(
        display_candidates[: max(qualify.display_slots(ctx) * 4, 4)],
        qualify.display_slots(ctx),
        diversity=0.6,
    )
    taken = formal_ids | {e["compound_id"] for e in preliminary}
    exploratory_pool = [
        e for e in by_tier if e["compound_id"] not in taken and e["tier"] in {"A", "B", "C"}
    ]
    exploratory = _select(exploratory_pool, limits["max_exploratory"], diversity=0.75)

    selection_id = selection_version(
        [e["compound_id"] for e in primary + preliminary + exploratory],
        sha256_text(json.dumps({"config": ctx.config.digest(), "limits": limits, "rules": rules,
                               "policy": "screening-evidence/v2", "strategy": strategy,
                               "evidence": [{k: v for k, v in e.items() if k not in {"selection_role", "selection_order", "selection_version", "qualification"}} for e in entries]}, sort_keys=True, allow_nan=False)),
    )
    for role, group in (("primary", primary), ("preliminary", preliminary), ("exploratory", exploratory)):
        for order, entry in enumerate(group, start=1):
            entry["selection_role"] = role
            entry["selection_order"] = order
            entry["selection_version"] = selection_id
            entry["qualification"] = verdicts[entry["compound_id"]]["qualification"]
            ctx.db.set_pose_status(entry["pose_id"], select_status=SelectionStatus.SELECTED)
            ctx.db.record_assessment("compound", entry["compound_id"], "selection", "selected", 1.0,
                                     status=None, note=role)

    for compound_id, verdict in verdicts.items():
        ctx.db.record_assessment(
            "compound", compound_id, "qualification", "formal_qualified",
            1.0 if verdict["formal_qualified"] else 0.0,
            AssessmentStatus.PASS if verdict["formal_qualified"] else AssessmentStatus.FLAG,
            None if verdict["formal_qualified"] else Reason.EVIDENCE_INCOMPLETE,
            note=verdict["statement"],
        )
        for requirement in verdict["requirements"]:
            ctx.db.record_assessment(
                "compound", compound_id, "qualification", f"requirement_{requirement['requirement']}",
                None,
                None,
                requirement["reason_code"],
                value_text=requirement["status"],
                note=requirement["detail"],
            )

    pareto = _pareto(eligible)
    return {
        "selection_version": selection_id,
        "policy_version": "screening-evidence/v2",
        "strategy": strategy,
        "pareto_objectives": {"affinity": "minimize kcal/mol", "mechanism": "maximize region contacts", "strain": "minimize kcal/mol"},
        "pareto_scope": "eligible candidates with complete finite objectives; not a qualification verdict",
        "pareto_unassessed": {e["compound_id"]: _missing_objectives(e) for e in eligible if _missing_objectives(e)},
        "formal_pareto_front": [e["compound_id"] for e in _pareto(formal_pool)],
        "limits": limits,
        "qualification_rules": rules,
        "display_primary_slots": qualify.display_slots(ctx),
        "candidates_considered": len(entries),
        "eligible": len(eligible),
        "excluded": len(exclusions),
        "formal_qualified": len(formal_pool),
        "primary": primary,
        "preliminary": preliminary,
        "exploratory": exploratory,
        "qualification": verdicts,
        "exclusions": exclusions,
        "pareto_front": [e["compound_id"] for e in pareto],
        "diversity": _diversity_report(primary + preliminary + exploratory),
        "notes": [
            "primary = formally qualified candidates only (every required evidence item passed)",
            "preliminary = display shortlist for inspection; it is NOT a formal candidate",
            "one pose per parent compound is selected; extra poses are kept as evidence, not as extra hits",
            f"formal qualification currently requires: {', '.join(rules)}",
        ],
    }


def _exclusion(entry, limits, ctx):
    if entry["tier"] == "D":
        return Reason.NO_CONTACTS if "no_contacts" in entry["tier_reasons"] else Reason.POSE_GEOMETRY_INVALID, "tier D"
    if entry["flags"].get("excessive_strain"):
        return Reason.EXCESSIVE_STRAIN, "strain above the failure threshold"
    if entry.get("components_missing"):
        missing = entry["components_missing"]
        if "affinity" in missing and "cnn_pose" in missing and "cnn_affinity" in missing:
            return Reason.SCORE_MISSING, "no usable docking score"
    return None


def _select(pool: list[dict[str, Any]], limit: int, *, diversity: float, prior=None) -> list[dict[str, Any]]:
    chosen: list[dict[str, Any]] = []
    prior = list(prior or [])
    remaining = sorted(pool, key=lambda e: (e["tier"], -e["evidence_score"], e["compound_id"]))
    while remaining and len(chosen) < limit:
        if not chosen and not prior:
            best = remaining[0]
        else:
            best = None
            best_key = None
            for entry in remaining:
                similarity = max((_similarity(entry, c) for c in prior + chosen), default=0.0)
                penalty = diversity * similarity
                key = (entry["tier"], -(entry["evidence_score"] - penalty), entry["compound_id"])
                if best_key is None or key < best_key:
                    best_key = key
                    best = entry
        chosen.append(best)
        remaining = [e for e in remaining if e["compound_id"] != best["compound_id"]]
    return chosen


_cache: dict[str, Any] = {}


def _similarity(a: dict[str, Any], b: dict[str, Any]) -> float:
    key = tuple(sorted((a["compound_id"], b["compound_id"])))
    if key in _cache:
        return _cache[key]
    fa = _fingerprint(a["compound_id"])
    fb = _fingerprint(b["compound_id"])
    if fa is None or fb is None:
        return 0.0
    value = rdio.morgan_similarity(fa, fb)
    _cache[key] = value
    return value


def _fingerprint(compound_id: str):
    return _cache.get(("fp", compound_id))


def prime_fingerprints(ctx, entries: list[dict[str, Any]]) -> None:
    """Pre-compute fingerprints from the stored SMILES so scoring never touches disk."""
    _cache.clear()  # Run-scoped identities/states must not reuse previous fingerprints.
    for entry in entries:
        row = ctx.db.one("SELECT canonical_smiles FROM chemical_state WHERE state_id=?", (entry["state_id"],))
        if not row or not row.get("canonical_smiles"):
            continue
        mol = rdio.mol_from_smiles(row["canonical_smiles"])
        if mol is not None:
            _cache[("fp", entry["compound_id"])] = rdio.morgan_fingerprint(mol)


def _missing_objectives(entry):
    return [k for k in ("affinity", "mechanism", "strain")
            if safe_float(entry.get("raw", {}).get(k)) is None]


def _objectives(entry):
    raw = entry["raw"]
    return (float(raw["affinity"]), -float(raw["mechanism"]), float(raw["strain"]))


def _pareto(entries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    # A missing, nonfinite or ineligible value must never masquerade as optimal.
    scored = [(e, _objectives(e)) for e in entries
              if e["tier"] != "D" and not e.get("flags", {}).get("excessive_strain")
              and not _missing_objectives(e)]
    front = [e for e, obj in scored if not any(
        all(a <= b for a, b in zip(other, obj)) and any(a < b for a, b in zip(other, obj))
        for candidate, other in scored if candidate is not e)]
    return sorted(front, key=lambda e: (e["tier"], -e["evidence_score"], e["compound_id"]))


def _select_policy(pool, limit, *, diversity, strategy):
    if strategy == "evidence_diversity":
        return _select(pool, limit, diversity=diversity)
    # Respect evidence tiers, then exhaust each raw-objective front before a
    # dominated layer. Apply diversity within each layer against all chosen rows.
    chosen = []
    for tier in sorted({e["tier"] for e in pool}):
        remaining = [e for e in pool if e["tier"] == tier and not _missing_objectives(e)]
        while remaining and len(chosen) < limit:
            front = _pareto(remaining)
            if not front:
                break
            selected = _select(front, limit - len(chosen), diversity=diversity, prior=chosen)
            chosen.extend(selected)
            ids = {e["compound_id"] for e in front}
            remaining = [e for e in remaining if e["compound_id"] not in ids]
        unknown = [e for e in pool if e["tier"] == tier and _missing_objectives(e)]
        chosen.extend(_select(unknown, max(0, limit-len(chosen)), diversity=diversity, prior=chosen))
    return chosen


def _diversity_report(selected: list[dict[str, Any]]) -> dict[str, Any]:
    if len(selected) < 2:
        return {"pairs": 0, "max_similarity": None, "mean_similarity": None}
    values = []
    for i in range(len(selected)):
        for j in range(i + 1, len(selected)):
            values.append(_similarity(selected[i], selected[j]))
    return {
        "pairs": len(values),
        "max_similarity": round(max(values), 4) if values else None,
        "mean_similarity": round(sum(values) / len(values), 4) if values else None,
        "note": "Morgan fingerprint Tanimoto over selected parents",
    }
