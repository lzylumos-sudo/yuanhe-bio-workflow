"""Y-MAVS membrane-aware evidence ranking.

This is an evidence-ordering layer, not a new binding-energy model.  GNINA's raw
affinity/CNN columns are preserved verbatim.  The ranking is lexicographic:

1. the existing evidence tier (pose QA / score / contact gate);
2. the existing within-pool evidence score;
3. completeness/quality of the membrane evidence;
4. steric and accessibility findings;
5. ensemble support.

The first two keys are intentionally the original screening decision.  Membrane
evidence can resolve equal baseline evidence, but it cannot turn a bad pose or a
missing GNINA score into a good one.  ``membrane_aware_score`` is a normalised
ordinal position in this final order, not an affinity, probability, or ΔG.
"""

from __future__ import annotations

import csv
import html
import math
import os
from typing import Any

from ..membrane.regions import REGION_UNKNOWN
from ..membrane.settings import is_enabled, scoring_rules
from ..util import atomic_write_json, now_iso, read_json, safe_float

ALGORITHM_VERSION = "ymavs-membrane-aware-ranking/1.0.0"
SCHEMA = "ymavs-membrane-aware-ranking/1"
SUPPORTED_STRATEGIES = {"evidence_first"}
TIER_ORDER = {"A": 0, "B": 1, "C": 2, "D": 3}


def _metric(ctx, pose_id: str, stage: str) -> dict[str, Any]:
    if not pose_id:
        return {}
    try:
        return ctx.db.metric_map("pose", pose_id, stage)
    except Exception:  # noqa: BLE001 - ranking must not crash on a missing assessment table
        return {}


def _value(metrics: dict[str, Any], key: str) -> Any:
    row = metrics.get(key) or {}
    if row.get("value") is not None:
        return row.get("value")
    return row.get("value_text")


def _item(
    value: Any,
    status: str,
    *,
    method: str,
    source: str,
    confidence_level: str,
    unit: str | None = None,
    detail: str = "",
) -> dict[str, Any]:
    return {
        "value": value,
        "status": status,
        "method": method,
        "source": source,
        "confidence_level": confidence_level,
        "unit": unit,
        "detail": detail,
    }


class MembraneAwareScorer:
    def __init__(self, ctx, *, settings: dict[str, Any]):
        self.ctx = ctx
        self.settings = settings
        self.rules = scoring_rules(settings)
        if self.rules["strategy"] not in SUPPORTED_STRATEGIES:
            raise ValueError(
                f"scoring strategy '{self.rules['strategy']}' is not implemented; "
                f"supported: {sorted(SUPPORTED_STRATEGIES)} (weighted_experimental remains a reserved interface)"
            )
        self.stage = read_json(os.path.join(ctx.stage_dir("membrane"), "summary.json"), {}) or {}
        ensemble_path = os.path.join(ctx.stage_dir("membrane"), "ensemble_analysis.json")
        self.ensemble = read_json(ensemble_path, {}) or {}

    def run(self, entries: list[dict[str, Any]]) -> dict[str, Any]:
        if not is_enabled(self.ctx) or not self.rules.get("enabled", True):
            return {
                "schema": SCHEMA,
                "algorithm_version": ALGORITHM_VERSION,
                "status": "disabled",
                "enabled": False,
                "strategy": self.rules["strategy"],
                "entries": entries,
                "notes": ["membrane-aware ranking is disabled; baseline ranking is unchanged"],
            }
        baseline = sorted(
            entries,
            key=lambda entry: (
                TIER_ORDER.get(str(entry.get("tier")), 9),
                -float(entry.get("evidence_score") or 0.0),
                str(entry.get("compound_id") or ""),
            ),
        )
        for index, entry in enumerate(baseline, start=1):
            entry["baseline_rank"] = index
            self._enrich(entry)
        ordered = sorted(
            baseline,
            key=lambda entry: (
                TIER_ORDER.get(str(entry.get("tier")), 9),
                -float(entry.get("evidence_score") or 0.0),
                int(entry.get("_membrane_support_rank", 9)),
                int(entry.get("_steric_rank", 1)),
                int(entry.get("_accessibility_rank", 1)),
                int(entry.get("_ensemble_rank", 1)),
                entry.get("baseline_rank", 0),
                str(entry.get("compound_id") or ""),
            ),
        )
        total = len(ordered)
        for index, entry in enumerate(ordered, start=1):
            entry["membrane_aware_rank"] = index
            entry["rank"] = index
            entry["rank_shift"] = int(entry.get("baseline_rank", index)) - index
            entry["membrane_aware_score"] = round((total - index + 1) / total, 6) if total else None
            entry["membrane_aware_score_scope"] = "normalised ordinal priority within this pool; not affinity or probability"
            entry["ranking_reasons"] = self._reasons(entry)
            for private in ("_membrane_support_rank", "_steric_rank", "_accessibility_rank", "_ensemble_rank"):
                entry.pop(private, None)
        counts: dict[str, int] = {}
        for entry in ordered:
            status = str(entry.get("membrane_evidence_status") or "incomplete")
            counts[status] = counts.get(status, 0) + 1
        summary = {
            "schema": SCHEMA,
            "algorithm_version": ALGORITHM_VERSION,
            "generated_at": now_iso(),
            "status": "evaluated",
            "enabled": True,
            "strategy": self.rules["strategy"],
            "policy": "evidence_first_lexicographic",
            "preserve_baseline": self.rules["preserve_baseline"],
            "allow_proxy_for_exploration": self.rules["allow_proxy_for_exploration"],
            "require_verified_orientation_for_formal": self.rules["require_verified_orientation_for_formal"],
            "pose_reselection": {
                "enabled": self.rules["allow_pose_reselection"],
                "status": "disabled_by_config" if not self.rules["allow_pose_reselection"] else "not_implemented",
                "baseline_pose_preserved": True,
                "note": (
                    "the baseline representative pose is retained; per-pose membrane analyses are stored, "
                    "and strain/redocking/qualification evidence is never silently transferred between poses"
                ),
            },
            "orientation": {
                "status": self.stage.get("status"),
                "orientation_verified": self.stage.get("orientation_verified"),
                "source": (self.stage.get("geometry") or {}).get("source")
                or (self.stage.get("plane") or {}).get("source"),
                "coordinate_frame": (self.stage.get("geometry") or {}).get("coordinate_frame"),
            },
            "evidence_status_counts": counts,
            "count": total,
            "entries": ordered,
            "limitations": [
                "membrane_aware_score is an ordinal position, not a binding free energy",
                "membrane evidence resolves ties after the existing tier and evidence score",
                "missing accessibility evidence is not proof that a ligand cannot bind",
                "the original GNINA scores and baseline rank are preserved",
            ],
        }
        return summary

    def _enrich(self, entry: dict[str, Any]) -> None:
        pose_id = str(entry.get("pose_id") or "")
        compound_id = str(entry.get("compound_id") or "")
        metrics = _metric(self.ctx, pose_id, "membrane")
        orientation_status = str(_value(metrics, "orientation_status") or self.stage.get("status") or "not_evaluated")
        geometry_status = str(_value(metrics, "membrane_geometry_status") or "")
        pocket_region = str(_value(metrics, "pocket_region") or REGION_UNKNOWN)
        pocket_contacts = safe_float(_value(metrics, "pocket_contact_count"))
        depth = safe_float(_value(metrics, "ligand_center_depth_A"))
        core_fraction = safe_float(_value(metrics, "ligand_core_fraction"))
        aqueous = str(_value(metrics, "aqueous_accessibility") or "not_evaluated")
        lateral = str(_value(metrics, "membrane_lateral_accessibility") or "not_evaluated")
        steric = str(_value(metrics, "steric_compatibility_status") or "not_evaluated")
        pocket_source = str(_value(metrics, "pocket_source") or "Y-MAVS membrane stage")
        verified = bool(
            self.stage.get("orientation_verified")
            and self.stage.get("status") == "sourced"
            and geometry_status == "evaluated"
            and orientation_status not in {"proxy_only", "needs_input", "not_evaluated"}
        )
        proxy = orientation_status == "proxy_only" or geometry_status == "proxy_only"
        orientation_item = _item(
            orientation_status,
            "pass" if verified else ("flag" if proxy else "not_evaluated"),
            method="PPM3 rigid-construct transfer and receptor-frame check"
            if verified
            else "membrane placement proxy or no verified source",
            source="membrane orientation evidence",
            confidence_level="high" if verified else ("low" if proxy else "none"),
            detail="computational placement, not an experimental inside/outside assignment",
        )
        pocket_item = _item(
            {
                "pocket_region": pocket_region,
                "pocket_contact_count": None if pocket_contacts is None else int(pocket_contacts),
                "ligand_center_depth_A": depth,
                "ligand_core_fraction": core_fraction,
            },
            "pass" if pocket_region != REGION_UNKNOWN and pocket_contacts not in (None, 0) else "not_evaluated",
            method="ligand-defined receptor contact set plus signed membrane depth",
            source=pocket_source,
            confidence_level="medium" if verified else ("low" if proxy else "none"),
            detail="geometric localization only; not a druggability or affinity prediction",
        )
        ligand_item = _item(
            {"ligand_center_depth_A": depth, "ligand_core_fraction": core_fraction},
            "pass" if depth is not None else "not_evaluated",
            method="signed membrane-plane depth",
            source="Y-MAVS membrane region analysis",
            confidence_level="high" if verified else ("low" if proxy else "none"),
            unit="angstrom",
            detail="position relative to the supplied membrane plane; no favorable direction is implied",
        )
        steric_item = _item(
            steric,
            "pass" if steric == "compatible" else ("flag" if steric == "minor_overlap" else ("fail" if steric == "clash" else "not_evaluated")),
            method="element van der Waals overlap proxy",
            source="Y-MAVS compatibility stage",
            confidence_level="medium",
            detail="not a replacement for PoseBusters and not an energy",
        )
        accessibility_item = _item(
            {
                "aqueous": aqueous,
                "membrane_lateral": lateral,
                "geometric_path_length_A": safe_float(_value(metrics, "geometric_path_length_A")),
            },
            "pass" if aqueous == "accessible" or lateral == "accessible" else ("flag" if aqueous != "not_evaluated" or lateral != "not_evaluated" else "not_evaluated"),
            method="1 A grid connectivity with 1.4 A probe",
            source="Y-MAVS accessibility stage",
            confidence_level="low",
            unit="angstrom",
            detail="geometric proxy only; a static closed path is not proof of no binding",
        )
        ensemble = self._ensemble_for(compound_id)
        ensemble_item = _item(
            ensemble.get("ensemble_status"),
            "pass" if ensemble.get("ensemble_status") == "complete" else ("flag" if ensemble.get("ensemble_status") == "partial" else "not_evaluated"),
            method="per-state evidence integration without cross-state score averaging",
            source="Y-MAVS ensemble stage",
            confidence_level="low",
            detail="single-state evidence can be biologically meaningful; coverage is not potency",
        )
        items = {
            "orientation_validity": orientation_item,
            "pocket_membrane_localization": pocket_item,
            "ligand_membrane_geometry": ligand_item,
            "protein_ligand_steric_compatibility": steric_item,
            "geometric_pocket_accessibility": accessibility_item,
            "multi_conformation_evidence": ensemble_item,
        }
        required_present = [depth is not None, pocket_region != REGION_UNKNOWN and pocket_contacts not in (None, 0)]
        accessibility_enabled = bool((self.settings.get("accessibility") or {}).get("enabled", True))
        if accessibility_enabled:
            required_present.append(aqueous not in ("", "not_evaluated") or lateral not in ("", "not_evaluated"))
        if not all(required_present):
            evidence_status = "incomplete"
        elif not verified and not proxy:
            evidence_status = "incomplete"
        elif proxy:
            evidence_status = "proxy_only"
        else:
            evidence_status = "evaluated"
        entry["membrane_evidence"] = items
        entry["membrane_evidence_status"] = evidence_status
        entry["pocket_region"] = pocket_region
        entry["ligand_center_depth_A"] = depth
        entry["ligand_core_fraction"] = core_fraction
        entry["geometric_accessibility_status"] = (
            "accessible" if (aqueous == "accessible" or lateral == "accessible")
            else ("uncertain" if (aqueous == "uncertain" or lateral == "uncertain") else "not_evaluated")
        )
        entry["aqueous_accessibility"] = aqueous
        entry["membrane_lateral_accessibility"] = lateral
        entry["steric_compatibility_status"] = steric
        entry["ensemble_status"] = ensemble.get("ensemble_status") or "not_evaluated"
        entry["ensemble_state_coverage"] = ensemble.get("state_coverage")
        entry["pose_reselection_status"] = (
            "disabled_by_config" if not self.rules["allow_pose_reselection"] else "not_implemented"
        )
        entry["_membrane_support_rank"] = 0 if verified and evidence_status == "evaluated" else (1 if proxy else 2)
        entry["_steric_rank"] = 2 if steric == "clash" else (1 if steric == "minor_overlap" else 0)
        entry["_accessibility_rank"] = 0 if entry["geometric_accessibility_status"] == "accessible" else 1
        entry["_ensemble_rank"] = 0 if ensemble.get("ensemble_status") == "complete" else 1

    def _ensemble_for(self, compound_id: str) -> dict[str, Any]:
        compounds = self.ensemble.get("compounds") or {}
        item = compounds.get(compound_id) or {}
        if not self.ensemble or self.ensemble.get("enabled") is False:
            return {"ensemble_status": "disabled", "state_coverage": None}
        return {
            "ensemble_status": item.get("ensemble_status") or "not_evaluated",
            "state_coverage": item.get("state_coverage"),
            "number_of_evaluated_states": item.get("number_of_evaluated_states"),
            "number_of_missing_states": item.get("number_of_missing_states"),
        }

    def _reasons(self, entry: dict[str, Any]) -> list[str]:
        reasons = [
            f"baseline tier {entry.get('tier')}",
            f"baseline GNINA evidence score {round(float(entry.get('evidence_score') or 0.0), 4)}",
            f"membrane evidence {entry.get('membrane_evidence_status')}",
            f"pocket region {entry.get('pocket_region')}",
        ]
        if entry.get("ligand_center_depth_A") is not None:
            reasons.append(f"ligand center depth {entry['ligand_center_depth_A']} A")
        reasons.append(f"steric {entry.get('steric_compatibility_status')}")
        reasons.append(f"accessibility {entry.get('geometric_accessibility_status')}")
        reasons.append(f"ensemble {entry.get('ensemble_status')}")
        shift = int(entry.get("rank_shift") or 0)
        if shift > 0:
            reasons.append(f"moved up {shift} position(s) among equal-or-better baseline evidence")
        elif shift < 0:
            reasons.append(f"moved down {abs(shift)} position(s) among equal-or-better baseline evidence")
        else:
            reasons.append("baseline order preserved")
        return reasons


def write_evidence(directory: str, entries: list[dict[str, Any]], summary: dict[str, Any]) -> str:
    payload = {
        "schema": "ymavs-membrane-evidence/1",
        "algorithm_version": ALGORITHM_VERSION,
        "generated_at": now_iso(),
        "orientation": summary.get("orientation"),
        "evidence_status_counts": summary.get("evidence_status_counts"),
        "compounds": {
            str(entry.get("compound_id")): {
                "compound_id": entry.get("compound_id"),
                "state_id": entry.get("state_id"),
                "pose_id": entry.get("pose_id"),
                "baseline_rank": entry.get("baseline_rank"),
                "membrane_aware_rank": entry.get("membrane_aware_rank"),
                "rank_shift": entry.get("rank_shift"),
                "membrane_evidence_status": entry.get("membrane_evidence_status"),
                "evidence_items": entry.get("membrane_evidence") or {},
                "ranking_reasons": entry.get("ranking_reasons") or [],
            }
            for entry in entries
        },
        "limitations": summary.get("limitations") or [],
    }
    path = os.path.join(directory, "membrane_evidence.json")
    atomic_write_json(path, payload)
    return path


def write_ranking_files(ctx, ranking: dict[str, Any], selection: dict[str, Any]) -> dict[str, str]:
    """Write final Y-MAVS ranking after formal qualification is known."""
    directory = ctx.reports_dir()
    membrane_summary = ranking.get("membrane_aware") or {}
    entries = ranking.get("ranking") or []
    if not entries or membrane_summary.get("status") != "evaluated":
        return {}
    verdicts = selection.get("qualification") or {}
    verdicts = selection.get("qualification") or {}
    exclusions = {
        str(item.get("compound_id")): item for item in (selection.get("exclusions") or [])
    }
    rows = []
    for entry in entries:
        verdict = verdicts.get(entry.get("compound_id")) or {}
        exclusion = exclusions.get(str(entry.get("compound_id"))) or {}
        if verdict:
            formal_qualified = verdict.get("formal_qualified")
            qualification_status = verdict.get("qualification") or verdict.get("qualification_status") or "not_evaluated"
        elif exclusion:
            formal_qualified = False
            qualification_status = "excluded"
        else:
            formal_qualified = None
            qualification_status = "not_evaluated"
        ranking_reasons = list(entry.get("ranking_reasons") or [])
        if exclusion:
            ranking_reasons.append(
                f"excluded before qualification: {exclusion.get('reason_code')} ({exclusion.get('detail')})"
            )
        rows.append(
            {
                "compound_id": entry.get("compound_id"),
                "state_id": entry.get("state_id"),
                "pose_id": entry.get("pose_id"),
                "baseline_rank": entry.get("baseline_rank"),
                "membrane_aware_rank": entry.get("membrane_aware_rank"),
                "rank_shift": entry.get("rank_shift"),
                "membrane_aware_score": entry.get("membrane_aware_score"),
                "gnina_affinity": (entry.get("raw") or {}).get("affinity"),
                "cnn_score": (entry.get("raw") or {}).get("cnn_pose"),
                "cnn_affinity": (entry.get("raw") or {}).get("cnn_affinity"),
                "pocket_region": entry.get("pocket_region"),
                "ligand_center_depth_A": entry.get("ligand_center_depth_A"),
                "ligand_core_fraction": entry.get("ligand_core_fraction"),
                "geometric_accessibility_status": entry.get("geometric_accessibility_status"),
                "steric_compatibility_status": entry.get("steric_compatibility_status"),
                "ensemble_status": entry.get("ensemble_status"),
                "membrane_evidence_status": entry.get("membrane_evidence_status"),
                "formal_qualified": formal_qualified,
                "qualification_status": qualification_status,
                "exclusion_reason_code": exclusion.get("reason_code"),
                "ranking_reasons": "; ".join(ranking_reasons),
            }
        )
    fields = [
        "compound_id", "state_id", "pose_id", "baseline_rank", "membrane_aware_rank", "rank_shift",
        "membrane_aware_score", "gnina_affinity", "cnn_score", "cnn_affinity", "pocket_region",
        "ligand_center_depth_A", "ligand_core_fraction", "geometric_accessibility_status",
        "steric_compatibility_status", "ensemble_status", "membrane_evidence_status",
        "formal_qualified", "qualification_status", "exclusion_reason_code", "ranking_reasons",
    ]
    csv_path = os.path.join(directory, "membrane_aware_ranking.csv")
    _write_csv(csv_path, fields, rows)
    json_path = os.path.join(directory, "membrane_aware_ranking.json")
    atomic_write_json(
        json_path,
        {
            "schema": SCHEMA,
            "algorithm_version": ALGORITHM_VERSION,
            "generated_at": now_iso(),
            "run_id": ctx.run_id,
            "config_digest": ctx.config.digest(),
            "strategy": membrane_summary.get("strategy"),
            "score_scope": "within_candidate_pool_only; membrane_aware_score is ordinal and not affinity",
            "orientation": membrane_summary.get("orientation"),
            "evidence_status_counts": membrane_summary.get("evidence_status_counts"),
            "pose_reselection": membrane_summary.get("pose_reselection"),
            "ranking": rows,
            "full_ranking": entries,
            "limitations": membrane_summary.get("limitations") or [],
        },
    )
    evidence_path = write_evidence(directory, entries, membrane_summary)
    html_path = os.path.join(directory, "membrane_screening_report.html")
    _write_html(html_path, ctx, rows, membrane_summary)
    return {
        "membrane_aware_ranking": csv_path,
        "membrane_aware_ranking_json": json_path,
        "membrane_evidence": evidence_path,
        "membrane_screening_report": html_path,
    }


def _write_csv(path: str, fields: list[str], rows: list[dict[str, Any]]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_html(path: str, ctx, rows: list[dict[str, Any]], summary: dict[str, Any]) -> None:
    esc = html.escape
    orientation = summary.get("orientation") or {}
    status_counts = summary.get("evidence_status_counts") or {}
    table_rows = []
    for row in rows[:100]:
        warnings = []
        if row.get("membrane_evidence_status") == "proxy_only":
            warnings.append("membrane orientation is proxy-only")
        if row.get("membrane_evidence_status") == "incomplete":
            warnings.append("membrane evidence incomplete")
        if row.get("steric_compatibility_status") == "clash":
            warnings.append("van der Waals clash proxy")
        if row.get("ensemble_status") in {"partial", "missing", "not_evaluated"}:
            warnings.append(f"ensemble {row.get('ensemble_status')}")
        table_rows.append(
            "<tr>"
            f"<td>{row.get('membrane_aware_rank')}</td>"
            f"<td>{row.get('baseline_rank')}</td>"
            f"<td>{row.get('rank_shift')}</td>"
            f"<td><code>{esc(str(row.get('compound_id')))}</code></td>"
            f"<td>{_fmt(row.get('gnina_affinity'))}</td>"
            f"<td>{_fmt(row.get('cnn_score'))}</td>"
            f"<td>{esc(str(row.get('pocket_region')))}</td>"
            f"<td>{_fmt(row.get('ligand_center_depth_A'))}</td>"
            f"<td>{_fmt(row.get('ligand_core_fraction'))}</td>"
            f"<td>{esc(str(row.get('geometric_accessibility_status')))}</td>"
            f"<td>{esc(str(row.get('steric_compatibility_status')))}</td>"
            f"<td>{esc(str(row.get('ensemble_status')))}</td>"
            f"<td>{esc(str(row.get('membrane_evidence_status')))}</td>"
            f"<td>{esc(str(row.get('formal_qualified')))}</td>"
            f"<td>{esc(str(row.get('qualification_status')))}</td>"
            f"<td>{esc('; '.join(warnings) or '-')}</td>"
            "</tr>"
        )
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Y-MAVS membrane-aware screening - {esc(ctx.run_id)}</title>
<style>
body {{ font-family: -apple-system, Segoe UI, Arial, sans-serif; margin: 2rem; color: #17202a; }}
table {{ border-collapse: collapse; width: 100%; font-size: 0.82rem; }}
th, td {{ border: 1px solid #d5d8dc; padding: 0.35rem 0.45rem; text-align: left; vertical-align: top; }}
th {{ background: #f4f6f7; position: sticky; top: 0; }}
.notice {{ background: #fef9e7; border-left: 4px solid #d4ac0d; padding: 0.8rem 1rem; margin: 1rem 0; }}
.meta {{ color: #566573; font-size: 0.9rem; }}
code {{ font-family: Consolas, monospace; }}
</style>
</head>
<body>
<h1>Y-MAVS membrane-aware screening report</h1>
<p class="meta">Run <code>{esc(ctx.run_id)}</code> · target <code>{esc(ctx.config.target.target_id)}</code> · generated {esc(now_iso())}</p>
<div class="notice">
<strong>Interpretation boundary.</strong> This report ranks candidates using membrane geometry evidence.
<code>membrane_aware_score</code> is an ordinal priority within this candidate pool, not a binding free energy,
activity probability, or validated hit-rate improvement. Membrane orientation:
<code>{esc(str(orientation.get('status')))}</code>, source <code>{esc(str(orientation.get('source')))}</code>.
Accessibility is a static geometric proxy; a missing static path is not proof that binding is impossible.
</div>
<p>Evidence status counts: {esc(str(status_counts))}. Formal qualification remains governed by the existing
qualification rules and is not upgraded by this ranking. Pose reselection:
<code>{esc(str((summary.get('pose_reselection') or {}).get('status', 'unknown')))}</code>.</p>
<table>
<thead><tr>
<th>Y-MAVS rank</th><th>Baseline rank</th><th>Shift</th><th>Compound</th><th>GNINA affinity</th>
<th>CNNscore</th><th>Pocket region</th><th>Ligand center depth (Å)</th><th>Core fraction</th>
<th>Accessibility</th><th>Steric</th><th>Ensemble</th><th>Membrane evidence</th>
<th>Formal</th><th>Qualification</th><th>Warnings</th>
</tr></thead>
<tbody>
{''.join(table_rows)}
</tbody></table>
<p class="meta">Generated by Y-MAVS {esc(ALGORITHM_VERSION)}. Files: membrane_aware_ranking.csv/json,
membrane_evidence.json, membrane_geometry.json, pocket_analysis.json, accessibility_analysis.json,
membrane_compatibility.json, ensemble_analysis.json.</p>
</body></html>
"""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(document)


def _fmt(value: Any) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(number):
        return ""
    return f"{number:.3f}"
