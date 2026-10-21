"""Report rendering: one honest, reproducible record of what was and was not done."""

from __future__ import annotations

import csv
import os
from typing import Any

from ..util import atomic_write_json, atomic_write_text, now_iso, read_json, relpath

DONE = "Implemented and validated with real data"
IMPL = "Implemented; real-data validation pending"
BLOCKED = "Blocked by external prerequisites"
NOTDONE = "Not implemented"


def write_all(ctx, state: dict[str, Any]) -> dict[str, Any]:
    directory = ctx.reports_dir()
    stages = _load_stages(ctx, state)
    ranking = stages.get("scoring") or {}
    selection = stages.get("selection") or {}
    md = stages.get("md") or {}
    target = stages.get("target") or {}
    imported = stages.get("import") or {}
    qa = stages.get("pose_qa") or {}
    mech = stages.get("mechanism") or {}
    redock = stages.get("redock") or {}
    membrane = stages.get("membrane") or {}
    membrane_ranking = read_json(os.path.join(directory, "membrane_aware_ranking.json"), {}) or {}
    membrane_files = {
        "membrane_aware_ranking": os.path.join(directory, "membrane_aware_ranking.csv"),
        "membrane_aware_ranking_json": os.path.join(directory, "membrane_aware_ranking.json"),
        "membrane_evidence": os.path.join(directory, "membrane_evidence.json"),
        "membrane_screening_report": os.path.join(directory, "membrane_screening_report.html"),
        "membrane_geometry": os.path.join(ctx.stage_dir("membrane"), "membrane_geometry.json"),
        "pocket_analysis": os.path.join(ctx.stage_dir("membrane"), "pocket_analysis.json"),
        "accessibility_analysis": os.path.join(ctx.stage_dir("membrane"), "accessibility_analysis.json"),
        "accessibility_paths": os.path.join(ctx.stage_dir("membrane"), "accessibility_paths.json"),
        "membrane_compatibility": os.path.join(ctx.stage_dir("membrane"), "membrane_compatibility.json"),
        "ensemble_analysis": os.path.join(ctx.stage_dir("membrane"), "ensemble_analysis.json"),
        "pose_membrane_evidence": os.path.join(ctx.stage_dir("membrane"), "pose_membrane_evidence.json"),
    }
    membrane_files = {name: path for name, path in membrane_files.items() if os.path.isfile(path)}

    candidates = os.path.join(directory, "candidates.csv")
    _write_rows(
        candidates,
        ["role", "order", "compound_id", "state_id", "pose_id", "tier", "qualification", "evidence_score",
         "components_used", "components_missing", "selection_version"],
        [
            {
                "role": entry.get("selection_role"),
                "order": entry.get("selection_order"),
                "compound_id": entry["compound_id"],
                "state_id": entry["state_id"],
                "pose_id": entry["pose_id"],
                "tier": entry["tier"],
                "qualification": entry.get("qualification"),
                "evidence_score": entry.get("evidence_score"),
                "components_used": ";".join(entry.get("components_used") or []),
                "components_missing": ";".join(entry.get("components_missing") or []),
                "selection_version": entry.get("selection_version"),
            }
            for entry in (selection.get("primary") or [])
            + (selection.get("preliminary") or [])
            + (selection.get("exploratory") or [])
        ],
    )
    qualification_csv = os.path.join(directory, "qualification.csv")
    qualification_rows: list[dict[str, Any]] = []
    for compound_id, verdict in sorted((selection.get("qualification") or {}).items()):
        for requirement in verdict.get("requirements") or []:
            qualification_rows.append(
                {
                    "compound_id": compound_id,
                    "pose_id": verdict.get("pose_id"),
                    "tier": verdict.get("tier"),
                    "formal_qualified": verdict.get("formal_qualified"),
                    "qualification": verdict.get("qualification"),
                    "requirement": requirement.get("requirement"),
                    "status": requirement.get("status"),
                    "satisfied": requirement.get("satisfied"),
                    "observed": _compact(requirement.get("observed")),
                    "reason_code": requirement.get("reason_code") or "",
                    "detail": requirement.get("detail"),
                }
            )
    _write_rows(
        qualification_csv,
        ["compound_id", "pose_id", "tier", "formal_qualified", "qualification",
         "requirement", "status", "satisfied", "observed", "reason_code", "detail"],
        qualification_rows,
    )
    exclusions = os.path.join(directory, "exclusions.csv")
    _write_rows(exclusions, ["compound_id", "pose_id", "tier", "reason_code", "detail"], selection.get("exclusions") or [])
    ranking_csv = os.path.join(directory, "ranking.csv")
    _write_rows(
        ranking_csv,
        ["rank", "compound_id", "tier", "evidence_score", "engine", "score_unit", "cnn_status", "vina_affinity_kcal_mol", "affinity", "cnn_pose", "cnn_affinity",
         "mechanism", "pose_quality", "missing"],
        [
            {
                "rank": entry.get("rank"),
                "engine":entry.get("engine"), "score_unit":entry.get("score_unit"), "cnn_status":entry.get("cnn_status"),
                "vina_affinity_kcal_mol":(entry.get("raw") or {}).get("affinity") if entry.get("engine")=="vina" else None,
                "compound_id": entry["compound_id"],
                "tier": entry["tier"],
                "evidence_score": entry.get("evidence_score"),
                "affinity": (entry.get("raw") or {}).get("affinity"),
                "cnn_pose": (entry.get("raw") or {}).get("cnn_pose"),
                "cnn_affinity": (entry.get("raw") or {}).get("cnn_affinity"),
                "mechanism": (entry.get("raw") or {}).get("mechanism"),
                "pose_quality": (entry.get("raw") or {}).get("pose_quality"),
                "missing": ";".join(entry.get("components_missing") or []),
            }
            for entry in ranking.get("ranking") or []
        ],
    )

    decision = {
        "schema": "screening-decision/2",
        "selection_version": selection.get("selection_version"),
        "strategy": selection.get("strategy", "legacy"),
        "score_scope": ranking.get("score_scope", "legacy/unspecified"),
        "pareto_front": selection.get("pareto_front", []),
        "formal_pareto_front": selection.get("formal_pareto_front", []),
        "pareto_unassessed": selection.get("pareto_unassessed", {}),
        "weight_diagnostics": ranking.get("diagnostics", {}),
    }
    decision_path = os.path.join(directory, "decision.json")
    atomic_write_json(decision_path, decision)
    evidence_path = os.path.join(directory, "candidate_evidence.csv")
    _write_rows(evidence_path,
        ["compound_id", "pose_id", "representative_policy", "qa", "affinity_kcal_mol", "strain_kcal_mol",
         "region_contacts", "evidence_coverage", "pareto", "formal", "best_rank", "worst_rank", "missing"],
        [{"compound_id": e["compound_id"], "pose_id": e["pose_id"],
          "representative_policy": e.get("representative_policy"), "qa": e.get("representative_qa"),
          "affinity_kcal_mol": e.get("raw", {}).get("affinity"),
          "strain_kcal_mol": e.get("raw", {}).get("strain"),
          "region_contacts": e.get("raw", {}).get("mechanism"),
          "evidence_coverage": e.get("evidence_coverage"),
          "pareto": e["compound_id"] in decision["pareto_front"],
          "formal": (selection.get("qualification", {}).get(e["compound_id"], {})).get("formal_qualified", False),
          "best_rank": (ranking.get("diagnostics", {}).get("rank_ranges", {}).get(e["compound_id"], {})).get("best"),
          "worst_rank": (ranking.get("diagnostics", {}).get("rank_ranges", {}).get(e["compound_id"], {})).get("worst"),
          "missing": ";".join(e.get("components_missing", []))}
         for e in ranking.get("ranking", [])])

    markdown = _markdown(ctx, state, ranking, selection, md, target, imported, qa, mech, redock, membrane, membrane_ranking)
    report_path = os.path.join(directory, "report.md")
    atomic_write_text(report_path, markdown)
    payload = {
        "run_id": ctx.run_id,
        "created": now_iso(),
        "config_digest": ctx.config.digest(),
        "code_version": ctx.config.code_version,
        "selection_version": selection.get("selection_version"),
        "candidates": len((selection.get("primary") or []) + (selection.get("preliminary") or []) + (selection.get("exploratory") or [])),
        "formal_candidates": len(selection.get("primary") or []),
        "preliminary_candidates": len(selection.get("preliminary") or []),
        "blocked": _blockers(md, redock, membrane),
        "files": {
            "report": relpath(report_path, ctx.run_dir),
            "candidates": relpath(candidates, ctx.run_dir),
            "qualification": relpath(qualification_csv, ctx.run_dir),
            "exclusions": relpath(exclusions, ctx.run_dir),
            "ranking": relpath(ranking_csv, ctx.run_dir),
            "decision": relpath(decision_path, ctx.run_dir),
            "candidate_evidence": relpath(evidence_path, ctx.run_dir),
            **{name: relpath(path, ctx.run_dir) for name, path in membrane_files.items()},
        },
    }
    atomic_write_json(os.path.join(directory, "report.json"), payload)
    return payload


def _load_stages(ctx, state: dict[str, Any]) -> dict[str, Any]:
    """Merge in-memory stage results with summaries already on disk.

    The report must describe the whole run even when it is regenerated by a
    later invocation that only re-ran a subset of stages.
    """
    from ..util import read_json

    merged: dict[str, Any] = {}
    directory = os.path.join(ctx.run_dir, "stages")
    if os.path.isdir(directory):
        for name in sorted(os.listdir(directory)):
            path = os.path.join(directory, name, "summary.json")
            if os.path.isfile(path):
                merged[name] = read_json(path, {}) or {}
    merged.update({k: v for k, v in (state.get("stages") or {}).items() if v})
    return merged


def _blockers(md, redock, membrane) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    if redock.get("status") == "blocked_external_tool":
        out.append({"stage": "redock", "code": redock.get("reason_code"), "detail": redock.get("note")})
    capability = (redock.get("capability") or {})
    if capability and not capability.get("available"):
        out.append({"stage": "redock", "code": "vina_unavailable", "detail": capability.get("reason")})
    md_cap = (membrane.get("md_capability") or {})
    if md_cap.get("status"):
        out.append({"stage": "md", "code": md_cap.get("status"), "detail": md_cap.get("note")})
    if md.get("blocked"):
        out.append({"stage": "md", "code": "packages_blocked", "detail": f"{md.get('blocked')} package(s) not executable"})
    return out


def _markdown(ctx, state, ranking, selection, md, target, imported, qa, mech, redock, membrane, membrane_ranking=None) -> str:
    lines: list[str] = []
    profile = ctx.config.target
    lines += [
        f"# ZVSS Screening report {ctx.run_id}",
        "",
        f"- Target: `{profile.target_id}` ({profile.display_name or profile.target_id})",
        f"- Reference structure: `{profile.primary.pdb_id}` chain `{profile.primary.chain}`",
        f"- Configuration digest: `{ctx.config.digest()}` | Code version `{ctx.config.code_version}`",
        f"- Budget: `{ctx.config.budget.budget_id}` | Selection version: `{(selection or {}).get('selection_version')}`",
        "",
        "Current docking: official Vina; affinity in kcal/mol, lower is better. CNN fields are N/A for Vina. Historical engines are ranked separately.",
        "MD research protocol: GAFF2 + HF/RESP; one 500 ns trajectory; core RMSD and PB window 400–500 ns. This screening report does not certify MD completion.",
        "## 1. Inputs and identity",
        "",
        f"- Files discovered: {len(imported.get('files') or [])} items",
        f"- Pose: Read {imported.get('poses_seen', 0)} / Import {imported.get('poses_imported', 0)} / Failed {imported.get('poses_failed', 0)}",
        f"- Compounds {imported.get('parents', 0)} | Chemical states {imported.get('states', 0)} | Docking jobs {imported.get('jobs', 0)}",
        f"- Receptor: {(imported.get('receptor') or {}).get('status')} | Box source: {(imported.get('box') or {}).get('source', 'Not found')}",
    ]
    if imported.get("needs_input"):
        lines.append(f"- Chemical structure requires manual completion: {len(imported['needs_input'])} entries")
    if imported.get("conflicts"):
        lines.append(f"- Identity conflicts: {len(imported['conflicts'])} entries, retained separately without merging")
    lines.append("")

    lines += ["## 2. Validation and local mechanism", ""]
    lines.append(f"- Pose QC: Passed {qa.get('passed', 0)} / Flagged {qa.get('flagged', 0)} / Failed {qa.get('failed', 0)}")
    pb = qa.get("posebusters") or {}
    lines.append(f"- PoseBusters: {'Available' if pb.get('available') else 'Unavailable'} {pb.get('reason') or ''}".rstrip())
    lines.append(f"- Contact analysis: Evaluated {mech.get('evaluated', 0)} poses, truncated {mech.get('cutoff')} Å")
    ref = mech.get("reference") or {}
    lines.append(f"- Reference-ligand contact baseline: {ref.get('status')} ({ref.get('contact_count', 'n/a')})")
    lines.append("")

    lines += ["## 3. Evidence tiers and ranking", ""]
    tiers = ranking.get("tiers") or {}
    lines.append(f"- Tier counts: A={tiers.get('A', 0)} B={tiers.get('B', 0)} C={tiers.get('C', 0)} D={tiers.get('D', 0)}")
    weights = ranking.get("weights_used") or {}
    lines.append(f"- Weights used: {weights}")
    skipped = ranking.get("weights_skipped") or []
    if skipped:
        lines.append(f"- Weights skipped for missing data: {skipped}")
    lines.append("")
    top = (ranking.get("ranking") or [])[:15]
    if top:
        lines += ["| Rank | Compounds | Tier | Composite score | Missing evidence |", "| --- | --- | --- | --- | --- |"]
        for entry in top:
            lines.append(
                f"| {entry.get('rank')} | `{entry['compound_id']}` | {entry['tier']} | "
                f"{entry.get('evidence_score'):.4f} | {', '.join(entry.get('components_missing') or []) or '-'} |"
            )
        lines.append("")

    lines += ["### Decision robustness", "",
              f"- Selection strategy: {selection.get('strategy', 'legacy')}",
              "- Scores are relative to this candidate pool, not comparable across runs and not activity probabilities.",
              "- decision.json contains fixed-pool weight perturbations and component ablations; rank ranges are not confidence intervals.",
              f"- Complete-evidence eligible Pareto front: {len(selection.get('pareto_front', []))}; "
              f"unassessed because objectives are missing: {len(selection.get('pareto_unassessed', {}))}.",
              "- candidate_evidence.csv preserves raw units, pose identity, missing values and evidence coverage.", ""]
    membrane_rows = membrane_ranking.get("ranking") or []
    if membrane_rows:
        orientation = membrane_ranking.get("orientation") or {}
        lines += [
            "### Y-MAVS membrane-aware evidence",
            "",
            f"- Membrane orientation: {orientation.get('status', 'not_evaluated')} "
            f"(source: {orientation.get('source') or 'none'}; verified: {orientation.get('orientation_verified')})",
            f"- Evidence-status counts: {membrane_ranking.get('evidence_status_counts') or {}}",
            "- membrane_aware_score is an ordinal priority within this pool, not a binding free energy or activity probability.",
            "- Membrane evidence is applied after the existing tier and evidence score; formal qualification is unchanged except for explicitly configured additive requirements.",
            "",
            "| Y-MAVS rank | Baseline | Shift | Compound | Affinity kcal/mol | Pocket region | Depth (Å) | Accessibility | Steric | Ensemble | Evidence | Formal |",
            "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
        ]
        for row in membrane_rows[:15]:
            lines.append(
                f"| {row.get('membrane_aware_rank')} | {row.get('baseline_rank')} | {row.get('rank_shift')} | "
                f"`{row.get('compound_id')}` | {row.get('docking_affinity_kcal_mol',row.get('gnina_affinity'))} | {row.get('pocket_region')} | "
                f"{row.get('ligand_center_depth_A')} | {row.get('geometric_accessibility_status')} | "
                f"{row.get('steric_compatibility_status')} | {row.get('ensemble_status')} | "
                f"{row.get('membrane_evidence_status')} | {row.get('formal_qualified')} |"
            )
        lines.append("")
        lines.append("Full tables: `membrane_aware_ranking.csv`, `membrane_evidence.json`, `membrane_screening_report.html`.")
        lines.append("")

    lines += ["## 4. Candidates and exclusions", ""]
    primary = (selection or {}).get("primary") or []
    preliminary = (selection or {}).get("preliminary") or []
    exploratory = (selection or {}).get("exploratory") or []
    verdicts = (selection or {}).get("qualification") or {}
    formal_count = (selection or {}).get("formal_qualified")
    lines.append(
        f"- Qualified primary candidates: {len(primary)} | Preliminary display candidates: {len(preliminary)} | "
        f"Exploratory candidates: {len(exploratory)} | Excluded: {(selection or {}).get('excluded', 0)}"
    )
    lines.append(
        f"- Compounds meeting all formal requirements: {formal_count if formal_count is not None else 'Not computed'} "
        f"(Formal requirements: {', '.join((selection or {}).get('qualification_rules') or []) or 'Not configured'})"
    )
    if not primary:
        lines.append(
            "- **No qualified primary candidates in this run**: Required evidence missing; details in `qualification.csv`), "
            "Rules were not relaxed; display lists do not establish scientific qualification."
        )
    if primary or preliminary or exploratory:
        lines += [
            "",
            "| Role | Index | Compounds | Tier | Qualification | Composite score | Unmet requirements |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
        for entry in primary + preliminary + exploratory:
            verdict = verdicts.get(entry["compound_id"]) or {}
            unmet = ", ".join(verdict.get("unmet_requirements") or []) or "-"
            lines.append(
                f"| {entry.get('selection_role')} | {entry.get('selection_order')} | `{entry['compound_id']}` | "
                f"{entry['tier']} | {verdict.get('qualification', '-')} | {entry.get('evidence_score'):.4f} | {unmet} |"
            )
        lines.append("")
    # per-candidate qualification detail for the shortlisted candidates
    shortlist = primary + preliminary + exploratory
    if shortlist:
        lines += ["Candidate qualification details:", "", "| Compounds | Requirement | Status | Observed value | Description |", "| --- | --- | --- | --- | --- |"]
        for entry in shortlist:
            verdict = verdicts.get(entry["compound_id"]) or {}
            for requirement in verdict.get("requirements") or []:
                lines.append(
                    f"| `{entry['compound_id']}` | {requirement['requirement']} | {requirement['status']} | "
                    f"{_compact(requirement.get('observed'))} | {requirement.get('detail', '')} |"
                )
        lines.append("")
    exclusions = (selection or {}).get("exclusions") or []
    if exclusions:
        lines += ["Selected exclusion reasons:", ""]
        for item in exclusions[:10]:
            lines.append(f"- `{item['compound_id']}`: {item['reason_code']} ({item['detail']})")
    lines.append("")

    lines += ["## 5. MD Task bundles", ""]
    lines.append(f"- Task bundles generated: {md.get('count', 0)} | Executable: {md.get('ready', 0)} | Blocked: {md.get('blocked', 0)}")
    for package in (md.get("packages") or [])[:6]:
        lines.append(f"- `{os.path.basename(package.get('packages', ''))}` → {package.get('execution_status')}")
    lines.append("")

    lines += ["## 6. Explicitly incomplete / Blocked", ""]
    blockers = _blockers(md, redock, membrane)
    if blockers:
        for item in blockers:
            lines.append(f"- `{item['stage']}` / `{item['code']}`: {item.get('detail')}")
    else:
        lines.append("- No recorded blockers; this only means none were triggered in this run")
    lines += [
        f"- Membrane orientation: {membrane.get('status', 'not_evaluated')}; source: {(membrane.get('plane') or {}).get('source', 'none')}. MD execution is assessed separately.",
        "- multiple seed Redocking: "
        + ("Completed" if (redock.get("status") == "completed") else "Blocked by Vina availability; no redocking evidence produced"),
        "",
        "## 7. Reproduction",
        "",
        "```",
        f"zvss run --input <Input directory> --target {profile.source_path} --budget {ctx.config.budget.budget_id} --label <Labels>",
        "```",
        "",
        "Same configuration digest + Same input + Same code version should produce identical selection_version. ",
        "",
    ]
    return "\n".join(lines)


def _compact(value: Any) -> str:
    """Render an observed value on one CSV cell without losing meaning."""
    if value is None:
        return ""
    if isinstance(value, dict):
        return ";".join(f"{key}={_compact(val)}" for key, val in sorted(value.items()))
    if isinstance(value, (list, tuple)):
        return ";".join(_compact(item) for item in value)
    return str(value)


def _write_rows(path: str, fields: list[str], rows: list[dict[str, Any]]) -> str:
    with open(path, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path
