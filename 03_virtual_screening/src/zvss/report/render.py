"""Report rendering: one honest, reproducible record of what was and was not done."""

from __future__ import annotations

import csv
import os
from typing import Any

from ..util import atomic_write_json, atomic_write_text, now_iso, relpath

DONE = "已实现并用真实数据验证"
IMPL = "已实现未用真实数据验证"
BLOCKED = "被外部条件阻塞"
NOTDONE = "未实现"


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
        ["rank", "compound_id", "tier", "evidence_score", "affinity", "cnn_pose", "cnn_affinity",
         "mechanism", "pose_quality", "missing"],
        [
            {
                "rank": entry.get("rank"),
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

    markdown = _markdown(ctx, state, ranking, selection, md, target, imported, qa, mech, redock, membrane)
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
        out.append({"stage": "redock", "code": "gnina_unavailable", "detail": capability.get("reason")})
    md_cap = (membrane.get("md_capability") or {})
    if md_cap.get("status"):
        out.append({"stage": "md", "code": md_cap.get("status"), "detail": md_cap.get("note")})
    if md.get("blocked"):
        out.append({"stage": "md", "code": "packages_blocked", "detail": f"{md.get('blocked')} package(s) not executable"})
    return out


def _markdown(ctx, state, ranking, selection, md, target, imported, qa, mech, redock, membrane) -> str:
    lines: list[str] = []
    profile = ctx.config.target
    lines += [
        f"# ZVSS 筛选报告 {ctx.run_id}",
        "",
        f"- 靶点: `{profile.target_id}` ({profile.display_name or profile.target_id})",
        f"- 参考结构: `{profile.primary.pdb_id}` chain `{profile.primary.chain}`",
        f"- 配置摘要: `{ctx.config.digest()}` | 代码版本 `{ctx.config.code_version}`",
        f"- 预算: `{ctx.config.budget.budget_id}` | 选择版本: `{(selection or {}).get('selection_version')}`",
        "",
        "## 1. 输入与身份",
        "",
        f"- 发现文件: {len(imported.get('files') or [])} 个",
        f"- 位姿: 读取 {imported.get('poses_seen', 0)} / 导入 {imported.get('poses_imported', 0)} / 失败 {imported.get('poses_failed', 0)}",
        f"- 化合物 {imported.get('parents', 0)} | 化学状态 {imported.get('states', 0)} | 对接作业 {imported.get('jobs', 0)}",
        f"- 受体: {(imported.get('receptor') or {}).get('status')} | 盒子来源: {(imported.get('box') or {}).get('source', '未找到')}",
    ]
    if imported.get("needs_input"):
        lines.append(f"- 需人工补充化学结构: {len(imported['needs_input'])} 条")
    if imported.get("conflicts"):
        lines.append(f"- 身份冲突: {len(imported['conflicts'])} 条（未合并，逐条保留）")
    lines.append("")

    lines += ["## 2. 校验与局部机制", ""]
    lines.append(f"- 位姿质检: 通过 {qa.get('passed', 0)} / 标记 {qa.get('flagged', 0)} / 失败 {qa.get('failed', 0)}")
    pb = qa.get("posebusters") or {}
    lines.append(f"- PoseBusters: {'可用' if pb.get('available') else '不可用'} {pb.get('reason') or ''}".rstrip())
    lines.append(f"- 接触分析: 已评估 {mech.get('evaluated', 0)} 个位姿, 截断 {mech.get('cutoff')} Å")
    ref = mech.get("reference") or {}
    lines.append(f"- 参考配体接触基线: {ref.get('status')} ({ref.get('contact_count', 'n/a')})")
    lines.append("")

    lines += ["## 3. 证据分层与排名", ""]
    tiers = ranking.get("tiers") or {}
    lines.append(f"- 分层计数: A={tiers.get('A', 0)} B={tiers.get('B', 0)} C={tiers.get('C', 0)} D={tiers.get('D', 0)}")
    weights = ranking.get("weights_used") or {}
    lines.append(f"- 使用权重: {weights}")
    skipped = ranking.get("weights_skipped") or []
    if skipped:
        lines.append(f"- 因缺少数据而跳过的权重: {skipped}")
    lines.append("")
    top = (ranking.get("ranking") or [])[:15]
    if top:
        lines += ["| 排名 | 化合物 | 层 | 综合分 | 缺失证据 |", "| --- | --- | --- | --- | --- |"]
        for entry in top:
            lines.append(
                f"| {entry.get('rank')} | `{entry['compound_id']}` | {entry['tier']} | "
                f"{entry.get('evidence_score'):.4f} | {', '.join(entry.get('components_missing') or []) or '-'} |"
            )
        lines.append("")

    lines += ["## 4. 候选与排除", ""]
    primary = (selection or {}).get("primary") or []
    preliminary = (selection or {}).get("preliminary") or []
    exploratory = (selection or {}).get("exploratory") or []
    verdicts = (selection or {}).get("qualification") or {}
    formal_count = (selection or {}).get("formal_qualified")
    lines.append(
        f"- 正式主候选: {len(primary)} | 展示性初步候选: {len(preliminary)} | "
        f"探索性候选: {len(exploratory)} | 排除: {(selection or {}).get('excluded', 0)}"
    )
    lines.append(
        f"- 满足全部正式要求的化合物: {formal_count if formal_count is not None else '未计算'} "
        f"(正式要求: {', '.join((selection or {}).get('qualification_rules') or []) or '未配置'})"
    )
    if not primary:
        lines.append(
            "- **本轮没有正式主候选**：缺少必需证据（逐项见 `qualification.csv`），"
            "规则未放宽，展示清单不作为科学资格。"
        )
    if primary or preliminary or exploratory:
        lines += [
            "",
            "| 角色 | 序号 | 化合物 | 层 | 资格 | 综合分 | 未满足要求 |",
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
        lines += ["逐候选资格明细:", "", "| 化合物 | 要求 | 状态 | 观察值 | 说明 |", "| --- | --- | --- | --- | --- |"]
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
        lines += ["排除理由（节选）:", ""]
        for item in exclusions[:10]:
            lines.append(f"- `{item['compound_id']}`: {item['reason_code']} ({item['detail']})")
    lines.append("")

    lines += ["## 5. MD 任务包", ""]
    lines.append(f"- 生成任务包: {md.get('count', 0)} | 可执行: {md.get('ready', 0)} | 阻塞: {md.get('blocked', 0)}")
    for package in (md.get("packages") or [])[:6]:
        lines.append(f"- `{os.path.basename(package.get('packages', ''))}` → {package.get('execution_status')}")
    lines.append("")

    lines += ["## 6. 明确未完成 / 阻塞", ""]
    blockers = _blockers(md, redock, membrane)
    if blockers:
        for item in blockers:
            lines.append(f"- `{item['stage']}` / `{item['code']}`: {item.get('detail')}")
    else:
        lines.append("- 无记录的阻塞项（仅代表本轮未触发）")
    lines += [
        "- 膜环境: 仅给出跨膜螺旋主轴几何代理，未构建膜体系",
        "- 多 seed 再对接: "
        + ("已完成" if (redock.get("status") == "completed") else "被 GNINA 可用性阻塞，未产生再对接证据"),
        "",
        "## 7. 复现",
        "",
        "```",
        f"zvss run --input <输入目录> --target {profile.source_path} --budget {ctx.config.budget.budget_id} --label <标签>",
        "```",
        "",
        "同一配置摘要 + 同一输入 + 同一代码版本应产生相同的 selection_version。",
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
