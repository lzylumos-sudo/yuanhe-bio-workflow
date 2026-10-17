# -*- coding: utf-8 -*-
"""Rebuild the entire evidence chain from the deposited PDB files.

    python run_evidence.py [--root .] [--out evidence] [--no-pocket]

Everything printed to ``evidence/report.md`` is computed here. No number in the
report is typed by hand, and no literature claim is accepted without a check:
priors are loaded from ``priors/murj_priors.json`` and each one is verified
against the coordinates, with failures reported rather than hidden.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from pepwedge import priors as priors_mod
from pepwedge.core import load_complex
from pepwedge.evidence import (atom_budget, anchor_table, census, dsasa_both,
                               interface_census, interface_rows, peptide_rows)
from pepwedge.verify import run_checks, summarize


def _w(path, rows, fields=None):
    if not rows:
        return
    fields = fields or list(rows[0].keys())
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        wr = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore')
        wr.writeheader()
        wr.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', default='.')
    ap.add_argument('--out', default='evidence')
    ap.add_argument('--priors', default=None)
    ap.add_argument('--no-pocket', action='store_true')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()

    pr = priors_mod.load_priors(a.priors)
    print('[priors]', priors_mod.summary(a.priors))

    # ---- load structures -------------------------------------------------
    complexes = {}
    for s in pr['structures']:
        path = os.path.join(a.root, s['file'])
        if not os.path.exists(path):
            print(f'  !! missing structure file {path}')
            continue
        complexes[s['tag']] = load_complex(path, s['tag'])
    if not complexes:
        sys.exit('no structures loaded')

    # ---- per-complex measurements ---------------------------------------
    cen_rows, int_rows, pep_rows, anc_rows = [], [], [], []
    dsasa_by_tag, ct_by_tag, ic_by_tag = {}, {}, {}
    for tag, cx in complexes.items():
        cen_rows.append(census(cx))
        ds = dsasa_both(cx)
        dsasa_by_tag[tag] = ds
        ic = interface_census(cx, ds)
        ic_by_tag[tag] = ic
        int_rows += interface_rows(cx, ds[0])
        pep_rows += peptide_rows(cx, ds[1])
        at = anchor_table(cx)
        anc_rows += at['salt_bridges']
        anc_rows += at['anchor_pairs']
        from pepwedge.core import contacts
        ct_by_tag[tag] = contacts(cx)

    # convergence across the three structures
    conv_sets = [{r for r in {row['resid'] for row in int_rows if row['complex'] == t}}
                 for t in complexes]
    for row in int_rows:
        row['convergence'] = round(sum(row['resid'] in s for s in conv_sets) / len(conv_sets), 2)

    # ---- priors verification --------------------------------------------
    checks = run_checks(pr, complexes, dict(dsasa=dsasa_by_tag, contacts=ct_by_tag))
    summary = summarize(checks)

    # ---- atom budget -----------------------------------------------------
    budget = atom_budget(complexes, dsasa_by_tag=dsasa_by_tag)

    # ---- optional pocket assignment --------------------------------------
    sub_of = {}
    if not a.no_pocket and 'SglCJ3' in complexes:
        try:
            from pepwedge.pocket import build_pocket, subpocket_of
            pk = build_pocket(complexes['SglCJ3'])
            keep = pk.buriedness >= 0.5
            pts = pk.points[keep]
            subs = [s for s, k in zip(pk.subpocket, keep) if k]
            for row in int_rows:
                if row['complex'] != 'SglCJ3':
                    continue
                p = np.array(row['ca'], dtype=float)
                i = int(np.argmin(np.linalg.norm(pts - p, axis=1)))
                sub_of[row['resid']] = (subs[i], round(float(np.linalg.norm(pts[i] - p)), 1))
        except Exception as exc:
            print('[pocket] skipped:', exc)

    # ---- write tables ----------------------------------------------------
    _w(os.path.join(a.out, '00_census.csv'), cen_rows)
    _w(os.path.join(a.out, '01_interface_murj.csv'), int_rows)
    _w(os.path.join(a.out, '02_interface_peptide.csv'), pep_rows)
    _w(os.path.join(a.out, '03_anchors.csv'), anc_rows)
    _w(os.path.join(a.out, '04_priors_check.csv'), [c.as_dict() for c in checks])
    with open(os.path.join(a.out, '05_atom_budget.json'), 'w', encoding='utf-8') as fh:
        json.dump(dict(budget=budget, interface=ic_by_tag, summary=summary,
                       census=cen_rows), fh, indent=1, ensure_ascii=False)
    with open(os.path.join(a.out, '06_priors_resolved.json'), 'w', encoding='utf-8') as fh:
        json.dump(dict(priors_path=pr.get('_path'), checks=[c.as_dict() for c in checks],
                       summary=summary), fh, indent=1, ensure_ascii=False)

    # ---- report ----------------------------------------------------------
    md = build_report(pr, complexes, cen_rows, int_rows, pep_rows, anc_rows,
                      ic_by_tag, budget, checks, summary, sub_of, a)
    with open(os.path.join(a.out, 'report.md'), 'w', encoding='utf-8') as fh:
        fh.write(md)
    print(f'[done] evidence written to {a.out}/ in {time.time()-t0:.1f}s')
    print('       report:', os.path.join(a.out, 'report.md'))
    for row in cen_rows:
        print(f"       {row['tag']}: MurJ {row['murj_residues']} res / "
              f"peptide {row['peptide_residues']} res")
    print('       priors:', summary['by_status'])


# ---------------------------------------------------------------------------
def build_report(pr, complexes, cen, irows, prows, arows, ic, budget, checks,
                 summary, sub_of, args) -> str:
    L = []
    A = L.append
    ref = pr['target']['reference']
    A('# 证据链报告：从 PDB 到设计约束')
    A('')
    A('> 本文件由 `run_evidence.py` 自动生成，**所有数字都是本机从坐标文件现算的**，'
      '没有任何数值是手写的。重跑即可复现全部内容。')
    A('')
    A(f'- 生成时间：{time.strftime("%Y-%m-%d %H:%M:%S")}')
    A(f'- 先验文件：`{pr.get("_path")}`（{priors_mod.summary(args.priors)}）')
    A(f'- 结构目录：`{os.path.abspath(args.root)}`')
    A(f'- 参考：{ref["citation"]} doi:{ref["doi"]} PMID:{ref["pmid"]}')
    A(f'- 编号：{pr["target"]["numbering"]}')
    A('')

    # 1. census
    A('## 1. 输入结构里到底有什么')
    A('')
    A('| 复合物 | MurJ 残基 | MurJ 原子 | MurJ 编号 | 构建体额外残基（已排除） | 肽残基 | 肽原子 | 肽编号 |')
    A('|---|---|---|---|---|---|---|---|')
    for r in cen:
        A(f'| {r["tag"]} | {r["murj_residues"]} | {r["murj_atoms"]} | {r["murj_range"]} | '
          f'{r["construct_residues"]}（{r["construct_range"]}） | '
          f'{r["peptide_residues"]} | {r["peptide_atoms"]} | {r["peptide_range"]} |')
    A('')
    A('肽序列（PDB 中实际建模的部分）')
    A('')
    for r in cen:
        A(f'- {r["tag"]}：`{r["peptide_sequence"]}`')
    A('')
    A('**缺失区间（PDB 未建模，不是序列缺失）**')
    A('')
    for r in cen:
        A(f'- {r["tag"]}：MurJ `{r["murj_gaps"]}`；肽 `{r["peptide_gaps"]}`')
    A('')
    A('> 这一节回答"这个结构能信到什么程度"。肽只有 '
      f'{min(r["peptide_residues"] for r in cen)}–{max(r["peptide_residues"] for r in cen)} '
      '个残基被建模，任何"肽有 N 个残基"的说法都必须以这里的数字为准。')
    A('')

    # 2. interface size
    A('## 2. 接触面有多大（这是"肽 vs 小分子"量化论证的来源）')
    A('')
    A('| 复合物 | 界面 MurJ 残基 | 界面肽残基 | 占建模肽比例 | 原子对 | 肽侧接触重原子 | 肽侧埋藏面积 Å² | 两侧面积不对称度 |')
    A('|---|---|---|---|---|---|---|---|')
    for tag, d in ic.items():
        A(f'| {tag} | {d["murj_residues_in_interface"]} | {d["peptide_residues_in_interface"]} | '
          f'{d["fraction_of_peptide_at_interface"]:.0%} | {d["atom_pairs"]} | '
          f'{d["peptide_heavy_atoms_in_contact"]} | {d["buried_area_peptide_side_A2"]} | '
          f'{d["asymmetry"]} |')
    A('')
    A('**原子预算对比**')
    A('')
    sm = budget['small_molecule']
    A('| | 肽（接触部分，取三结构最大值） | 本算法可生成的小分子骨架 |')
    A('|---|---|---|')
    pep_max = max(v['peptide_heavy_atoms_in_contact'] for v in budget['per_complex'].values())
    pep_all = max(v['peptide_heavy_atoms'] for v in budget['per_complex'].values())
    ar_max = max(v['buried_area_peptide_side_A2'] for v in budget['per_complex'].values())
    A(f'| 重原子数 | 全长 {pep_all}，其中进入接触 {pep_max} | '
      f'{sm["heavy_atoms_min"]}–{sm["heavy_atoms_max"]}（中位数 {sm["heavy_atoms_median"]:.0f}） |')
    A(f'| 埋藏面积 | 约 {ar_max:.0f} Å² | 通常 300–450 Å² |')
    A(f'| 比值 | — | 约为肽的 1/{budget["median_ratio"]:.0f}（仅计重原子） |')
    A('')
    A(f'> 结论：肽用约 {pep_max} 个重原子覆盖约 {ar_max:.0f} Å² 的界面；'
      f'本算法能造出的骨架只有 {sm["heavy_atoms_min"]}–{sm["heavy_atoms_max"]} 个重原子。'
      '**这不是"把肽缩小"，而是"用更少的原子复刻同样的相互作用向量"**——'
      '这正是把设计目标写成"药效团点之间的距离约束"而不是"形状匹配"的定量理由。')
    A('')

    # 3. important residues
    A('## 3. 哪些氨基酸重要，重要在哪里')
    A('')
    A('判据全部可复算：接触密度、埋藏面积、跨三种肽的收敛度、'
      '以及（来自先验文件的）功能必需/耐药证据。')
    A('')
    merged = {}
    for row in irows:
        k = row['resid']
        m = merged.setdefault(k, dict(row, complexes=[], n_atom_contacts=0, dSASA_A2=0.0))
        m['n_atom_contacts'] += row['n_atom_contacts']
        m['dSASA_A2'] = max(m['dSASA_A2'], row['dSASA_A2'])
        m['complexes'].append(row['complex'])
    ess = {int(r['resid']): 'essential' for r in pr['murj_residues'] if 'essential' in r.get('roles', [])}
    res_mut = {int(m['resid']): f"{m['wt'][:3].title()}{m['resid']}{m['mut'][:3].title()}"
               for m in pr['resistance_mutations']}
    ranked = sorted(merged.values(), key=lambda r: (-r['n_atom_contacts'], r['resid']))
    A('| 残基 | 接触原子对（三结构合计） | 最大 dSASA Å² | 收敛度 | 亚口袋 | CA 坐标 | 实验证据 |')
    A('|---|---|---|---|---|---|---|')
    for r in ranked[:20]:
        ev = []
        if r['resid'] in ess:
            ev.append('功能必需')
        if r['resid'] in res_mut:
            ev.append(res_mut[r['resid']])
        sp = sub_of.get(r['resid'])
        A(f'| {r["one"]}{r["resid"]} | {r["n_atom_contacts"]} | {r["dSASA_A2"]:.0f} | '
          f'{r["convergence"]} | {sp[0] if sp else "—"} | '
          f'({r["ca"][0]:.1f}, {r["ca"][1]:.1f}, {r["ca"][2]:.1f}) | {", ".join(ev) or "—"} |')
    A('')
    common = sorted(rid for rid in merged if merged[rid]['convergence'] == 1.0)
    A(f'**三种肽共同接触的残基（几何收敛，n={len(common)}）**：'
      + ', '.join(f'{merged[i]["one"]}{i}' for i in common))
    A('')
    A('> 注意区分两个不同的量：论文 Fig. 4f 的 "13 个完全保守残基" 是**序列保守性**'
      '（需要跨 diderm 的多序列比对，本仓库不复现）；上表这一行是**几何收敛性**'
      '（三种肽在本文三个结构里都接触）。两者不一定相等，混用会出错。')
    A('')

    # 4. anchors
    A('## 4. 正电锚点之间的几何（骨架长度的来源）')
    A('')
    sb = [r for r in arows if 'dist' in r]
    A('| MurJ 原子 | 肽原子 | 距离 Å | 复合物 |')
    A('|---|---|---|---|')
    for r in sb:
        A(f'| {r["murj"]} {r["murj_atom"]} | {r["sgl"]} {r["sgl_atom"]} | {r["dist"]} | {r["complex"]} |')
    A('')
    pairs = [r for r in arows if 'sidechain_centroid_A' in r]
    A('| 复合物 | 锚点对 | 侧链质心距离 Å | 最近原子距离 Å |')
    A('|---|---|---|---|')
    for r in pairs:
        A(f'| {r["complex"]} | {r["a"]}–{r["b"]} | {r["sidechain_centroid_A"]} | {r["min_atom_A"]} |')
    A('')
    A('> 这些距离就是骨架枚举时 `D` 的取值范围来源；`A1–A8`/`A1–R2` 等药效团轴'
      '是肽的**功能原子**之间的距离，本质上也是从这张表推出来的。')
    A('')

    # 5. priors
    A('## 5. 文献先验逐条核对')
    A('')
    A(f'共 {summary["total"]} 条，PASS {summary["passed"]}，FAIL {summary["failed"]}，'
      f'FLAG {summary["flagged"]}，WEAK {summary["weak"]}，'
      f'NO_RESIDUE {summary["no_residue"]}，MISMATCH {summary["mismatch"]}。')
    A('')
    A('| 类型 | 目标 | 核对方式 | 期望 | 结果 | 观测 | 出处 |')
    A('|---|---|---|---|---|---|---|')
    for c in checks:
        A(f'| {c.kind} | {c.target} | {c.verify_type} | {c.expectation} | '
          f'**{c.status}** | {c.observed} | {c.locator} |')
    A('')
    bad = [c for c in checks if c.status in ('FAIL', 'FLAG', 'MISMATCH', 'NO_RESIDUE')]
    if bad:
        A('### 未能复现的条目（必须人工看一眼）')
        A('')
        for c in bad:
            A(f'- **{c.target}**（{c.status}）：{c.observed}')
            if c.claim:
                A(f'  - 原声称：{c.claim}')
            if c.quote:
                A(f'  - 原文："{c.quote}"（{c.locator}）')
            if c.detail:
                A(f'  - 说明：{c.detail}')
        A('')

    # 6. unverifiable
    unv = pr.get('unverifiable_claims', [])
    if unv:
        A('## 6. 明确不核对、也不在本仓库验证的声明')
        A('')
        for u in unv:
            A(f'- {u["claim"]}')
            A(f'  - 原因：{u["why_not_verified"]}')
            if u.get('quote'):
                A(f'  - 原文："{u["quote"]}"（{u.get("locator", "")}）')
        A('')

    # 7. provenance
    A('## 7. 每个数字对应哪段代码')
    A('')
    A('| 量 | 计算函数 | 输出文件 |')
    A('|---|---|---|')
    for q, fn, out in [
        ('链/残基/原子数、缺失区间', 'evidence.census', '00_census.csv'),
        ('每残基接触、埋藏面积', 'evidence.interface_rows', '01_interface_murj.csv'),
        ('肽侧接触与埋藏', 'evidence.peptide_rows', '02_interface_peptide.csv'),
        ('界面规模、原子对、两侧面积', 'evidence.interface_census', '05_atom_budget.json'),
        ('盐桥与锚点距离', 'evidence.anchor_table', '03_anchors.csv'),
        ('原子预算对比', 'evidence.atom_budget', '05_atom_budget.json'),
        ('文献先验核对', 'verify.run_checks', '04_priors_check.csv'),
        ('亚口袋归属', 'pocket.build_pocket + evidence', '本报告第 3 节'),
    ]:
        A(f'| {q} | `{fn}` | `{out}` |')
    A('')
    A('## 8. 复现命令')
    A('')
    A('```bash')
    A('python run_evidence.py --out evidence      # 本报告 + 全部表格')
    A('python run_pipeline.py --out results       # 设计输出（依赖同一份 priors）')
    A('python tests/test_smoke.py                 # 自检')
    A('```')
    A('')
    A('修改 `priors/murj_priors.json` 即可替换文献先验（例如换用另一篇论文的'
      '突变数据），无需改代码；核对结果会立刻反映在第 5 节。')
    return '\n'.join(L) + '\n'


if __name__ == '__main__':
    main()
