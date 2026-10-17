# -*- coding: utf-8 -*-
"""run_filters.py -- apply structure-derived filters to the scored scaffolds.

Unlike the property envelope, these filters ask whether the molecule does what
the peptide does, using only the three deposited complexes:

  F1 anchor specificity  anionic head must reach a validated anchor
  F2 like-charge guard   no carboxylate...carboxylate contact
  F3 sub-pocket bridge   touches both a polar (S1/S2) and a hydrophobic (S3/S4) pocket
  F4 ensemble            re-place in all three structures; scores must agree

Usage:
    python run_filters.py                 # all scored scaffolds
    python run_filters.py --targets 538   # only the triple-pass subset
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pepwedge import filters as F                      # noqa: E402
from pepwedge import pipeline as PL                    # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', default='.')
    ap.add_argument('--out', default='results')
    ap.add_argument('--targets', type=int, default=0,
                    help='0 = every scored scaffold; N = keep only the first N '
                         '(the drug_like subset is already sorted by total)')
    ap.add_argument('--max-spread', type=float, default=5.0,
                    help='max allowed score spread across the three structures')
    a = ap.parse_args()

    t0 = time.time()
    cxs = PL.load_all(a.root)
    feats, _ = PL.pharmacophore_stage(cxs, a.out, verbose=False)
    pk, _ = PL.pocket_stage(cxs, a.out, verbose=False)
    results, _ = PL.design_stage(cxs, feats, a.out, verbose=False)

    pocket = json.load(open(os.path.join(a.out, 'pocket.json'), encoding='utf-8'))
    drug_like = [s for s in results
                 if s.components.get('property', 9) <= 0.35
                 and s.components.get('peptide_overlap', 0) >= 0.5
                 and s.components.get('clash', 9) < 0.25]
    cand = drug_like if a.targets == 0 else drug_like[:a.targets]
    print(f'[input] {len(results)} scored scaffolds, {len(drug_like)} triple-pass, '
          f'evaluating {len(cand)}')

    ref = cxs['SglCJ3']
    recs = {t: F.receptor_in_frame(cx, ref) for t, cx in cxs.items()}
    base = F.evaluate(results, cxs, pocket)          # F1-F3 over everything
    by_smiles = {}
    for row in base:
        by_smiles.setdefault(row['smiles'], row)

    fd = {f.fid: f for f in feats}
    rows = []
    for sc in cand:
        r1 = by_smiles.get(sc.smiles)
        F1, F2 = fd[sc.axis[0]].pos, fd[sc.axis[1]].pos
        sub = {t: F.pocket_subset(r, (F1 + F2) / 2) for t, r in recs.items()}
        r4 = F.repose_ensemble(sc, F1, F2, sub)
        rows.append(dict(
            smiles=sc.smiles, axis=f'{sc.axis[0]}-{sc.axis[1]}',
            total=sc.components['total'], heavy_atoms=len(sc.elements),
            n_anchor=r1['n_anchor'], anchor_ids=r1['anchors'],
            n_repulsion=r1['n_repulsion'],
            n_polar=r1['n_polar'], n_hydro=r1['n_hydro'], bridge=r1['bridge'],
            ens_rmsd=round(r4['rmsd'], 2), ens_spread=round(r4['score_spread'], 2),
            ens_min_score=round(r4['min_score'], 2), ens_max_clash=round(r4['max_clash'], 2),
            **{f'score_{t}': round(v, 2) for t, v in r4['scores'].items()},
            keep=bool(r1['n_anchor'] >= 1 and r1['n_repulsion'] == 0 and r1['bridge']
                      and r4['score_spread'] <= a.max_spread)))

    path = os.path.join(a.out, 'structure_filters.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(a.out, 'structure_filters_summary.json'), 'w',
              encoding='utf-8') as fh:
        json.dump(dict(n_input=len(cand), max_spread=a.max_spread,
                       survivors=sum(1 for r in rows if r['keep']),
                       by_stage=dict(
                           anchor=sum(1 for r in rows if r['n_anchor'] >= 1),
                           anchor2=sum(1 for r in rows if r['n_anchor'] >= 2),
                           no_repulsion=sum(1 for r in rows if r['n_repulsion'] == 0),
                           bridge=sum(1 for r in rows if r['bridge']),
                           ensemble=sum(1 for r in rows if r['ens_spread'] <= a.max_spread)),
                       calibration=F.calibration_report(cxs, pocket)),
                  fh, ensure_ascii=False, indent=2)

    n = len(rows)
    print(f'[funnel] {n} -> '
          f"anchor {sum(1 for r in rows if r['n_anchor'] >= 1)} -> "
          f"bridge {sum(1 for r in rows if r['bridge'])} -> "
          f"ensemble {sum(1 for r in rows if r['ens_spread'] <= a.max_spread)} -> "
          f"keep {sum(1 for r in rows if r['keep'])}")
    print(f'[done] {path} in {time.time() - t0:.0f}s')
    for r in rows[:5]:
        print(f"   {r['total']:6.2f} spread={r['ens_spread']:5.2f} {r['smiles']}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
