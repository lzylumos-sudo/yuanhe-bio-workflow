# -*- coding: utf-8 -*-
"""run_export_poses.py -- dump the 3D poses of the 538 drug-like scaffolds.

The web console needs real coordinates for the scaffolds that pass the
structure-derived filters (only the Top-25 poses were written to SDF before).
This re-runs the design stage (about 50 s) and writes

    results/drug_like_poses.json   {key: {smiles, axis, el, xyz, bonds, total}}

key = "<smiles>|<axis>".  Coordinates are in the SglCJ3 (9NU8) frame, the same
frame as /api/scene.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pepwedge import pipeline as PL                    # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', default='.')
    ap.add_argument('--out', default='results')
    a = ap.parse_args()

    cxs = PL.load_all(a.root)
    feats, _ = PL.pharmacophore_stage(cxs, a.out, verbose=False)
    PL.pocket_stage(cxs, a.out, verbose=False)
    results, _ = PL.design_stage(cxs, feats, a.out, verbose=False)
    drug_like = [s for s in results
                 if s.components.get('property', 9) <= 0.35
                 and s.components.get('peptide_overlap', 0) >= 0.5
                 and s.components.get('clash', 9) < 0.25]

    out = {}
    for sc in drug_like:
        key = f'{sc.smiles}|{sc.axis[0]}-{sc.axis[1]}'
        if key in out:
            continue
        xyz = np.asarray(sc.pose)
        bonds = [[i, j] for i in range(len(xyz)) for j in range(i + 1, len(xyz))
                 if float(np.linalg.norm(xyz[i] - xyz[j])) < 1.75]
        out[key] = dict(smiles=sc.smiles, axis=f'{sc.axis[0]}-{sc.axis[1]}',
                        el=list(sc.elements),
                        xyz=[[round(float(c), 3) for c in p] for p in xyz],
                        bonds=bonds, total=float(sc.components['total']))
    path = os.path.join(a.out, 'drug_like_poses.json')
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(out, fh, ensure_ascii=False, separators=(',', ':'))
    print(f'[done] {len(out)} poses -> {path}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
