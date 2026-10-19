# -*- coding: utf-8 -*-
"""pepwedge.pipeline -- end-to-end run: interface -> pharmacophore -> pocket
-> constrained scaffold assembly -> scoring -> ranked result set."""
from __future__ import annotations

import csv
import json
import os
import time
from typing import Dict, List

import numpy as np

from .core import (COMPLEXES, SGL_FUNCTIONAL, contacts, hotspot_table, load_complex,
                   membrane_frame, residue_sasa, salt_bridges)
from .design import enumerate_scaffolds, property_penalty
from .pharmacophore import anchor_pairs, consensus_features, superpose_complex
from .pocket import build_pocket, pocket_report, subpocket_centroids
from .pose import Receptor, place, receptor_from_complex, score_pose

# contact string used for a fast "which residues does the ligand touch" readout
POLAR = {'N', 'O'}


def _ensure(d):
    os.makedirs(d, exist_ok=True)
    return d


def load_all(root):
    cxs = {}
    for tag, f in COMPLEXES.items():
        cxs[tag] = load_complex(os.path.join(root, f), tag)
    return cxs


def interface_stage(cxs, outdir, verbose=True):
    dsasa_by_tag = {}
    conv_sets = []
    rows = []
    for tag, cx in cxs.items():
        free, bound, dsasa = residue_sasa(cx)
        dsasa_by_tag[tag] = dsasa
        ct = contacts(cx)
        conv_sets.append(set(rid for rid, _ in ct))
    for tag, cx in cxs.items():
        for row in hotspot_table(cx, dsasa=dsasa_by_tag[tag], conv_sets=conv_sets):
            row['complex'] = tag
            rows.append(row)
    # consensus table: best evidence per residue across complexes
    merged = {}
    for r in rows:
        k = r['resid']
        if k not in merged or r['score'] > merged[k]['score']:
            merged[k] = dict(r)
        merged[k]['complexes'] = sorted(set(merged[k].get('complexes', []) + [r['complex']]))
    cons = sorted(merged.values(), key=lambda r: -r['score'])
    with open(os.path.join(outdir, 'interface_hotspots.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['murj_resid', 'resname', 'score', 'convergence', 'dSASA_A2',
                    'n_contacts', 'n_peptide_residues', 'dmin_A', 'complexes', 'evidence'])
        for r in cons:
            w.writerow([r['resid'], r['one'], r['score'], r['convergence'], r['dSASA'],
                        r['n_contacts'], r['n_peptide_residues'], r['dmin'],
                        ';'.join(r['complexes']), r['evidence']])
    if verbose:
        print(f'[interface] {len(cons)} MurJ residues at the peptide interface')
        for r in cons[:8]:
            print(f'   {r["one"]}{r["resid"]:<4d} score={r["score"]:.2f}  dSASA={r["dSASA"]:.0f}  {r["evidence"]}')
    return cons


def pharmacophore_stage(cxs, outdir, verbose=True):
    feats = consensus_features(cxs, min_structures=1, max_features=12)
    pairs = anchor_pairs(feats)
    with open(os.path.join(outdir, 'pharmacophore.json'), 'w', encoding='utf-8') as fh:
        json.dump(dict(features=[f.as_dict() for f in feats],
                       pairs=[p for p in pairs if p['weight'] >= 0.45]), fh, indent=1, ensure_ascii=False)
    if verbose:
        print(f'[pharmacophore] {len(feats)} consensus features')
        for f in feats:
            print(f'   {f.fid:<4s} {f.ftype:<10s} conv={f.convergence:.2f} w={f.weight:.2f} '
                  f'-> {" ".join(f.murj_contacts[:4])}')
    return feats, pairs


def pocket_stage(cxs, outdir, verbose=True):
    cx = cxs['SglCJ3']
    pk = build_pocket(cx)
    keep = pk.buriedness >= 0.5
    pk.points = pk.points[keep]
    pk.buriedness = pk.buriedness[keep]
    pk.chemistry = [c for c, k in zip(pk.chemistry, keep) if k]
    pk.subpocket = [s for s, k in zip(pk.subpocket, keep) if k]
    rep = pocket_report(pk, cx)
    with open(os.path.join(outdir, 'pocket.json'), 'w', encoding='utf-8') as fh:
        json.dump(dict(total_volume_A3=round(float(pk.volume), 1),
                       subpockets=rep,
                       centroids={k: [round(float(x), 2) for x in v]
                                  for k, v in subpocket_centroids(pk).items()}), fh,
                  indent=1, ensure_ascii=False)
    if verbose:
        print(f'[pocket] {pk.volume:.0f} A^3 buried cavity, {len(rep)} sub-pockets')
        for r in rep:
            print(f'   {r["subpocket"]:<22s} {r["volume_A3"]:>7.0f} A^3  buriedness={r["mean_buriedness"]}')
    return pk, rep


def peptide_atoms_in_frame(cxs, ref='SglCJ3'):
    ref_cx = cxs[ref]
    out = []
    for tag, cx in cxs.items():
        if tag == ref:
            out.append(cx.l_xyz)
            continue
        shared = [r for r in ref_cx.receptor
                  if 'CA' in ref_cx.receptor[r].atoms and r in cx.receptor
                  and 'CA' in cx.receptor[r].atoms]
        mob = np.array([cx.receptor[r].atoms['CA'] for r in shared])
        tgt = np.array([ref_cx.receptor[r].atoms['CA'] for r in shared])
        from .core import kabsch
        rt = kabsch(mob, tgt)
        from .core import apply_rt
        out.append(apply_rt(cx.l_xyz, rt))
    return np.vstack(out)


def peptide_overlap(sc, pep):
    d = np.linalg.norm(sc.pose[:, None, :] - pep[None, :, :], axis=2)
    return float((d.min(1) < 2.3).mean())


def pick_axes(feats, dmin=9.0, dmax=18.0, min_weight=0.45):
    by_id = {f.fid: f for f in feats}
    axes = []
    for a in feats:
        if a.ftype != 'anion':
            continue
        for b in feats:
            if b.fid == a.fid or b.ftype == 'cation':
                continue
            d = float(np.linalg.norm(a.pos - b.pos))
            if not (dmin <= d <= dmax):
                continue
            w = a.weight * b.weight
            if w < min_weight:
                continue
            axes.append(dict(F1=a, F2=b, distance=d, weight=w))
    axes.sort(key=lambda x: -x['weight'])
    # keep the best axis per (F1, F2) unordered pair
    seen, keep = set(), []
    for ax in axes:
        k = tuple(sorted((ax['F1'].fid, ax['F2'].fid)))
        if k in seen:
            continue
        seen.add(k)
        keep.append(ax)
    return keep


def design_stage(cxs, feats, outdir, max_per_axis=2500, verbose=True):
    cx = cxs['SglCJ3']
    rec_full = receptor_from_complex(cx)
    pep = peptide_atoms_in_frame(cxs)
    axes = pick_axes(feats)
    results: List = []
    axis_rows = []
    for ai, ax in enumerate(axes):
        F1, F2 = ax['F1'], ax['F2']
        D = float(np.linalg.norm(F2.pos - F1.pos))
        cands = enumerate_scaffolds(D, tol=0.55, ring_counts=(2, 3), max_mw=500)
        # property pre-filter
        pre = []
        for sc in cands:
            pen = property_penalty(sc)
            if pen > 0.75:
                continue
            pre.append((pen, sc))
        pre.sort(key=lambda t: t[0])
        pre = pre[:max_per_axis]
        rec = rec_full.subset((F1.pos + F2.pos) / 2.0, r=15.0)
        t0 = time.time()
        for pen, sc in pre:
            place(sc, F1.pos, F2.pos, rec, n_rot=12)
            sc.axis = (F1.fid, F2.fid)
            sc.axis_length = round(D, 2)
            sc.components['pharm'] = round(0.5 * (F1.weight + F2.weight), 3)
            sc.components['property'] = round(pen, 3)
            rigidity = min(1.0, D / (3.0 * (1 + sc.rot)))
            overlap = peptide_overlap(sc, pep)
            total = (sc.pose_score - 1.5 * pen + 0.6 * rigidity
                     + 0.4 * (F1.weight + F2.weight) + 1.6 * overlap)
            sc.components['peptide_overlap'] = round(overlap, 2)
            sc.components['total'] = round(total, 3)
            sc.components['rigidity'] = round(rigidity, 2)
            results.append(sc)
        axis_rows.append(dict(axis=f'{F1.fid}-{F2.fid}', types=f'{F1.ftype}-{F2.ftype}',
                              distance=round(D, 2), weight=round(ax['weight'], 2),
                              candidates=len(cands), scored=len(pre),
                              seconds=round(time.time() - t0, 1)))
        if verbose:
            print(f'[design] axis {F1.fid}({F1.ftype}) -> {F2.fid}({F2.ftype})  D={D:.1f} A  '
                  f'{len(cands)} skeletons -> {len(pre)} scored  ({time.time()-t0:.1f}s)')
    results.sort(key=lambda s: -s.components['total'])
    return results, axis_rows


def annotate_contacts(sc, rec: 'Receptor', cx):
    """Which MurJ residues the top pose touches (polar + hydrophobic)."""
    d = np.linalg.norm(sc.pose[:, None, :] - rec.xyz[None, :, :], axis=2)
    ids = set()
    for i in range(len(sc.pose)):
        near = np.where(d[i] < 4.2)[0]
        for j in near:
            ids.add(int(rec.resid[j]))
    return sorted(ids)


def write_results(results, axes, outdir, cx, top=200, rec=None):
    if rec is None:
        rec = receptor_from_complex(cx)
    rows = []
    for i, sc in enumerate(results[:top]):
        ids = annotate_contacts(sc, rec, cx)
        contacts = [f'{cx.r(j).one}{j}' for j in ids]
        rows.append(dict(
            rank=i + 1, smiles=sc.smiles, axis=f'{sc.axis[0]}-{sc.axis[1]}',
            axis_length_A=sc.axis_length,
            family=f'{sc.head}|{sc.n_ring}ring|{sc.tail}',
            head=sc.head, rings='/'.join(sc.rings), links='/'.join(sc.links), tail=sc.tail,
            MW=round(sc.mw, 1), cLogP=round(sc.clogp, 2), TPSA=round(sc.tpsa, 1),
            HBD=sc.hbd, HBA=sc.hba, RotB=sc.rot, Rings=sc.n_ring,
            heavy_atoms=len(sc.elements),
            pose_score=sc.pose_score, total=sc.components['total'],
            clash=sc.components.get('clash'), salt=sc.components.get('salt'),
            hb=sc.components.get('hb'), hyd=sc.components.get('hyd'),
            peptide_overlap=sc.components.get('peptide_overlap'),
            buried=sc.components.get('buried'), property_penalty=sc.components.get('property'),
            rigidity=sc.components.get('rigidity'),
            murj_contacts=' '.join(contacts)))
    with open(os.path.join(outdir, 'scaffolds_ranked.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    write_sdf(results, os.path.join(outdir, 'top_scaffold_poses.sdf'), n=25)
    return rows


def write_sdf(results, path, n=25):
    with open(path, 'w', encoding='utf-8') as fh:
        for sc in results[:n]:
            xyz = sc.pose
            fh.write(f'{sc.smiles}\n  pepwedge  {sc.axis[0]}-{sc.axis[1]} score={sc.components["total"]}\n')
            fh.write(f'{len(xyz):3d}{0:3d}  0  0  0  0  0  0  0  0  0999 V2000\n')
            for a, x in zip(sc.elements, xyz):
                fh.write(f'{x[0]:10.4f}{x[1]:10.4f}{x[2]:10.4f} {a:<3s} 0  0  0  0  0  0  0  0  0  0  0  0\n')
            bonds = []
            for a in range(len(xyz)):
                for b in range(a + 1, len(xyz)):
                    if float(np.linalg.norm(xyz[a] - xyz[b])) < 1.75:
                        bonds.append((a, b))
            for a, b in bonds:
                fh.write(f'{a+1:3d}{b+1:3d}  1  0  0  0  0\n')
            fh.write('M  END\n$$$$\n')


def select_diverse(results, per_family=3, n=40):
    """Stratified pick: at most `per_family` skeletons of the same composition."""
    out, count = [], {}
    for sc in results:
        fam = f'{sc.head}|{sc.n_ring}|{sc.tail}'
        if count.get(fam, 0) >= per_family:
            continue
        count[fam] = count.get(fam, 0) + 1
        out.append(sc)
        if len(out) >= n:
            break
    return out


def run(root, outdir='results', verbose=True, figdir='figures'):
    t0 = time.time()
    _ensure(outdir)
    cxs = load_all(root)
    if verbose:
        for tag, cx in cxs.items():
            print(f'[load] {tag}: {len(cx.receptor)} MurJ residues, {len(cx.ligand)} peptide residues')
    hotspots = interface_stage(cxs, outdir, verbose)
    feats, pairs = pharmacophore_stage(cxs, outdir, verbose)
    pk, rep = pocket_stage(cxs, outdir, verbose)
    results, axis_rows = design_stage(cxs, feats, outdir, verbose=verbose)
    rec = receptor_from_complex(cxs['SglCJ3'])
    rows = write_results(results, axis_rows, outdir, cxs['SglCJ3'], rec=rec)
    diverse = select_diverse(results, per_family=2, n=24)
    write_sdf(diverse, os.path.join(outdir, 'diverse_scaffolds.sdf'), n=24)
    stats = dict(
        enumerated=int(sum(a['candidates'] for a in axis_rows)),
        property_pass=int(sum(a['scored'] for a in axis_rows)),
        posed=int(len(results)),
        no_clash=int(sum(1 for s in results if s.components.get('clash', 1) < 0.25)),
        double_anchor=int(sum(1 for s in results if s.components.get('salt', 0) >= 2)),
        peptide_mimic=int(sum(1 for s in results if s.components.get('peptide_overlap', 0) >= 0.5)),
        drug_like=int(sum(1 for s in results
                          if s.components.get('property', 9) <= 0.35
                          and s.components.get('peptide_overlap', 0) >= 0.5
                          and s.components.get('clash', 9) < 0.25)),
    )
    summary = dict(funnel=stats,
                   hotspots=[dict(resid=h['resid'], resname=h['one'], score=h['score'],
                                   convergence=h['convergence'], dSASA=h['dSASA'],
                                   evidence=h['evidence']) for h in hotspots[:20]],
                   features=[f.as_dict() for f in feats],
                   axes=axis_rows,
                   n_scaffolds=len(results),
                   best=[dict(smiles=r['smiles'], total=r['total'], axis=r['axis'],
                              props={k: r[k] for k in ('MW', 'cLogP', 'TPSA', 'HBD',
                                                       'HBA', 'RotB', 'Rings',
                                                       'heavy_atoms')},
                              contacts=r['murj_contacts'].split())
                         for r in rows[:10]])
    with open(os.path.join(outdir, 'summary.json'), 'w', encoding='utf-8') as fh:
        json.dump(summary, fh, indent=1, ensure_ascii=False)
    if verbose:
        print('[funnel]', stats)
    with open(os.path.join(outdir, 'axes.json'), 'w', encoding='utf-8') as fh:
        json.dump(axis_rows, fh, indent=1)
    # ---- figures
    try:
        from .figures import fig_funnel, fig_hotspots, fig_pocket_map, fig_scaffold
        _ensure(figdir)
        fig_hotspots(hotspots, os.path.join(figdir, 'fig1_hotspots.png'))
        fig_pocket_map(pk, cxs, feats, os.path.join(figdir, 'fig2_pocket_map.png'))
        fig_scaffold(diverse[0], cxs, rec, os.path.join(figdir, 'fig3_top_scaffold.png'),
                     label=diverse[0].smiles)
        fig_funnel(stats, os.path.join(figdir, 'fig4_funnel.png'))
        if verbose:
            print(f'[figures] written to {figdir}')
    except Exception as exc:                      # figures are optional
        print('[figures] skipped:', exc)

    if verbose:
        print(f'[done] {len(results)} scored scaffolds in {time.time()-t0:.1f}s')
        print('top 10:')
        for r in rows[:10]:
            print(f'   {r["rank"]:>2d} {r["total"]:>6.2f}  {r["smiles"]}  '
                  f'MW={r["MW"]} cLogP={r["cLogP"]} TPSA={r["TPSA"]}')
    return dict(complexes=cxs, hotspots=hotspots, features=feats, pocket=pk,
                pocket_report=rep, results=results, rows=rows, axes=axis_rows)
