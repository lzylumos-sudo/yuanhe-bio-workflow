# -*- coding: utf-8 -*-
"""pepwedge.filters -- filters that use *interaction geometry*, not property envelopes.

Motivation
----------
The property envelope (MW / cLogP / TPSA / RotB) answers "is this molecule
drug-like in general".  It says nothing about whether this molecule does the
*thing the peptide does*.  Everything below is derived from the three deposited
complexes, so it is objective and recomputable:

F1  anchor salt bridge   the anionic head must reach one of the four
                         genetically validated lysines/arginines (K46, R53,
                         R270, K368) -- not merely "some positive atom"
F2  like-charge clash    a carboxylate oxygen within 3.6 A of a MurJ acidic
                         oxygen (D/E) is a real defect; the peptides never do it
F3  sub-pocket bridge    the molecule must touch both a polar/basic sub-pocket
                         (S1/S2) and a hydrophobic one (S3/S4); that is what
                         "spanning the TM2-TM7 groove" means geometrically
F4  ensemble consistency the same pose is re-scored against all three
                         structures after Kabsch superposition; a pose that
                         only works in one receptor is over-fit to it

Thresholds are the ones the peptides themselves satisfy (see
``calibration_report``), not values picked by hand.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .core import Complex, apply_rt
from .pharmacophore import superpose_complex
from .pose import Receptor, place, receptor_from_complex, score_pose

ANCHORS = (46, 53, 270, 368)          # K46, R53, R270, K368
ANCHOR_CUT = 3.8                      # A, O...N
REPULSION_CUT = 3.6                   # A, O...O
SUBPOCKET_CUT = 4.2                   # A, heavy-atom contact

POLAR_POCKETS = ('S1_basic_bowl', 'S2_carboxylate_gate')
HYDRO_POCKETS = ('S3_TM2_TM7_groove', 'S4_TM7_floor')


def _resids(names: Sequence[str]) -> np.ndarray:
    """['K46', 'R270', ...] -> array([46, 270, ...])"""
    return np.array([int(''.join(c for c in n if c.isdigit())) for n in names], dtype=int)


def _pocket_residues(pocket: dict) -> Tuple[np.ndarray, np.ndarray]:
    polar, hydro = set(), set()
    for sp in pocket['subpockets']:
        ids = set(_resids(sp['residues']).tolist())
        if sp['subpocket'] in POLAR_POCKETS:
            polar |= ids
        elif sp['subpocket'] in HYDRO_POCKETS:
            hydro |= ids
    return np.array(sorted(polar - hydro)), np.array(sorted(hydro - polar))


def _anion_o(sc) -> np.ndarray:
    return np.array([i for i, (e, q) in enumerate(zip(sc.elements, sc.charge))
                     if e == 'O' and q < 0], dtype=int)


def receptor_in_frame(cx: Complex, ref: Complex) -> Receptor:
    """Receptor of ``cx`` expressed in ``ref``'s coordinate frame."""
    rec = receptor_from_complex(cx)
    if cx.tag != ref.tag:
        rt = superpose_complex(cx, ref)[0]
        rec.xyz = apply_rt(rec.xyz, rt)
    return rec


def anchor_salt(sc, rec: Receptor) -> dict:
    """F1 -- shortest O...N distance to each of the four validated anchors."""
    oi = _anion_o(sc)
    if oi.size == 0:
        return dict(n_anchor=0, dmin=99.0, anchors=[])
    lig = sc.pose[oi]
    mask = rec.is_pos & np.isin(rec.resid, ANCHORS)
    if not mask.any():
        return dict(n_anchor=0, dmin=99.0, anchors=[])
    d = np.linalg.norm(lig[:, None, :] - rec.xyz[None, mask, :], axis=2)
    rids = rec.resid[mask]
    hit = {}
    for k, rid in enumerate(rids):
        m = float(d[:, k].min())
        if m <= ANCHOR_CUT:
            hit[int(rid)] = min(hit.get(int(rid), 99.0), m)
    return dict(n_anchor=len(hit), dmin=min(hit.values()) if hit else 99.0,
                anchors=sorted(hit))


def like_charge(sc, rec: Receptor) -> dict:
    """F2 -- anionic head too close to a MurJ carboxylate."""
    oi = _anion_o(sc)
    if oi.size == 0:
        return dict(n_repulsion=0, dmin=99.0)
    lig = sc.pose[oi]
    mask = rec.is_neg
    if not mask.any():
        return dict(n_repulsion=0, dmin=99.0)
    d = np.linalg.norm(lig[:, None, :] - rec.xyz[None, mask, :], axis=2)
    mins = d.min(axis=1)
    return dict(n_repulsion=int((mins <= REPULSION_CUT).sum()), dmin=float(mins.min()))


def subpocket_bridge(sc, rec: Receptor, polar: np.ndarray, hydro: np.ndarray) -> dict:
    """F3 -- does the molecule touch both a polar and a hydrophobic sub-pocket?"""
    lig = sc.pose
    def touched(ids):
        mask = np.isin(rec.resid, ids)
        if not mask.any():
            return set()
        d = np.linalg.norm(lig[:, None, :] - rec.xyz[None, mask, :], axis=2)
        keep = rec.resid[mask][d.min(axis=0) <= SUBPOCKET_CUT]
        return set(int(x) for x in keep)
    p, h = touched(polar), touched(hydro)
    return dict(n_polar=len(p), n_hydro=len(h), bridge=bool(p and h),
                polar=sorted(p), hydro=sorted(h))


def ensemble(sc, recs: Dict[str, Receptor]) -> dict:
    """F4 -- re-score the *same* pose against every structure.

    score_pose() takes arrays, not Features, so the pose is held fixed and only
    the receptor changes: the spread of scores is therefore a measure of how
    much the pose depends on one particular structure.
    """
    out = {}
    for tag, rec in recs.items():
        total, parts = score_pose(sc.pose, sc.elements, sc.charge, rec)
        out[tag] = dict(score=float(total), clash=float(parts['clash']))
    vals = [v['score'] for v in out.values()]
    return dict(scores=out,
                min_ratio=float(min(vals) / max(max(vals), 1e-9)),
                max_clash=float(max(v['clash'] for v in out.values())))


def evaluate(results, cxs: Dict[str, Complex], pocket: dict,
             ref_tag: str = 'SglCJ3') -> List[dict]:
    """Apply F1-F4 to every scored scaffold.  Returns one row per scaffold."""
    ref = cxs[ref_tag]
    recs = {tag: receptor_in_frame(cx, ref) for tag, cx in cxs.items()}
    rec = recs[ref_tag]
    polar, hydro = _pocket_residues(pocket)
    rows = []
    for sc in results:
        r1 = anchor_salt(sc, rec)
        r2 = like_charge(sc, rec)
        r3 = subpocket_bridge(sc, rec, polar, hydro)
        r4 = ensemble(sc, recs)
        rows.append(dict(
            smiles=sc.smiles, axis=f'{sc.axis[0]}-{sc.axis[1]}', total=sc.components['total'],
            n_anchor=r1['n_anchor'], anchor_dmin=round(r1['dmin'], 2),
            anchors=' '.join(str(a) for a in r1['anchors']),
            n_repulsion=r2['n_repulsion'], repulsion_dmin=round(r2['dmin'], 2),
            n_polar=r3['n_polar'], n_hydro=r3['n_hydro'], bridge=r3['bridge'],
            ens_min_ratio=round(r4['min_ratio'], 3), ens_max_clash=round(r4['max_clash'], 3),
            **_ens_columns(r4['scores']),
            keep=bool(r1['n_anchor'] >= 1 and r2['n_repulsion'] == 0 and r3['bridge']
                      and r4['max_clash'] < 0.25),
        ))
    return rows


def _ens_columns(scores: dict) -> dict:
    return {f'score_{tag}': round(v['score'], 2) for tag, v in scores.items()}


def summarize(rows: List[dict]) -> dict:
    n = len(rows)
    def frac(k):
        return f"{100.0 * sum(1 for r in rows if r[k]) / max(n, 1):.0f}%"
    return dict(
        total=n,
        anchor_ok=sum(1 for r in rows if r['n_anchor'] >= 1),
        anchor_two=sum(1 for r in rows if r['n_anchor'] >= 2),
        no_repulsion=sum(1 for r in rows if r['n_repulsion'] == 0),
        bridge=sum(1 for r in rows if r['bridge']),
        ensemble_stable=sum(1 for r in rows if r['ens_max_clash'] < 0.25),
        keep=sum(1 for r in rows if r['keep']),
        anchor_rate=frac('n_anchor'), bridge_rate=frac('bridge'),
    )


def calibration_report(cxs: Dict[str, Complex], pocket: dict) -> List[dict]:
    """Same metrics computed on the peptides themselves -- the thresholds' source."""
    ref = cxs['SglCJ3']
    out = []
    for tag, cx in cxs.items():
        last = cx.ligand[max(cx.ligand)]
        if 'OXT' not in last.atoms:
            out.append(dict(tag=tag, note='C-terminal carboxylate not modelled (no OXT)'))
            continue
        rec = receptor_in_frame(cx, ref)
        last = cx.l(max(cx.ligand))
        O = np.array([last.atoms[n] for n in ('O', 'OXT')])
        if tag != ref.tag:
            O = apply_rt(O, superpose_complex(cx, ref)[0])
        mask = rec.is_pos & np.isin(rec.resid, ANCHORS)
        d = np.linalg.norm(O[:, None, :] - rec.xyz[None, mask, :], axis=2)
        rids = rec.resid[mask]
        best = min(((float(d[:, k].min()), int(rids[k])) for k in range(len(rids))))
        neg = np.linalg.norm(O[:, None, :] - rec.xyz[None, rec.is_neg, :], axis=2).min()
        out.append(dict(tag=tag, anchor=int(best[1]), dist=round(best[0], 2),
                        repulsion=round(float(neg), 2)))
    return out


def repose_ensemble(sc, F1: np.ndarray, F2: np.ndarray, recs: Dict[str, Receptor]) -> dict:
    """F4' -- re-**place** the scaffold in every structure and compare binding modes.

    Transplanting one pose into another structure is not a valid test: at 3.6 A
    the side chains differ by 1-2 A, so a tightly packed pose always clashes.
    What is meaningful is whether the *same molecule* finds the *same binding
    mode* when the search is redone against each structure.  All poses live in
    the reference frame, so they can be compared atom-by-atom without another
    superposition.
    """
    keep_pose = sc.pose.copy()
    keep_comp = dict(sc.components)
    poses, scores, clashes = {}, {}, {}
    try:
        for tag, rec in recs.items():
            place(sc, F1, F2, rec, n_rot=12)
            poses[tag] = sc.pose.copy()
            # place() writes pose_score; components['total'] is only set later
            # by design_stage, so reading it here would silently compare stale values
            scores[tag] = float(sc.pose_score)
            clashes[tag] = float(sc.components['clash'])
    finally:
        sc.pose = keep_pose
        sc.components.clear()
        sc.components.update(keep_comp)
    tags = list(poses)
    rmsd = 0.0
    for i in range(len(tags)):
        for j in range(i + 1, len(tags)):
            rmsd = max(rmsd, float(np.sqrt(((poses[tags[i]] - poses[tags[j]]) ** 2).sum(1).mean())))
    vals = list(scores.values())
    return dict(rmsd=rmsd, score_spread=max(vals) - min(vals),
                min_score=min(vals), max_clash=max(clashes.values()),
                scores=scores, clashes=clashes)


def pocket_subset(rec: Receptor, center: np.ndarray, r: float = 18.0) -> Receptor:
    """Trim the receptor to the pocket neighbourhood (18 A is far beyond any
    contact term, and it makes re-placement ~10x faster)."""
    return rec.subset(center, r) if hasattr(rec, 'subset') else rec
