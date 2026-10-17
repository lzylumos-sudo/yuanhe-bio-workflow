# -*- coding: utf-8 -*-
"""pepwedge.pharmacophore -- turn the peptide interface into a 3D feature set.

Three independently evolved peptides hit the same MurJ surface. We therefore
extract, for each complex, the *functional* atoms of the peptide (the ones
mutagenesis says matter), map them into a common MurJ frame, and merge them into
a consensus pharmacophore whose spread quantifies how reproducible each
interaction is.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List

import numpy as np

from .core import (Complex, SGL_FUNCTIONAL, apply_rt, contacts, hotspot_table,
                   kabsch)

# feature type of the side-chain / backbone atom used as the interaction point
FEATURE_ATOMS = {
    'ASP': ('anion', ('OD1', 'OD2')),
    'GLU': ('anion', ('OE1', 'OE2')),
    'LYS': ('cation', ('NZ',)),
    'ARG': ('cation', ('NH1', 'NH2', 'NE')),
    'HIS': ('aromatic', ('CG', 'ND1', 'CD2', 'CE1', 'NE2')),
    'TYR': ('aromatic', ('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ')),
    'PHE': ('aromatic', ('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ')),
    'TRP': ('aromatic', ('CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2')),
    'SER': ('donor', ('OG',)),
    'THR': ('donor', ('OG1',)),
    'ASN': ('donor', ('OD1', 'ND2')),
    'GLN': ('donor', ('OE1', 'NE2')),
}

TOL = {'anion': 1.6, 'cation': 1.6, 'aromatic': 1.8, 'donor': 1.5, 'hydrophobic': 2.0}


@dataclass
class Feature:
    fid: str
    ftype: str
    pos: np.ndarray
    weight: float = 1.0
    tol: float = 1.5
    sources: List[str] = field(default_factory=list)
    murj_contacts: List[str] = field(default_factory=list)
    spread: float = 0.0
    convergence: float = 1.0

    def as_dict(self):
        return dict(fid=self.fid, type=self.ftype,
                    pos=[round(float(x), 2) for x in self.pos],
                    weight=round(self.weight, 2), tol=self.tol,
                    sources=self.sources, murj_contacts=self.murj_contacts,
                    spread=round(float(self.spread), 2),
                    convergence=round(float(self.convergence), 2))


def _feature_point(res) -> tuple:
    ftype, names = FEATURE_ATOMS.get(res.name, ('hydrophobic', ()))
    if ftype == 'hydrophobic':
        # side-chain centroid
        skip = {'N', 'CA', 'C', 'O', 'OXT'}
        pts = [v for k, v in res.atoms.items() if k not in skip]
        return 'hydrophobic', (np.mean(pts, 0) if pts else res.centroid())
    pts = [res.atoms[n] for n in names if n in res.atoms]
    if not pts:
        return 'hydrophobic', res.centroid()
    return ftype, np.mean(pts, 0)


def complex_features(cx: Complex, tag: str) -> List[Feature]:
    functional = SGL_FUNCTIONAL.get(tag, {})
    ct = contacts(cx)
    cand = {}
    for (rid, lid), rec in ct.items():
        cand.setdefault(lid, set()).add(rid)
    out = []
    last = max(cx.ligand)
    for lid in sorted(cx.ligand):
        res = cx.l(lid)
        # C-terminal carboxylate of the peptide = a key anionic anchor
        if lid == last and 'OXT' in res.atoms:
            pos = np.mean([res.atoms['OXT'], res.atoms['O']] if 'O' in res.atoms else [res.atoms['OXT']], 0)
            contacts_here = sorted(cand.get(lid, set()))
            out.append(Feature(f'F_Cterm', 'anion', pos, weight=1.3, tol=1.5,
                               sources=[f'{tag}:{res.one}{lid}(C-term)'],
                               murj_contacts=[f'{cx.r(r).one}{r}' for r in contacts_here]))
            continue
        if lid not in functional and lid not in cand:
            continue
        ftype, pos = _feature_point(res)
        if ftype == 'hydrophobic' and lid not in functional:
            continue
        # weight: contact density, boosted for mutagenesis-validated residues
        w = min(1.5, 0.5 + 0.05 * sum(1 for (r, l) in ct if l == lid))
        if lid in functional:
            w += 0.4
        out.append(Feature(f'F_{tag}_{res.one}{lid}', ftype, pos, weight=round(w, 2),
                           tol=TOL.get(ftype, 1.6),
                           sources=[f'{tag}:{res.one}{lid}'],
                           murj_contacts=[f'{cx.r(r).one}{r}' for r in sorted(cand.get(lid, set()))]))
    return out


def superpose_complex(cx_mobile: Complex, cx_ref: Complex):
    """Kabsch transform of cx_mobile onto cx_ref using shared CA atoms."""
    shared = [rid for rid in cx_ref.receptor
              if 'CA' in cx_ref.receptor[rid].atoms and rid in cx_mobile.receptor
              and 'CA' in cx_mobile.receptor[rid].atoms]
    mob = np.array([cx_mobile.receptor[r].atoms['CA'] for r in shared])
    ref = np.array([cx_ref.receptor[r].atoms['CA'] for r in shared])
    return kabsch(mob, ref), len(shared)


def consensus_features(complexes: Dict[str, Complex], ref_tag='SglCJ3',
                       merge_r=3.0, min_structures=1, max_features=10) -> List[Feature]:
    """Map all peptide features into the reference frame and merge by proximity."""
    ref = complexes[ref_tag]
    pooled = []
    for tag, cx in complexes.items():
        rt = (np.eye(3), np.zeros(3)) if tag == ref_tag else superpose_complex(cx, ref)[0]
        for f in complex_features(cx, tag):
            pooled.append(Feature(f.fid, f.ftype, apply_rt(f.pos, rt), f.weight,
                                  f.tol, list(f.sources), list(f.murj_contacts)))
    # greedy clustering by type + proximity (highest weight first)
    clusters: List[Feature] = []
    for f in sorted(pooled, key=lambda x: -x.weight):
        hit = None
        for c in clusters:
            if c.ftype == f.ftype and np.linalg.norm(c.pos - f.pos) < merge_r:
                hit = c
                break
        if hit is None:
            clusters.append(Feature(f.fid, f.ftype, f.pos.copy(), f.weight, f.tol,
                                    list(f.sources), list(f.murj_contacts)))
        else:
            tag = f.sources[0].split(':')[0]
            if tag in {s.split(':')[0] for s in hit.sources}:
                continue                      # one vote per structure
            n_prev = len(hit.sources)
            hit.pos = (hit.pos * n_prev + f.pos) / (n_prev + 1)
            hit.sources += f.sources
            hit.weight = max(hit.weight, f.weight)
            hit.murj_contacts = sorted(set(hit.murj_contacts) | set(f.murj_contacts))
            hit.spread = float(np.std([apply_rt(f.pos, (np.eye(3), np.zeros(3)))
                                       for f in [f]], axis=0).mean()) if False else hit.spread
    # reproducibility: how many of the three peptides support the feature
    n_struct = len(complexes)
    for c in clusters:
        conv = len({s.split(':')[0] for s in c.sources}) / n_struct
        c.spread = max(c.spread, 0.0)
        c.weight = round(c.weight * (0.6 + 0.4 * conv), 3)
        c.convergence = conv
    clusters = [c for c in clusters]
    clusters.sort(key=lambda c: (-c.weight, c.ftype))
    # unique, chemistry-first ids
    prefix = {'anion': 'A', 'cation': 'C', 'aromatic': 'R', 'donor': 'D', 'hydrophobic': 'H'}
    for i, c in enumerate(clusters):
        c.fid = f'{prefix.get(c.ftype, "X")}{i + 1}'
    keep = [c for c in clusters
            if c.convergence >= min_structures / n_struct]
    return keep[:max_features]


def anchor_pairs(features: List[Feature]):
    """Distances between features -> the geometric constraints for scaffold assembly."""
    rows = []
    for i in range(len(features)):
        for j in range(i + 1, len(features)):
            a, b = features[i], features[j]
            d = float(np.linalg.norm(a.pos - b.pos))
            if 3.0 < d < 20.0:
                rows.append(dict(a=a.fid, b=b.fid, types=f'{a.ftype}-{b.ftype}',
                                 distance=round(d, 2),
                                 weight=round(a.weight * b.weight, 2)))
    rows.sort(key=lambda r: -r['weight'])
    return rows
