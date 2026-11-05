# -*- coding: utf-8 -*-
"""pepwedge.pocket -- cavity voxelisation, druggability and sub-pocket partition.

The peptide defines a *negative image* of the cavity. We rasterise that cavity,
score how enclosed it is (ray-cast buriedness, a lightweight SiteMap analogue),
label its chemistry, and split it into sub-pockets that the scaffold assembler
can address one fragment at a time.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List

import numpy as np

from .core import Complex, _cell_lists, _neighbors

# sub-pocket definitions: MurJ residues that line each region (numbering = EcMurJ)
SUBPOCKETS = {
    'S1_basic_bowl':   [46, 270, 269, 244, 243, 267, 266],
    'S2_carboxylate_gate': [53, 368, 273, 367, 57, 364],
    'S3_TM2_TM7_groove':   [42, 47, 50, 51, 54, 112, 113, 115, 117],
    'S4_TM7_floor':        [237, 238, 240, 241, 245, 247, 248, 251, 252, 255],
}

# Execution budget only: keep the same distances, cutoffs and point ordering.
# Account conservatively for vector temporaries and the previous reduction.
_DISTANCE_WORK_BYTES = 64 * 1024 * 1024


def _distance_rows(reference_atoms, requested=200000):
    if requested <= 0:
        raise ValueError('Distance chunk size must be positive')
    bytes_per_row = max(1, reference_atoms) * np.dtype(np.float64).itemsize * 10
    return min(requested, max(1, _DISTANCE_WORK_BYTES // bytes_per_row))


@dataclass
class Pocket:
    points: np.ndarray
    step: float
    buriedness: np.ndarray
    chemistry: List[str]
    subpocket: List[str]
    extra: Dict = field(default_factory=dict)

    @property
    def volume(self):
        return len(self.points) * self.step ** 3


def pocket_voxels(cx: Complex, step=0.8, clearance=3.0, envelope=7.0, chunk=200000):
    """Grid points inside the peptide envelope that are at least `clearance` A
    away from every MurJ heavy atom (i.e. sterically accessible to a ligand)."""
    ra = cx.r_xyz
    la = cx.l_xyz
    lo = la.min(0) - envelope
    hi = la.max(0) + envelope
    axes = [np.arange(lo[i], hi[i] + step / 2, step) for i in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, 3)

    # keep only voxels inside the peptide envelope
    keep = np.zeros(len(grid), bool)
    for s in range(0, len(la), 200):
        atoms = la[s:s + 200]
        rows = _distance_rows(len(atoms))
        for start in range(0, len(grid), rows):
            g = grid[start:start + rows]
            d = np.linalg.norm(g[:, None, :] - atoms[None, :, :], axis=2)
            keep[start:start + rows] |= d.min(1) <= envelope
        if keep.all():
            break
    grid = grid[keep]

    # remove voxels clashing with the receptor
    clash = np.zeros(len(grid), bool)
    chunk = _distance_rows(len(ra), chunk)
    for s in range(0, len(grid), chunk):
        g = grid[s:s + chunk]
        d2 = ((g[:, None, :] - ra[None, :, :]) ** 2).sum(-1)
        clash[s:s + chunk] = d2.min(1) < clearance ** 2
    grid = grid[~clash]
    return grid, step


def _dirs(n=50):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


def buriedness(points: np.ndarray, receptor_xyz: np.ndarray, n_dir=50, probe=1.8, rmax=13.0):
    """Fraction of rays that hit a receptor atom within rmax -> enclosure score."""
    dirs = _dirs(n_dir)
    grid, lo, cell = _cell_lists(receptor_xyz, rmax)
    out = np.zeros(len(points))
    for i, p in enumerate(points):
        k = tuple(np.floor((p - lo) / cell).astype(int))
        nb = _neighbors(receptor_xyz, grid, lo, cell, k, rmax)
        if len(nb) == 0:
            continue
        v = receptor_xyz[nb] - p
        t = v @ dirs.T                       # (natoms, ndir) projection
        d2 = (v ** 2).sum(1)[:, None] - t ** 2
        hit = (t > 0) & (d2 < probe ** 2)
        out[i] = hit.any(0).mean()
    return out


def _frag(rid):
    return rid


CHEM_FROM_RES = {
    'LYS': 'pos', 'ARG': 'pos', 'HIS': 'pos',
    'ASP': 'neg', 'GLU': 'neg',
    'SER': 'polar', 'THR': 'polar', 'ASN': 'polar', 'GLN': 'polar',
    'TYR': 'aromatic', 'TRP': 'aromatic', 'PHE': 'aromatic',
}


def label_chemistry(cx: Complex, points):
    """Chemistry of the receptor wall nearest to each voxel."""
    ra = cx.r_xyz
    labels = []
    owner_res = []
    for rid, r in cx.receptor.items():
        for _ in r.atoms:
            owner_res.append(rid)
    owner_res = np.array(owner_res)
    rows = _distance_rows(len(ra), 20000)
    for s in range(0, len(points), rows):
        g = points[s:s + rows]
        d2 = ((g[:, None, :] - ra[None, :, :]) ** 2).sum(-1)
        idx = d2.argmin(1)
        for j in idx:
            r = cx.r(owner_res[j])
            labels.append(CHEM_FROM_RES.get(r.name, 'hydrophobic'))
    return labels


def subpocket_of(cx: Complex, points):
    """Assign each voxel to the nearest sub-pocket centroid (chemical partition)."""
    cents = {}
    for name, rids in SUBPOCKETS.items():
        pts = [cx.r(rid).centroid() for rid in rids if rid in cx.receptor]
        if pts:
            cents[name] = np.mean(pts, 0)
    names = list(cents)
    C = np.array([cents[n] for n in names])
    d = np.linalg.norm(points[:, None, :] - C[None, :, :], axis=2)
    return [names[i] for i in d.argmin(1)]


def build_pocket(cx: Complex, step=0.8, clearance=3.0, envelope=7.0, n_dir=50):
    pts, step = pocket_voxels(cx, step=step, clearance=clearance, envelope=envelope)
    bur = buriedness(pts, cx.r_xyz, n_dir=n_dir)
    chem = label_chemistry(cx, pts)
    sp = subpocket_of(cx, pts)
    return Pocket(pts, step, bur, chem, sp)


def pocket_report(pk: Pocket, cx: Complex):
    rows = []
    for name in SUBPOCKETS:
        sel = np.array([s == name for s in pk.subpocket])
        if not sel.any():
            continue
        p = pk.points[sel]
        chem = [c for c, s in zip(pk.chemistry, pk.subpocket) if s == name]
        rows.append(dict(subpocket=name,
                         n_voxels=int(sel.sum()),
                         volume_A3=round(float(sel.sum() * pk.step ** 3), 1),
                         centroid=np.round(p.mean(0), 2).tolist(),
                         extent_A=np.round(p.max(0) - p.min(0), 1).tolist(),
                         mean_buriedness=round(float(pk.buriedness[sel].mean()), 3),
                         chemistry={k: chem.count(k) for k in ('pos', 'neg', 'polar', 'aromatic', 'hydrophobic')},
                         residues=[f'{cx.r(r).one}{r}' for r in SUBPOCKETS[name] if r in cx.receptor]))
    rows.sort(key=lambda r: -r['volume_A3'])
    return rows


def subpocket_centroids(pk: Pocket):
    out = {}
    for name in SUBPOCKETS:
        sel = np.array([s == name for s in pk.subpocket])
        if sel.any():
            out[name] = pk.points[sel].mean(0)
    return out
