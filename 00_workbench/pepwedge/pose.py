# -*- coding: utf-8 -*-
"""pepwedge.pose -- rigid placement of an assembled skeleton on a design axis
and a contact-based scoring function.

The score deliberately mirrors the interactions that make the Sgl peptides
inhibitory, so a good number means "reproduces the peptide's chemical logic",
not merely "fits the pocket".
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np

from .core import VDW
from .design import Scaffold

RADII = np.array([VDW.get('C', 1.7), VDW.get('N', 1.55), VDW.get('O', 1.52),
                  VDW.get('S', 1.8)])


@dataclass
class Receptor:
    xyz: np.ndarray
    radius: np.ndarray
    is_pos: np.ndarray          # Lys NZ / Arg NH1,NH2
    is_neg: np.ndarray          # Asp/Glu carboxylate O
    is_polar: np.ndarray        # N / O
    is_carbon: np.ndarray       # C / S
    resname: List[str]
    resid: np.ndarray

    def subset(self, center: np.ndarray, r: float = 14.0):
        d = np.linalg.norm(self.xyz - center, axis=1)
        m = d < r
        return Receptor(self.xyz[m], self.radius[m], self.is_pos[m], self.is_neg[m],
                        self.is_polar[m], self.is_carbon[m],
                        [n for n, k in zip(self.resname, m) if k], self.resid[m])


def receptor_from_complex(cx) -> Receptor:
    xyz, radius, pos, neg, polar, carbon, names, ids = [], [], [], [], [], [], [], []
    for rid, r in cx.receptor.items():
        for an, coord in r.atoms.items():
            el = an[0] if an[0].isalpha() else 'C'
            xyz.append(coord)
            radius.append(VDW.get(el.upper(), 1.8))
            pos.append(an in ('NZ', 'NH1', 'NH2') and r.name in ('LYS', 'ARG'))
            neg.append(r.name in ('ASP', 'GLU') and el.upper() == 'O')
            polar.append(el.upper() in ('N', 'O'))
            carbon.append(el.upper() in ('C', 'S'))
            names.append(r.name)
            ids.append(rid)
    return Receptor(np.array(xyz), np.array(radius), np.array(pos), np.array(neg),
                    np.array(polar), np.array(carbon), names, np.array(ids))


def _rot_to_axis(u: np.ndarray) -> np.ndarray:
    """Rotation matrix mapping +x onto unit vector u."""
    x = np.array([1.0, 0.0, 0.0])
    v = np.cross(x, u)
    c = float(np.dot(x, u))
    if np.linalg.norm(v) < 1e-8:
        return np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))


def score_pose(lig: np.ndarray, els: List[str], chg: np.ndarray, rec: Receptor,
               clash_k=6.0, salt_w=3.0, hb_w=0.5, hyd_w=0.10, exp_w=0.40):
    d = np.linalg.norm(lig[:, None, :] - rec.xyz[None, :, :], axis=2)
    el = np.array([e.upper() for e in els])
    r_lig = np.where(el == 'O', VDW['O'], np.where(el == 'N', VDW['N'],
                    np.where(el == 'S', VDW['S'], VDW['C'])))
    thr = (r_lig[:, None] + rec.radius[None, :]) * 0.78
    clash = float(np.clip(thr - d, 0, None).__pow__(2).sum())

    is_anion = chg < 0                       # carboxylate oxygens only
    is_lig_no = np.isin(el, ['N', 'O'])
    salt = 0
    if is_anion.any() and rec.is_pos.any():
        sub = d[np.ix_(is_anion, rec.is_pos)]
        salt = int((sub < 4.2).any(1).sum())
    hb = 0
    if is_lig_no.any() and rec.is_polar.any():
        sub = d[np.ix_(is_lig_no, rec.is_polar)]
        # exclude pairs that are already counted as salt bridges
        hb = int(((sub < 3.65) & (sub > 2.2)).any(1).sum())
    is_lig_c = np.isin(el, ['C', 'S'])
    hyd = int((d[np.ix_(is_lig_c, rec.is_carbon)] < 4.6).any(1).sum())
    exposed = int((d.min(1) > 5.5).sum())

    n = len(lig)
    buried = 1.0 - exposed / max(1, n)
    total = (-clash_k * clash + salt_w * salt + hb_w * hb + hyd_w * hyd
             - exp_w * exposed + 0.8 * buried)
    return total, dict(clash=round(clash, 2), salt=salt, hb=hb, hyd=hyd,
                       exposed=exposed, buried=round(buried, 2))


def place(sc: Scaffold, F1: np.ndarray, F2: np.ndarray, rec: Receptor, n_rot=12):
    """Rigid-body placement: head feature on F1, axis pointing to F2."""
    axis = F2 - F1
    L = float(np.linalg.norm(axis))
    if L < 1e-3:
        return None
    u = axis / L
    R0 = _rot_to_axis(u)
    best = None
    for mirror in (1.0, -1.0):
        for k in range(n_rot):
            phi = 2 * math.pi * k / n_rot
            c, s = math.cos(phi), math.sin(phi)
            Rz = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
            M = R0 @ Rz
            loc = sc.atoms.copy()
            loc[:, 1] *= mirror
            xyz = loc @ M.T + F1
            tot, comp = score_pose(xyz, sc.elements, sc.charge, rec)
            if best is None or tot > best[0]:
                best = (tot, comp, xyz, mirror, phi)
    tot, comp, xyz, mirror, phi = best
    sc.pose = xyz
    sc.pose_score = round(tot, 3)
    sc.components.update(comp)
    sc.components['mirror'] = mirror
    sc.components['phi'] = round(phi, 2)
    return sc


def pharmacophore_match(sc: Scaffold, feats, F1, F2, dtol=1.8):
    """How close the realised feature positions are to the 3D pharmacophore."""
    if sc.pose is None:
        return 0.0
    # head feature = first atom of the head cap (x=0)
    p_head = sc.pose[0]
    # tail feature: nearest ligand atom to F2
    d = np.linalg.norm(sc.pose - F2, axis=1)
    p_tail = sc.pose[int(d.argmin())]
    e1 = np.linalg.norm(p_head - F1)
    e2 = np.linalg.norm(p_tail - F2)
    return float(max(0.0, 1.5 - (e1 + e2) / 2.0) / 1.5)
