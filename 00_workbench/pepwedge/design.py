# -*- coding: utf-8 -*-
"""pepwedge.design -- axis-constrained fragment assembly (the scaffold engine).

Idea
----
A pharmacophore pair (F1, F2) with separation D defines a *design axis*.  A
scaffold is a linear string of rigid units laid along that axis:

    [head cap]==[link]==[ring]==[link]==[ring]==...==[link]==[tail cap]

Every unit contributes a fixed length along the axis, so "does this skeleton
reach both anchors" becomes a length-constraint problem that is solved by DFS /
interval lookup -- no stochastic search, fully enumerable and explainable.

A surviving skeleton is then placed in 3D (rigid body, rotation sampled about
the axis) and scored with a contact function that rewards exactly the
interactions the peptides use: salt bridges to K46/R53/R270/K368, H-bonds to
Q244/S241/T251, hydrophobic burial in the TM2/TM7 groove.

The output is a *scaffold* (skeleton + attachment points), which a downstream
generative model (Pocket2Mol / TargetDiff / REINVENT) elaborates into full
structures.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

MASS = {'C': 12.011, 'N': 14.007, 'O': 15.999, 'S': 32.06, 'H': 1.008}


# ---------------------------------------------------------------------------
@dataclass
class Atom:
    el: str
    xyz: Tuple[float, float, float]
    nH: int = 0
    charge: int = 0


@dataclass
class Unit:
    uid: str
    kind: str                 # 'ring' | 'link' | 'cap'
    span: float               # axis length from the unit's left port to its right port
    atoms: List[Atom] = field(default_factory=list)
    feature: Optional[str] = None      # pharmacophore type carried by a cap
    feat_x: float = 0.0                # x of the feature point in unit frame
    smi: str = ''
    tpsa: float = 0.0
    clogp: float = 0.0
    rot: int = 0
    rings: int = 0
    label: str = ''


def _hexagon(radius=1.40):
    """Six vertices; vertex 0 at the origin, vertex 3 at (2r, 0)."""
    c = np.array([radius, 0.0])
    return [c + radius * np.array([math.cos(math.pi + k * math.pi / 3),
                                   math.sin(math.pi + k * math.pi / 3)]) for k in range(6)]


def _pentagon_ports(radius=1.19):
    """Pentagon with port vertices 0 and 2 on the x axis (2,5-relationship)."""
    pts = [np.array([radius * math.cos(math.pi / 2 + 2 * math.pi * k / 5),
                     radius * math.sin(math.pi / 2 + 2 * math.pi * k / 5)]) for k in range(5)]
    p, q = pts[0], pts[2]
    d = q - p
    ang = math.atan2(d[1], d[0])
    R = np.array([[math.cos(-ang), -math.sin(-ang)], [math.sin(-ang), math.cos(-ang)]])
    out = [R @ (v - p) for v in pts]
    return out, float(np.linalg.norm(d))


# --- ring units (ports at vertex 0 and n//2) -------------------------------
def _ring(uid, verts, span, hetero, hmap, smi, tpsa, clogp, label):
    atoms = [Atom(hetero.get(k, 'C'), (float(v[0]), float(v[1]), 0.0), hmap.get(k, 0))
             for k, v in enumerate(verts)]
    return Unit(uid, 'ring', span, atoms, None, 0.0, smi, tpsa, clogp, 0, 1, label)


PHENYLENE = _ring('Ph', _hexagon(1.40), 2.80, {}, {0: 0, 1: 1, 2: 1, 3: 0, 4: 1, 5: 1},
                  '{prev}c{d}ccc({next})cc{d}', 0.0, 1.7, '1,4-phenylene')
_cy = _hexagon(1.40)
CYCLOHEXYLENE = _ring('Cx', _cy, 2.80, {}, {0: 1, 1: 2, 2: 2, 3: 1, 4: 2, 5: 2},
                      '{prev}C{d}CCC({next})CC{d}', 0.0, 2.5, '1,4-cyclohexylene')
_th_v, _th_L = _pentagon_ports(1.19)
THIOPHENE = _ring('Th', _th_v, _th_L, {1: 'S'}, {0: 0, 1: 0, 2: 0, 3: 1, 4: 1},
                  '{prev}c{d}sc({next})cc{d}', 28.2, 1.8, 'thiophene-2,5-diyl')

RINGS = {u.uid: u for u in (PHENYLENE, CYCLOHEXYLENE, THIOPHENE)}

# --- link units ------------------------------------------------------------
LINKS = {
    'bond':   Unit('bond', 'link', 1.48, [], None, 0, '', 0.0, 0.0, 1, 0, 'biaryl bond'),
    'alkyne': Unit('alkyne', 'link', 4.06,
                   [Atom('C', (1.43, 0, 0), 0), Atom('C', (2.63, 0, 0), 0)],
                   None, 0, 'C#C', 0.0, 0.4, 0, 0, 'C#C'),
    'ch2':    Unit('ch2', 'link', 2.96, [Atom('C', (1.48, 0, 0), 2)],
                   None, 0, 'C', 0.0, 0.5, 1, 0, 'CH2'),
    'ether':  Unit('ether', 'link', 2.72, [Atom('O', (1.36, 0, 0), 0)],
                   None, 0, 'O', 9.2, -0.3, 1, 0, 'ether'),
    'amide':  Unit('amide', 'link', 4.25,
                   [Atom('C', (1.50, 0, 0), 0), Atom('O', (1.98, 1.02, 0.0), 0),
                    Atom('N', (2.85, 0.0, 0.0), 1)],
                   None, 0, 'C(=O)N', 29.1, -0.8, 0, 0, 'amide'),
    'ch2ch2': Unit('ch2ch2', 'link', 4.49,
                   [Atom('C', (1.48, 0, 0), 2), Atom('C', (3.01, 0, 0), 2)],
                   None, 0, 'CC', 0.0, 1.0, 2, 0, 'CH2CH2'),
}


# --- caps ------------------------------------------------------------------
def _cap_acid():
    return Unit('COOH', 'cap', 0.0,
                [Atom('C', (0, 0, 0), 0, 0), Atom('O', (-0.62, 0.86, 0), 1, -1),
                 Atom('O', (-0.62, -0.86, 0), 0, -1)],
                'anion', 0.0, 'OC(=O)', 37.3, -0.3, 0, 0, 'carboxylic acid')


def _cap_acid_ch2():
    return Unit('CH2COOH', 'cap', 1.52,
                [Atom('C', (0, 0, 0), 0, 0), Atom('O', (-0.62, 0.86, 0), 1, -1),
                 Atom('O', (-0.62, -0.86, 0), 0, -1), Atom('C', (1.52, 0, 0), 2)],
                'anion', 0.0, 'OC(=O)C', 37.3, 0.2, 1, 0, 'phenylacetic acid')


def _cap_phenyl_like(uid, hetero, hmap, smi, tpsa, clogp, label):
    v = _hexagon(1.40)
    atoms = [Atom(hetero.get(k, 'C'), (float(v[k][0]), float(v[k][1]), 0.0), hmap.get(k, 0))
             for k in range(6)]
    return Unit(uid, 'cap', 1.40, atoms, 'aromatic', 1.40, smi, tpsa, clogp, 0, 1, label)


CAP_PH = _cap_phenyl_like('Ph', {}, {0: 0, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1},
                          'c{d}ccccc{d}', 0.0, 1.9, 'phenyl')
CAP_PY = _cap_phenyl_like('Py', {2: 'N'}, {0: 0, 1: 1, 2: 0, 3: 1, 4: 1, 5: 1},
                          'c{d}cnccc{d}', 12.9, 0.9, 'pyridin-3-yl')
CAP_PYR = _cap_phenyl_like('Pyr', {2: 'N', 4: 'N'}, {0: 0, 1: 1, 2: 0, 3: 1, 4: 0, 5: 1},
                           'c{d}nccnc{d}', 25.8, 0.3, 'pyrimidin-5-yl')
CAP_CY = _cap_phenyl_like('Cy', {}, {0: 1, 1: 2, 2: 2, 3: 2, 4: 2, 5: 2},
                          'C{d}CCCCC{d}', 0.0, 2.5, 'cyclohexyl')
CAP_CY.feature = 'hydrophobic'

CAP_TBU = Unit('tBu', 'cap', 0.0,
               [Atom('C', (0, 0, 0), 0), Atom('C', (1.45, 0.6, 0.5), 3),
                Atom('C', (1.45, -0.6, 0.6), 3), Atom('C', (1.30, 0.0, -1.0), 3)],
               'hydrophobic', 0.0, 'C(C)(C)C', 0.0, 2.3, 0, 0, 'tert-butyl')


def _cap_benzoic():
    v = _hexagon(1.40)
    atoms = [Atom('C', (float(v[k][0]), float(v[k][1]), 0.0),
                  0 if k in (0, 3) else 1) for k in range(6)]
    atoms += [Atom('C', (4.28, 0, 0), 0, 0), Atom('O', (4.90, 0.86, 0), 1, -1),
              Atom('O', (4.90, -0.86, 0), 0, -1)]
    return Unit('PhCOOH', 'cap', 4.28, atoms, 'anion', 4.28,
                'c{d}ccc(C(=O)O)cc{d}', 37.3, 1.4, 0, 1, 'benzoic acid')


HEADS = {u.uid: u for u in (_cap_acid(), _cap_acid_ch2())}
TAILS = {u.uid: u for u in (CAP_PH, CAP_PY, CAP_PYR, CAP_CY, CAP_TBU, _cap_benzoic())}


# ---------------------------------------------------------------------------
@dataclass
class Scaffold:
    head: str
    tail: str
    rings: Tuple[str, ...]
    links: Tuple[str, ...]
    span: float
    smiles: str
    atoms: np.ndarray
    elements: List[str]
    nH: np.ndarray
    charge: np.ndarray
    mw: float
    clogp: float
    tpsa: float
    hbd: int
    hba: int
    rot: int
    n_ring: int
    pose: Optional[np.ndarray] = None
    pose_score: float = 0.0
    components: Dict[str, float] = field(default_factory=dict)
    axis: Tuple[str, str] = ('', '')
    axis_length: float = 0.0
    sid: str = ''

    def props(self):
        return dict(MW=round(self.mw, 1), cLogP=round(self.clogp, 2),
                    TPSA=round(self.tpsa, 1), HBD=self.hbd, HBA=self.hba,
                    RotB=self.rot, Rings=self.n_ring,
                    heavy_atoms=len(self.elements))


def compose(head: Unit, rings: Sequence[Unit], links: Sequence[Unit], tail: Unit) -> Scaffold:
    """Lay the units along +x; build coordinates and SMILES together."""
    # ---- SMILES: recursive, each ring gets its own ring-closure digit
    n = len(rings)
    digits = list(range(1, n + 3))
    rest = links[n].smi + tail.smi.format(d=digits[n], next='')
    for i in reversed(range(n)):
        rest = links[i].smi + rings[i].smi.format(prev='', d=digits[i], next=rest)
    smiles = head.smi + rest

    # ---- 3D
    atoms: List[Atom] = [Atom(a.el, a.xyz, a.nH, a.charge) for a in head.atoms]
    x = head.span

    def push(unit):
        nonlocal x
        for a in unit.atoms:
            atoms.append(Atom(a.el, (x + a.xyz[0], a.xyz[1], a.xyz[2]), a.nH, a.charge))
        x += unit.span

    for i, ring in enumerate(rings):
        push(links[i])
        push(ring)
    push(links[n])
    push(tail)

    xyz = np.array([a.xyz for a in atoms], dtype=float)
    els = [a.el for a in atoms]
    nH = np.array([a.nH for a in atoms], dtype=float)
    chg = np.array([a.charge for a in atoms], dtype=float)
    mw = float(sum(MASS[e] for e in els) + nH.sum() * MASS['H'])
    tpsa = sum(u.tpsa for u in [head, tail, *rings, *links])
    clogp = sum(u.clogp for u in [head, tail, *rings, *links])
    hbd = int(sum(1 for e, h in zip(els, nH) if e in 'NO' and h > 0))
    hba = int(sum(1 for e in els if e in 'NO'))
    rot = int(sum(u.rot for u in [head, tail, *rings, *links]))
    nring = int(sum(u.rings for u in [head, tail, *rings]))
    return Scaffold(head.uid, tail.uid, tuple(r.uid for r in rings),
                    tuple(l.uid for l in links), float(x), smiles, xyz, els, nH, chg,
                    mw, clogp, tpsa, hbd, hba, rot, nring)


# ---------------------------------------------------------------------------
_PATH_CACHE: Dict[int, Tuple[np.ndarray, list]] = {}


def _link_paths(k: int):
    if k in _PATH_CACHE:
        return _PATH_CACHE[k]
    units = list(LINKS.values())
    combos = list(itertools.product(units, repeat=k))
    sums = np.array([sum(u.span for u in c) for c in combos])
    order = np.argsort(sums)
    _PATH_CACHE[k] = (sums[order], [combos[i] for i in order])
    return _PATH_CACHE[k]


FIRST_LINKS = {'bond', 'alkyne', 'ch2', 'ch2ch2'}


def enumerate_scaffolds(D: float, tol=0.55, ring_counts=(1, 2, 3), max_mw=520.0):
    """All unit sequences whose axis length matches D within tol."""
    out: List[Scaffold] = []
    ring_units = list(RINGS.values())
    for hname, head in HEADS.items():
        for tname, tail in TAILS.items():
            for nr in ring_counts:
                base = head.span + tail.span
                # cheap span screen
                if base + nr * 2.80 + (nr + 1) * 1.48 > D + tol + 6:
                    pass
                for rcombo in itertools.product(ring_units, repeat=nr):
                    ring_span = sum(r.span for r in rcombo)
                    need = D - head.span - tail.span - ring_span
                    sums, combos = _link_paths(nr + 1)
                    i0 = int(np.searchsorted(sums, need - tol))
                    i1 = int(np.searchsorted(sums, need + tol))
                    for lcombo in combos[i0:i1]:
                        if lcombo[0].uid not in FIRST_LINKS:
                            continue
                        uids = [l.uid for l in lcombo]
                        if max(uids.count(u) for u in set(uids)) > 2:
                            continue
                        if uids.count('amide') > 1:
                            continue
                        approx_mw = (sum(MASS[a.el] for u in (head, tail, *lcombo, *rcombo) for a in u.atoms)
                                     + sum(a.nH for u in (head, tail, *lcombo, *rcombo) for a in u.atoms) * 1.008)
                        if approx_mw > max_mw:
                            continue
                        out.append(compose(head, rcombo, lcombo, tail))
    uniq = {}
    for sc in out:
        uniq.setdefault(sc.smiles, sc)
    return list(uniq.values())


# ---------------------------------------------------------------------------
# property envelope tuned for a periplasm-facing Gram-negative target
# ---------------------------------------------------------------------------
PROP_TARGET = dict(mw=(260, 470), clogp=(0.8, 4.5), tpsa_max=95.0,
                   hbd_max=2, hba_max=9, rot_max=6, ring_max=4)


def property_penalty(sc: Scaffold):
    pen, why = 0.0, []
    lo, hi = PROP_TARGET['mw']
    if sc.mw > hi:
        pen += (sc.mw - hi) / 40.0
        why.append('MW')
    elif sc.mw < lo:
        pen += (lo - sc.mw) / 150.0
    lo, hi = PROP_TARGET['clogp']
    if sc.clogp > hi:
        pen += 0.6 * (sc.clogp - hi)
        why.append('cLogP')
    elif sc.clogp < lo:
        pen += 0.4 * (lo - sc.clogp)
        why.append('cLogP')
    if sc.tpsa > PROP_TARGET['tpsa_max']:
        pen += (sc.tpsa - 90) / 12.0
        why.append('TPSA')
    if sc.hbd > PROP_TARGET['hbd_max']:
        pen += 0.5 * (sc.hbd - 2)
        why.append('HBD')
    if sc.hba > PROP_TARGET['hba_max']:
        pen += 0.3 * (sc.hba - 9)
        why.append('HBA')
    if sc.rot > PROP_TARGET['rot_max']:
        pen += 0.5 * (sc.rot - 6)
        why.append('RotB')
    if sc.n_ring > PROP_TARGET['ring_max']:
        pen += 0.5 * (sc.n_ring - 4)
        why.append('Rings')
    sc.components['property'] = round(pen, 3)
    sc.components['property_flags'] = ','.join(sorted(set(why)))
    return pen
