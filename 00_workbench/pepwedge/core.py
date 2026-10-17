# -*- coding: utf-8 -*-
"""pepwedge.core -- dependency-free structure layer (numpy only).

* minimal PDB reader
* Shrake-Rupley solvent accessible surface (cell-list accelerated)
* residue contact / salt-bridge analysis
* genetics-anchored hotspot scoring

Keeping this layer free of Bio.PDB / RDKit / PyMOL means the whole design
pipeline runs on a laptop CPU, supporting local development and reproducible
structure analysis.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np

AA3 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
       'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
       'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V'}

POS_ATOMS = {'LYS': ('NZ',), 'ARG': ('NE', 'NH1', 'NH2'), 'HIS': ('ND1', 'NE2')}
NEG_ATOMS = {'ASP': ('OD1', 'OD2'), 'GLU': ('OE1', 'OE2')}

VDW = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80, 'F': 1.47,
       'CL': 1.75, 'BR': 1.85, 'I': 1.98, 'H': 1.20}
PROBE = 1.4

# ---------------------------------------------------------------------------
# Literature priors are NOT hard-coded here any more: they live in
# priors/murj_priors.json, carry their citation, and are checked against the
# coordinates by verify.py. The names below are kept for backward compatibility
# with the rest of the package.
# ---------------------------------------------------------------------------
from . import priors as _priors_mod            # noqa: E402

COMPLEXES = _priors_mod.structure_map()
ESSENTIAL = _priors_mod.essential_map()
RESISTANCE = _priors_mod.resistance_map()
SGL_FUNCTIONAL = _priors_mod.peptide_functional_map()
CONTACT_CUTOFF = float(_priors_mod.thresholds().get('contact_cutoff_A', 4.5))
SALT_CUTOFF = float(_priors_mod.thresholds().get('salt_bridge_cutoff_A', 4.0))


# ---------------------------------------------------------------------------
@dataclass
class Residue:
    num: int
    name: str
    atoms: Dict[str, np.ndarray] = field(default_factory=dict)

    @property
    def one(self) -> str:
        return AA3.get(self.name, 'X')

    def xyz(self, names=None) -> np.ndarray:
        keys = list(self.atoms) if names is None else [n for n in names if n in self.atoms]
        if not keys:
            return np.zeros((0, 3))
        return np.array([self.atoms[k] for k in keys])

    def centroid(self, names=None) -> np.ndarray:
        v = self.xyz(names)
        return v.mean(0) if len(v) else np.zeros(3)


@dataclass
class Complex:
    tag: str
    receptor: Dict[int, Residue]
    ligand: Dict[int, Residue]
    path: str = ''
    # residues of the same chain that fall OUTSIDE the target's numbering range
    # (expression tag / fusion partner). They are kept for bookkeeping but are
    # never treated as MurJ by any analysis.
    extra: Dict[int, Residue] = field(default_factory=dict)

    @property
    def r_xyz(self):
        return np.concatenate([r.xyz() for r in self.receptor.values()])

    @property
    def l_xyz(self):
        return np.concatenate([r.xyz() for r in self.ligand.values()])

    def r(self, num) -> Residue:
        return self.receptor[int(num)]

    def l(self, num) -> Residue:
        return self.ligand[int(num)]

    def seq(self, which='ligand'):
        d = self.ligand if which == 'ligand' else self.receptor
        return ''.join(d[k].one for k in sorted(d))


def read_pdb(path) -> Dict[str, Dict[int, Residue]]:
    chains: Dict[str, Dict[int, Residue]] = {}
    seen = set()
    with open(path, 'r', errors='ignore') as fh:
        for line in fh:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            alt = line[16]
            if alt not in (' ', 'A'):
                continue
            name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21]
            if resname in ('HOH', 'WAT'):
                continue
            try:
                num = int(line[22:26])
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            elem = line[76:78].strip().upper() or name[0]
            if elem == 'H' or name.startswith('H'):
                continue
            key = (chain, num, resname)
            if key in seen and name in chains.setdefault(chain, {}).setdefault(num, Residue(num, resname)).atoms:
                continue
            seen.add(key)
            chains.setdefault(chain, {}).setdefault(num, Residue(num, resname))
            chains[chain][num].atoms[name] = np.array([x, y, z], dtype=float)
    return chains


def load_complex(path, tag='complex', residue_range=None) -> Complex:
    """MurJ = the larger chain, Sgl = the other protein chain.

    Construct-only residues (negative numbering, fusion partners) are split off
    into ``extra`` so that no downstream analysis mistakes them for MurJ.
    """
    chains = read_pdb(path)
    order = sorted(chains.items(), key=lambda kv: -len(kv[1]))
    rec_raw = order[0][1]
    lig = order[1][1]
    lo, hi = residue_range or _priors_mod.residue_range()
    rec = {k: v for k, v in rec_raw.items() if lo <= k <= hi}
    extra = {k: v for k, v in rec_raw.items() if not (lo <= k <= hi)}
    return Complex(tag, rec, lig, str(path), extra)


# ---------------------------------------------------------------------------
def _cell_lists(xyz: np.ndarray, cell: float):
    lo = xyz.min(0) - cell
    idx = np.floor((xyz - lo) / cell).astype(int)
    grid = {}
    for i, k in enumerate(map(tuple, idx)):
        grid.setdefault(k, []).append(i)
    return grid, lo, cell


def _neighbors(xyz, grid, lo, cell, k, radius):
    rng = int(np.ceil(radius / cell))
    out = []
    kx, ky, kz = k
    for dx in range(-rng, rng + 1):
        for dy in range(-rng, rng + 1):
            for dz in range(-rng, rng + 1):
                lst = grid.get((kx + dx, ky + dy, kz + dz))
                if lst:
                    out.extend(lst)
    return np.array(out, dtype=int)


_SPHERE = None


def _sphere_points(n=92):
    """Golden-spiral sphere."""
    global _SPHERE
    if _SPHERE is not None and len(_SPHERE) == n:
        return _SPHERE
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    pts = np.stack([np.cos(theta) * np.sin(phi),
                    np.sin(theta) * np.sin(phi),
                    np.cos(phi)], axis=1)
    _SPHERE = pts
    return pts


def sasa_per_atom(xyz: np.ndarray, elements: List[str], n_points=92):
    """Shrake-Rupley SASA (A^2) per atom."""
    rad = np.array([VDW.get(e.upper(), 1.8) + PROBE for e in elements])
    cell = float(2 * rad.max() + 1.0)
    grid, lo, cell = _cell_lists(xyz, cell)
    pts = _sphere_points(n_points)
    out = np.zeros(len(xyz))
    for i in range(len(xyz)):
        k = tuple(np.floor((xyz[i] - lo) / cell).astype(int))
        nb = _neighbors(xyz, grid, lo, cell, k, 2 * rad.max() + 1.0)
        nb = nb[nb != i]
        if len(nb) == 0:
            out[i] = 4 * np.pi * rad[i] ** 2
            continue
        test = xyz[i] + rad[i] * pts
        d = np.linalg.norm(test[:, None, :] - xyz[nb][None, :, :], axis=2)
        acc = np.all(d > rad[nb][None, :], axis=1)
        out[i] = 4 * np.pi * rad[i] ** 2 * acc.mean()
    return out


def _flatten(res_dict, prefix='R'):
    atoms, elems, owner = [], [], []
    for rid, r in res_dict.items():
        for an, xyz in r.atoms.items():
            atoms.append(xyz)
            elems.append(an[0] if an[0].isalpha() else 'C')
            owner.append((prefix, rid))
    return np.array(atoms), elems, owner


def residue_sasa(cx: Complex):
    """Per-residue SASA of the isolated MurJ and of MurJ in the complex."""
    xyz, elems, owner = _flatten(cx.receptor, 'R')
    sasa = sasa_per_atom(xyz, elems)
    free = {}
    for (_, rid), s in zip(owner, sasa):
        free[rid] = free.get(rid, 0.0) + s

    xa, ea, oa = _flatten(cx.receptor, 'R')
    xb, eb, ob = _flatten(cx.ligand, 'L')
    xyz2 = np.vstack([xa, xb])
    elems2 = ea + eb
    owner2 = oa + ob
    sasa2 = sasa_per_atom(xyz2, elems2)
    bound = {}
    for (tag, rid), s in zip(owner2, sasa2):
        if tag == 'R':
            bound[rid] = bound.get(rid, 0.0) + s
    dsasa = {k: max(0.0, free.get(k, 0) - bound.get(k, 0)) for k in free}
    return free, bound, dsasa


# ---------------------------------------------------------------------------
def contacts(cx: Complex, cutoff=4.5):
    """{(murj_resid, sgl_resid): {'dmin', 'n', 'pairs'}}"""
    ra = cx.r_xyz
    rlab = [(rid, an) for rid, r in cx.receptor.items() for an in r.atoms]
    out = {}
    grid, lo, cell = _cell_lists(ra, 8.0)
    for lid, l in cx.ligand.items():
        for an, xyz in l.atoms.items():
            nb = _neighbors(ra, grid, lo, cell, tuple(np.floor((xyz - lo) / cell).astype(int)), cutoff)
            for j in nb:
                d = float(np.linalg.norm(ra[j] - xyz))
                if d >= cutoff:
                    continue
                key = (int(rlab[j][0]), lid)
                rec = out.setdefault(key, {'dmin': 9.9, 'n': 0, 'pairs': []})
                rec['dmin'] = min(rec['dmin'], d)
                rec['n'] += 1
                rec['pairs'].append((rlab[j][1], an, round(d, 2)))
    return out


def salt_bridges(cx: Complex, cutoff=4.0):
    out = []
    for rid, r in cx.receptor.items():
        for an in POS_ATOMS.get(r.name, ()):
            if an not in r.atoms:
                continue
            for lid, l in cx.ligand.items():
                for bn in NEG_ATOMS.get(l.name, ()):
                    if bn not in l.atoms:
                        continue
                    d = float(np.linalg.norm(r.atoms[an] - l.atoms[bn]))
                    if d < cutoff:
                        out.append(dict(murj=f'{r.one}{rid}', murj_atom=an,
                                        sgl=f'{l.one}{lid}', sgl_atom=bn, dist=round(d, 2)))
    last = cx.ligand[max(cx.ligand)]
    term_o = [n for n in ('O', 'OXT') if n in last.atoms]
    if 'OXT' in last.atoms:
        for rid, r in cx.receptor.items():
            for an in POS_ATOMS.get(r.name, ()):
                if an not in r.atoms:
                    continue
                for on in term_o:
                    d = float(np.linalg.norm(r.atoms[an] - last.atoms[on]))
                    if d < cutoff:
                        out.append(dict(murj=f'{r.one}{rid}', murj_atom=an,
                                        sgl=f'{last.one}{last.num}(C-term)', sgl_atom=on,
                                        dist=round(d, 2)))
    return out


def membrane_frame(cx: Complex):
    ca = np.array([r.atoms['CA'] for r in cx.receptor.values() if 'CA' in r.atoms])
    c0 = ca.mean(0)
    _, _, vt = np.linalg.svd(ca - c0, full_matrices=False)
    n = vt[0] / np.linalg.norm(vt[0])
    last = cx.ligand[max(cx.ligand)]
    cterm = last.centroid([k for k in ('C', 'O', 'OXT') if k in last.atoms])
    if np.dot(cterm - c0, n) < 0:
        n = -n
    return n, c0


def kabsch(mobile: np.ndarray, target: np.ndarray):
    """Return (R, t) mapping mobile -> target."""
    mc, tc = mobile.mean(0), target.mean(0)
    H = (mobile - mc).T @ (target - tc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = tc - R @ mc
    return R, t


def apply_rt(xyz, rt):
    R, t = rt
    return np.asarray(xyz) @ R.T + t


# ---------------------------------------------------------------------------
# genetics-anchored hotspot scoring
# ---------------------------------------------------------------------------
def hotspot_table(cx: Complex, dsasa=None, conv_sets=None, weights=(0.30, 0.25, 0.15, 0.20, 0.10)):
    """score = w1*n_peptides(convergence) + w2*burial + w3*contact_density
             + w4*essentiality + w5*resistance_evidence"""
    ct = contacts(cx)
    if dsasa is None:
        dsasa = residue_sasa(cx)[2]
    per = {}
    for (rid, lid), rec in ct.items():
        d = per.setdefault(rid, {'n': 0, 'dmin': 9.9, 'sgl': set()})
        d['n'] += rec['n']
        d['dmin'] = min(d['dmin'], rec['dmin'])
        d['sgl'].add(lid)
    maxb = max(dsasa.values()) if dsasa else 1.0
    maxn = max((v['n'] for v in per.values()), default=1)
    w1, w2, w3, w4, w5 = weights
    rows = []
    for rid, d in per.items():
        r = cx.r(rid)
        conv = 1.0 if conv_sets is None else sum(rid in s for s in conv_sets) / len(conv_sets)
        bur = dsasa.get(rid, 0.0) / maxb if maxb else 0.0
        dens = d['n'] / maxn
        ess = 1.0 if rid in ESSENTIAL else 0.0
        res = 1.0 if rid in RESISTANCE else 0.0
        score = w1 * conv + w2 * bur + w3 * dens + w4 * ess + w5 * res
        rows.append(dict(resid=rid, resname=r.name, one=r.one, score=round(score, 3),
                         convergence=round(conv, 2), dSASA=round(dsasa.get(rid, 0.0), 1),
                         n_contacts=d['n'], n_peptide_residues=len(d['sgl']),
                         dmin=round(d['dmin'], 2),
                         evidence='; '.join(x for x in (ESSENTIAL.get(rid), RESISTANCE.get(rid)) if x)))
    rows.sort(key=lambda x: (-x['score'], x['resid']))
    return rows
