# -*- coding: utf-8 -*-
"""pepwedge.evidence -- every structural number the design rests on, measured.

Nothing in this module comes from the literature. Each function recomputes a
quantity from the deposited coordinates so that a reviewer (or a collaborator
who has never seen this repository) can re-derive it:

    census()             chains, residues, atoms, gaps, element composition
    dsasa_both()         buried surface on BOTH sides of the interface
    interface_rows()     per MurJ residue: contacts, dSASA, convergence, anchors
    peptide_rows()       per peptide residue: contacts, buried area, role
    interface_census()   the headline numbers of the peptide-protein interface
    anchor_table()       salt bridges + distances between the four basic anchors
    atom_budget()        peptide heavy atoms vs small-molecule heavy atoms

All angles/distances are in Angstrom, all areas in Angstrom^2.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

from . import priors as _priors
from .core import (NEG_ATOMS, POS_ATOMS, VDW, Complex, _flatten, contacts,
                   sasa_per_atom, salt_bridges)

CONTACT_CUTOFF = 4.5
SALT_CUTOFF = 4.0


# ---------------------------------------------------------------------------
# 1. census
# ---------------------------------------------------------------------------
def _gaps(nums: List[int]) -> List[str]:
    out = []
    s = sorted(nums)
    for a, b in zip(s, s[1:]):
        if b - a > 1:
            out.append(f'{a + 1}-{b - 1}' if b - a > 2 else f'{a + 1}')
    return out


def census(cx: Complex) -> dict:
    rec, lig = cx.receptor, cx.ligand
    r_atoms = sum(len(r.atoms) for r in rec.values())
    l_atoms = sum(len(l.atoms) for l in lig.values())
    elems: Dict[str, int] = {}
    for r in rec.values():
        for an in r.atoms:
            e = an[0] if an[0].isalpha() else 'C'
            elems[e] = elems.get(e, 0) + 1
    l_elems: Dict[str, int] = {}
    for r in lig.values():
        for an in r.atoms:
            e = an[0] if an[0].isalpha() else 'C'
            l_elems[e] = l_elems.get(e, 0) + 1
    ex = cx.extra
    return dict(
        tag=cx.tag,
        construct_residues=len(ex),
        construct_atoms=sum(len(r.atoms) for r in ex.values()),
        construct_range=(f'{min(ex)}…{max(ex)}' if ex else 'none'),
        murj_residues=len(rec), murj_atoms=r_atoms,
        murj_range=f'{min(rec)}-{max(rec)}',
        murj_gaps=';'.join(_gaps(list(rec))) or 'none',
        peptide_residues=len(lig), peptide_atoms=l_atoms,
        peptide_range=f'{min(lig)}-{max(lig)}',
        peptide_gaps=';'.join(_gaps(list(lig))) or 'none',
        peptide_sequence=cx.seq('ligand'),
        contacts_heavy_atoms=l_atoms,
    )


# ---------------------------------------------------------------------------
# 2. buried area on both sides
# ---------------------------------------------------------------------------
def dsasa_both(cx: Complex):
    """(rec_dsasa, lig_dsasa, rec_free_sasa, lig_free_sasa).

    dSASA = SASA(isolated chain) - SASA(chain inside the complex), i.e. the
    interface area contributed by each residue. The two totals should agree to
    within a few percent; a large disagreement means the probe/grid is too
    coarse for the structure.
    """
    xa, ea, oa = _flatten(cx.receptor, 'R')
    xb, eb, ob = _flatten(cx.ligand, 'L')
    sa = sasa_per_atom(xa, ea)
    sb = sasa_per_atom(xb, eb)
    sab = sasa_per_atom(np.vstack([xa, xb]), ea + eb)
    free_r, free_l, bound_r, bound_l = {}, {}, {}, {}
    for (_, rid), v in zip(oa, sa):
        free_r[rid] = free_r.get(rid, 0.0) + v
    for (_, lid), v in zip(ob, sb):
        free_l[lid] = free_l.get(lid, 0.0) + v
    for (tag, rid), v in zip(oa + ob, sab):
        d = bound_r if tag == 'R' else bound_l
        d[rid] = d.get(rid, 0.0) + v
    dsasa_r = {k: max(0.0, free_r.get(k, 0.0) - bound_r.get(k, 0.0)) for k in free_r}
    dsasa_l = {k: max(0.0, free_l.get(k, 0.0) - bound_l.get(k, 0.0)) for k in free_l}
    return dsasa_r, dsasa_l, sum(free_r.values()), sum(free_l.values())


# ---------------------------------------------------------------------------
# 3. per-residue interface tables
# ---------------------------------------------------------------------------
def interface_rows(cx: Complex, dsasa_r=None, conv_sets=None) -> List[dict]:
    """One row per MurJ residue that touches the peptide."""
    if dsasa_r is None:
        dsasa_r = dsasa_both(cx)[0]
    ct = contacts(cx, CONTACT_CUTOFF)
    per: Dict[int, dict] = {}
    for (rid, lid), rec in ct.items():
        d = per.setdefault(rid, {'n': 0, 'dmin': 9.9, 'pep': set()})
        d['n'] += rec['n']
        d['dmin'] = min(d['dmin'], rec['dmin'])
        d['pep'].add(lid)
    rows = []
    for rid, d in per.items():
        r = cx.r(rid)
        conv = (sum(rid in s for s in conv_sets) / len(conv_sets)) if conv_sets else None
        rows.append(dict(
            complex=cx.tag, resid=rid, resname=r.name, one=r.one,
            n_atom_contacts=d['n'], n_peptide_residues=len(d['pep']),
            dmin_A=round(d['dmin'], 2), dSASA_A2=round(dsasa_r.get(rid, 0.0), 1),
            convergence=None if conv is None else round(conv, 2),
            ca=[round(float(x), 2) for x in r.atoms.get('CA', r.centroid())],
        ))
    rows.sort(key=lambda x: (-x['n_atom_contacts'], x['resid']))
    return rows


def peptide_rows(cx: Complex, dsasa_l=None) -> List[dict]:
    if dsasa_l is None:
        dsasa_l = dsasa_both(cx)[1]
    ct = contacts(cx, CONTACT_CUTOFF)
    per: Dict[int, dict] = {}
    for (rid, lid), rec in ct.items():
        d = per.setdefault(lid, {'n': 0, 'dmin': 9.9, 'murj': set()})
        d['n'] += rec['n']
        d['dmin'] = min(d['dmin'], rec['dmin'])
        d['murj'].add(rid)
    rows = []
    for lid in sorted(cx.ligand):
        res = cx.l(lid)
        d = per.get(lid, {'n': 0, 'dmin': None, 'murj': set()})
        rows.append(dict(
            complex=cx.tag, pep_resid=lid, resname=res.name, one=res.one,
            n_atom_contacts=d['n'], n_murj_residues=len(d['murj']),
            dmin_A=None if d['dmin'] is None else round(d['dmin'], 2),
            dSASA_A2=round(dsasa_l.get(lid, 0.0), 1),
            heavy_atoms=len(res.atoms),
            buried=bool(dsasa_l.get(lid, 0.0) > 1.0),
        ))
    return rows


# ---------------------------------------------------------------------------
# 4. headline interface numbers
# ---------------------------------------------------------------------------
def interface_census(cx: Complex, dsasa=None) -> dict:
    if dsasa is None:
        dsasa = dsasa_both(cx)
    dsasa_r, dsasa_l, free_r, free_l = dsasa
    pr = peptide_rows(cx, dsasa_l)
    ct = contacts(cx, CONTACT_CUTOFF)
    touched_r = {rid for (rid, _) in ct}
    touched_l = {lid for (_, lid) in ct}

    # heavy atoms that are actually within the cutoff of the partner
    xa, _, _ = _flatten(cx.receptor, 'R')
    xb, _, _ = _flatten(cx.ligand, 'L')
    d = np.linalg.norm(xb[:, None, :] - xa[None, :, :], axis=2)
    near = d < CONTACT_CUTOFF
    pep_atoms_in_contact = int(near.any(1).sum())
    murj_atoms_in_contact = int(near.any(0).sum())

    return dict(
        tag=cx.tag,
        murj_residues_in_interface=len(touched_r),
        peptide_residues_in_interface=len(touched_l),
        peptide_residues_modelled=len(cx.ligand),
        fraction_of_peptide_at_interface=round(len(touched_l) / len(cx.ligand), 2),
        atom_pairs=int(near.sum()),
        peptide_heavy_atoms_in_contact=pep_atoms_in_contact,
        murj_heavy_atoms_in_contact=murj_atoms_in_contact,
        buried_area_murj_side_A2=round(sum(dsasa_r.values()), 1),
        buried_area_peptide_side_A2=round(sum(dsasa_l.values()), 1),
        asymmetry=round(abs(sum(dsasa_r.values()) - sum(dsasa_l.values()))
                        / max(1.0, sum(dsasa_r.values())), 3),
    )


# ---------------------------------------------------------------------------
# 5. anchors
# ---------------------------------------------------------------------------
_ANCHORS = [46, 53, 270, 368]


def _basic_atoms(res):
    return [res.atoms[a] for a in POS_ATOMS.get(res.name, ()) if a in res.atoms]


def _acidic_atoms(cx: Complex, lid: int):
    """Carboxylate oxygens of a peptide residue.

    For a C-terminal residue the carboxylate is made of the backbone O and OXT,
    so BOTH must be considered; using only OXT under-measures the interaction
    by about 2 A and produces spurious "no salt bridge" verdicts.
    """
    l = cx.l(lid)
    out = [l.atoms[a] for a in NEG_ATOMS.get(l.name, ()) if a in l.atoms]
    if lid == max(cx.ligand) and 'OXT' in l.atoms:
        for n in ('O', 'OXT'):
            if n in l.atoms:
                out.append(l.atoms[n])
    return out


def anchor_table(cx: Complex) -> dict:
    """Salt bridges + pairwise geometry of the four basic anchors."""
    sb = [dict(r, complex=cx.tag) for r in salt_bridges(cx, SALT_CUTOFF)]
    present = [i for i in _ANCHORS if i in cx.receptor]
    pairs = []
    for i, a in enumerate(present):
        for b in present[i + 1:]:
            ra, rb = cx.r(a), cx.r(b)
            pa, pb = _basic_atoms(ra), _basic_atoms(rb)
            if not pa or not pb:
                continue
            dmin = min(float(np.linalg.norm(x - y)) for x in pa for y in pb)
            dsc = float(np.linalg.norm(np.mean(pa, 0) - np.mean(pb, 0)))
            pairs.append(dict(complex=cx.tag, a=f'{ra.one}{a}', b=f'{rb.one}{b}',
                              sidechain_centroid_A=round(dsc, 2),
                              min_atom_A=round(dmin, 2)))
    return dict(salt_bridges=sb, anchor_pairs=pairs,
                anchors_present=[f'{cx.r(i).one}{i}' for i in present])


# ---------------------------------------------------------------------------
# 6. atom budget: why a peptide is not a scaled-up small molecule
# ---------------------------------------------------------------------------
def atom_budget(complexes: Dict[str, Complex], dsasa_by_tag=None,
                n_design: int = 4000) -> dict:
    """Compare the heavy-atom budget of the peptide interface with that of the
    scaffolds this pipeline can actually build."""
    per_complex = {}
    for tag, cx in complexes.items():
        dsasa = dsasa_by_tag[tag] if dsasa_by_tag else dsasa_both(cx)
        ic = interface_census(cx, dsasa)
        per_complex[tag] = dict(
            peptide_heavy_atoms=sum(len(l.atoms) for l in cx.ligand.values()),
            peptide_heavy_atoms_in_contact=ic['peptide_heavy_atoms_in_contact'],
            peptide_residues_in_interface=ic['peptide_residues_in_interface'],
            buried_area_peptide_side_A2=ic['buried_area_peptide_side_A2'],
        )

    # small-molecule side: enumerate what the design stage can produce
    from .design import enumerate_scaffolds
    counts = []
    for D in (12.6, 13.3, 14.1, 16.5):
        for sc in enumerate_scaffolds(D, tol=0.55, ring_counts=(2, 3)):
            counts.append(len(sc.atoms))
            if len(counts) >= n_design:
                break
    counts = np.array(counts) if counts else np.array([0])
    sm = dict(n_scaffolds=int(len(counts)),
              heavy_atoms_min=int(counts.min()),
              heavy_atoms_median=float(np.median(counts)),
              heavy_atoms_max=int(counts.max()))

    pep = max(v['peptide_heavy_atoms_in_contact'] for v in per_complex.values())
    ratio = pep / max(1.0, sm['heavy_atoms_median'])
    return dict(per_complex=per_complex, small_molecule=sm,
                median_ratio=round(ratio, 1))
