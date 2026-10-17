"""Symmetry-corrected heavy-atom RMSD in a shared frame, no re-alignment (plan §19 L-class)."""
from __future__ import annotations
import numpy as np
from rdkit import Chem


def _skeleton(mol):
    """Heavy-atom graph with bond orders/charges/aromaticity erased -> matches across protonation states."""
    m = Chem.RemoveHs(mol, sanitize=False)
    rw = Chem.RWMol(m)
    for a in rw.GetAtoms():
        a.SetFormalCharge(0); a.SetIsAromatic(False); a.SetNoImplicit(True); a.SetNumExplicitHs(0)
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
    return rw.GetMol()


def symmetric_rmsd(pose, ref, conf_pose=-1, conf_ref=-1) -> float:
    sp, sr = _skeleton(pose), _skeleton(ref)
    if sp.GetNumAtoms() != sr.GetNumAtoms():
        raise ValueError("HEAVY_ATOM_COUNT_MISMATCH")
    P = Chem.RemoveHs(pose, sanitize=False).GetConformer(conf_pose).GetPositions()
    R = Chem.RemoveHs(ref, sanitize=False).GetConformer(conf_ref).GetPositions()
    matches = sr.GetSubstructMatches(sp, uniquify=False, maxMatches=100000)
    if not matches:
        raise ValueError("NO_ATOM_MAPPING")
    el_p = [a.GetSymbol() for a in sp.GetAtoms()]; el_r = [a.GetSymbol() for a in sr.GetAtoms()]
    best = None
    for mt in matches:                       # mt[i] = ref index of pose atom i
        if any(el_p[i] != el_r[j] for i, j in enumerate(mt)):
            continue
        d = float(np.sqrt(np.mean(np.sum((P - R[list(mt)]) ** 2, axis=1))))
        best = d if best is None or d < best else best
    if best is None:
        raise ValueError("NO_ELEMENT_CONSISTENT_MAPPING")
    return best
