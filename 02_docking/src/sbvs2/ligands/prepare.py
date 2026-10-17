"""Ligand states, 3D conformers and charges (plan §13-14).

State generation: MolScrub is used when importable; otherwise a *conservative rule fallback* is used and
every state is tagged method='rdkit_rules_v0' so reviewers see it is not MolScrub output.
"""
from __future__ import annotations
import math
from dataclasses import dataclass, field
from rdkit import Chem
from rdkit.Chem import AllChem, rdPartialCharges
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

from ..hashing import hash_obj

try:                                         # optional; the importable module is `molscrub` (pip name molscrub)
    from molscrub import Scrub as _Scrub
    import molscrub as _molscrub
    HAVE_MOLSCRUB = True
    MOLSCRUB_VERSION = getattr(_molscrub, "__version__", "0.2.2")
except Exception:
    HAVE_MOLSCRUB = False
    MOLSCRUB_VERSION = None

_RULES = [  # (name, SMARTS of neutral group, pKa-ish midpoint, builder)
    ("carboxylic_acid", "[CX3](=O)[OX2H1]", 4.5, "acid"),
    ("aliphatic_amine", "[NX3;H2,H1,H0;!$(NC=[O,S,N]);!$(N-a);!$(N-[SX4]);!$(N#*);!$(N=*);!$(N-N);!$(N-O)]", 9.5, "base"),
]


def enumerate_stereo(mol, max_n=8):
    opts = StereoEnumerationOptions(onlyUnassigned=True, unique=True, maxIsomers=max_n)   # known centers kept
    isos = list(EnumerateStereoisomers(mol, options=opts))
    truncated = len(isos) >= max_n
    return isos, truncated


def _states_molscrub(mol, ph_lo, ph_hi, max_states):
    """pKa-model based protonation/tautomer states (MolScrub, GPL-3 tool, run as a separate pip dependency).
    Returns [(mol_no_Hs, why, method, smiles)] or None when MolScrub cannot process the molecule."""
    try:
        sc = _Scrub(ph_low=ph_lo, ph_high=ph_hi, skip_gen3d=True, skip_ringfix=True)
        outs = [Chem.RemoveHs(x) for x in sc(Chem.Mol(mol))]
    except Exception:
        return None
    uniq, seen = [], set()
    for m in outs:
        smi = Chem.MolToSmiles(m)
        if smi not in seen:
            seen.add(smi); uniq.append((m, f"molscrub_pH[{ph_lo},{ph_hi}]", f"molscrub_{MOLSCRUB_VERSION}", smi))
    return (uniq[:max_states], len(uniq) > max_states) if uniq else None


def _states_rules(mol, ph_lo, ph_hi, max_states):
    """Always include the neutral form; add ionised form when pKa midpoint lies below/above the pH window."""
    states = [(mol, "neutral_input", "rdkit_rules_v0")]
    rw = Chem.RWMol(mol)
    changed = False
    for name, sma, pka, kind in _RULES:
        patt = Chem.MolFromSmarts(sma)
        for match in rw.GetSubstructMatches(patt):
            if kind == "acid" and ph_hi > pka - 1.0:
                at = rw.GetAtomWithIdx(match[2]); at.SetFormalCharge(-1); at.SetNumExplicitHs(0); changed = True
            elif kind == "base" and ph_lo < pka + 1.0:
                at = rw.GetAtomWithIdx(match[0])
                if at.GetFormalCharge() == 0 and at.GetTotalNumHs() + at.GetDegree() < 4:
                    at.SetFormalCharge(1); at.SetNumExplicitHs(at.GetTotalNumHs() + 1); changed = True
    if changed:
        try:
            m2 = rw.GetMol(); Chem.SanitizeMol(m2)
            states.insert(0, (m2, f"rule_ionised_pH[{ph_lo},{ph_hi}]", "rdkit_rules_v0"))
        except Exception:
            pass
    uniq, seen = [], set()
    for m, why, meth in states:
        s = Chem.MolToSmiles(m)
        if s not in seen:
            seen.add(s); uniq.append((m, why, meth, s))
    return uniq[:max_states], len(uniq) > max_states


@dataclass
class PrepSettings:
    ph_lo: float = 6.4
    ph_hi: float = 8.4
    ph_confirmed: bool = False
    max_stereo: int = 8
    max_states: int = 8
    n_conf_generate: int = 20
    n_conf_keep: int = 3
    charge_model: str = "gasteiger"
    seed: int = 20261012
    unspecified_stereo: str = "review_or_enumerate"
    state_method: str = "auto"        # auto (MolScrub if importable, else rules_v0) | rules

    def policy_hash(self) -> str:
        return hash_obj({**self.__dict__, "molscrub": MOLSCRUB_VERSION})


def make_conformers(mol, n_gen, n_keep, seed, rms=0.5):
    mh = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3(); ps.randomSeed = seed
    cids = list(AllChem.EmbedMultipleConfs(mh, n_gen, ps))
    if not cids:
        return None, "EMBED_FAILED", {}
    ff = "MMFF"
    if AllChem.MMFFHasAllMoleculeParams(mh):
        res = AllChem.MMFFOptimizeMoleculeConfs(mh, maxIters=2000)
    else:
        ff = "UFF_FALLBACK"; res = AllChem.UFFOptimizeMoleculeConfs(mh, maxIters=2000)
    conv = {c: (r[0] == 0) for c, r in zip(cids, res)}
    energies = {c: r[1] for c, r in zip(cids, res)}
    order = sorted(cids, key=lambda c: energies[c])
    keep = []
    for c in order:
        if not conv[c]:
            continue
        if all(AllChem.GetConformerRMS(mh, k, c, prealigned=False) > rms for k in keep):
            keep.append(c)
        if len(keep) >= n_keep:
            break
    info = {"forcefield": ff, "generated": len(cids), "kept": keep, "energies": {c: energies[c] for c in keep},
            "converged": {c: conv[c] for c in keep}}
    return (mh, keep), ("OK" if keep else "NOT_CONVERGED"), info


def gasteiger(mh):
    rdPartialCharges.ComputeGasteigerCharges(mh)
    q = [float(a.GetDoubleProp("_GasteigerCharge")) for a in mh.GetAtoms()]
    return q


def charge_check(mh, q, tol=1e-3):
    """Plan §14: finite, complete map, sum equals integer formal charge within 1e-3 e."""
    target = sum(a.GetFormalCharge() for a in mh.GetAtoms())
    ok_fin = len(q) == mh.GetNumAtoms() and all(math.isfinite(x) for x in q)
    diff = abs(sum(q) - target) if ok_fin else float("nan")
    return {"finite_complete": ok_fin, "sum": sum(q) if ok_fin else None, "target": target, "abs_diff": diff,
            "pass": bool(ok_fin and diff <= tol)}


def prepare_parent(parent: dict, s: PrepSettings):
    """Yield state records for one parent. Never silently truncates; truncation is recorded."""
    base = Chem.MolFromSmiles(parent["isomeric_smiles"])
    out = []
    if parent["unspecified_stereocenters"] and s.unspecified_stereo == "review_or_enumerate":
        stereos, trunc_st = enumerate_stereo(base, s.max_stereo)
        stereo_note = "enumerated_unknown_centers"
    else:
        stereos, trunc_st, stereo_note = [base], False, "kept_as_given"
    for si, sm in enumerate(stereos, 1):
        stereo_id = f"{parent['parent_id']}_S{si}"
        got = _states_molscrub(sm, s.ph_lo, s.ph_hi, s.max_states) if (HAVE_MOLSCRUB and s.state_method != "rules") else None
        states, trunc_state = got if got else _states_rules(sm, s.ph_lo, s.ph_hi, s.max_states)
        for ti, (m, why, meth, smi) in enumerate(states, 1):
            sid = f"{stereo_id}_T{ti}"
            rec = {"parent_id": parent["parent_id"], "stereo_id": stereo_id, "state_id": sid, "smiles": smi,
                   "formal_charge": Chem.GetFormalCharge(m), "evidence": why, "method": meth, "stereo_note": stereo_note,
                   "truncated": bool(trunc_st or trunc_state),
                   "truncation_reason": "enumeration_cap_reached" if (trunc_st or trunc_state) else "",
                   "ph_range": [s.ph_lo, s.ph_hi], "ph_confirmed": s.ph_confirmed, "status": "PENDING_REVIEW"}
            built, code, info = make_conformers(m, s.n_conf_generate, s.n_conf_keep, s.seed)
            rec["conformer_status"] = code; rec["conformer_info"] = info
            if built:
                mh, keep = built
                q = gasteiger(mh)
                rec["charges"] = {"method": "gasteiger", "values": q, "check": charge_check(mh, q)}
                rec["_mol"] = mh; rec["_keep"] = keep
                if not rec["charges"]["check"]["pass"]:
                    rec["status"] = "BLOCKED"; rec["block_reason"] = "CHARGE_ERROR"
            else:
                rec["status"] = "BLOCKED"; rec["block_reason"] = code
            rec["state_hash"] = hash_obj({k: rec[k] for k in ("smiles", "formal_charge", "method")} | {"seed": s.seed, "policy": s.policy_hash()})
            out.append(rec)
    return out


def estimate_jobs(n_states: int, n_conf: int, n_receptors: int, n_seeds: int, max_poses=9):
    jobs = n_states * n_conf * n_receptors * n_seeds
    return {"n_jobs": jobs, "max_poses": jobs * max_poses}
