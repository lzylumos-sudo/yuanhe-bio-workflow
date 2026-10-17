"""Build an RDKit protein molecule with chemically consistent formal charges.

PDB proximity bonding gives connectivity and (for standard residues) bond
orders but not reliable charges. Charges are re-derived here from valence and
the hydrogens that are really present, so carboxylates, ammonium, guanidinium
and imidazolium are perceived correctly by the interaction SMARTS.
"""
from __future__ import annotations

from pathlib import Path

from rdkit import Chem

BACKBONE_NAMES = {"N", "CA", "C", "O", "OXT", "H", "HA", "HA2", "HA3", "H1", "H2", "H3"}


def build_protein_mol(pdb_path: Path) -> tuple[Chem.Mol, dict]:
    m = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=False, proximityBonding=True)
    if m is None:
        raise ValueError(f"RDKit could not read prepared protein {pdb_path}")
    for a in m.GetAtoms():
        a.SetFormalCharge(0)
        a.SetNoImplicit(True)
        a.SetNumRadicalElectrons(0)

    stats = {"cationic_N": 0, "anionic_O": 0, "his_states": {}}
    m = _fix_histidines(m, stats)
    for a in m.GetAtoms():
        sym = a.GetSymbol()
        val = sum(b.GetBondTypeAsDouble() for b in a.GetBonds())
        if sym == "N" and round(val) == 4:
            a.SetFormalCharge(1)
            stats["cationic_N"] += 1
        elif sym == "O" and round(val) == 1 and a.GetDegree() == 1:
            nb = a.GetNeighbors()[0]
            # carboxylate O- : single bond to a carbon that also carries a double-bonded O
            b = m.GetBondBetweenAtoms(a.GetIdx(), nb.GetIdx())
            if nb.GetSymbol() == "C" and b.GetBondType() == Chem.BondType.SINGLE and any(
                n.GetSymbol() == "O" and n.GetIdx() != a.GetIdx() and
                m.GetBondBetweenAtoms(nb.GetIdx(), n.GetIdx()).GetBondType() == Chem.BondType.DOUBLE
                for n in nb.GetNeighbors()
            ):
                a.SetFormalCharge(-1)
                stats["anionic_O"] += 1
    Chem.SanitizeMol(m, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES, catchErrors=True)
    m.UpdatePropertyCache(strict=False)
    Chem.SetAromaticity(m)
    stats["net_charge"] = Chem.GetFormalCharge(m)
    return m, stats


def _fix_histidines(m: Chem.Mol, stats: dict) -> Chem.Mol:
    """Set imidazole bond orders from the hydrogens that are really present.

    RDKit's PDB template always uses one Kekule form, which wrongly charges the
    N that carries a hydrogen in neutral HID/HIE. States: HID (HD1 only),
    HIE (HE2 only), HIP (both, +1 on the ND1=CE1 side), no H on either (deprotonated, flagged).
    """
    rw = Chem.RWMol(m)
    residues: dict[tuple, dict[str, int]] = {}
    for a in rw.GetAtoms():
        ri = a.GetPDBResidueInfo()
        if ri is not None and ri.GetResidueName().strip() in ("HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP"):
            key = (ri.GetChainId(), ri.GetResidueNumber(), ri.GetInsertionCode(), ri.GetResidueName())
            residues.setdefault(key, {})[ri.GetName().strip()] = a.GetIdx()
    S, D = Chem.BondType.SINGLE, Chem.BondType.DOUBLE
    for key, at in residues.items():
        need = ("CG", "ND1", "CD2", "CE1", "NE2")
        if not all(n in at for n in need):
            continue
        hd1 = "HD1" in at
        he2 = "HE2" in at
        state = "HIP" if hd1 and he2 else "HID" if hd1 else "HIE" if he2 else "deprotonated"
        orders = {("CG", "CD2"): D, ("CG", "ND1"): S, ("CD2", "NE2"): S, ("ND1", "CE1"): S, ("CE1", "NE2"): S}
        if state in ("HID", "deprotonated"):
            orders[("CE1", "NE2")] = D
        else:                       # HIE / HIP
            orders[("ND1", "CE1")] = D
        for (x, y), bt in orders.items():
            b = rw.GetBondBetweenAtoms(at[x], at[y])
            if b is not None:
                b.SetBondType(bt)
                b.SetIsAromatic(False)
        for n in need:
            rw.GetAtomWithIdx(at[n]).SetIsAromatic(False)
        stats["his_states"][f"{key[0]}:{key[1]}{key[2].strip()}"] = state
    return rw.GetMol()


def is_backbone(atom_name: str) -> bool:
    return atom_name.strip() in BACKBONE_NAMES
