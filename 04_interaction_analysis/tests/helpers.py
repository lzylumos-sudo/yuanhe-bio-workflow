"""Synthetic geometry for per-interaction threshold tests.

Molecules are built from SMILES, embedded once, then rigidly moved so that the
tested distance/angle is exact. Everything goes through zlis.interactions.run_analysis,
the same code path as real structures.
"""
from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from zlis.chemistry import LigandChem
from zlis.interactions import ProteinIndex, run_analysis
from zlis.structure import Atom, ResidueGroup


def embed(smiles: str, seed: int = 7) -> Chem.Mol:
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(m, randomSeed=seed) == 0
    AllChem.MMFFOptimizeMolecule(m)
    return m


def xyz(m: Chem.Mol) -> np.ndarray:
    return np.array(m.GetConformer().GetPositions())


def set_xyz(m: Chem.Mol, pts: np.ndarray) -> None:
    c = m.GetConformer()
    for i, p in enumerate(pts):
        c.SetAtomPosition(i, [float(v) for v in p])


def rot_between(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a, b = a / np.linalg.norm(a), b / np.linalg.norm(b)
    v = np.cross(a, b)
    s, c = np.linalg.norm(v), float(np.dot(a, b))
    if s < 1e-9:
        if c > 0:
            return np.eye(3)
        ortho = np.cross(a, [1, 0, 0]) if abs(a[0]) < 0.9 else np.cross(a, [0, 1, 0])
        ortho /= np.linalg.norm(ortho)
        return 2 * np.outer(ortho, ortho) - np.eye(3)
    k = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + k + k @ k * ((1 - c) / s ** 2)


def orient(m: Chem.Mol, origin_idx: int, axis_idx: int, axis_dir=(1, 0, 0)) -> None:
    """Put atom origin_idx at (0,0,0) with bond origin->axis along axis_dir."""
    p = xyz(m) - xyz(m)[origin_idx]
    R = rot_between(p[axis_idx], np.array(axis_dir, float))
    set_xyz(m, p @ R.T)


def place(m: Chem.Mol, idx: int, target) -> None:
    set_xyz(m, xyz(m) - xyz(m)[idx] + np.array(target, float))


def name_atoms(m: Chem.Mol, resname: str, chain: str, resseq: int, names: dict[int, str], hetero: bool) -> None:
    counters: dict[str, int] = {}
    for a in m.GetAtoms():
        nm = names.get(a.GetIdx())
        if nm is None:
            el = a.GetSymbol()
            counters[el] = counters.get(el, 0) + 1
            nm = f"{el}{counters[el]}X"
        ri = Chem.AtomPDBResidueInfo()
        ri.SetName(f"{nm:<4s}")
        ri.SetResidueName(resname)
        ri.SetResidueNumber(resseq)
        ri.SetChainId(chain)
        ri.SetIsHeteroAtom(hetero)
        a.SetMonomerInfo(ri)


def make_ligand(m: Chem.Mol, names: dict[int, str]) -> LigandChem:
    name_atoms(m, "LIG", "B", 1, names, True)
    infos = []
    for a in m.GetAtoms():
        nm = a.GetPDBResidueInfo().GetName().strip()
        infos.append({"uid": f"M1:B:LIG1:{nm}", "serial": a.GetIdx() + 1, "name": nm, "element": a.GetSymbol(), "source": "input"})
    return LigandChem(mol=m, atoms=infos, chem_source="template", status="ok",
                      formal_charge=Chem.GetFormalCharge(m), smiles=Chem.MolToSmiles(Chem.RemoveHs(m)))


def make_protein(m: Chem.Mol, resname: str, names: dict[int, str]) -> ProteinIndex:
    name_atoms(m, resname, "A", 1, names, False)
    uid, nm, res = {}, {}, {}
    for a in m.GetAtoms():
        n = a.GetPDBResidueInfo().GetName().strip()
        uid[a.GetIdx()] = f"M1:A:{resname}1:{n}"
        nm[a.GetIdx()] = n
        res[a.GetIdx()] = {"chain": "A", "resname": resname, "resseq": 1, "icode": "", "label": f"A:{resname}1"}
    return ProteinIndex(m, uid, nm, res)


def atom(name, el, x, y, z, resname="HOH", chain="A", resseq=500):
    return Atom(serial=1, name=name, altloc="", resname=resname, chain=chain, resseq=resseq, icode="", x=x, y=y, z=z,
                occupancy=1.0, bfactor=0.0, element=el, record="HETATM", model=1)


def analyse(lig: LigandChem, prot: ProteinIndex, waters=(), ions=()) -> list[dict]:
    w = [ResidueGroup(key=a.res_key, atoms=[a]) for a in waters]
    i = [ResidueGroup(key=a.res_key, atoms=[a]) for a in ions]
    return run_analysis(lig, prot, [], w, i, 1, "B:LIG1")


def types(events):
    return sorted({e["type"] for e in events})


# ---- scenario builders -------------------------------------------------

def hbond_ligand_donor(dist: float, dha: float = 180.0) -> list[dict]:
    """Ligand N-H donates to protein acetamide O (neutral)."""
    lig = embed("CN")
    n = [a.GetIdx() for a in lig.GetAtoms() if a.GetSymbol() == "N"][0]
    h = [x.GetIdx() for x in lig.GetAtomWithIdx(n).GetNeighbors() if x.GetSymbol() == "H"][0]
    orient(lig, n, h)
    hh = float(np.linalg.norm(xyz(lig)[h]))
    th = np.radians(dha)
    r = hh * np.cos(th) + np.sqrt((hh * np.cos(th)) ** 2 - hh ** 2 + dist ** 2)
    target = np.array([hh - r * np.cos(th), r * np.sin(th), 0.0])
    pro = embed("CC(C)=O")
    o = [a.GetIdx() for a in pro.GetAtoms() if a.GetSymbol() == "O"][0]
    place(pro, o, target)
    # swing the protein molecule away from the ligand so nothing else is in contact
    c = [a.GetIdx() for a in pro.GetAtoms() if a.GetSymbol() == "C"][0]
    L = make_ligand(lig, {n: "N1", h: "H1"})
    P = make_protein(pro, "ASN", {o: "OD1", c: "CG"})
    return analyse(L, P)


def hbond_protein_donor(dist: float, dha: float = 180.0) -> list[dict]:
    """Protein side-chain N-H donates to ligand carbonyl O."""
    pro = embed("CN")
    n = [a.GetIdx() for a in pro.GetAtoms() if a.GetSymbol() == "N"][0]
    h = [x.GetIdx() for x in pro.GetAtomWithIdx(n).GetNeighbors() if x.GetSymbol() == "H"][0]
    orient(pro, n, h)
    hh = float(np.linalg.norm(xyz(pro)[h]))
    th = np.radians(dha)
    r = hh * np.cos(th) + np.sqrt((hh * np.cos(th)) ** 2 - hh ** 2 + dist ** 2)
    target = np.array([hh - r * np.cos(th), r * np.sin(th), 0.0])
    lig = embed("CC(C)=O")
    o = [a.GetIdx() for a in lig.GetAtoms() if a.GetSymbol() == "O"][0]
    place(lig, o, target)
    c = [a.GetIdx() for a in lig.GetAtoms() if a.GetSymbol() == "C"][0]
    L = make_ligand(lig, {o: "O1", c: "C1"})
    P = make_protein(pro, "LYS", {n: "NZ", h: "HZ1"})
    return analyse(L, P)


def ionic(dist: float, ligand_cation=True, same_sign=False) -> list[dict]:
    cat, an = embed("C[NH3+]"), embed("CC(=O)[O-]")
    lig, pro = (cat, an) if ligand_cation else (an, cat)
    if same_sign:
        pro = embed("C[NH3+]") if ligand_cation else embed("CC(=O)[O-]")
    ln = [a.GetIdx() for a in lig.GetAtoms() if abs(a.GetFormalCharge()) == 1][0]
    pn = [a.GetIdx() for a in pro.GetAtoms() if abs(a.GetFormalCharge()) == 1][0]
    place(lig, ln, (0, 0, 0))
    place(pro, pn, (dist, 0, 0))
    # point both molecules' tails away from each other (+/-x) to avoid accidental H contacts
    lc = [x.GetIdx() for x in lig.GetAtomWithIdx(ln).GetNeighbors() if x.GetSymbol() == "C"][0]
    pc = [x.GetIdx() for x in pro.GetAtomWithIdx(pn).GetNeighbors() if x.GetSymbol() == "C"][0]
    for m, i, c, sign in ((lig, ln, lc, -1), (pro, pn, pc, 1)):
        p = xyz(m)
        R = rot_between(p[c] - p[i], np.array([sign, 0.0, 0.0]))
        set_xyz(m, (p - p[i]) @ R.T + p[i])
    L = make_ligand(lig, {ln: "X1"})
    P = make_protein(pro, "ASP", {pn: "OD1"})
    return analyse(L, P)


def hydrophobic(dist: float) -> list[dict]:
    # isobutane central carbon (CH, matches the hydrophobe pattern); one methyl points away along +/-x
    lig, pro = embed("CC(C)C"), embed("CC(C)C")
    place(lig, 1, (0, 0, 0))
    place(pro, 1, (dist, 0, 0))
    for m, sign in ((lig, -1), (pro, 1)):
        p = xyz(m)
        R = rot_between(p[0] - p[1], np.array([sign, 0.0, 0.0]))
        set_xyz(m, (p - p[1]) @ R.T + p[1])
    L = make_ligand(lig, {1: "C1"})
    P = make_protein(pro, "LEU", {1: "CG"})
    return analyse(L, P)


def pi_stack(centroid_dist: float) -> list[dict]:
    lig, pro = embed("c1ccccc1"), embed("c1ccccc1")
    for m in (lig, pro):
        p = xyz(m)
        p = p - p[:6].mean(axis=0)
        # plane normal -> z
        u, s, vt = np.linalg.svd(p[:6])
        R = rot_between(vt[2], np.array([0, 0, 1.0]))
        set_xyz(m, p @ R.T)
    p = xyz(pro)
    set_xyz(pro, p + np.array([0, 0, centroid_dist]))
    L = make_ligand(lig, {})
    P = make_protein(pro, "PHE", {})
    return analyse(L, P)


def cation_pi(dist: float) -> list[dict]:
    pro = embed("c1ccccc1")
    p = xyz(pro) - xyz(pro)[:6].mean(axis=0)
    _, _, vt = np.linalg.svd(p[:6])
    set_xyz(pro, p @ rot_between(vt[2], np.array([0, 0, 1.0])).T)
    lig = embed("C[NH3+]")
    n = [a.GetIdx() for a in lig.GetAtoms() if a.GetSymbol() == "N"][0]
    c = [x.GetIdx() for x in lig.GetAtomWithIdx(n).GetNeighbors() if x.GetSymbol() == "C"][0]
    place(lig, n, (0, 0, dist))
    pp = xyz(lig)
    set_xyz(lig, (pp - pp[n]) @ rot_between(pp[c] - pp[n], np.array([0, 0, 1.0])).T + pp[n])
    L = make_ligand(lig, {n: "N1"})
    P = make_protein(pro, "PHE", {})
    return analyse(L, P)


def halogen(dist: float) -> list[dict]:
    lig = embed("CBr")
    br = [a.GetIdx() for a in lig.GetAtoms() if a.GetSymbol() == "Br"][0]
    c = [a.GetIdx() for a in lig.GetAtoms() if a.GetSymbol() == "C"][0]
    orient(lig, br, c, axis_dir=(-1, 0, 0))   # C at -x of Br, so sigma-hole points to +x
    pro = embed("CC(C)=O")
    o = [a.GetIdx() for a in pro.GetAtoms() if a.GetSymbol() == "O"][0]
    cc = [a.GetIdx() for a in pro.GetAtoms() if a.GetSymbol() == "C" and a.GetDegree() == 3][0]
    place(pro, o, (dist, 0, 0))
    pp = xyz(pro)
    # carbonyl C placed so that Br...O=C is 120 degrees (inside the 80-140 window)
    want = np.array([np.cos(np.radians(60)), np.sin(np.radians(60)), 0.0])
    set_xyz(pro, (pp - pp[o]) @ rot_between(pp[cc] - pp[o], want).T + pp[o])
    L = make_ligand(lig, {br: "BR1", c: "C1"})
    P = make_protein(pro, "ASN", {o: "OD1", cc: "CG"})
    return analyse(L, P)


def metal(dist: float) -> list[dict]:
    lig = embed("CC(=O)[O-]")
    o = [a.GetIdx() for a in lig.GetAtoms() if a.GetSymbol() == "O"][0]
    place(lig, o, (0, 0, 0))
    pro = embed("C")
    place(pro, 0, (30, 0, 0))
    L = make_ligand(lig, {o: "O1"})
    P = make_protein(pro, "ALA", {0: "CB"})
    zn = atom("ZN", "Zn", dist, 0.0, 0.0, resname="ZN", resseq=300)
    return analyse(L, P, ions=[zn])


def water_bridge(d_lig: float, d_prot: float) -> list[dict]:
    lig = embed("CC(=O)C")
    o = [a.GetIdx() for a in lig.GetAtoms() if a.GetSymbol() == "O"][0]
    place(lig, o, (-d_lig, 0, 0))
    pro = embed("CC(C)=O")
    po = [a.GetIdx() for a in pro.GetAtoms() if a.GetSymbol() == "O"][0]
    place(pro, po, (d_prot, 0.0, 0.0))
    L = make_ligand(lig, {o: "O1"})
    P = make_protein(pro, "ASN", {po: "OD1"})
    w = atom("O", "O", 0.0, 0.0, 0.0)
    return analyse(L, P, waters=[w])
