"""Atom-level interaction detection.

ProLIF (count=True) is the base engine for H-bonds, halogen bonds, ionic,
cation-pi, pi-stacking, hydrophobic and van der Waals contacts. ZLIS adds
independent geometry for things ProLIF does not cover in this setup:
water-mediated bridges, ligand-metal coordination and steric clashes.

Every event carries the original-PDB atom identity (uid), roles, geometry and
Cartesian anchors, so tables, 2D diagrams and PyMOL all draw the same event.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
from rdkit import Chem

from .chemistry import LigandChem, atom_uid
from .protein_mol import is_backbone
from .structure import Atom, ResidueGroup, METALS

warnings.filterwarnings("ignore")

PROLIF_TYPES = ["HBDonor", "HBAcceptor", "XBDonor", "XBAcceptor", "Cationic", "Anionic",
                "CationPi", "PiCation", "PiStacking", "Hydrophobic", "VdWContact"]

TYPE_MAP = {
    "HBDonor": "hydrogen_bond", "HBAcceptor": "hydrogen_bond",
    "XBDonor": "halogen_bond", "XBAcceptor": "halogen_bond",
    "Cationic": "ionic", "Anionic": "ionic",
    "CationPi": "cation_pi", "PiCation": "cation_pi",
    "FaceToFace": "pi_pi", "EdgeToFace": "pi_pi", "PiStacking": "pi_pi",
    "Hydrophobic": "hydrophobic", "VdWContact": "vdw_contact",
}

# Distance thresholds used for ZLIS-native geometry (documented in the report).
WATER_BRIDGE_MAX = 3.5
METAL_COORD_MAX = 2.8
CLASH_OVERLAP_MIN = 0.5   # Angstrom of vdW overlap beyond which a pair is flagged
POLAR = {"N", "O"}


def _vdw(sym: str) -> float:
    return Chem.GetPeriodicTable().GetRvdw(sym)


@dataclass
class ProteinIndex:
    mol: Chem.Mol
    uid_by_idx: dict[int, str]
    name_by_idx: dict[int, str]
    res_by_idx: dict[int, dict]


def index_protein(mol: Chem.Mol, orig_atoms: list[Atom], model: int) -> ProteinIndex:
    """Map every prepared-protein atom to its ORIGINAL PDB atom identity."""
    lookup = {(a.chain, a.resseq, a.icode, a.resname, a.name.strip()): a for a in orig_atoms}
    uid_by, name_by, res_by = {}, {}, {}
    for at in mol.GetAtoms():
        ri = at.GetPDBResidueInfo()
        chain = ri.GetChainId().strip() or "_"
        key = (chain, ri.GetResidueNumber(), ri.GetInsertionCode().strip(), ri.GetResidueName().strip(), ri.GetName().strip())
        o = lookup.get(key)
        name = ri.GetName().strip()
        if o is not None:
            uid = atom_uid(o)
        else:  # hydrogen added in preparation: inherit parent identity + suffix
            nb = [n for n in at.GetNeighbors()]
            parent = nb[0].GetPDBResidueInfo().GetName().strip() if nb else "?"
            uid = f"M{model}:{chain}:{key[3]}{key[1]}{key[2]}:{name}(H on {parent})"
        idx = at.GetIdx()
        uid_by[idx], name_by[idx] = uid, name
        res_by[idx] = {"chain": chain, "resname": key[3], "resseq": key[1], "icode": key[2],
                       "label": f"{chain}:{key[3]}{key[1]}{key[2]}"}
    return ProteinIndex(mol, uid_by, name_by, res_by)


def _xyz(mol: Chem.Mol, idx: int) -> list[float]:
    p = mol.GetConformer().GetAtomPosition(idx)
    return [round(p.x, 3), round(p.y, 3), round(p.z, 3)]


def _centroid(mol: Chem.Mol, idxs) -> list[float]:
    pts = np.array([_xyz(mol, i) for i in idxs])
    return [round(float(v), 3) for v in pts.mean(axis=0)]


def _heavy_parent(mol: Chem.Mol, idx: int) -> int:
    a = mol.GetAtomWithIdx(idx)
    if a.GetSymbol() == "H" and a.GetDegree() == 1:
        return a.GetNeighbors()[0].GetIdx()
    return idx


def run_analysis(lig: LigandChem, prot: ProteinIndex, polymer_atoms: list[Atom],
                 waters: list[ResidueGroup], ions: list[ResidueGroup], model: int,
                 lig_label: str) -> list[dict]:
    import prolif as plf

    lmol, pmol = lig.mol, prot.mol
    lig_info = lig.atoms

    def L(idx: int) -> dict:
        i = _heavy_parent(lmol, idx)
        info = lig_info[i]
        return {"uid": info["uid"], "name": info["name"], "element": info["element"],
                "rdkit_idx": i, "xyz": _xyz(lmol, i)}

    def LH(idx: int) -> dict:
        info = lig_info[idx]
        return {"uid": info["uid"], "name": info["name"], "element": "H", "rdkit_idx": idx, "xyz": _xyz(lmol, idx)}

    def P(idx: int) -> dict:
        i = _heavy_parent(pmol, idx)
        return {"uid": prot.uid_by_idx[i], "name": prot.name_by_idx[i],
                "element": pmol.GetAtomWithIdx(i).GetSymbol(), "rdkit_idx": i, "xyz": _xyz(pmol, i)}

    def PH(idx: int) -> dict:
        return {"uid": prot.uid_by_idx[idx], "name": prot.name_by_idx[idx], "element": "H",
                "rdkit_idx": idx, "xyz": _xyz(pmol, idx)}

    pm = plf.Molecule.from_rdkit(pmol)
    lm = plf.Molecule.from_rdkit(lmol)
    fp = plf.Fingerprint(PROLIF_TYPES, count=True)
    fp.run_from_iterable([lm], pm, progress=False)

    events: list[dict] = []

    def add(ev: dict) -> None:
        events.append(ev)

    for (_lres, pres), per_type in fp.ifp[0].items():
        for itype, items in per_type.items():
            cat = TYPE_MAP.get(itype, itype)
            for meta in items:
                lidx = list(meta["parent_indices"]["ligand"])
                pidx = list(meta["parent_indices"]["protein"])
                ev: dict = {"engine": "prolif", "prolif_type": itype, "type": cat, "model": model,
                            "distance": round(float(meta.get("distance", float("nan"))), 2)}
                res = prot.res_by_idx[pidx[0]]
                ev["residue"] = res

                if cat == "hydrogen_bond":
                    if itype == "HBDonor":   # ligand donor -> protein acceptor
                        d, h, a = L(lidx[0]), LH(lidx[1]), P(pidx[0])
                        ev.update(donor_side="ligand", donor=d, hydrogen=h, acceptor=a)
                    else:                    # protein donor -> ligand acceptor
                        d, h, a = P(pidx[0]), PH(pidx[1]), L(lidx[0])
                        ev.update(donor_side="protein", donor=d, hydrogen=h, acceptor=a)
                    ev["angle_DHA"] = round(float(meta["DHA_angle"]), 1)
                    ev["ligand_atoms"], ev["protein_atoms"] = [L(lidx[0])], [P(pidx[0])]
                    ev["distance_label"] = "D...A"
                    ev["h_a_distance"] = round(float(np.linalg.norm(
                        np.array(ev["hydrogen"]["xyz"]) - np.array(ev["acceptor"]["xyz"]))), 2)
                    ev["ligand_atoms"] = [d if itype == "HBDonor" else a]
                    ev["protein_atoms"] = [a if itype == "HBDonor" else d]
                elif cat == "halogen_bond":
                    # ProLIF XBDonor: ligand donates the halogen; indices (halogen-bearing C, X)
                    lat = [L(i) for i in lidx]
                    pat = [P(i) for i in pidx]
                    ev["ligand_atoms"], ev["protein_atoms"] = lat, pat
                    ev["angles"] = {k: round(float(v), 1) for k, v in meta.items() if k.endswith("angle")}
                elif cat in ("pi_pi", "cation_pi"):
                    ring_l = itype in ("PiCation", "FaceToFace", "EdgeToFace", "PiStacking")
                    lat = [L(i) for i in lidx]
                    pat = [P(i) for i in pidx]
                    ev["ligand_atoms"], ev["protein_atoms"] = lat, pat
                    ev["ligand_anchor"] = _centroid(lmol, lidx) if len(lidx) > 2 else lat[0]["xyz"]
                    ev["protein_anchor"] = _centroid(pmol, pidx) if len(pidx) > 2 else pat[0]["xyz"]
                    if cat == "pi_pi":
                        ev["geometry"] = "face_to_face" if itype == "FaceToFace" else (
                            "edge_to_face" if itype == "EdgeToFace" else "stacking")
                        if "plane_angle" in meta:
                            ev["plane_angle"] = round(float(meta["plane_angle"]), 1)
                    else:
                        ev["cation_side"] = "ligand" if itype == "CationPi" else "protein"
                        if "angle" in meta:
                            ev["angle"] = round(float(meta["angle"]), 1)
                else:
                    ev["ligand_atoms"] = [L(i) for i in lidx]
                    ev["protein_atoms"] = [P(i) for i in pidx]

                if cat == "vdw_contact":
                    # heavy-atom pairs only; H contacts are dominated by riding-H geometry
                    if any(a["element"] == "H" for a in ev["ligand_atoms"] + ev["protein_atoms"]):
                        continue
                    if any(lmol.GetAtomWithIdx(i).GetSymbol() == "H" for i in lidx) or \
                       any(pmol.GetAtomWithIdx(i).GetSymbol() == "H" for i in pidx):
                        continue
                    ev["radius_model"] = "RDKit periodic-table (Bondi) vdW radii, tolerance 0.0 A"

                names = {a["name"] for a in ev["protein_atoms"]}
                if cat == "pi_pi" or cat == "cation_pi":
                    ev["protein_part"] = "sidechain"
                else:
                    flags = [is_backbone(n) for n in names]
                    ev["protein_part"] = "backbone" if all(flags) else ("sidechain" if not any(flags) else "mixed")
                add(ev)

    _add_metal_coordination(events, lig, lmol, ions, model, L)
    _add_water_bridges(events, lig, lmol, pmol, prot, polymer_atoms, waters, model, L)
    _add_clashes(events, lmol, pmol, prot, L, P)
    _annotate(events, lig, prot)

    order = {"hydrogen_bond": 0, "ionic": 1, "metal_coordination": 2, "halogen_bond": 3, "pi_pi": 4,
             "cation_pi": 5, "water_bridge": 6, "hydrophobic": 7, "vdw_contact": 8, "steric_clash": 9}
    events.sort(key=lambda e: (order.get(e["type"], 99), e["residue"]["chain"], e["residue"]["resseq"],
                               e.get("distance", 0)))
    for n, ev in enumerate(events, 1):
        ev["event_id"] = f"{lig_label.replace(':', '_')}_E{n:03d}"
        ev["ligand_label"] = lig_label
    return events


def _add_metal_coordination(events, lig, lmol, ions, model, L):
    for g in ions:
        for ion in g.atoms:
            ip = np.array(ion.xyz)
            for at in lmol.GetAtoms():
                if at.GetSymbol() not in ("O", "N", "S") or at.GetIdx() >= len(lig.atoms):
                    continue
                p = lmol.GetConformer().GetAtomPosition(at.GetIdx())
                d = float(np.linalg.norm(ip - np.array([p.x, p.y, p.z])))
                if d <= METAL_COORD_MAX:
                    events.append({
                        "engine": "zlis", "type": "metal_coordination", "model": model,
                        "distance": round(d, 2), "residue": {"chain": ion.chain, "resname": ion.resname,
                        "resseq": ion.resseq, "icode": ion.icode, "label": f"{ion.chain}:{ion.resname}{ion.resseq}{ion.icode}"},
                        "ligand_atoms": [L(at.GetIdx())],
                        "protein_atoms": [{"uid": atom_uid(ion), "name": ion.name.strip(), "element": ion.element,
                                           "rdkit_idx": None, "xyz": [ion.x, ion.y, ion.z]}],
                        "protein_part": "metal",
                        "rule": f"ligand N/O/S within {METAL_COORD_MAX} A of metal ion",
                    })


def _add_water_bridges(events, lig, lmol, pmol, prot, polymer_atoms, waters, model, L):
    """Ligand polar atom - water O - protein polar atom, all within WATER_BRIDGE_MAX.

    Crystal waters rarely carry hydrogens, so only the O...X distances are
    evaluated; these are reported as candidates, never as confirmed H-bond paths.
    """
    if not waters:
        return
    lpolar = [a for a in lmol.GetAtoms() if a.GetSymbol() in POLAR]
    ppolar = [a for a in pmol.GetAtoms() if a.GetSymbol() in POLAR]
    if not lpolar or not ppolar:
        return
    lxyz = np.array([[*lmol.GetConformer().GetAtomPosition(a.GetIdx())] for a in lpolar])
    pxyz = np.array([[*pmol.GetConformer().GetAtomPosition(a.GetIdx())] for a in ppolar])
    for g in waters:
        o = next((a for a in g.atoms if a.element.upper() == "O"), None)
        if o is None:
            continue
        w = np.array(o.xyz)
        dl = np.linalg.norm(lxyz - w, axis=1)
        dp = np.linalg.norm(pxyz - w, axis=1)
        for i in np.where(dl <= WATER_BRIDGE_MAX)[0]:
            for j in np.where(dp <= WATER_BRIDGE_MAX)[0]:
                pa = ppolar[j]
                la = lpolar[i]
                res = prot.res_by_idx[pa.GetIdx()]
                events.append({
                    "engine": "zlis", "type": "water_bridge", "model": model,
                    "distance": round(float(max(dl[i], dp[j])), 2), "residue": res,
                    "ligand_atoms": [L(la.GetIdx())],
                    "protein_atoms": [{"uid": prot.uid_by_idx[pa.GetIdx()], "name": prot.name_by_idx[pa.GetIdx()],
                                       "element": pa.GetSymbol(), "rdkit_idx": pa.GetIdx(), "xyz": _xyz(pmol, pa.GetIdx())}],
                    "water": {"uid": atom_uid(o), "label": f"{o.chain}:{o.resname}{o.resseq}{o.icode}", "xyz": [o.x, o.y, o.z]},
                    "geometry": f"ligand...O {dl[i]:.2f} A / O...protein {dp[j]:.2f} A",
                    "segment_distances": {"ligand_to_water": round(float(dl[i]), 2), "water_to_protein": round(float(dp[j]), 2)},
                    "protein_part": "backbone" if is_backbone(prot.name_by_idx[pa.GetIdx()]) else "sidechain",
                    "confidence": "candidate",
                    "rule": f"both O...polar distances <= {WATER_BRIDGE_MAX} A (reported distance = the longer segment); water hydrogens absent so donor/acceptor roles are not assigned",
                })


def _add_clashes(events, lmol, pmol, prot, L, P):
    lh = [a.GetIdx() for a in lmol.GetAtoms() if a.GetSymbol() != "H"]
    ph = [a.GetIdx() for a in pmol.GetAtoms() if a.GetSymbol() != "H"]
    if not lh or not ph:
        return
    lx = np.array([[*lmol.GetConformer().GetAtomPosition(i)] for i in lh])
    px = np.array([[*pmol.GetConformer().GetAtomPosition(i)] for i in ph])
    d = np.linalg.norm(lx[:, None, :] - px[None, :, :], axis=2)
    lrad = np.array([_vdw(lmol.GetAtomWithIdx(i).GetSymbol()) for i in lh])
    prad = np.array([_vdw(pmol.GetAtomWithIdx(i).GetSymbol()) for i in ph])
    overlap = (lrad[:, None] + prad[None, :]) - d
    for a, b in zip(*np.where(overlap >= CLASH_OVERLAP_MIN)):
        li, pi = lh[a], ph[b]
        ls, ps = lmol.GetAtomWithIdx(li).GetSymbol(), pmol.GetAtomWithIdx(pi).GetSymbol()
        if ls in POLAR and ps in POLAR and d[a, b] >= 2.5:
            continue  # normal H-bond / ionic contact range
        if ls in POLAR and ps in POLAR and overlap[a, b] < 1.0:
            continue
        events.append({
            "engine": "zlis", "type": "steric_clash", "model": 0, "distance": round(float(d[a, b]), 2),
            "residue": prot.res_by_idx[pi], "ligand_atoms": [L(li)], "protein_atoms": [P(pi)],
            "overlap": round(float(overlap[a, b]), 2),
            "protein_part": "backbone" if is_backbone(prot.name_by_idx[pi]) else "sidechain",
            "rule": f"heavy-atom vdW overlap >= {CLASH_OVERLAP_MIN} A (Bondi radii)",
        })


def _annotate(events, lig: LigandChem, prot: ProteinIndex) -> None:
    for ev in events:
        flags = []
        if ev["type"] in ("hydrogen_bond", "ionic", "cation_pi") and lig.status != "ok":
            flags.append("ligand_chemistry_unconfirmed")
        if ev["type"] == "hydrogen_bond" and ev["residue"]["resname"] in ("HIS", "HID", "HIE", "HIP"):
            flags.append("his_protonation_state_assumed")
        if ev["type"] in ("hydrogen_bond",) and lig.chem_source == "unverified":
            flags.append("bond_orders_unverified")
        ev.setdefault("confidence", "provisional" if flags else "ok")
        ev["flags"] = flags
