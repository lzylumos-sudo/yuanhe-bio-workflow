"""PDB reading, model/altloc handling and ligand candidate detection.

Every atom keeps its original identity (model, chain, resname, resseq, icode,
atom name, serial). Nothing here moves coordinates.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

WATER_NAMES = {"HOH", "WAT", "H2O", "DOD", "TIP", "TIP3", "SOL"}
AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
    "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "SEC", "PYL",
    "HID", "HIE", "HIP", "CYX", "CYM", "ASH", "GLH", "LYN", "MSE",
}
NUCLEIC = {"A", "C", "G", "U", "DA", "DC", "DG", "DT", "DU", "RA", "RC", "RG", "RU"}
METALS = {
    "LI", "NA", "MG", "K", "CA", "MN", "FE", "CO", "NI", "CU", "ZN", "CD", "HG", "AL",
    "MO", "W", "PT", "AU", "AG", "SR", "BA", "CS", "RB", "FE2", "CU1", "NI2",
}
# Common crystallisation additives that should not be offered as a ligand.
BUFFER_NOISE = {"SO4", "PO4", "GOL", "EDO", "PEG", "ACT", "DMS", "CL", "BR", "IOD", "NO3", "FMT", "TRS", "MPD", "PG4"}


@dataclass
class Atom:
    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    resseq: int
    icode: str
    x: float
    y: float
    z: float
    occupancy: float
    bfactor: float
    element: str
    record: str  # ATOM / HETATM
    model: int
    line: str = ""

    @property
    def res_key(self) -> tuple:
        return (self.chain, self.resseq, self.icode, self.resname)

    @property
    def is_h(self) -> bool:
        return self.element.upper() in ("H", "D")

    @property
    def xyz(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)


@dataclass
class ResidueGroup:
    key: tuple
    atoms: list[Atom] = field(default_factory=list)

    @property
    def label(self) -> str:
        c, n, i, r = self.key
        return f"{c}:{r}{n}{i}".strip()


@dataclass
class Structure:
    path: str
    models: dict[int, list[Atom]]
    warnings: list[str] = field(default_factory=list)


def _element_from_line(line: str, name: str) -> str:
    el = line[76:78].strip()
    if el:
        return el.capitalize() if len(el) > 1 else el.upper()
    stripped = name.strip()
    if len(name) >= 2 and name[0] == " ":
        return stripped[:1].upper()
    letters = "".join(ch for ch in stripped if ch.isalpha())
    return letters[:1].upper() if letters else "X"


def read_pdb(path: str | Path) -> Structure:
    path = Path(path)
    models: dict[int, list[Atom]] = {}
    model = 1
    warnings: list[str] = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec == "MODEL":
                try:
                    model = int(line[10:14])
                except ValueError:
                    model += 1
                continue
            if rec not in ("ATOM", "HETATM"):
                continue
            line = line.rstrip("\n").ljust(80)
            try:
                a = Atom(
                    serial=int(line[6:11]),
                    name=line[12:16],
                    altloc=line[16].strip(),
                    resname=line[17:20].strip(),
                    chain=line[21].strip() or "_",
                    resseq=int(line[22:26]),
                    icode=line[26].strip(),
                    x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54]),
                    occupancy=float(line[54:60] or 1.0),
                    bfactor=float(line[60:66] or 0.0),
                    element=_element_from_line(line, line[12:16]),
                    record=rec, model=model, line=line,
                )
            except ValueError:
                warnings.append(f"unparseable line skipped: {line[:30]!r}")
                continue
            models.setdefault(model, []).append(a)
    if not models:
        raise ValueError(f"No ATOM/HETATM records in {path}")
    if len(models) > 1:
        warnings.append(
            f"{len(models)} MODELs found; each is analysed separately, atoms are never mixed between models."
        )
    return Structure(str(path), models, warnings)


def resolve_altlocs(atoms: list[Atom]) -> tuple[list[Atom], list[str]]:
    """Keep the highest-occupancy (first on ties) altloc for each atom."""
    best: dict[tuple, Atom] = {}
    order: list[tuple] = []
    dropped = 0
    for a in atoms:
        k = (a.chain, a.resseq, a.icode, a.resname, a.name)
        if k not in best:
            best[k] = a
            order.append(k)
        elif a.altloc and a.occupancy > best[k].occupancy + 1e-6:
            best[k] = a
            dropped += 1
        else:
            dropped += 1
    notes = []
    if dropped:
        notes.append(f"{dropped} alternate-location atoms dropped (kept highest occupancy, first on ties).")
    return [best[k] for k in order], notes


def group_residues(atoms: list[Atom]) -> list[ResidueGroup]:
    groups: dict[tuple, ResidueGroup] = {}
    for a in atoms:
        groups.setdefault(a.res_key, ResidueGroup(a.res_key)).atoms.append(a)
    return list(groups.values())


@dataclass
class LigandCandidate:
    group: ResidueGroup
    heavy_atoms: int
    min_dist_to_protein: float
    kind: str

    @property
    def label(self) -> str:
        return self.group.label

    def as_dict(self) -> dict:
        c, n, i, r = self.group.key
        return {
            "label": self.label, "resname": r, "chain": c, "resseq": n, "icode": i,
            "heavy_atoms": self.heavy_atoms,
            "min_dist_to_protein": round(self.min_dist_to_protein, 2),
            "kind": self.kind,
        }


def classify_components(atoms: list[Atom]):
    """Split into polymer, waters, ions and hetero (ligand-like) residues."""
    polymer, waters, ions, hets = [], [], [], []
    for g in group_residues(atoms):
        rn = g.key[3]
        if rn in WATER_NAMES:
            waters.append(g)
        elif rn in AMINO_ACIDS or rn in NUCLEIC:
            polymer.append(g)
        elif len(g.atoms) == 1 and (g.atoms[0].element.upper() in METALS or rn.upper() in METALS):
            ions.append(g)
        else:
            hets.append(g)
    return polymer, waters, ions, hets


def _min_dist(g: ResidueGroup, prot_xyz) -> float:
    import numpy as np

    pts = np.array([a.xyz for a in g.atoms if not a.is_h])
    if len(pts) == 0 or len(prot_xyz) == 0:
        return float("inf")
    d = np.linalg.norm(pts[:, None, :] - prot_xyz[None, :, :], axis=2)
    return float(d.min())


def find_ligand_candidates(atoms: list[Atom], min_heavy: int = 6) -> list[LigandCandidate]:
    import numpy as np

    polymer, _w, _i, hets = classify_components(atoms)
    prot_xyz = np.array([a.xyz for g in polymer for a in g.atoms if not a.is_h])
    out: list[LigandCandidate] = []
    for g in hets:
        if g.key[3].upper() in BUFFER_NOISE:
            continue
        heavy = sum(1 for a in g.atoms if not a.is_h)
        if heavy < 4:
            continue
        md = _min_dist(g, prot_xyz)
        if md > 8.0:
            continue
        kind = "ligand" if heavy >= min_heavy else "small-fragment"
        out.append(LigandCandidate(g, heavy, md, kind))
    out.sort(key=lambda c: (-(c.heavy_atoms >= 10), c.min_dist_to_protein, -c.heavy_atoms))
    return out
