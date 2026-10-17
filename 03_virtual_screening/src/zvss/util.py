"""Small shared helpers: atomic file writes, timestamps, PDB I/O, geometry.

The PDB reader here is deliberately explicit about models, altlocs and insertion
codes.  Docking outputs and experimental structures are full of those, and any
pipeline that reduces a residue to a single integer loses information that later
becomes an unexplainable mismatch.
"""

from __future__ import annotations

import datetime as _dt
import json
import os
import tempfile
from dataclasses import dataclass, field
from typing import Any, Iterable, Iterator, Sequence

import numpy as np

WATER_NAMES = {"HOH", "WAT", "DOD", "H2O", "TIP", "SOL"}
ION_NAMES = {
    "NA", "K", "CL", "MG", "CA", "ZN", "MN", "FE", "CU", "CO", "NI", "CD", "SO4", "PO4", "IOD", "BR",
}
# Small-molecule-looking HETATM names we must NOT silently treat as protein.
LIPID_NAMES = {"POP", "DOP", "DPP", "POC", "PEA", "CHL", "CHL1", "CLR", "POPC", "POPE", "DPPC", "OLC"}
AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE", "SEC", "PYL", "HSD", "HSE", "HSP", "CYX", "ASH", "GLH", "LYN",
}


def now_iso() -> str:
    return _dt.datetime.now(_dt.timezone.utc).astimezone().isoformat(timespec="seconds")


def utc_stamp() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def read_text(path: str, encoding: str = "utf-8") -> str:
    with open(path, "r", encoding=encoding, errors="replace") as handle:
        return handle.read()


def read_bytes(path: str) -> bytes:
    with open(path, "rb") as handle:
        return handle.read()


def atomic_write_text(path: str, text: str, encoding: str = "utf-8") -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    directory = os.path.dirname(os.path.abspath(path))
    fd, tmp = tempfile.mkstemp(prefix=".zvss-", suffix=".tmp", dir=directory)
    try:
        with os.fdopen(fd, "w", encoding=encoding, newline="\n") as handle:
            handle.write(text)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    return path


def atomic_write_bytes(path: str, data: bytes) -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    directory = os.path.dirname(os.path.abspath(path))
    fd, tmp = tempfile.mkstemp(prefix=".zvss-", suffix=".tmp", dir=directory)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    return path


def atomic_write_json(path: str, payload: Any, indent: int = 2) -> str:
    return atomic_write_text(path, json.dumps(payload, ensure_ascii=False, indent=indent, default=str) + "\n")


def read_json(path: str, default: Any = None) -> Any:
    if not os.path.exists(path):
        return default
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path


def relpath(path: str, start: str) -> str:
    try:
        return os.path.relpath(path, start).replace(os.sep, "/")
    except ValueError:
        return path.replace(os.sep, "/")


# --------------------------------------------------------------------------- #
# PDB parsing
# --------------------------------------------------------------------------- #


@dataclass
class PdbAtom:
    record: str
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
    occupancy: float = 1.0
    bfactor: float = 0.0
    element: str = ""
    charge: str = ""
    model: int = 1
    line_no: int = 0

    @property
    def coord(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z], dtype=float)

    @property
    def key(self) -> tuple[str, int, str, str, str, int]:
        return (self.chain, self.resseq, self.icode, self.resname, self.name, self.model)

    def is_hydrogen(self) -> bool:
        element = (self.element or self.name).strip().upper()
        if element in {"H", "D"}:
            return True
        stripped = self.name.strip()
        return bool(stripped) and stripped[0] in {"H", "D"} and not stripped[:2].strip().isdigit()

    def to_line(self) -> str:
        record = self.record if len(self.record) == 6 else "ATOM  "
        name = self.name
        if len(name) < 4:
            name = f" {name:<3}" if len(name) <= 3 else name
        element = (self.element or self.name.strip()[:1]).upper()
        # PDB resSeq has exactly four columns. Wrapping with ``% 10000`` silently
        # rewrote negative author numbers (fusion constructs) into 9xxx values and
        # destroyed residue identity, so out-of-range numbers now fail loudly.
        if not -999 <= int(self.resseq) <= 9999:
            raise ValueError(
                f"resSeq {self.resseq} for {self.chain}:{self.resname} cannot be written in PDB columns 23-26; "
                "choose a numbering that fits, or keep the original mmCIF"
            )
        return (
            f"{record:<6}{self.serial % 100000:>5} {name:<4}{self.altloc or ' ':1}{self.resname:>3} "
            f"{(self.chain or ' '):1}{int(self.resseq):>4}{(self.icode or ' '):1}   "
            f"{self.x:8.3f}{self.y:8.3f}{self.z:8.3f}{self.occupancy:6.2f}{self.bfactor:6.2f}          "
            f"{element:>2}{self.charge or '':>2}"
        )


@dataclass
class PdbStructure:
    atoms: list[PdbAtom] = field(default_factory=list)
    models: list[int] = field(default_factory=list)
    header_lines: list[str] = field(default_factory=list)

    def chains(self) -> list[str]:
        seen: list[str] = []
        for atom in self.atoms:
            if atom.chain not in seen:
                seen.append(atom.chain)
        return seen

    def residues(self, chain: str | None = None) -> list[tuple[str, int, str, str]]:
        out: list[tuple[str, int, str, str]] = []
        for atom in self.atoms:
            if chain is not None and atom.chain != chain:
                continue
            key = (atom.chain, atom.resseq, atom.icode, atom.resname)
            if key not in out:
                out.append(key)
        return out

    def by_model(self, model: int) -> "PdbStructure":
        return PdbStructure([a for a in self.atoms if a.model == model], [model], list(self.header_lines))

    def heavy_atoms(self) -> list[PdbAtom]:
        return [a for a in self.atoms if not a.is_hydrogen()]

    def coords(self) -> np.ndarray:
        if not self.atoms:
            return np.zeros((0, 3))
        return np.array([[a.x, a.y, a.z] for a in self.atoms], dtype=float)


def parse_pdb_text(text: str, *, keep_altloc: str = "first", keep_hydrogens: bool = True) -> PdbStructure:
    structure = PdbStructure()
    current_model = 1
    seen_altloc: set[tuple[str, int, str, str, str]] = set()
    for line_no, line in enumerate(text.splitlines(), start=1):
        record = line[:6]
        if record.startswith("MODEL"):
            try:
                current_model = int(line[10:14].strip() or 1)
            except ValueError:
                current_model = len(structure.models) + 1
            if current_model not in structure.models:
                structure.models.append(current_model)
            continue
        if record.startswith("ENDMDL"):
            continue
        if not (record.startswith("ATOM") or record.startswith("HETATM")):
            if not structure.atoms and record.startswith(("HEADER", "TITLE", "REMARK", "CRYST1", "SEQRES", "DBREF", "SEQADV")):
                structure.header_lines.append(line.rstrip("\n"))
            continue
        altloc = line[16:17]
        if altloc.strip() and keep_altloc == "first":
            tag = (line[21:22], int(line[22:26] or 0), line[26:27], line[17:20], line[12:16])
            if tag in seen_altloc:
                continue
            seen_altloc.add(tag)
        elif altloc.strip() and keep_altloc == "none":
            continue
        atom = PdbAtom(
            record="HETATM" if line.startswith("HETATM") else "ATOM  ",
            serial=int(line[6:11].strip() or 0) if line[6:11].strip().isdigit() else len(structure.atoms) + 1,
            name=line[12:16].strip(),
            altloc=altloc.strip(),
            resname=line[17:20].strip(),
            chain=line[21:22].strip(),
            resseq=int(line[22:26].strip() or 0) if line[22:26].strip().lstrip("-").isdigit() else 0,
            icode=line[26:27].strip(),
            x=float(line[30:38]),
            y=float(line[38:46]),
            z=float(line[46:54]),
            occupancy=float(line[54:60] or 1.0),
            bfactor=float(line[60:66] or 0.0),
            element=line[76:78].strip(),
            charge=line[78:80].strip(),
            model=current_model,
            line_no=line_no,
        )
        if not keep_hydrogens and atom.is_hydrogen():
            continue
        structure.atoms.append(atom)
    if not structure.models:
        structure.models = [1]
    return structure


def read_pdb(path: str, **kwargs: Any) -> PdbStructure:
    return parse_pdb_text(read_text(path), **kwargs)


def write_pdb(path: str, atoms: Iterable[PdbAtom], extra_header: Sequence[str] = ()) -> str:
    lines: list[str] = []
    for header in extra_header:
        lines.append(header if len(header) <= 80 else header[:80])
    for index, atom in enumerate(atoms, start=1):
        atom.serial = index
        lines.append(atom.to_line())
    lines.append("END")
    return atomic_write_text(path, "\n".join(lines) + "\n")


# --------------------------------------------------------------------------- #
# geometry
# --------------------------------------------------------------------------- #


def pairwise_distances(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    if a.size == 0 or b.size == 0:
        return np.zeros((len(a), len(b)))
    diff = a[:, None, :] - b[None, :, :]
    return np.sqrt((diff * diff).sum(axis=-1))


def min_distance(a: np.ndarray, b: np.ndarray) -> float | None:
    if a.size == 0 or b.size == 0:
        return None
    return float(pairwise_distances(a, b).min())


def centroid(coords: np.ndarray) -> np.ndarray:
    if coords.size == 0:
        return np.zeros(3)
    return coords.mean(axis=0)


def unit(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm < 1e-12:
        return vector * 0.0
    return vector / norm


def kabsch_transform(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Return (rotation, translation, rmsd) so that mobile @ R + t ~= target."""
    if len(mobile) < 3 or len(mobile) != len(target):
        raise ValueError("kabsch needs >=3 paired points of equal length")
    mobile_center = mobile.mean(axis=0)
    target_center = target.mean(axis=0)
    a = mobile - mobile_center
    b = target - target_center
    covariance = a.T @ b
    v, _s, wt = np.linalg.svd(covariance)
    d = np.sign(np.linalg.det(v @ wt))
    correction = np.diag([1.0, 1.0, d])
    rotation = (v @ correction @ wt).T
    translation = target_center - rotation @ mobile_center
    aligned = mobile @ rotation.T + translation
    rmsd = float(np.sqrt(((aligned - target) ** 2).sum(axis=1).mean()))
    return rotation, translation, rmsd


def apply_transform(coords: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    return coords @ rotation.T + translation


def as_matrix3(payload: Any) -> np.ndarray:
    return np.array(payload, dtype=float).reshape(3, 3)


def chunked(items: Sequence[Any], size: int) -> Iterator[Sequence[Any]]:
    for start in range(0, len(items), size):
        yield items[start : start + size]


def median_or_none(values: Iterable[float | None]) -> float | None:
    clean = [float(v) for v in values if v is not None]
    if not clean:
        return None
    return float(np.median(clean))


def safe_int(value: Any) -> int | None:
    try:
        if value is None or (isinstance(value, str) and not value.strip()):
            return None
        return int(float(value))
    except (TypeError, ValueError):
        return None


def safe_float(value: Any) -> float | None:
    try:
        if value is None or (isinstance(value, str) and not value.strip()):
            return None
        out = float(value)
    except (TypeError, ValueError):
        return None
    if out != out or out in (float("inf"), float("-inf")):
        return None
    return out
