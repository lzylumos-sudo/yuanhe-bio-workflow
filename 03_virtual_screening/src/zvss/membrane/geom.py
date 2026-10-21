"""Small, dependency-light geometry helpers shared by the Y-MAVS modules.

The functions here deliberately do one thing each: convert parsed PDB atoms to
coordinates, look up a van der Waals radius, classify a residue for composition
reporting, and return receptor contacts grouped by residue.  None of them make a
binding-affinity claim.
"""

from __future__ import annotations

import math
from typing import Any, Iterable

import numpy as np

# Bondi/van der Waals radii in Angstrom.  Values are used only to build geometric
# occupancy grids and to distinguish close contact from atomic overlap.  They are
# not a force-field parameter set.
VDW_RADII: dict[str, float] = {
    "H": 1.20,
    "D": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
    "SE": 1.90,
    "B": 1.92,
    "SI": 2.10,
    "NA": 2.27,
    "K": 2.75,
    "MG": 1.73,
    "CA": 2.31,
    "ZN": 1.39,
    "FE": 2.00,
    "MN": 2.00,
    "CU": 1.40,
    "CO": 2.00,
    "NI": 1.63,
    "CD": 1.58,
}

HYDROPHOBIC_RESIDUES = frozenset(
    {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO", "TYR"}
)
CHARGED_RESIDUES = frozenset({"ASP", "GLU", "LYS", "ARG", "HIS", "HSD", "HSE", "HSP"})
POLAR_RESIDUES = frozenset(
    {"SER", "THR", "ASN", "GLN", "CYS", "TYR", "HIS", "HSD", "HSE", "HSP", "TRP"}
)
# Glycine has no side chain and is reported separately from the three broad
# composition classes rather than being silently forced into "hydrophobic".
OTHER_RESIDUES = frozenset({"GLY", "MSE", "SEC", "PYL"})


def element_radius(element: str, default: float = 1.70) -> float:
    return float(VDW_RADII.get((element or "").strip().upper(), default))


def as_coords(atoms: Iterable[Any]) -> np.ndarray:
    values = []
    for atom in atoms:
        if getattr(atom, "is_hydrogen", None) and atom.is_hydrogen():
            continue
        values.append([float(atom.x), float(atom.y), float(atom.z)])
    return np.asarray(values, dtype=float).reshape((-1, 3))


def residue_key(atom: Any) -> tuple[str, int, str, str]:
    """Stable residue identity, including insertion code and chain."""
    return (
        str(getattr(atom, "chain", "") or ""),
        int(getattr(atom, "resseq", 0) or 0),
        str(getattr(atom, "icode", "") or ""),
        str(getattr(atom, "resname", "") or "").upper(),
    )


def residue_label(key: tuple[str, int, str, str]) -> str:
    chain, resseq, icode, resname = key
    return f"{chain}{resseq}{icode}:{resname}" if chain else f"{resseq}{icode}:{resname}"


def residue_class(resname: str) -> str:
    """Broad composition class; not a solvent-exposure or energy descriptor."""
    name = (resname or "").strip().upper()
    if name in HYDROPHOBIC_RESIDUES:
        return "hydrophobic"
    if name in CHARGED_RESIDUES:
        return "charged"
    if name in POLAR_RESIDUES:
        return "polar"
    return "other"


def is_polar_atom(element: str) -> bool:
    return (element or "").strip().upper() in {"N", "O", "S"}


def _distance_blocks(a: np.ndarray, b: np.ndarray, block: int = 2048) -> np.ndarray:
    """Minimum distance from each row of ``a`` to any row of ``b``.

    cKDTree is preferred; this chunked implementation exists only as a fallback
    for installations where SciPy is unavailable.  It never expands a full
    pairwise matrix for a protein-sized receptor.
    """
    if len(a) == 0 or len(b) == 0:
        return np.full(len(a), math.inf, dtype=float)
    out = np.full(len(a), math.inf, dtype=float)
    for start in range(0, len(a), block):
        part = a[start : start + block]
        diff = part[:, None, :] - b[None, :, :]
        out[start : start + len(part)] = np.sqrt((diff * diff).sum(axis=2)).min(axis=1)
    return out


def nearest_distances(query: np.ndarray, reference: np.ndarray) -> np.ndarray:
    if len(query) == 0 or len(reference) == 0:
        return np.full(len(query), math.inf, dtype=float)
    try:
        from scipy.spatial import cKDTree

        distances, _ = cKDTree(reference).query(query, k=1, workers=-1)
        return np.asarray(distances, dtype=float)
    except ImportError:  # pragma: no cover - SciPy is a declared dependency
        return _distance_blocks(query, reference)


def receptor_contacts_by_residue(
    ligand_xyz: np.ndarray,
    receptor_atoms: Iterable[Any],
    cutoff_A: float,
) -> list[dict[str, Any]]:
    """Return per-residue receptor contacts within ``cutoff_A``.

    Only heavy atoms are considered.  The result is a geometric contact set;
    it is not a contact energy and says nothing about binding affinity.
    """
    ligand = np.asarray(ligand_xyz, dtype=float).reshape((-1, 3))
    atoms = [a for a in receptor_atoms if not (hasattr(a, "is_hydrogen") and a.is_hydrogen())]
    if len(ligand) == 0 or not atoms:
        return []
    rec_xyz = as_coords(atoms)
    try:
        from scipy.spatial import cKDTree

        tree = cKDTree(rec_xyz)
        # A ligand atom can have more than one receptor neighbour inside the
        # cutoff; take the union of indices so a residue is not lost because
        # only the nearest atom was queried.
        hits: set[int] = set()
        for point in ligand:
            hits.update(tree.query_ball_point(point, float(cutoff_A)))
        indices = sorted(hits)
    except ImportError:  # pragma: no cover
        distances = _distance_blocks(ligand, rec_xyz)
        indices = sorted(int(i) for i in np.nonzero((distances <= float(cutoff_A)).any(axis=0))[0])
    grouped: dict[tuple[str, int, str, str], dict[str, Any]] = {}
    for index in indices:
        atom = atoms[index]
        key = residue_key(atom)
        entry = grouped.setdefault(
            key,
            {
                "chain": key[0],
                "resseq": key[1],
                "icode": key[2],
                "resname": key[3],
                "residue": residue_label(key),
                "atom_count": 0,
                "polar_atom_count": 0,
                "min_distance_A": math.inf,
                "contact_atom_indices": [],
            },
        )
        point = np.asarray([atom.x, atom.y, atom.z], dtype=float)
        distance = float(np.sqrt(((ligand - point) ** 2).sum(axis=1)).min())
        if distance <= float(cutoff_A):
            entry["atom_count"] += 1
            if is_polar_atom(getattr(atom, "element", "")):
                entry["polar_atom_count"] += 1
            entry["min_distance_A"] = min(entry["min_distance_A"], distance)
            entry["contact_atom_indices"].append(int(getattr(atom, "serial", index + 1)))
    out = sorted(grouped.values(), key=lambda row: (row["min_distance_A"], row["chain"], row["resseq"], row["icode"]))
    for row in out:
        row["min_distance_A"] = round(float(row["min_distance_A"]), 3)
        row["residue_class"] = residue_class(row["resname"])
    return out


def residue_composition(contact_rows: list[dict[str, Any]]) -> dict[str, Any]:
    total = len(contact_rows)
    counts = {"hydrophobic": 0, "charged": 0, "polar": 0, "other": 0}
    for row in contact_rows:
        counts[str(row.get("residue_class") or "other")] = counts.get(str(row.get("residue_class") or "other"), 0) + 1
    denominator = total or 1
    return {
        "residue_count": total,
        "hydrophobic_residues": counts["hydrophobic"],
        "charged_residues": counts["charged"],
        "polar_residues": counts["polar"],
        "other_residues": counts["other"],
        "hydrophobic_fraction": round(counts["hydrophobic"] / denominator, 6) if total else None,
        "charged_fraction": round(counts["charged"] / denominator, 6) if total else None,
        "polar_fraction": round(counts["polar"] / denominator, 6) if total else None,
        "other_fraction": round(counts["other"] / denominator, 6) if total else None,
        "note": "composition is a geometric residue count, not solvent exposure or binding energy",
    }
