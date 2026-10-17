"""Receptor inventory, component stripping, versions and local boxes.

No residue number is hard-coded: everything comes from the target profile or
from what is actually present in the file.  A receptor change that could move a
score (atom types, protonation, chains, pocket coordinates) must produce a new
``receptor_version`` so cached results invalidate.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any, Iterable

import numpy as np

from ..ids import receptor_version, sha256_file, sha256_text
from ..util import (
    AMINO_ACIDS,
    ION_NAMES,
    LIPID_NAMES,
    WATER_NAMES,
    PdbAtom,
    PdbStructure,
    atomic_write_json,
    centroid,
    parse_pdb_text,
    read_text,
    write_pdb,
)


@dataclass
class ComponentGroup:
    kind: str
    resname: str
    chain: str
    resseq: int
    icode: str
    atoms: int
    heavy_atoms: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "resname": self.resname,
            "chain": self.chain,
            "resseq": self.resseq,
            "icode": self.icode,
            "atoms": self.atoms,
            "heavy_atoms": self.heavy_atoms,
        }


@dataclass
class ReceptorInventory:
    path: str
    sha256: str
    models: list[int]
    chains: list[str]
    residues_per_chain: dict[str, dict[str, Any]]
    hetero: list[ComponentGroup]
    water_count: int
    hydrogen_count: int
    heavy_atom_count: int
    altloc_residues: list[dict[str, Any]]
    insertion_codes: list[dict[str, Any]]
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "sha256": self.sha256,
            "models": self.models,
            "chains": self.chains,
            "residues_per_chain": self.residues_per_chain,
            "hetero": [group.to_dict() for group in self.hetero],
            "water_count": self.water_count,
            "hydrogen_count": self.hydrogen_count,
            "heavy_atom_count": self.heavy_atom_count,
            "altloc_residues": self.altloc_residues,
            "insertion_codes": self.insertion_codes,
            "notes": self.notes,
        }


def load_structure(path: str) -> PdbStructure:
    text = read_text(path)
    if path.lower().endswith(".gz"):
        import gzip

        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            text = handle.read()
    return parse_pdb_text(text, keep_hydrogens=True)


def inventory(path: str, *, model: int | None = None) -> ReceptorInventory:
    structure = load_structure(path)
    if model is not None:
        structure = structure.by_model(model)
    elif len(structure.models) > 1:
        structure = structure.by_model(structure.models[0])

    chains: list[str] = []
    per_chain: dict[str, dict[str, Any]] = {}
    hetero: dict[tuple[str, int, str, str], ComponentGroup] = {}
    water = 0
    hydrogens = 0
    heavy = 0
    altloc: dict[tuple[str, int, str, str], set[str]] = {}
    insertion: list[dict[str, Any]] = []

    for atom in structure.atoms:
        upper = atom.resname.upper()
        if atom.is_hydrogen():
            hydrogens += 1
        else:
            heavy += 1
        if atom.altloc:
            altloc.setdefault((atom.chain, atom.resseq, atom.icode, atom.resname), set()).add(atom.altloc)
        if atom.icode:
            insertion.append({"chain": atom.chain, "resseq": atom.resseq, "icode": atom.icode, "resname": atom.resname})
        if atom.record.strip().upper() == "ATOM" and upper in AMINO_ACIDS:
            if atom.chain not in chains:
                chains.append(atom.chain)
            entry = per_chain.setdefault(
                atom.chain,
                {"residues": [], "first": atom.resseq, "last": atom.resseq, "atoms": 0, "gaps": []},
            )
            entry["atoms"] += 1
            if atom.resseq not in entry["residues"]:
                entry["residues"].append(atom.resseq)
            entry["first"] = min(entry["first"], atom.resseq)
            entry["last"] = max(entry["last"], atom.resseq)
        elif atom.record.strip().upper() == "HETATM":
            if upper in WATER_NAMES:
                water += 1
                continue
            key = (atom.chain, atom.resseq, atom.icode, atom.resname)
            group = hetero.setdefault(
                key,
                ComponentGroup(
                    kind=_classify_hetero(upper), resname=atom.resname, chain=atom.chain,
                    resseq=atom.resseq, icode=atom.icode, atoms=0, heavy_atoms=0,
                ),
            )
            group.atoms += 1
            if not atom.is_hydrogen():
                group.heavy_atoms += 1

    for entry in per_chain.values():
        entry["residues"].sort()
        residues = entry["residues"]
        entry["gaps"] = [
            [residues[i], residues[i + 1]] for i in range(len(residues) - 1) if residues[i + 1] - residues[i] > 1
        ]
        entry["count"] = len(residues)
        entry.pop("residues")

    return ReceptorInventory(
        path=path,
        sha256=sha256_file(path),
        models=structure.models,
        chains=chains,
        residues_per_chain=per_chain,
        hetero=list(hetero.values()),
        water_count=water,
        hydrogen_count=hydrogens,
        heavy_atom_count=heavy,
        altloc_residues=[
            {"chain": k[0], "resseq": k[1], "icode": k[2], "resname": k[3], "altlocs": sorted(v)}
            for k, v in altloc.items()
        ],
        insertion_codes=insertion,
        notes=[],
    )


def _classify_hetero(resname: str) -> str:
    if resname in WATER_NAMES:
        return "water"
    if resname in ION_NAMES:
        return "ion"
    if resname in LIPID_NAMES:
        return "lipid"
    if resname in AMINO_ACIDS:
        return "modified_residue"
    return "ligand"


def strip_receptor(
    structure: PdbStructure,
    *,
    keep_chains: Iterable[str] | None = None,
    keep_residues: Iterable[tuple[str, int]] | None = None,
    drop_water: bool = True,
    drop_ions: bool = True,
    drop_lipids: bool = True,
    drop_ligands: bool = True,
    drop_hydrogens: bool = False,
) -> tuple[PdbStructure, dict[str, Any]]:
    keep = set(keep_chains) if keep_chains else None
    allowed = {(str(chain), int(resseq)) for chain, resseq in keep_residues} if keep_residues else None
    removed: dict[str, int] = {}
    removed_residues: dict[str, set[tuple[str, int]]] = {}
    kept: list[PdbAtom] = []
    for atom in structure.atoms:
        upper = atom.resname.upper()
        reason = None
        if keep is not None and atom.chain not in keep:
            reason = "chain_not_kept"
        elif allowed is not None and (atom.chain, int(atom.resseq)) not in allowed:
            reason = "outside_selected_residues"
        elif drop_hydrogens and atom.is_hydrogen():
            reason = "hydrogen"
        elif atom.record.strip().upper() == "HETATM":
            if upper in WATER_NAMES and drop_water:
                reason = "water"
            elif upper in ION_NAMES and drop_ions:
                reason = "ion"
            elif upper in LIPID_NAMES and drop_lipids:
                reason = "lipid"
            elif drop_ligands and upper not in AMINO_ACIDS:
                reason = "ligand_or_additive"
        if reason:
            removed[reason] = removed.get(reason, 0) + 1
            removed_residues.setdefault(reason, set()).add((atom.chain, int(atom.resseq)))
        else:
            kept.append(atom)
    log = {
        "kept_atoms": len(kept),
        "removed": removed,
        "removed_residue_counts": {key: len(value) for key, value in sorted(removed_residues.items())},
        "keep_chains": sorted(keep) if keep else None,
        "keep_residue_range": (
            [min(r for _, r in allowed), max(r for _, r in allowed)] if allowed else None
        ),
    }
    return PdbStructure(kept, structure.models, structure.header_lines), log


def derive_receptor(
    src_path: str,
    out_path: str,
    *,
    keep_chains: Iterable[str] | None = None,
    keep_residues: Iterable[tuple[str, int]] | None = None,
    drop_water: bool = True,
    drop_ions: bool = True,
    drop_lipids: bool = True,
    drop_ligands: bool = True,
    drop_hydrogens: bool = False,
    recipe: str = "strip_components",
    coordinate_frame: str = "source_frame",
) -> dict[str, Any]:
    """Write a derived receptor and return its version record."""
    structure = load_structure(src_path)
    stripped, log = strip_receptor(
        structure,
        keep_chains=keep_chains,
        keep_residues=keep_residues,
        drop_water=drop_water,
        drop_ions=drop_ions,
        drop_lipids=drop_lipids,
        drop_ligands=drop_ligands,
        drop_hydrogens=drop_hydrogens,
    )
    if not stripped.atoms:
        raise ValueError("receptor stripping removed every atom; refusing to write an empty receptor")
    write_pdb(out_path, stripped.atoms)
    digest = sha256_file(out_path)
    version = receptor_version(digest, recipe, coordinate_frame)
    record = {
        "source": src_path,
        "source_sha256": sha256_file(src_path),
        "output": out_path,
        "output_sha256": digest,
        "version": version,
        "recipe": recipe,
        "coordinate_frame": coordinate_frame,
        "strip_log": log,
        "atoms": len(stripped.atoms),
    }
    atomic_write_json(os.path.splitext(out_path)[0] + ".version.json", record)
    return record


def protein_atoms(structure: PdbStructure) -> list[PdbAtom]:
    return [
        atom for atom in structure.atoms
        if atom.record.strip().upper() == "ATOM" and atom.resname.upper() in AMINO_ACIDS and not atom.is_hydrogen()
    ]


def atoms_for_residues(
    structure: PdbStructure, residues: Iterable[tuple[str, int]], *, include_hydrogens: bool = False
) -> list[PdbAtom]:
    wanted = {(chain, int(resseq)) for chain, resseq in residues}
    return [
        atom
        for atom in structure.atoms
        if (atom.chain, atom.resseq) in wanted and (include_hydrogens or not atom.is_hydrogen())
    ]


def local_box(
    atoms: list[PdbAtom],
    *,
    padding: float = 5.0,
    min_size: float = 14.0,
    max_size: float = 30.0,
    method: str = "region_bbox_plus_padding",
) -> dict[str, Any]:
    """Axis-aligned box around a local region, clipped to sane docking sizes."""
    if not atoms:
        raise ValueError("cannot build a box from zero atoms")
    coords = np.array([[a.x, a.y, a.z] for a in atoms], dtype=float)
    low = coords.min(axis=0) - padding
    high = coords.max(axis=0) + padding
    size = np.clip(high - low, min_size, max_size)
    center = (low + high) / 2.0
    return {
        "center": [round(float(v), 3) for v in center],
        "size": [round(float(v), 3) for v in size],
        "raw_extent": [round(float(v), 3) for v in (high - low)],
        "clipped": bool(np.any((high - low) < min_size) or np.any((high - low) > max_size)),
        "padding_per_side": padding,
        "units": "angstrom",
        "method": method,
        "atoms_used": len(atoms),
    }


def box_contains(box: dict[str, Any], coords: np.ndarray) -> tuple[int, int]:
    center = np.asarray(box["center"], dtype=float)
    size = np.asarray(box["size"], dtype=float)
    low = center - size / 2.0
    high = center + size / 2.0
    inside = int(np.all((coords >= low) & (coords <= high), axis=1).sum())
    return inside, int(len(coords))


def box_signature(box: dict[str, Any]) -> str:
    return sha256_text(f"{box.get('center')}|{box.get('size')}|{box.get('units')}")[:16]
