"""Coordinate-derived local binding hypotheses for prepared MurJ complexes."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

try:  # Supports both ``python -m src.hypotheses`` and a direct script invocation.
    from .structure import Atom, group_residues, parse_pdb
except ImportError:  # pragma: no cover
    from structure import Atom, group_residues, parse_pdb


HYDROPHOBIC_RESIDUES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "TYR", "PRO"}
POLAR_ATOMS = {"N", "O", "S"}
BACKBONE = {"N", "CA", "C", "O", "OXT"}


def _element(atom: Atom) -> str:
    return (atom.element or next((c for c in atom.name if c.isalpha()), "")).upper()


def _locator(atom: Atom) -> str:
    return f"{atom.chain}:{atom.resname}{atom.resseq}{atom.icode}:{atom.name}"


def _center(atoms: Sequence[Atom]) -> list[float]:
    return np.mean(np.asarray([a.xyz for a in atoms], dtype=float), axis=0).round(3).tolist()


def _direction(source: Sequence[Atom], target: Sequence[Atom]) -> list[float] | None:
    vector = np.asarray(_center(target)) - np.asarray(_center(source))
    norm = float(np.linalg.norm(vector))
    return None if norm == 0 else (vector / norm).round(4).tolist()


def _contacts(peptide: Sequence[Atom], receptor: Sequence[Atom], cutoff_A: float) -> list[tuple[Atom, Atom, float]]:
    if not peptide or not receptor:
        return []
    pxyz = np.asarray([a.xyz for a in peptide])
    rxyz = np.asarray([a.xyz for a in receptor])
    distances = np.sqrt(np.sum((pxyz[:, None, :] - rxyz[None, :, :]) ** 2, axis=2))
    pi, ri = np.where(distances <= cutoff_A)
    return [(peptide[int(i)], receptor[int(j)], float(distances[i, j])) for i, j in zip(pi, ri)]


def _sidechain_atoms(atoms: Sequence[Atom]) -> list[Atom]:
    sidechain = [a for a in atoms if a.name not in BACKBONE and _element(a) != "H"]
    return sidechain or [a for a in atoms if _element(a) != "H"]


def _feature(feature_type: str, source_atoms: Sequence[Atom], receptor_atoms: Sequence[Atom], *,
             radius_A: float, required: bool, confidence: str) -> dict[str, object]:
    return {
        "type": feature_type,
        "center_xyz": _center(source_atoms),
        "radius_A": radius_A,
        "direction": _direction(source_atoms, receptor_atoms),
        "required": required,
        "evidence_locator": {
            "source_atoms": [_locator(a) for a in source_atoms],
            "contact_receptor_atoms": [_locator(a) for a in receptor_atoms],
            "coordinate_basis": "deposited_atom_coordinates",
        },
        "confidence": confidence,
    }


def _exclusion_spheres(receptor: Sequence[Atom], centers: Sequence[Sequence[float]], radius_A: float = 1.5,
                       shell_A: float = 6.0, max_spheres: int = 96) -> list[dict[str, object]]:
    if not centers:
        return []
    center_array = np.asarray(centers, dtype=float)
    candidates: list[tuple[float, Atom]] = []
    for atom in receptor:
        if _element(atom) == "H":
            continue
        distance = float(np.min(np.linalg.norm(center_array - np.asarray(atom.xyz), axis=1)))
        if distance <= shell_A:
            candidates.append((distance, atom))
    candidates.sort(key=lambda item: item[0])
    return [{"center_xyz": list(atom.xyz), "radius_A": radius_A, "evidence_locator": _locator(atom)}
            for _, atom in candidates[:max_spheres]]


def derive_hypotheses(
    pdb_path: str | Path,
    *,
    structure_id: str = "9NU4",
    model: int = 1,
    receptor_chain: str = "A",
    peptide_chain: str = "B",
    receptor_start: int = 5,
    receptor_end: int = 500,
    contact_cutoff_A: float = 4.5,
    feature_radius_A: float = 2.0,
    frame_id: str = "9NU4_deposited_model_1",
    residue_mapping_version: str = "auth-v1",
    threshold_version: str = "murj-local-v1",
) -> dict[str, object]:
    """Derive H1 and H2 from measured source-peptide/receptor contacts.

    H1 uses the best-contacting hydrophobic peptide side chain. H2 combines that
    observed region with the best-contacting polar peptide atom. No salt bridge or
    terminal charge is inferred.
    """
    structure = parse_pdb(pdb_path, model=model)
    receptor = [a for a in structure.atoms if a.chain == receptor_chain and receptor_start <= a.resseq <= receptor_end]
    peptide = [a for a in structure.atoms if a.chain == peptide_chain]
    contacts = _contacts(peptide, receptor, contact_cutoff_A)
    if not contacts:
        raise ValueError(f"no peptide/receptor contacts within {contact_cutoff_A} A")

    peptide_residues = group_residues(peptide)
    contact_by_residue: dict[tuple[str, int, str], list[tuple[Atom, Atom, float]]] = {}
    for pa, ra, distance in contacts:
        contact_by_residue.setdefault(pa.residue_id, []).append((pa, ra, distance))
    hydro_candidates: list[tuple[float, int, tuple[str, int, str]]] = []
    for rid, atoms in peptide_residues.items():
        if atoms[0].resname in HYDROPHOBIC_RESIDUES and rid in contact_by_residue:
            values = contact_by_residue[rid]
            hydro_candidates.append((min(v[2] for v in values), -len(values), rid))
    if not hydro_candidates:
        raise ValueError("source peptide has no contacting hydrophobic residue")
    hydro_rid = min(hydro_candidates)[2]
    hydro_source = _sidechain_atoms(peptide_residues[hydro_rid])
    hydro_contact_atoms = list(dict.fromkeys(v[1] for v in contact_by_residue[hydro_rid]))

    polar_pairs = [(pa, ra, d) for pa, ra, d in contacts
                   if _element(pa) in POLAR_ATOMS and pa.name not in {"C", "CA"}]
    if not polar_pairs:
        raise ValueError("source peptide has no coordinate-supported polar contact")
    polar_atom, polar_receptor_atom, polar_distance = min(polar_pairs, key=lambda item: item[2])

    hydro_feature = _feature("hydrophobic", hydro_source, hydro_contact_atoms,
                             radius_A=feature_radius_A, required=True, confidence="coordinate_supported")
    polar_kind = "negative_ionizable_or_acceptor" if polar_atom.resname in {"ASP", "GLU"} and polar_atom.name.startswith(("OD", "OE")) else "polar_donor_or_acceptor"
    polar_feature = _feature(polar_kind, [polar_atom], [polar_receptor_atom],
                             radius_A=feature_radius_A, required=True, confidence="coordinate_supported_type_ambiguous")
    centers = [hydro_feature["center_xyz"], polar_feature["center_xyz"]]
    exclusions = _exclusion_spheres(receptor, centers)
    common = {
        "source_structure_ids": [structure_id], "residue_mapping_version": residue_mapping_version,
        "frame_id": frame_id, "comparable_receptors": [structure_id],
        "threshold_version": threshold_version,
        "thresholds": {"contact_cutoff_A": contact_cutoff_A, "feature_radius_A": feature_radius_A,
                       "exclusion_radius_A": 1.5},
        "exclusion_spheres": exclusions,
    }
    h1 = {**common, "hypothesis_id": "H1", "name": "hydrophobic_cleft_occupancy",
          "features": [hydro_feature],
          "pass_rule": "all required features matched within radius and no severe exclusion-sphere overlap",
          "evidence_summary": {"source_residue": f"{hydro_rid[0]}:{peptide_residues[hydro_rid][0].resname}{hydro_rid[1]}{hydro_rid[2]}",
                               "minimum_contact_A": round(min(v[2] for v in contact_by_residue[hydro_rid]), 3)}}
    h2 = {**common, "hypothesis_id": "H2", "name": "polar_hydrophobic_cooperation",
          "features": [hydro_feature, polar_feature],
          "pass_rule": "hydrophobic and polar required features matched within radii and no severe exclusion-sphere overlap",
          "evidence_summary": {"hydrophobic_source_residue": h1["evidence_summary"]["source_residue"],
                               "polar_source_atom": _locator(polar_atom), "polar_contact_atom": _locator(polar_receptor_atom),
                               "polar_contact_A": round(polar_distance, 3),
                               "interpretation_limit": "coordinate contact; not asserted to be a salt bridge or validated pharmacophore"}}
    return {"schema_version": "hypotheses-v1", "structure_id": structure_id, "frame_id": frame_id,
            "threshold_version": threshold_version, "hypotheses": [h1, h2]}


def generate_hypotheses(pdb_path: str | Path, output_json: str | Path, **kwargs: object) -> dict[str, object]:
    result = derive_hypotheses(pdb_path, **kwargs)
    output_json = Path(output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Derive local hypotheses from peptide contacts")
    parser.add_argument("pdb")
    parser.add_argument("output_json")
    parser.add_argument("--structure-id", default="9NU4")
    parser.add_argument("--model", type=int, default=1)
    parser.add_argument("--receptor-chain", default="A")
    parser.add_argument("--peptide-chain", default="B")
    parser.add_argument("--contact-cutoff", type=float, default=4.5)
    parser.add_argument("--feature-radius", type=float, default=2.0)
    args = parser.parse_args(argv)
    generate_hypotheses(args.pdb, args.output_json, structure_id=args.structure_id, model=args.model,
                        receptor_chain=args.receptor_chain, peptide_chain=args.peptide_chain,
                        contact_cutoff_A=args.contact_cutoff, feature_radius_A=args.feature_radius,
                        frame_id=f"{args.structure_id}_deposited_model_{args.model}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
