"""PDB parsing and auditable structure preparation for the v2 pipeline.

The implementation deliberately has no chemistry-toolkit dependency.  It preserves
PDB author residue identifiers, chooses alternate conformers by occupancy, and uses
COMPND declarations (rather than chain length) to identify deposited entities.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import asdict, dataclass
from datetime import date
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class Atom:
    record: str
    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    resseq: int
    icode: str
    xyz: tuple[float, float, float]
    occupancy: float
    bfactor: float
    element: str
    charge: str
    model: int = 1

    @property
    def residue_id(self) -> tuple[str, int, str]:
        return self.chain, self.resseq, self.icode


@dataclass(frozen=True)
class Entity:
    mol_id: str
    molecule: str
    chains: tuple[str, ...]


@dataclass
class PDBStructure:
    atoms: list[Atom]
    entities: list[Entity]
    source: str | None = None
    selected_model: int = 1
    available_models: tuple[int, ...] = (1,)
    discarded_altloc_atoms: int = 0

    def residues(self, atoms: Iterable[Atom] | None = None) -> dict[tuple[str, int, str], list[Atom]]:
        grouped: dict[tuple[str, int, str], list[Atom]] = {}
        for atom in self.atoms if atoms is None else atoms:
            grouped.setdefault(atom.residue_id, []).append(atom)
        return grouped


def _safe_int(text: str, default: int = 0) -> int:
    try:
        return int(text.strip())
    except ValueError:
        return default


def _safe_float(text: str, default: float = 0.0) -> float:
    try:
        return float(text.strip())
    except ValueError:
        return default


def parse_compnd(lines: Iterable[str]) -> list[Entity]:
    """Parse COMPND MOL_ID/MOLECULE/CHAIN declarations, including continuations."""
    payload = " ".join(line[10:].strip() for line in lines if line.startswith("COMPND"))
    entities: list[Entity] = []
    current: dict[str, str] = {}
    for item in payload.split(";"):
        item = item.strip()
        if not item or ":" not in item:
            continue
        key, value = (part.strip() for part in item.split(":", 1))
        key = key.upper()
        if key == "MOL_ID" and current:
            entities.append(_entity_from_fields(current))
            current = {}
        current[key] = value
    if current:
        entities.append(_entity_from_fields(current))
    return entities


def _entity_from_fields(fields: Mapping[str, str]) -> Entity:
    chains = tuple(c.strip() for c in fields.get("CHAIN", "").split(",") if c.strip())
    return Entity(fields.get("MOL_ID", ""), fields.get("MOLECULE", ""), chains)


def entity_for_chain(entities: Sequence[Entity], chain: str) -> Entity:
    matches = [entity for entity in entities if chain in entity.chains]
    if len(matches) != 1:
        raise ValueError(f"chain {chain!r} must occur in exactly one COMPND entity; found {len(matches)}")
    return matches[0]


def _parse_atom(line: str, model: int) -> Atom:
    return Atom(
        record=line[0:6].strip(), serial=_safe_int(line[6:11]), name=line[12:16].strip(),
        altloc=line[16:17].strip(), resname=line[17:20].strip(), chain=line[21:22].strip(),
        resseq=_safe_int(line[22:26]), icode=line[26:27].strip(),
        xyz=(_safe_float(line[30:38]), _safe_float(line[38:46]), _safe_float(line[46:54])),
        occupancy=_safe_float(line[54:60], 1.0), bfactor=_safe_float(line[60:66]),
        element=(line[76:78].strip() if len(line) >= 78 else ""),
        charge=(line[78:80].strip() if len(line) >= 80 else ""), model=model,
    )


def parse_pdb(path: str | Path, model: int = 1) -> PDBStructure:
    """Read one PDB model and resolve altlocs by highest occupancy.

    Altloc selection is per atom name within an author residue.  Occupancy ties
    prefer blank, then ``A``, then lexical order, making selection deterministic.
    """
    path = Path(path)
    lines = path.read_text(encoding="utf-8").splitlines()
    entities = parse_compnd(lines)
    current_model = 1
    explicit_models: list[int] = []
    candidates: list[Atom] = []
    for line in lines:
        record = line[0:6].strip()
        if record == "MODEL":
            current_model = _safe_int(line[10:14], len(explicit_models) + 1)
            explicit_models.append(current_model)
        elif record in {"ATOM", "HETATM"} and current_model == model:
            candidates.append(_parse_atom(line, current_model))
    available = tuple(dict.fromkeys(explicit_models)) if explicit_models else (1,)
    if model not in available:
        raise ValueError(f"model {model} not present; available models: {available}")
    chosen: dict[tuple[str, int, str, str, str], Atom] = {}
    order: list[tuple[str, int, str, str, str]] = []
    for atom in candidates:
        key = (atom.chain, atom.resseq, atom.icode, atom.resname, atom.name)
        if key not in chosen:
            chosen[key] = atom
            order.append(key)
            continue
        old = chosen[key]
        def rank(a: Atom) -> tuple[float, int, int, str]:
            return (a.occupancy, int(a.altloc == ""), int(a.altloc == "A"), chr(255) if a.altloc == "" else "".join(chr(255-ord(c)) for c in a.altloc))
        if rank(atom) > rank(old):
            chosen[key] = atom
    atoms = [chosen[key] for key in order]
    return PDBStructure(atoms, entities, str(path), model, available, len(candidates) - len(atoms))


def group_residues(atoms: Iterable[Atom]) -> dict[tuple[str, int, str], list[Atom]]:
    grouped: dict[tuple[str, int, str], list[Atom]] = {}
    for atom in atoms:
        grouped.setdefault(atom.residue_id, []).append(atom)
    return grouped


def _selected_peptide_atoms(atoms: Sequence[Atom], selected: Sequence[int | str | tuple[int, str]] | None) -> list[Atom]:
    if selected is None:
        return list(atoms)
    wanted: set[tuple[int, str]] = set()
    for item in selected:
        if isinstance(item, tuple):
            wanted.add((int(item[0]), str(item[1]).strip()))
        elif isinstance(item, int):
            wanted.add((item, ""))
        else:
            text = str(item).strip()
            digits = text.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
            wanted.add((int(digits), text[len(digits):]))
    return [a for a in atoms if (a.resseq, a.icode) in wanted]


def crop_complete_residues(receptor_atoms: Sequence[Atom], peptide_atoms: Sequence[Atom], distance_A: float = 6.0) -> list[Atom]:
    """Return every atom of receptor residues having any atom within ``distance_A``."""
    if distance_A <= 0:
        raise ValueError("distance_A must be positive")
    if not peptide_atoms:
        raise ValueError("selected peptide region has no atoms")
    peptide_xyz = np.asarray([a.xyz for a in peptide_atoms], dtype=float)
    kept: set[tuple[str, int, str]] = set()
    for residue_id, atoms in group_residues(receptor_atoms).items():
        xyz = np.asarray([a.xyz for a in atoms], dtype=float)
        if np.any(np.sum((xyz[:, None, :] - peptide_xyz[None, :, :]) ** 2, axis=2) <= distance_A ** 2):
            kept.add(residue_id)
    return [a for a in receptor_atoms if a.residue_id in kept]


def terminal_carboxylate_status(residue_atoms: Sequence[Atom], explicit_terminal: bool) -> dict[str, object]:
    """Report terminal-carboxylate evidence without inferring it from backbone O."""
    names = {a.name for a in residue_atoms}
    has_oxt = "OXT" in names
    supported = bool(has_oxt and explicit_terminal)
    reason = "OXT_and_explicit_terminal" if supported else ("OXT_without_explicit_terminal" if has_oxt else "no_OXT")
    return {"has_OXT": has_oxt, "explicit_terminal": bool(explicit_terminal),
            "terminal_carboxylate": supported, "reason": reason}


def _format_atom(atom: Atom, serial: int) -> str:
    element = atom.element or next((c for c in atom.name if c.isalpha()), "")
    # Prepared files contain the chosen conformer only, so altloc is intentionally blank.
    return (f"{atom.record:<6}{serial:5d} {atom.name:>4s} {atom.resname:>3s} {atom.chain:1s}"
            f"{atom.resseq:4d}{atom.icode:1s}   {atom.xyz[0]:8.3f}{atom.xyz[1]:8.3f}{atom.xyz[2]:8.3f}"
            f"{atom.occupancy:6.2f}{atom.bfactor:6.2f}          {element:>2s}{atom.charge:>2s}\n")


def write_pdb(path: str | Path, atoms: Sequence[Atom], remarks: Sequence[str] = ()) -> Path:
    path = Path(path)
    with path.open("w", encoding="utf-8") as handle:
        for remark in remarks:
            handle.write(f"REMARK 950 {remark}\n")
        previous: tuple[str, int, str] | None = None
        serial = 1
        for atom in atoms:
            if previous is not None and atom.chain != previous[0]:
                handle.write("TER\n")
            handle.write(_format_atom(atom, serial))
            serial += 1
            previous = atom.residue_id
        handle.write("TER\nEND\n")
    return path


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def prepare_structure(
    pdb_path: str | Path,
    output_dir: str | Path,
    *,
    structure_id: str = "9NU4",
    model: int = 1,
    receptor_chain: str = "A",
    peptide_chain: str = "B",
    receptor_start: int = 5,
    receptor_end: int = 500,
    selected_peptide_residues: Sequence[int | str | tuple[int, str]] | None = None,
    crop_distance_A: float | None = 6.0,
    peptide_explicit_c_terminus: bool = False,
    source_url: str = "https://www.rcsb.org/structure/9NU4",
) -> dict[str, object]:
    """Prepare receptor/peptide PDBs plus auditable mapping, manifest and checksums."""
    pdb_path, output_dir = Path(pdb_path), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    structure = parse_pdb(pdb_path, model=model)
    receptor_entity = entity_for_chain(structure.entities, receptor_chain)
    peptide_entity = entity_for_chain(structure.entities, peptide_chain)
    receptor_range_atoms = [a for a in structure.atoms if a.chain == receptor_chain and receptor_start <= a.resseq <= receptor_end]
    peptide_atoms = [a for a in structure.atoms if a.chain == peptide_chain]
    if not receptor_range_atoms or not peptide_atoms:
        raise ValueError("configured receptor or peptide selection is empty")
    evidence_peptide = _selected_peptide_atoms(peptide_atoms, selected_peptide_residues)
    receptor_atoms = (crop_complete_residues(receptor_range_atoms, evidence_peptide, crop_distance_A)
                      if crop_distance_A is not None else receptor_range_atoms)
    full_receptor_path = write_pdb(output_dir / "receptor_full.pdb", receptor_range_atoms,
                                   [f"SOURCE {structure_id}; SELECTED RECEPTOR ONLY; NOT PROTONATED"])
    receptor_path = write_pdb(output_dir / "receptor_prepared.pdb", receptor_atoms,
                              [f"SOURCE {structure_id} MODEL {model}", f"ENTITY {receptor_entity.mol_id}: {receptor_entity.molecule}"])
    peptide_path = write_pdb(output_dir / "peptide_reference.pdb", peptide_atoms,
                             [f"SOURCE {structure_id} MODEL {model}", f"ENTITY {peptide_entity.mol_id}: {peptide_entity.molecule}",
                              "REFERENCE COORDINATES ONLY; NO DISTANCE-INFERRED BOND ORDERS"])

    map_rows: list[dict[str, object]] = []
    for role, atoms, entity in (("receptor", receptor_atoms, receptor_entity), ("peptide", peptide_atoms, peptide_entity)):
        for index, ((chain, resseq, icode), residue_atoms) in enumerate(group_residues(atoms).items(), 1):
            map_rows.append({"structure_id": structure_id, "model": model, "role": role,
                             "entity_id": entity.mol_id, "molecule": entity.molecule, "auth_chain": chain,
                             "auth_resseq": resseq, "auth_icode": icode, "resname": residue_atoms[0].resname,
                             "prepared_chain": chain, "prepared_resseq": resseq, "prepared_icode": icode,
                             "canonical_index": index})
    residue_map_path = output_dir / "residue_map.csv"
    _write_csv(residue_map_path, list(map_rows[0]), map_rows)

    manifest_path = output_dir / "structures_manifest.csv"
    manifest_fields = ["structure_id", "source_file", "source_url", "retrieval_date", "sha256", "model",
                       "receptor_entity_id", "receptor_chain", "receptor_range", "peptide_entity_id", "peptide_chain",
                       "prepared_receptor", "peptide_reference", "crop_distance_A"]
    _write_csv(manifest_path, manifest_fields, [{
        "structure_id": structure_id, "source_file": pdb_path.name, "source_url": source_url,
        "retrieval_date": date.today().isoformat(), "sha256": sha256_file(pdb_path), "model": model,
        "receptor_entity_id": receptor_entity.mol_id, "receptor_chain": receptor_chain,
        "receptor_range": f"{receptor_start}-{receptor_end}", "peptide_entity_id": peptide_entity.mol_id,
        "peptide_chain": peptide_chain, "prepared_receptor": receptor_path.name,
        "peptide_reference": peptide_path.name, "crop_distance_A": crop_distance_A if crop_distance_A is not None else "",
    }])

    peptide_residues = list(group_residues(peptide_atoms).values())
    terminal = terminal_carboxylate_status(peptide_residues[-1], peptide_explicit_c_terminus)
    report: dict[str, object] = {
        "structure_id": structure_id, "source_file": str(pdb_path), "selected_model": model,
        "available_models": list(structure.available_models), "entity_assignment_basis": "PDB_COMPND_CHAIN",
        "receptor": {"chain": receptor_chain, "entity": asdict(receptor_entity), "configured_range": [receptor_start, receptor_end],
                     "range_atom_count": len(receptor_range_atoms), "prepared_atom_count": len(receptor_atoms),
                     "prepared_residue_count": len(group_residues(receptor_atoms))},
        "peptide": {"chain": peptide_chain, "entity": asdict(peptide_entity), "atom_count": len(peptide_atoms),
                    "residue_count": len(peptide_residues), "selected_contact_residues": list(map(str, selected_peptide_residues or [])),
                    "c_terminus": terminal},
        "crop": {"distance_A": crop_distance_A, "whole_residue_selection": True},
        "alternate_locations": {"selection": "highest_occupancy_blank_then_A_tiebreak", "discarded_atom_records": structure.discarded_altloc_atoms},
        "warnings": ([] if terminal["terminal_carboxylate"] else ["Peptide terminal carboxylate not assigned: requires OXT and explicit biological-terminal status."]),
    }
    report_path = output_dir / "preparation_report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checksum_targets = [full_receptor_path, receptor_path, peptide_path, manifest_path, residue_map_path, report_path]
    checksums = {p.name: sha256_file(p) for p in checksum_targets}
    checksum_path = output_dir / "checksums.sha256"
    checksum_path.write_text("".join(f"{digest}  {name}\n" for name, digest in sorted(checksums.items())), encoding="utf-8")
    report["checksums"] = checksums
    return {"receptor_full": full_receptor_path, "receptor_prepared": receptor_path, "peptide_reference": peptide_path,
            "structures_manifest": manifest_path, "residue_map": residue_map_path,
            "preparation_report": report_path, "checksums": checksum_path, "report": report}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Prepare a PDB receptor/peptide complex")
    parser.add_argument("pdb")
    parser.add_argument("output_dir")
    parser.add_argument("--structure-id", default="9NU4")
    parser.add_argument("--model", type=int, default=1)
    parser.add_argument("--receptor-chain", default="A")
    parser.add_argument("--peptide-chain", default="B")
    parser.add_argument("--receptor-start", type=int, default=5)
    parser.add_argument("--receptor-end", type=int, default=500)
    parser.add_argument("--crop-distance", type=float, default=6.0)
    parser.add_argument("--peptide-residues", nargs="*")
    parser.add_argument("--peptide-explicit-c-terminus", action="store_true")
    args = parser.parse_args(argv)
    prepare_structure(args.pdb, args.output_dir, structure_id=args.structure_id, model=args.model,
                      receptor_chain=args.receptor_chain, peptide_chain=args.peptide_chain,
                      receptor_start=args.receptor_start, receptor_end=args.receptor_end,
                      selected_peptide_residues=args.peptide_residues,
                      crop_distance_A=args.crop_distance,
                      peptide_explicit_c_terminus=args.peptide_explicit_c_terminus)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
