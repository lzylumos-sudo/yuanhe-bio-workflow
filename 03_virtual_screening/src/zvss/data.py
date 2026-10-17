"""Shared accessors for stage implementations."""

from __future__ import annotations

import functools
import os
from typing import Any

import numpy as np
from rdkit import Chem

from .chem import rdio
from .util import PdbAtom, parse_pdb_text, read_json, read_text


def pose_rows(ctx, *, run_only: bool = True) -> list[dict[str, Any]]:
    sql = """
        SELECT p.*, s.canonical_smiles, s.formal_charge, c.inchikey, c.formula, c.heavy_atoms,
               j.receptor_version, j.box_json AS box, j.protocol_json AS protocol,
               j.seed, j.engine, j.engine_version,
               j.model_id, j.scoring_function, j.receptor_sha256
        FROM pose p
        LEFT JOIN chemical_state s ON s.state_id = p.state_id
        LEFT JOIN compound c ON c.compound_id = p.compound_id
        LEFT JOIN docking_job j ON j.job_id = p.job_id
        ORDER BY p.pose_id
    """
    rows = ctx.db.query(sql)
    out = []
    for row in rows:
        entry = dict(row)
        entry["scores"] = read_json_text(row["scores"]) if isinstance(row["scores"], str) else (row["scores"] or {})
        entry["box"] = read_json_text(row["box"]) if isinstance(row["box"], str) else (row["box"] or {})
        entry["protocol"] = read_json_text(row["protocol"]) if isinstance(row["protocol"], str) else (row["protocol"] or {})
        entry["meta"] = read_json_text(row["meta"]) if isinstance(row["meta"], str) else (row["meta"] or {})
        out.append(entry)
    return out


def read_json_text(text: str) -> Any:
    import json

    try:
        return json.loads(text) if text else {}
    except (TypeError, ValueError):
        return {}


@functools.lru_cache(maxsize=512)
def _mol_from_path(path: str) -> Chem.Mol | None:
    if not path or not os.path.exists(path):
        return None
    for mol, _props in rdio.iter_sdf(path):
        if mol is not None:
            return mol
    return None


def pose_mol(ctx, pose: dict[str, Any]) -> Chem.Mol | None:
    path = pose.get("coords_path")
    if not path:
        candidate = os.path.join(ctx.poses_dir(), f"{pose['pose_id']}.sdf")
        path = candidate if os.path.exists(candidate) else ""
    return _mol_from_path(path) if path else None


def receptor_record(ctx) -> dict[str, Any]:
    path = os.path.join(ctx.stage_dir("import"), "receptor.json")
    record = read_json(path, {}) or {}
    if not record:
        record = read_json(os.path.join(ctx.stage_dir("target"), "receptor.json"), {}) or {}
    return record


def receptor_atoms(ctx) -> tuple[list[PdbAtom], str]:
    record = receptor_record(ctx)
    path = record.get("archived") or record.get("output") or record.get("path") or ""
    if not path or not os.path.exists(path):
        return [], ""
    structure = parse_pdb_text(read_text(path), keep_hydrogens=False)
    atoms = [
        a for a in structure.atoms
        if a.record.upper().startswith("ATOM") and not a.is_hydrogen()
    ]
    return atoms, path


def reference_structure(ctx) -> tuple[list[PdbAtom], dict[str, Any]]:
    record = read_json(os.path.join(ctx.ref_dir(), "reference.json"), {}) or {}
    path = record.get("path") or ""
    if not path or not os.path.exists(path):
        return [], record
    structure = parse_pdb_text(read_text(path), keep_hydrogens=False)
    return structure.atoms, record


def coords_of(atoms: list[PdbAtom]) -> np.ndarray:
    return np.array([[a.x, a.y, a.z] for a in atoms], dtype=float)
