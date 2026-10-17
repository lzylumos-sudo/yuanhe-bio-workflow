"""Stable identifiers and content hashes.

Every identifier is derived from content with a documented recipe so that the
same input always produces the same id on any machine, and so that ids are never
derived from file order, directory order or ``hash()`` (which is salted per
process and therefore useless for reproducibility).
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from typing import Any, Iterable, Sequence

CHUNK = 1 << 20


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_text(text: str) -> str:
    return sha256_bytes(text.encode("utf-8"))


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def _digest(*parts: Any) -> str:
    payload = json.dumps(parts, ensure_ascii=False, sort_keys=True, default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def compound_id(canonical_isomeric_smiles: str) -> str:
    """Parent compound id.

    Definition: ``CMP`` + first 16 hex chars of sha256 over the *canonical
    isomeric* SMILES.  Stereochemistry is retained, so enantiomers never merge;
    protonation is handled one level down by :func:`chemical_state_id`.
    """
    return "CMP" + sha256_text(canonical_isomeric_smiles)[:16].upper()


def chemical_state_id(canonical_isomeric_smiles: str, formal_charge: int, protonation_note: str = "") -> str:
    return "STA" + _digest("state", canonical_isomeric_smiles, int(formal_charge), protonation_note)[:16].upper()


def docking_job_id(
    state_id: str,
    receptor_version: str,
    box_signature: str,
    protocol_signature: str,
    seed: int | None,
    engine: str,
    model: str | None,
) -> str:
    return "JOB" + _digest(
        "job", state_id, receptor_version, box_signature, protocol_signature, seed, engine, model
    )[:16].upper()


def pose_id(job_id: str, sdf_record_index: int | None, coords_digest: str, rank: int | None) -> str:
    """Pose id binds the job, the record position inside its source and the
    rounded coordinates that were actually scored."""
    return "POS" + _digest("pose", job_id, sdf_record_index, coords_digest, rank)[:16].upper()


def receptor_version(
    receptor_sha256: str,
    prep_recipe: str,
    coordinate_frame: str,
    mutant_recipe: str = "",
) -> str:
    """A receptor *version*, not just a file hash.

    Any change in protonation, chain content, pocket coordinates or geometry has
    to move this value so stale scores are re-evaluated.
    """
    return "REC" + _digest("rec", receptor_sha256, prep_recipe, coordinate_frame, mutant_recipe)[:16].upper()


def selection_version(members: Iterable[str], rules_digest: str) -> str:
    ordered = sorted(members)
    return "SEL" + _digest("sel", ordered, rules_digest)[:16].upper()


def coords_digest(coords: Sequence[Sequence[float]], ndigits: int = 3) -> str:
    """Hash rounded coordinates.

    Rounding to 1e-3 A is far below any docking precision change but stops
    float formatting noise from creating new pose ids for identical geometry.
    """
    flat: list[str] = []
    for xyz in coords:
        flat.append(",".join(f"{c:.{ndigits}f}" for c in xyz))
    return sha256_text("|".join(flat))


def stable_key(*values: Any) -> str:
    """Tie-break key: deterministic, independent of dict/list ordering."""
    return _digest("tie", *values)


def cache_key(stage: str, input_hash: str, config_hash: str, code_version: str, tool_versions: dict) -> str:
    return _digest("cache", stage, input_hash, config_hash, code_version, tool_versions)


def short_hash(text: str, length: int = 12) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:length].upper()


def quantize(value: float | None, digits: int = 4) -> float | None:
    if value is None:
        return None
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    return round(float(value), digits)


def file_fingerprint(path: str) -> dict[str, Any]:
    stat = os.stat(path)
    return {
        "path": path,
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": sha256_file(path),
    }
