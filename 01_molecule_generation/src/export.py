"""Portable, strictly validated formal library exports (no synthetic provenance)."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Mapping

from .chemistry import (STANDARDIZATION_VERSION, calculate_properties, compound_id,
                        state_id, pose_id, provenance_id, read_sdf, same_molecule, standardize_smiles)

LIBRARY_FILENAME = "murj_predocking_library.csv"

REQUIRED_FIELDS = (
    "compound_id", "canonical_smiles", "standardization_version", "state_id", "pose_id",
    "provenance_id", "source", "route_status", "geometry_status", "qc_status", "sdf_path",
    "MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "formal_charge",
)
OPTIONAL_FIELDS = ("pKa_source", "microstate_method", "synthesis_evidence_id", "hypothesis_id",
                   "diversity_group", "parent_seed_id", "reinvent_run_id", "selection_reason")


def _read_csv(path):
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f"missing/duplicate CSV columns: {path.name}")
        rows = list(reader)
        if any(None in row for row in rows):
            raise ValueError(f"extra CSV values: {path.name}")
        return reader.fieldnames, rows


def _write_csv(path, rows, required=()):
    rows = [dict(r) for r in rows]
    fields = list(dict.fromkeys([*required, *(k for row in rows for k in row)]))
    if not fields:
        fields = ["id"]
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def validate_export(output_dir: str | Path) -> dict:
    root = Path(output_dir).resolve()
    errors = []
    count = 0
    try:
        columns, rows = _read_csv(root / LIBRARY_FILENAME)
        missing = set(REQUIRED_FIELDS) - set(columns)
        if missing:
            errors.append(f"missing required columns: {sorted(missing)}")
        if not rows:
            errors.append("formal library is empty")
        count = len(rows)
        seen = set()
        for index, row in enumerate(rows, 2):
            prefix = f"library.csv:{index}"
            try:
                empty = [k for k in REQUIRED_FIELDS if not str(row.get(k) or "").strip()]
                if empty:
                    raise ValueError(f"empty required values: {empty}")
                if row["standardization_version"] != STANDARDIZATION_VERSION:
                    raise ValueError("unsupported standardization version")
                canonical = standardize_smiles(row["canonical_smiles"])
                if canonical != row["canonical_smiles"]:
                    raise ValueError("SMILES is not canonical")
                if row["compound_id"] != compound_id(canonical) or row["state_id"] != state_id(canonical):
                    raise ValueError("compound/state ID mismatch")
                identity = row["pose_id"]
                if identity in seen:
                    raise ValueError("duplicate pose ID")
                seen.add(identity)
                if row["source"] not in {"Peptide2Mol", "REINVENT4", "REINVENT4/LibInvent"}:
                    raise ValueError("formal export requires real model source")
                if row["qc_status"] != "pass" or row["geometry_status"] != "pass":
                    raise ValueError("QC/geometry must pass")
                if row["route_status"] not in {"pass", "validated", "supported"}:
                    raise ValueError("route must have synthesis support")
                rel = Path(row["sdf_path"])
                sdf = (root / rel).resolve()
                if rel.is_absolute() or not sdf.is_relative_to(root):
                    raise ValueError("SDF must be a relative path within export directory")
                molecules = read_sdf(sdf)
                if len(molecules) != 1 or not same_molecule(molecules[0], canonical):
                    raise ValueError("CSV/SDF chemistry mismatch")
                for key in ("compound_id", "state_id", "pose_id", "provenance_id"):
                    if not molecules[0].HasProp(key) or molecules[0].GetProp(key) != row[key]:
                        raise ValueError(f"SDF/CSV {key} mismatch")
                for key, kind in (("pose_id", "pos"), ("provenance_id", "prv")):
                    if not re.fullmatch(kind + r"_[0-9a-f]{24}", row[key]):
                        raise ValueError(f"invalid {key}")
                frame = row.get("coordinate_frame", "")
                if pose_id(molecules[0], row["state_id"], frame) != row["pose_id"]:
                    raise ValueError("pose ID does not match coordinates")
                for key, value in calculate_properties(canonical).items():
                    if not math.isfinite(float(row[key])) or abs(float(row[key]) - value) > 0.001:
                        raise ValueError(f"descriptor mismatch: {key}")
            except Exception as exc:
                errors.append(f"{prefix}: {exc}")
        # Provenance is mandatory evidence, not just a model-name string.
        _, provenance = _read_csv(root / "provenance.csv")
        ids = [p.get("provenance_id") for p in provenance]
        if len(ids) != len(set(ids)):
            errors.append("duplicate provenance IDs")
        by_id = {p.get("provenance_id"): p for p in provenance}
        for row in rows:
            evidence = by_id.get(row.get("provenance_id"), {})
            needed = ("model", "checkpoint_sha256", "run_id", "seed")
            if any(not str(evidence.get(k) or "").strip() for k in needed):
                errors.append(f"missing real-model provenance: {row.get('provenance_id')}")
            digest = evidence.get("checkpoint_sha256", "")
            if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest.lower()):
                errors.append("invalid checkpoint SHA-256")
            if evidence.get("model") != row.get("source"):
                errors.append("provenance model/source mismatch")
    except Exception as exc:
        errors.append(str(exc))
    return {"valid": not errors, "errors": errors, "row_count": count,
            "standardization_version": STANDARDIZATION_VERSION}


def export_library(rows, output_dir: str | Path, side_tables: Mapping | None = None,
                   routes=None) -> dict:
    """Write supplied rows/evidence, returning paths and strict validation.

    Caller supplies already selected rows and SDFs relative to output_dir. Missing
    QC, synthesis evidence or model provenance is not invented by this exporter.
    Side tables are sequences of mappings keyed by simple table names, including
    the mandatory ``provenance`` table. No successful-stage marker is created.
    """
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    rows = [dict(row) for row in rows]
    paths = {"library": root / LIBRARY_FILENAME}
    _write_csv(paths["library"], rows, REQUIRED_FIELDS)
    tables = dict(side_tables or {})
    if routes is not None:
        tables["routes"] = routes
    for name, records in tables.items():
        if not isinstance(name, str) or not name.replace("_", "").isalnum() or name == "library":
            raise ValueError(f"unsafe/reserved side table name: {name!r}")
        path = root / f"{name}.csv"
        _write_csv(path, records)
        paths[name] = path
    validation = validate_export(root)
    report = root / "validation.json"
    report.write_text(json.dumps(validation, indent=2), encoding="utf-8")
    paths["validation"] = report
    return {"paths": {k: str(p) for k, p in paths.items()}, "validation": validation}
