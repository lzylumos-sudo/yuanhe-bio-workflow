"""Membrane-protein pocket analysis for Y-MAVS.

The pocket is defined by an actual docked ligand and the receptor atoms around it,
not by a new blind pocket predictor.  The module reports membrane depth, residue
composition, annotated transmembrane-segment contact and a small set of geometric
descriptors.  It does not compute a cavity volume when no validated volume method
is available, and it never invents a transmembrane helix from distance alone.
"""

from __future__ import annotations

import csv
import math
import os
from typing import Any, Iterable

import numpy as np

from ..util import atomic_write_json, now_iso
from .geom import (
    receptor_contacts_by_residue,
    residue_composition,
    residue_key,
)
from .regions import MembraneGeometry, MembraneRegionAnalyzer, REGION_UNKNOWN

ALGORITHM_VERSION = "ymavs-pocket/1.0.0"
UNKNOWN = "unknown"


class MembranePocketAnalyzer:
    """Local, ligand-defined pocket analysis in a verified membrane frame."""

    def __init__(
        self,
        receptor_atoms: Iterable[Any],
        geometry: MembraneGeometry,
        region_analyzer: MembraneRegionAnalyzer,
        *,
        contact_cutoff_A: float = 5.0,
        embedded_core_fraction: float = 0.60,
        interface_fraction: float = 0.50,
        aqueous_outside_fraction: float = 0.60,
        mixed_min_fraction: float = 0.20,
        annotated_segments: dict[str, list[tuple[str, int, str]]] | None = None,
    ):
        self.receptor_atoms = list(receptor_atoms)
        self.geometry = geometry
        self.regions = region_analyzer
        self.contact_cutoff_A = float(contact_cutoff_A)
        if self.contact_cutoff_A <= 0 or not math.isfinite(self.contact_cutoff_A):
            raise ValueError("pocket contact_cutoff_A must be finite and positive")
        self.embedded_core_fraction = float(embedded_core_fraction)
        self.interface_fraction = float(interface_fraction)
        self.aqueous_outside_fraction = float(aqueous_outside_fraction)
        self.mixed_min_fraction = float(mixed_min_fraction)
        for name, value in (
            ("embedded_core_fraction", self.embedded_core_fraction),
            ("interface_fraction", self.interface_fraction),
            ("aqueous_outside_fraction", self.aqueous_outside_fraction),
            ("mixed_min_fraction", self.mixed_min_fraction),
        ):
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be in [0, 1]")
        self.annotated_segments = annotated_segments or {}
        self._residue_coords: dict[tuple[str, int, str, str], np.ndarray] = {}
        for atom in self.receptor_atoms:
            if hasattr(atom, "is_hydrogen") and atom.is_hydrogen():
                continue
            key = residue_key(atom)
            self._residue_coords.setdefault(key, []).append([atom.x, atom.y, atom.z])  # type: ignore[arg-type]
        self._residue_coords = {key: np.asarray(value, dtype=float) for key, value in self._residue_coords.items()}

    def analyze(
        self,
        ligand_xyz: np.ndarray,
        *,
        pose_id: str = "",
        compound_id: str = "",
        box_center: Iterable[float] | None = None,
        pocket_id: str | None = None,
    ) -> dict[str, Any]:
        ligand = np.asarray(ligand_xyz, dtype=float).reshape((-1, 3))
        if len(ligand) == 0:
            return self._unavailable(pose_id=pose_id, compound_id=compound_id, reason="ligand has no coordinates")
        contacts = receptor_contacts_by_residue(ligand, self.receptor_atoms, self.contact_cutoff_A)
        pocket_center = ligand.mean(axis=0)
        center_method = "ligand_heavy_atom_centroid"
        box_descriptor = None
        if box_center is not None:
            try:
                box_point = np.asarray(list(box_center), dtype=float).reshape(3)
            except (TypeError, ValueError):
                box_point = None
            if box_point is not None:
                box_descriptor = {
                    "box_center_A": [round(float(v), 4) for v in box_point],
                    "box_center_depth_A": (
                        round(float(self.geometry.signed_depth(box_point.reshape(1, 3))[0]), 4)
                        if self.geometry.usable
                        else None
                    ),
                }
        segment_hits = self._segment_hits(contacts)
        residue_rows = self._residue_rows(contacts)
        core_values = [row["core_fraction"] for row in residue_rows if row.get("core_fraction") is not None]
        interface_values = [row["interface_fraction"] for row in residue_rows if row.get("interface_fraction") is not None]
        outside_values = [row["outside_core_fraction"] for row in residue_rows if row.get("outside_core_fraction") is not None]
        core_fraction = float(np.mean(core_values)) if core_values else None
        interface_fraction = float(np.mean(interface_values)) if interface_values else None
        outside_fraction = float(np.mean(outside_values)) if outside_values else None
        contact_xyz = self._contact_coords(contacts)
        contact_centroid_depth = None
        if contact_xyz.size and self.geometry.usable:
            contact_centroid = contact_xyz.mean(axis=0)
            contact_centroid_depth = float(self.geometry.signed_depth(contact_centroid.reshape(1, 3))[0])
        if contact_xyz.size:
            center_distances = np.sqrt(((contact_xyz - pocket_center) ** 2).sum(axis=1))
            extent = float(2.0 * center_distances.max())
            rg = float(np.sqrt(((contact_xyz - contact_xyz.mean(axis=0)) ** 2).sum(axis=1).mean()))
        else:
            extent = 0.0
            rg = 0.0
        depth = self.geometry.signed_depth(pocket_center.reshape(1, 3))[0] if self.geometry.usable else float("nan")
        region = self._classify(core_fraction, interface_fraction, outside_fraction, len(contacts))
        composition = residue_composition(residue_rows)
        status = "evaluated" if self.geometry.usable and contacts else ("needs_input" if not self.geometry.usable else "not_evaluated")
        result = {
            "schema": "ymavs-pocket-analysis/1",
            "algorithm_version": ALGORITHM_VERSION,
            "generated_at": now_iso(),
            "status": status,
            "pose_id": pose_id,
            "compound_id": compound_id,
            "pocket_id": pocket_id or (f"{compound_id}:{pose_id}" if (compound_id or pose_id) else "pocket"),
            "pocket_definition": "receptor residues within contact_cutoff_A of the docked ligand heavy atoms",
            "contact_cutoff_A": self.contact_cutoff_A,
            "pocket_center_A": [round(float(v), 4) for v in pocket_center],
            "pocket_center_method": center_method,
            "pocket_contact_centroid_depth_A": (
                None if contact_centroid_depth is None or not math.isfinite(float(contact_centroid_depth))
                else round(float(contact_centroid_depth), 4)
            ),
            "docking_box_geometry": box_descriptor,
            "pocket_depth_A": None if not math.isfinite(float(depth)) else round(float(depth), 4),
            "pocket_region": region if self.geometry.usable else REGION_UNKNOWN,
            "pocket_core_fraction": None if core_fraction is None else round(core_fraction, 6),
            "pocket_interface_fraction": None if interface_fraction is None else round(interface_fraction, 6),
            "pocket_outside_core_fraction": None if outside_fraction is None else round(outside_fraction, 6),
            "pocket_region_counts": {
                "membrane_core": sum(1 for row in residue_rows if row.get("region") == "membrane_core"),
                "membrane_interface": sum(1 for row in residue_rows if row.get("region") == "membrane_interface"),
                "outside_membrane_core": sum(1 for row in residue_rows if row.get("region") == "outside_membrane_core"),
                "unknown": sum(1 for row in residue_rows if row.get("region") == REGION_UNKNOWN),
            },
            "pocket_region_thresholds": {
                "embedded_core_fraction": self.embedded_core_fraction,
                "interface_fraction": self.interface_fraction,
                "aqueous_outside_fraction": self.aqueous_outside_fraction,
                "mixed_min_fraction": self.mixed_min_fraction,
            },
            "pocket_contact_residues": [row["residue"] for row in residue_rows],
            "pocket_contact_count": len(residue_rows),
            "pocket_atom_contact_count": int(sum(int(row.get("atom_count") or 0) for row in residue_rows)),
            "pocket_residue_count": len({(row["chain"], row["resseq"], row["icode"]) for row in residue_rows}),
            "pocket_extent_A": round(extent, 4),
            "pocket_radius_of_gyration_A": round(rg, 4),
            "pocket_volume_A3": None,
            "pocket_volume_status": "not_evaluated",
            "pocket_volume_note": "no validated cavity-volume method is applied; a bounding box is not reported as a cavity volume",
            "residue_composition": composition,
            "pocket_tm_segments": [name for name, payload in segment_hits.items() if payload["contacted"]],
            "pocket_tm_segment_evidence": segment_hits,
            "tm_annotation_status": (
                "configured_annotation_present" if self.annotated_segments else "no_configured_annotation"
            ),
            "tm_annotation_scope": "configured target annotation residues only; no helix is invented from distances",
            "orientation_status": self.geometry.orientation_status,
            "source": self.geometry.source or None,
            "coordinate_frame": self.geometry.coordinate_frame,
            "evidence_level": "sourced_geometry" if self.geometry.verified else ("geometric_proxy" if self.geometry.usable else "not_evaluated"),
            "limitations": [
                "contact residues are selected by a geometric cutoff, not by energetic decomposition",
                "residue composition is not a measurement of solvent exposure or pocket druggability",
                "pocket volume is not evaluated in this version",
            ],
            "residue_rows": residue_rows,
        }
        return result

    def _unavailable(self, *, pose_id: str, compound_id: str, reason: str) -> dict[str, Any]:
        return {
            "schema": "ymavs-pocket-analysis/1",
            "algorithm_version": ALGORITHM_VERSION,
            "generated_at": now_iso(),
            "status": "not_evaluated",
            "pose_id": pose_id,
            "compound_id": compound_id,
            "pocket_id": f"{compound_id}:{pose_id}" if (compound_id or pose_id) else "pocket",
            "reason": reason,
            "pocket_region": REGION_UNKNOWN,
            "pocket_contact_count": 0,
            "pocket_residue_count": 0,
            "pocket_tm_segments": [],
            "pocket_volume_A3": None,
            "pocket_volume_status": "not_evaluated",
            "orientation_status": self.geometry.orientation_status,
            "coordinate_frame": self.geometry.coordinate_frame,
            "residue_rows": [],
        }

    def _contact_coords(self, contacts: list[dict[str, Any]]) -> np.ndarray:
        values = []
        for row in contacts:
            key = (row["chain"], int(row["resseq"]), row.get("icode") or "", str(row["resname"]).upper())
            coords = self._residue_coords.get(key)
            if coords is not None and len(coords):
                values.append(coords)
        return np.concatenate(values, axis=0) if values else np.zeros((0, 3))

    def _residue_rows(self, contacts: list[dict[str, Any]]) -> list[dict[str, Any]]:
        rows = []
        for contact in contacts:
            key = (contact["chain"], int(contact["resseq"]), contact.get("icode") or "", str(contact["resname"]).upper())
            coords = self._residue_coords.get(key)
            if coords is None or len(coords) == 0:
                continue
            descriptor = self.regions.analyze_points(coords)
            row = dict(contact)
            row.update(
                {
                    "key": key,
                    "region": descriptor["region"] if self.geometry.usable else REGION_UNKNOWN,
                    "signed_depth_mean_A": descriptor["center_depth_A"],
                    "signed_depth_min_A": descriptor["depth_min_A"],
                    "signed_depth_max_A": descriptor["depth_max_A"],
                    "core_fraction": descriptor["core_fraction"],
                    "interface_fraction": descriptor["interface_fraction"],
                    "outside_core_fraction": descriptor["outside_core_fraction"],
                    "orientation_status": self.geometry.orientation_status,
                }
            )
            rows.append(row)
        return rows

    def _segment_hits(self, contacts: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
        contact_keys = {(row["chain"], int(row["resseq"]), row.get("icode") or "") for row in contacts}
        result: dict[str, dict[str, Any]] = {}
        for name, residues in self.annotated_segments.items():
            normalized = {(str(chain), int(resseq), str(icode or "")) for chain, resseq, icode in residues}
            hit_keys = sorted(contact_keys & normalized)
            result[name] = {
                "contacted": bool(hit_keys),
                "annotated_residue_count": len(normalized),
                "contacting_residues": [f"{chain}{resseq}{icode}" for chain, resseq, icode in hit_keys],
                "scope": "configured annotation residue set only",
            }
        return result

    def _classify(
        self,
        core: float | None,
        interface: float | None,
        outside: float | None,
        contact_count: int,
    ) -> str:
        if not self.geometry.usable or contact_count == 0:
            return UNKNOWN
        values = [value for value in (core, interface, outside) if value is not None]
        if not values:
            return UNKNOWN
        if core is not None and core >= self.embedded_core_fraction:
            return "membrane_embedded"
        if outside is not None and outside >= self.aqueous_outside_fraction:
            return "aqueous_side"
        if interface is not None and interface >= self.interface_fraction:
            return "membrane_interface"
        present = [value for value in (core, interface, outside) if value is not None and value >= self.mixed_min_fraction]
        if len(present) >= 2:
            return "mixed"
        return UNKNOWN

    def write_outputs(self, directory: str, analyses: list[dict[str, Any]]) -> dict[str, str]:
        clean = [_strip_private(item) for item in analyses]
        path = os.path.join(directory, "pocket_analysis.json")
        atomic_write_json(
            path,
            {
                "schema": "ymavs-pocket-analysis-collection/1",
                "algorithm_version": ALGORITHM_VERSION,
                "generated_at": now_iso(),
                "count": len(clean),
                "status": "evaluated" if any(item.get("status") == "evaluated" for item in clean) else "needs_input",
                "analyses": clean,
            },
        )
        rows: list[dict[str, Any]] = []
        for item in analyses:
            for row in item.get("residue_rows") or []:
                rows.append(
                    {
                        "pose_id": item.get("pose_id"),
                        "compound_id": item.get("compound_id"),
                        "pocket_region": item.get("pocket_region"),
                        "chain": row.get("chain"),
                        "resseq": row.get("resseq"),
                        "icode": row.get("icode"),
                        "resname": row.get("resname"),
                        "residue": row.get("residue"),
                        "min_distance_A": row.get("min_distance_A"),
                        "atom_count": row.get("atom_count"),
                        "residue_class": row.get("residue_class"),
                        "region": row.get("region"),
                        "signed_depth_mean_A": row.get("signed_depth_mean_A"),
                        "core_fraction": row.get("core_fraction"),
                        "interface_fraction": row.get("interface_fraction"),
                        "outside_core_fraction": row.get("outside_core_fraction"),
                    }
                )
        csv_path = os.path.join(directory, "pocket_residues.csv")
        fields = [
            "pose_id", "compound_id", "pocket_region", "chain", "resseq", "icode", "resname",
            "residue", "min_distance_A", "atom_count", "residue_class", "region",
            "signed_depth_mean_A", "core_fraction", "interface_fraction", "outside_core_fraction",
        ]
        _write_csv(csv_path, fields, rows)
        return {"pocket_analysis": path, "pocket_residues": csv_path}


def _strip_private(item: dict[str, Any]) -> dict[str, Any]:
    payload = dict(item)
    payload.pop("residue_rows", None)
    return payload


def _write_csv(path: str, fields: list[str], rows: list[dict[str, Any]]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
