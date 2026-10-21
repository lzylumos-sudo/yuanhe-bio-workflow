"""Membrane geometry and region classification for Y-MAVS.

This is a coordinate-aware extension of the existing PPM3 orientation transfer.
It never re-derives a membrane plane by PCA when a sourced placement exists, and
it never converts a PPM3 placement into an experimental inside/outside claim.

Definitions written by this module:

``z = (r - center) . normal``
    Signed depth along the membrane normal.  The sign is only a side of the
    chosen normal; it is not cytoplasmic/extracellular without topology evidence.

``membrane_core``
    ``|z| < half_core_thickness``.  This is a hydrophobic-core *depth range*,
    not the lipid-occupied volume and not the full bilayer thickness.

``membrane_interface``
    ``half_core_thickness <= |z| <= half_core_thickness + interface_buffer``.
    The buffer is a configurable geometric partition parameter.

``outside_membrane_core``
    ``|z| > half_core_thickness + interface_buffer``.  It is not automatically
    aqueous: a membrane protein may have a lateral opening at that depth.

``unknown``
    No usable membrane geometry, or a non-finite coordinate.
"""

from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass, field
from typing import Any, Iterable

import numpy as np

from ..util import atomic_write_json, now_iso
from .geom import as_coords

ALGORITHM_VERSION = "ymavs-membrane-regions/1.0.0"
REGION_CORE = "membrane_core"
REGION_INTERFACE = "membrane_interface"
REGION_OUTSIDE = "outside_membrane_core"
REGION_UNKNOWN = "unknown"
REGION_LABELS = (REGION_CORE, REGION_INTERFACE, REGION_OUTSIDE, REGION_UNKNOWN)


def _unit(vector: Iterable[float]) -> np.ndarray:
    value = np.asarray(list(vector), dtype=float).reshape(3)
    norm = float(np.linalg.norm(value))
    if not math.isfinite(norm) or norm <= 1e-12:
        raise ValueError("membrane normal must be finite and non-zero")
    return value / norm


@dataclass
class MembraneGeometry:
    """A verified or proxy membrane plane in a declared coordinate frame."""

    center: np.ndarray
    normal: np.ndarray
    half_thickness_A: float
    source: str
    method: str
    orientation_status: str = "needs_input"
    coordinate_frame: str = ""
    receptor_sha256: str = ""
    uncertainty: str = ""
    evidence_hashes: dict[str, str] = field(default_factory=dict)
    limitations: list[str] = field(default_factory=list)
    algorithm_version: str = ALGORITHM_VERSION

    @classmethod
    def unavailable(cls, *, reason: str = "membrane orientation unavailable") -> "MembraneGeometry":
        return cls(
            center=np.zeros(3),
            normal=np.zeros(3),
            half_thickness_A=0.0,
            source="",
            method="not_computed",
            orientation_status="needs_input",
            coordinate_frame="",
            uncertainty=reason,
            limitations=[reason],
        )

    @classmethod
    def from_plane(
        cls,
        plane: dict[str, Any],
        *,
        orientation_status: str,
        uncertainty: str = "",
    ) -> "MembraneGeometry":
        """Create geometry from an existing orientation-transfer ``plane`` record."""
        center = np.asarray(plane.get("center"), dtype=float).reshape(3)
        normal = _unit(plane.get("normal"))
        thickness = plane.get("thickness_A")
        if thickness is None:
            half = plane.get("membrane_half_thickness_A")
            if half is None:
                raise ValueError("membrane plane has no thickness_A or half thickness")
            thickness = 2.0 * float(half)
        thickness = float(thickness)
        if not math.isfinite(thickness) or thickness <= 0:
            raise ValueError("membrane thickness must be finite and positive")
        if orientation_status not in {"verified_computational", "proxy_only"}:
            raise ValueError("unsupported orientation status: " + str(orientation_status))
        limitations = list(plane.get("limitations") or [])
        limitations.extend(
            [
                "signed depth is relative to the supplied normal; no cytoplasmic/extracellular assignment is made",
                "hydrophobic-core depth is not a lipid-occupancy map or a full bilayer extent",
            ]
        )
        return cls(
            center=center,
            normal=normal,
            half_thickness_A=thickness / 2.0,
            source=str(plane.get("source") or ""),
            method=str(plane.get("method") or ""),
            orientation_status=orientation_status,
            coordinate_frame=str(plane.get("coordinate_frame") or ""),
            receptor_sha256=str(plane.get("receptor_sha256") or ""),
            uncertainty=uncertainty,
            evidence_hashes=dict(plane.get("evidence_hashes") or {}),
            limitations=sorted(set(limitations)),
        )

    @property
    def usable(self) -> bool:
        return bool(
            self.orientation_status in {"verified_computational", "proxy_only"}
            and self.half_thickness_A > 0
            and np.isfinite(self.center).all()
            and np.isfinite(self.normal).all()
            and abs(float(np.linalg.norm(self.normal)) - 1.0) <= 1e-6
        )

    @property
    def thickness_A(self) -> float | None:
        return 2.0 * self.half_thickness_A if self.usable else None

    @property
    def verified(self) -> bool:
        return self.orientation_status == "verified_computational"

    def signed_depth(self, points: np.ndarray) -> np.ndarray:
        values = np.asarray(points, dtype=float).reshape((-1, 3))
        if not self.usable or len(values) == 0:
            return np.full(len(values), np.nan, dtype=float)
        return (values - self.center) @ self.normal

    def region_for_depth(self, depth: float, buffer_A: float) -> str:
        if not math.isfinite(float(depth)):
            return REGION_UNKNOWN
        absolute = abs(float(depth))
        half = float(self.half_thickness_A)
        buffer = max(0.0, float(buffer_A))
        if absolute < half:
            return REGION_CORE
        if absolute <= half + buffer:
            return REGION_INTERFACE
        return REGION_OUTSIDE

    def side_for_depth(self, depth: float) -> str:
        if not math.isfinite(float(depth)):
            return "unknown"
        if abs(float(depth)) < 1e-9:
            return "at_midplane"
        # A core atom lies on the midplane's positive or negative side, but that
        # side is not a membrane side.  The distinction is kept explicit.
        return "positive_normal_side" if float(depth) > 0 else "negative_normal_side"

    def point_descriptor(self, point: Iterable[float], buffer_A: float) -> dict[str, Any]:
        xyz = np.asarray(list(point), dtype=float).reshape(3)
        depth = float(self.signed_depth(xyz.reshape(1, 3))[0]) if self.usable else float("nan")
        return {
            "region": self.region_for_depth(depth, buffer_A) if self.usable else REGION_UNKNOWN,
            "signed_depth_A": None if not math.isfinite(depth) else round(depth, 4),
            "absolute_depth_A": None if not math.isfinite(depth) else round(abs(depth), 4),
            "side": self.side_for_depth(depth) if self.usable else "unknown",
            "orientation_status": self.orientation_status,
            "source": self.source,
            "coordinate_frame": self.coordinate_frame,
        }

    def to_dict(self, *, buffer_A: float) -> dict[str, Any]:
        status = "evaluated" if self.verified else ("proxy_only" if self.usable else "needs_input")
        return {
            "schema": "ymavs-membrane-geometry/1",
            "algorithm_version": self.algorithm_version,
            "generated_at": now_iso(),
            "status": status,
            "orientation_status": self.orientation_status,
            "usable": self.usable,
            "verified_for_formal_use": self.verified,
            "coordinate_frame": self.coordinate_frame,
            "receptor_sha256": self.receptor_sha256,
            "source": self.source,
            "method": self.method,
            "uncertainty": self.uncertainty,
            "geometry": {
                "center_A": None if not self.usable else [round(float(v), 6) for v in self.center],
                "normal_unit": None if not self.usable else [round(float(v), 6) for v in self.normal],
                "half_core_thickness_A": round(float(self.half_thickness_A), 6) if self.usable else None,
                "full_thickness_A": round(float(self.thickness_A), 6) if self.thickness_A else None,
                "interface_buffer_A": round(float(buffer_A), 6),
            },
            "region_definition": {
                "membrane_core": "abs(signed_depth_A) < half_core_thickness_A",
                "membrane_interface": "half_core_thickness_A <= abs(signed_depth_A) <= half_core_thickness_A + interface_buffer_A",
                "outside_membrane_core": "abs(signed_depth_A) > half_core_thickness_A + interface_buffer_A",
                "unknown": "no usable membrane geometry or non-finite coordinate",
                "side_warning": "positive/negative refers only to the supplied normal; no cytoplasmic/extracellular assignment is made",
            },
            "topology_assignment": "not_attempted",
            "evidence_hashes": dict(self.evidence_hashes),
            "limitations": list(self.limitations),
        }


class MembraneRegionAnalyzer:
    """Classify atoms, residues, pockets and ligands against one membrane plane."""

    def __init__(
        self,
        geometry: MembraneGeometry,
        *,
        interface_buffer_A: float = 5.0,
    ):
        if not math.isfinite(float(interface_buffer_A)) or float(interface_buffer_A) < 0:
            raise ValueError("interface_buffer_A must be finite and non-negative")
        self.geometry = geometry
        self.interface_buffer_A = float(interface_buffer_A)

    def analyze_points(self, points: np.ndarray) -> dict[str, Any]:
        values = np.asarray(points, dtype=float).reshape((-1, 3))
        if not self.geometry.usable or len(values) == 0:
            counts = {label: 0 for label in REGION_LABELS}
            counts[REGION_UNKNOWN] = int(len(values))
            return {
                "count": int(len(values)),
                "region": REGION_UNKNOWN,
                "region_counts": counts,
                "core_fraction": None,
                "interface_fraction": None,
                "outside_core_fraction": None,
                "center_depth_A": None,
                "depth_min_A": None,
                "depth_max_A": None,
                "depth_median_A": None,
                "orientation_status": self.geometry.orientation_status,
                "coordinate_frame": self.geometry.coordinate_frame,
                "method": "signed_plane_depth",
                "source": self.geometry.source or None,
                "evidence_level": "not_evaluated",
            }
        depths = self.geometry.signed_depth(values)
        finite = np.isfinite(depths)
        if not finite.any():
            counts = {label: 0 for label in REGION_LABELS}
            counts[REGION_UNKNOWN] = int(len(values))
            return {
                "count": int(len(values)),
                "region": REGION_UNKNOWN,
                "region_counts": counts,
                "core_fraction": None,
                "interface_fraction": None,
                "outside_core_fraction": None,
                "center_depth_A": None,
                "depth_min_A": None,
                "depth_max_A": None,
                "depth_median_A": None,
                "orientation_status": self.geometry.orientation_status,
                "coordinate_frame": self.geometry.coordinate_frame,
                "method": "signed_plane_depth",
                "source": self.geometry.source or None,
                "evidence_level": "not_evaluated",
            }
        valid = depths[finite]
        half = float(self.geometry.half_thickness_A)
        buffer = self.interface_buffer_A
        absolute = np.abs(valid)
        core = absolute < half
        interface = (absolute >= half) & (absolute <= half + buffer)
        outside = absolute > half + buffer
        count = int(len(valid))
        counts = {
            REGION_CORE: int(core.sum()),
            REGION_INTERFACE: int(interface.sum()),
            REGION_OUTSIDE: int(outside.sum()),
            REGION_UNKNOWN: int((~finite).sum()),
        }
        centroid_depth = None
        if len(values[finite]):
            centroid_depth = float(self.geometry.signed_depth(values[finite].mean(axis=0).reshape(1, 3))[0])
        return {
            "count": count,
            "region_counts": counts,
            "core_fraction": round(float(core.mean()), 6),
            "interface_fraction": round(float(interface.mean()), 6),
            "outside_core_fraction": round(float(outside.mean()), 6),
            "center_depth_A": None if centroid_depth is None else round(centroid_depth, 4),
            "depth_min_A": round(float(valid.min()), 4),
            "depth_max_A": round(float(valid.max()), 4),
            "depth_median_A": round(float(np.median(valid)), 4),
            "region": self._dominant_region(counts, count),
            "orientation_status": self.geometry.orientation_status,
            "coordinate_frame": self.geometry.coordinate_frame,
            "method": "signed_plane_depth",
            "source": self.geometry.source or None,
            "evidence_level": "sourced_geometry" if self.geometry.verified else "geometric_proxy",
        }

    def _dominant_region(self, counts: dict[str, int], total: int) -> str:
        if total <= 0:
            return REGION_UNKNOWN
        core = counts.get(REGION_CORE, 0) / total
        interface = counts.get(REGION_INTERFACE, 0) / total
        outside = counts.get(REGION_OUTSIDE, 0) / total
        largest = max((core, REGION_CORE), (interface, REGION_INTERFACE), (outside, REGION_OUTSIDE))
        if largest[0] >= 0.5:
            return largest[1]
        return "mixed"

    def analyze_residues(self, atoms: Iterable[Any]) -> list[dict[str, Any]]:
        grouped: dict[tuple[str, int, str, str], list[Any]] = {}
        for atom in atoms:
            if hasattr(atom, "is_hydrogen") and atom.is_hydrogen():
                continue
            key = (
                str(getattr(atom, "chain", "") or ""),
                int(getattr(atom, "resseq", 0) or 0),
                str(getattr(atom, "icode", "") or ""),
                str(getattr(atom, "resname", "") or "").upper(),
            )
            grouped.setdefault(key, []).append(atom)
        rows: list[dict[str, Any]] = []
        for key in sorted(grouped):
            atoms_here = grouped[key]
            coords = as_coords(atoms_here)
            descriptor = self.analyze_points(coords)
            row = {
                "chain": key[0],
                "resseq": key[1],
                "icode": key[2],
                "resname": key[3],
                "residue": f"{key[0]}{key[1]}{key[2]}:{key[3]}" if key[0] else f"{key[1]}{key[2]}:{key[3]}",
                "atom_count": len(atoms_here),
                "signed_depth_mean_A": descriptor["center_depth_A"],
                "signed_depth_min_A": descriptor["depth_min_A"],
                "signed_depth_max_A": descriptor["depth_max_A"],
                "region": descriptor["region"] if self.geometry.usable else REGION_UNKNOWN,
                "core_fraction": descriptor["core_fraction"],
                "interface_fraction": descriptor["interface_fraction"],
                "outside_core_fraction": descriptor["outside_core_fraction"],
                "side": self.geometry.side_for_depth(descriptor["center_depth_A"]) if descriptor["center_depth_A"] is not None else "unknown",
                "orientation_status": self.geometry.orientation_status,
                "source": self.geometry.source or None,
            }
            rows.append(row)
        return rows

    def analyze_atoms(self, atoms: Iterable[Any]) -> list[dict[str, Any]]:
        rows = []
        for index, atom in enumerate(atoms):
            if hasattr(atom, "is_hydrogen") and atom.is_hydrogen():
                continue
            point = np.asarray([[float(atom.x), float(atom.y), float(atom.z)]], dtype=float)
            depth = float(self.geometry.signed_depth(point)[0]) if self.geometry.usable else float("nan")
            rows.append(
                {
                    "atom_index": index,
                    "chain": str(getattr(atom, "chain", "") or ""),
                    "resseq": int(getattr(atom, "resseq", 0) or 0),
                    "icode": str(getattr(atom, "icode", "") or ""),
                    "resname": str(getattr(atom, "resname", "") or "").upper(),
                    "atom_name": str(getattr(atom, "name", "") or ""),
                    "element": str(getattr(atom, "element", "") or ""),
                    "signed_depth_A": None if not math.isfinite(depth) else round(depth, 4),
                    "absolute_depth_A": None if not math.isfinite(depth) else round(abs(depth), 4),
                    "region": self.geometry.region_for_depth(depth, self.interface_buffer_A) if self.geometry.usable else REGION_UNKNOWN,
                }
            )
        return rows

    def analyze_ligand(
        self,
        coordinates: np.ndarray,
        *,
        pose_id: str = "",
        compound_id: str = "",
        elements: Iterable[str] | None = None,
    ) -> dict[str, Any]:
        coords = np.asarray(coordinates, dtype=float).reshape((-1, 3))
        descriptor = self.analyze_points(coords)
        depth_values = self.geometry.signed_depth(coords) if self.geometry.usable else np.full(len(coords), np.nan)
        elements_list = list(elements or [])
        atom_rows = []
        for index, point in enumerate(coords):
            depth = float(depth_values[index]) if self.geometry.usable else float("nan")
            atom_rows.append(
                {
                    "pose_id": pose_id,
                    "compound_id": compound_id,
                    "atom_index": index,
                    "element": elements_list[index] if index < len(elements_list) else "",
                    "x": round(float(point[0]), 4),
                    "y": round(float(point[1]), 4),
                    "z": round(float(point[2]), 4),
                    "signed_depth_A": None if not math.isfinite(depth) else round(depth, 4),
                    "absolute_depth_A": None if not math.isfinite(depth) else round(abs(depth), 4),
                    "region": self.geometry.region_for_depth(depth, self.interface_buffer_A) if self.geometry.usable else REGION_UNKNOWN,
                    "side": self.geometry.side_for_depth(depth) if self.geometry.usable else "unknown",
                }
            )
        return {
            "pose_id": pose_id,
            "compound_id": compound_id,
            "atom_count": len(coords),
            "ligand_center_depth_A": descriptor["center_depth_A"],
            "ligand_depth_min_A": descriptor["depth_min_A"],
            "ligand_depth_max_A": descriptor["depth_max_A"],
            "ligand_depth_median_A": descriptor["depth_median_A"],
            "ligand_core_fraction": descriptor["core_fraction"],
            "ligand_interface_fraction": descriptor["interface_fraction"],
            "ligand_outside_core_fraction": descriptor["outside_core_fraction"],
            "ligand_region": descriptor["region"],
            "dominant_region": descriptor["region"],
            "region_counts": descriptor["region_counts"],
            "orientation_status": self.geometry.orientation_status,
            "coordinate_frame": self.geometry.coordinate_frame,
            "source": self.geometry.source or None,
            "method": "signed_plane_depth",
            "evidence_level": descriptor["evidence_level"],
            "atom_rows": atom_rows,
        }

    def analyze_pocket(self, pocket_coordinates: np.ndarray, *, pocket_id: str = "") -> dict[str, Any]:
        descriptor = self.analyze_points(pocket_coordinates)
        return {
            "pocket_id": pocket_id,
            "pocket_depth_A": descriptor["center_depth_A"],
            "signed_depth_min_A": descriptor["depth_min_A"],
            "signed_depth_max_A": descriptor["depth_max_A"],
            "pocket_core_fraction": descriptor["core_fraction"],
            "pocket_interface_fraction": descriptor["interface_fraction"],
            "pocket_outside_core_fraction": descriptor["outside_core_fraction"],
            "pocket_region": descriptor["region"],
            "region_counts": descriptor["region_counts"],
            "orientation_status": self.geometry.orientation_status,
            "coordinate_frame": self.geometry.coordinate_frame,
            "source": self.geometry.source or None,
            "method": "signed_plane_depth",
            "evidence_level": descriptor["evidence_level"],
        }

    def write_geometry(
        self,
        directory: str,
        *,
        receptor_path: str = "",
        input_hashes: dict[str, str] | None = None,
        atom_rows: list[dict[str, Any]] | None = None,
    ) -> str:
        payload = self.geometry.to_dict(buffer_A=self.interface_buffer_A)
        payload["algorithm_version"] = ALGORITHM_VERSION
        payload["input_hashes"] = dict(input_hashes or {})
        payload["receptor_path"] = receptor_path or None
        if atom_rows:
            counts = {label: 0 for label in REGION_LABELS}
            for row in atom_rows:
                counts[str(row.get("region") or REGION_UNKNOWN)] = counts.get(str(row.get("region") or REGION_UNKNOWN), 0) + 1
            payload["receptor_atom_region_counts"] = counts
        path = os.path.join(directory, "membrane_geometry.json")
        return atomic_write_json(path, payload)

    def write_residue_csv(self, directory: str, rows: list[dict[str, Any]]) -> str:
        path = os.path.join(directory, "residue_membrane_regions.csv")
        fields = [
            "chain", "resseq", "icode", "resname", "residue", "atom_count",
            "signed_depth_mean_A", "signed_depth_min_A", "signed_depth_max_A",
            "region", "core_fraction", "interface_fraction", "outside_core_fraction",
            "side", "orientation_status", "source",
        ]
        _write_csv(path, fields, rows)
        return path

    def write_ligand_csv(self, directory: str, ligand_rows: list[dict[str, Any]]) -> str:
        path = os.path.join(directory, "ligand_membrane_regions.csv")
        fields = [
            "pose_id", "compound_id", "atom_index", "element", "x", "y", "z",
            "signed_depth_A", "absolute_depth_A", "region", "side",
        ]
        _write_csv(path, fields, ligand_rows)
        return path


def _write_csv(path: str, fields: list[str], rows: list[dict[str, Any]]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
