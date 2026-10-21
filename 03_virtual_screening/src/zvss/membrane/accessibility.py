"""Static geometric accessibility for a membrane-protein binding site.

Method: build a coarse occupancy grid from heavy receptor atoms expanded by a
probe radius, seed free cells around the docked ligand, and test six-connected
free-space paths to (a) the aqueous side beyond the membrane interface and (b)
the lateral membrane boundary.  This is deliberately a *geometric proxy*.

It does not establish that a drug molecule can permeate a membrane, that a lipid
path exists in a dynamic bilayer, or that a static closed conformation is
unbindable.  A missing path is therefore reported as ``uncertain``, never as a
biological exclusion.
"""

from __future__ import annotations

import csv
import math
import os
from collections import deque
from typing import Any, Iterable

import numpy as np

from ..util import atomic_write_json, now_iso
from .geom import as_coords, element_radius
from .regions import MembraneGeometry

ALGORITHM_VERSION = "ymavs-accessibility/1.0.0"
ACCESSIBLE = "accessible"
UNCERTAIN = "uncertain"
NOT_EVALUATED = "not_evaluated"


class MembraneAccessibilityAnalyzer:
    def __init__(
        self,
        geometry: MembraneGeometry,
        *,
        grid_spacing_A: float = 1.0,
        probe_radius_A: float = 1.4,
        margin_A: float = 8.0,
        seed_radius_A: float = 3.0,
        max_grid_points: int = 1_500_000,
        interface_buffer_A: float = 5.0,
    ):
        self.geometry = geometry
        self.grid_spacing_A = float(grid_spacing_A)
        self.probe_radius_A = float(probe_radius_A)
        self.margin_A = float(margin_A)
        self.seed_radius_A = float(seed_radius_A)
        self.max_grid_points = int(max_grid_points)
        self.interface_buffer_A = float(interface_buffer_A)
        if self.grid_spacing_A <= 0 or not math.isfinite(self.grid_spacing_A):
            raise ValueError("grid_spacing_A must be finite and positive")
        if self.probe_radius_A < 0 or not math.isfinite(self.probe_radius_A):
            raise ValueError("probe_radius_A must be finite and non-negative")
        self._prepared = False
        self._prepare_error: str | None = None
        self._free: np.ndarray | None = None
        self._shape: tuple[int, int, int] = (0, 0, 0)
        self._origin = np.zeros(3)
        self._basis: np.ndarray | None = None
        self._z_values = np.zeros(0)
        self._receptor_atom_count = 0
        self._grid_points = 0
        self._occupied_fraction: float | None = None

    # ------------------------------------------------------------------ #
    # coordinate frame and occupancy grid
    # ------------------------------------------------------------------ #
    def _frame(self) -> tuple[np.ndarray, np.ndarray]:
        """Return (center, orthonormal basis with the membrane normal third)."""
        normal = np.asarray(self.geometry.normal, dtype=float)
        axis = np.array([1.0, 0.0, 0.0])
        if abs(float(axis @ normal)) > 0.9:
            axis = np.array([0.0, 1.0, 0.0])
        u = np.cross(normal, axis)
        norm = float(np.linalg.norm(u))
        if norm < 1e-9:
            raise ValueError("could not construct a membrane-plane basis")
        u = u / norm
        v = np.cross(normal, u)
        v = v / float(np.linalg.norm(v))
        basis = np.column_stack([u, v, normal])
        return np.asarray(self.geometry.center, dtype=float), basis

    def _to_local(self, xyz: np.ndarray) -> np.ndarray:
        assert self._basis is not None
        return (np.asarray(xyz, dtype=float) - self.geometry.center) @ self._basis

    def _to_world(self, local: np.ndarray) -> np.ndarray:
        assert self._basis is not None
        return self.geometry.center + np.asarray(local, dtype=float) @ self._basis.T

    def prepare(self, receptor_atoms: Iterable[Any]) -> dict[str, Any]:
        atoms = [a for a in receptor_atoms if not (hasattr(a, "is_hydrogen") and a.is_hydrogen())]
        self._receptor_atom_count = len(atoms)
        if not self.geometry.usable:
            self._prepare_error = "membrane_geometry_needs_input"
            return self._prepare_summary()
        coords = as_coords(atoms)
        if len(coords) == 0:
            self._prepare_error = "receptor_has_no_heavy_atoms"
            return self._prepare_summary()
        try:
            self.geometry.center, self._basis = self._frame()
        except ValueError as exc:
            self._prepare_error = str(exc)
            return self._prepare_summary()
        local = self._to_local(coords)
        half = float(self.geometry.half_thickness_A)
        required = half + self.interface_buffer_A + max(2.0, self.margin_A * 0.5)
        z_min = min(float(local[:, 2].min()) - self.margin_A, -required)
        z_max = max(float(local[:, 2].max()) + self.margin_A, required)
        x_min, y_min = local[:, 0].min() - self.margin_A, local[:, 1].min() - self.margin_A
        x_max, y_max = local[:, 0].max() + self.margin_A, local[:, 1].max() + self.margin_A
        spacing = self.grid_spacing_A
        nx = int(math.ceil((x_max - x_min) / spacing)) + 1
        ny = int(math.ceil((y_max - y_min) / spacing)) + 1
        nz = int(math.ceil((z_max - z_min) / spacing)) + 1
        points = nx * ny * nz
        self._shape = (nx, ny, nz)
        self._grid_points = points
        if points > self.max_grid_points:
            self._prepare_error = (
                f"grid has {points} points, above max_grid_points={self.max_grid_points}; "
                "coarsen grid_spacing_A or reduce the receptor region explicitly"
            )
            return self._prepare_summary()
        self._origin = np.array([x_min, y_min, z_min], dtype=float)
        self._z_values = self._origin[2] + np.arange(nz) * spacing
        free = np.ones((nx, ny, nz), dtype=bool)
        # Rasterise the expanded van der Waals volume into the local grid.  Only
        # each atom's local neighbourhood is visited, not the full grid.
        for atom, point in zip(atoms, local):
            radius = element_radius(getattr(atom, "element", "")) + self.probe_radius_A
            low = np.floor((point - radius - self._origin) / spacing).astype(int)
            high = np.ceil((point + radius - self._origin) / spacing).astype(int)
            low = np.maximum(low, 0)
            high = np.minimum(high, np.asarray(self._shape) - 1)
            if np.any(high < low):
                continue
            i0, j0, k0 = (int(v) for v in low)
            i1, j1, k1 = (int(v) for v in high)
            xs = self._origin[0] + np.arange(i0, i1 + 1) * spacing
            ys = self._origin[1] + np.arange(j0, j1 + 1) * spacing
            zs = self._origin[2] + np.arange(k0, k1 + 1) * spacing
            dx = (xs - point[0]) ** 2
            dy = (ys - point[1]) ** 2
            dz = (zs - point[2]) ** 2
            mask = dx[:, None, None] + dy[None, :, None] + dz[None, None, :] <= radius * radius
            free[i0 : i1 + 1, j0 : j1 + 1, k0 : k1 + 1] &= ~mask
        self._free = free
        self._occupied_fraction = round(float(1.0 - free.mean()), 6)
        self._prepared = True
        self._prepare_error = None
        return self._prepare_summary()

    def _prepare_summary(self) -> dict[str, Any]:
        return {
            "prepared": self._prepared,
            "status": "ready" if self._prepared else "not_evaluated",
            "reason": self._prepare_error,
            "receptor_atoms": self._receptor_atom_count,
            "grid_shape": list(self._shape),
            "grid_points": self._grid_points,
            "grid_spacing_A": self.grid_spacing_A,
            "probe_radius_A": self.probe_radius_A,
            "occupied_fraction": self._occupied_fraction,
            "max_grid_points": self.max_grid_points,
            "coordinate_frame": self.geometry.coordinate_frame,
            "orientation_status": self.geometry.orientation_status,
        }

    # ------------------------------------------------------------------ #
    # per-pose path analysis
    # ------------------------------------------------------------------ #
    def analyze(
        self,
        ligand_xyz: np.ndarray,
        *,
        pose_id: str = "",
        compound_id: str = "",
        pocket_center: Iterable[float] | None = None,
    ) -> dict[str, Any]:
        ligand = np.asarray(ligand_xyz, dtype=float).reshape((-1, 3))
        if not self._prepared or self._free is None:
            return self._not_evaluated(
                pose_id,
                compound_id,
                reason=self._prepare_error or "accessibility grid was not prepared",
                ligand=ligand,
            )
        if len(ligand) == 0:
            return self._not_evaluated(pose_id, compound_id, reason="ligand has no coordinates")
        assert self._basis is not None
        local = self._to_local(ligand)
        origin = self._origin
        spacing = self.grid_spacing_A
        shape = self._shape
        if (
            local[:, 0].min() < origin[0] - spacing
            or local[:, 0].max() > origin[0] + (shape[0] - 1) * spacing + spacing
            or local[:, 1].min() < origin[1] - spacing
            or local[:, 1].max() > origin[1] + (shape[1] - 1) * spacing + spacing
            or local[:, 2].min() < origin[2] - spacing
            or local[:, 2].max() > origin[2] + (shape[2] - 1) * spacing + spacing
        ):
            return self._not_evaluated(
                pose_id,
                compound_id,
                reason="ligand lies outside the prepared receptor grid",
                ligand=ligand,
            )
        seeds = self._seed_indices(local)
        if not seeds:
            return self._uncertain(
                pose_id,
                compound_id,
                reason="no free grid cell was found within seed_radius_A of the ligand",
                seeds=0,
            )
        search = self._bfs(seeds, local)
        aqueous = search.get("aqueous")
        lateral = search.get("membrane_lateral")
        aqueous_status = ACCESSIBLE if aqueous else UNCERTAIN
        lateral_status = ACCESSIBLE if lateral else UNCERTAIN
        paths = {}
        if aqueous:
            paths["aqueous"] = aqueous
        if lateral:
            paths["membrane_lateral"] = lateral
        lengths = [entry["length_A"] for entry in paths.values()]
        result = {
            "schema": "ymavs-accessibility-analysis/1",
            "algorithm_version": ALGORITHM_VERSION,
            "generated_at": now_iso(),
            "pose_id": pose_id,
            "compound_id": compound_id,
            "status": "evaluated",
            "aqueous_accessibility": aqueous_status,
            "membrane_lateral_accessibility": lateral_status,
            "geometric_path_length_A": round(min(lengths), 3) if lengths else None,
            "aqueous_path_length_A": paths.get("aqueous", {}).get("length_A"),
            "membrane_lateral_path_length_A": paths.get("membrane_lateral", {}).get("length_A"),
            "probe_radius_A": self.probe_radius_A,
            "grid_spacing_A": self.grid_spacing_A,
            "method": "grid_connectivity",
            "evidence_level": "geometric_proxy",
            "geometric_accessibility_only": True,
            "static_connectivity": {
                "aqueous": bool(aqueous),
                "membrane_lateral": bool(lateral),
            },
            "seed_cell_count": len(seeds),
            "visited_free_cells": int(search["visited"]),
            "ligand_local_bounds": {
                "x": [round(float(local[:, 0].min()), 3), round(float(local[:, 0].max()), 3)],
                "y": [round(float(local[:, 1].min()), 3), round(float(local[:, 1].max()), 3)],
                "z": [round(float(local[:, 2].min()), 3), round(float(local[:, 2].max()), 3)],
            },
            "orientation_status": self.geometry.orientation_status,
            "coordinate_frame": self.geometry.coordinate_frame,
            "limitations": [
                "a 1.4 A probe describes a small solvent-sized probe, not a whole drug molecule",
                "the static crystal/docking structure may close a channel that opens dynamically",
                "lipid-side entry requires lipid rearrangement and is not proven by an empty grid path",
                "a missing static path is reported as uncertain, not as proof of no binding",
                "grid_connectivity is a geometric proxy, not a free-energy or kinetic calculation",
            ],
            "paths": paths,
        }
        if aqueous is None and lateral is None:
            result["status"] = "uncertain"
            result["reason"] = (
                "no connected free-space path reached either the aqueous or lateral membrane boundary "
                "of this static grid"
            )
        return result

    def _seed_indices(self, ligand_local: np.ndarray) -> list[int]:
        assert self._free is not None
        nx, ny, nz = self._shape
        spacing = self.grid_spacing_A
        indices: set[int] = set()
        radius = self.seed_radius_A
        for point in ligand_local:
            low = np.floor((point - radius - self._origin) / spacing).astype(int)
            high = np.ceil((point + radius - self._origin) / spacing).astype(int)
            low = np.maximum(low, 0)
            high = np.minimum(high, np.asarray(self._shape) - 1)
            if np.any(high < low):
                continue
            i0, j0, k0 = (int(v) for v in low)
            i1, j1, k1 = (int(v) for v in high)
            xs = self._origin[0] + np.arange(i0, i1 + 1) * spacing
            ys = self._origin[1] + np.arange(j0, j1 + 1) * spacing
            zs = self._origin[2] + np.arange(k0, k1 + 1) * spacing
            dx = (xs - point[0]) ** 2
            dy = (ys - point[1]) ** 2
            dz = (zs - point[2]) ** 2
            near = dx[:, None, None] + dy[None, :, None] + dz[None, None, :] <= radius * radius
            block = self._free[i0 : i1 + 1, j0 : j1 + 1, k0 : k1 + 1] & near
            for i, j, k in np.argwhere(block):
                indices.add((int(i + i0) * ny + int(j + j0)) * nz + int(k + k0))
        return sorted(indices)

    def _bfs(self, seeds: list[int], ligand_local: np.ndarray) -> dict[str, Any]:
        assert self._free is not None
        nx, ny, nz = self._shape
        total = nx * ny * nz
        plane = ny * nz
        visited = bytearray(total)
        previous = np.full(total, -1, dtype=np.int64)
        queue: deque[int] = deque()
        for seed in seeds:
            if self._free.flat[seed] and not visited[seed]:
                visited[seed] = 1
                previous[seed] = seed
                queue.append(seed)
        aqueous = None
        lateral = None
        half = float(self.geometry.half_thickness_A)
        interface = half + self.interface_buffer_A
        z_values = self._z_values
        visited_count = 0
        while queue:
            current = queue.popleft()
            visited_count += 1
            i = current // plane
            remainder = current - i * plane
            j = remainder // nz
            k = remainder - j * nz
            boundary = i in (0, nx - 1) or j in (0, ny - 1) or k in (0, nz - 1)
            if boundary:
                abs_depth = abs(float(z_values[k]))
                if aqueous is None and abs_depth > interface:
                    aqueous = current
                lateral_boundary = i in (0, nx - 1) or j in (0, ny - 1)
                if lateral is None and lateral_boundary and abs_depth <= interface:
                    lateral = current
                if aqueous is not None and lateral is not None:
                    break
            neighbours = []
            if i > 0:
                neighbours.append(current - plane)
            if i + 1 < nx:
                neighbours.append(current + plane)
            if j > 0:
                neighbours.append(current - nz)
            if j + 1 < ny:
                neighbours.append(current + nz)
            if k > 0:
                neighbours.append(current - 1)
            if k + 1 < nz:
                neighbours.append(current + 1)
            for neighbour in neighbours:
                if not visited[neighbour] and self._free.flat[neighbour]:
                    visited[neighbour] = 1
                    previous[neighbour] = current
                    queue.append(neighbour)
        return {
            "aqueous": None if aqueous is None else self._trace(aqueous, previous),
            "membrane_lateral": None if lateral is None else self._trace(lateral, previous),
            "visited": visited_count,
        }

    def _trace(self, end: int, previous: np.ndarray) -> dict[str, Any]:
        nx, ny, nz = self._shape
        plane = ny * nz
        path: list[int] = []
        current = int(end)
        guard = 0
        while guard < previous.size + 1:
            path.append(current)
            parent = int(previous[current])
            if parent == current or parent < 0:
                break
            current = parent
            guard += 1
        path.reverse()
        local = np.array(
            [
                [
                    self._origin[0] + (index // plane) * self.grid_spacing_A,
                    self._origin[1] + ((index - (index // plane) * plane) // nz) * self.grid_spacing_A,
                    self._origin[2] + ((index - (index // plane) * plane) - ((index - (index // plane) * plane) // nz) * nz) * self.grid_spacing_A,
                ]
                for index in path
            ],
            dtype=float,
        )
        world = self._to_world(local)
        length = float((len(path) - 1) * self.grid_spacing_A)
        return {
            "length_A": round(length, 3),
            "grid_point_count": len(path),
            "endpoint": [round(float(v), 3) for v in world[-1]],
            "path_xyz_A": [[round(float(v), 3) for v in point] for point in world],
        }

    def _not_evaluated(self, pose_id: str, compound_id: str, *, reason: str, ligand: np.ndarray | None = None) -> dict[str, Any]:
        return {
            "schema": "ymavs-accessibility-analysis/1",
            "algorithm_version": ALGORITHM_VERSION,
            "generated_at": now_iso(),
            "pose_id": pose_id,
            "compound_id": compound_id,
            "status": NOT_EVALUATED,
            "aqueous_accessibility": NOT_EVALUATED,
            "membrane_lateral_accessibility": NOT_EVALUATED,
            "geometric_path_length_A": None,
            "method": "grid_connectivity",
            "evidence_level": "not_evaluated",
            "geometric_accessibility_only": True,
            "reason": reason,
            "probe_radius_A": self.probe_radius_A,
            "grid_spacing_A": self.grid_spacing_A,
            "orientation_status": self.geometry.orientation_status,
            "coordinate_frame": self.geometry.coordinate_frame,
            "paths": {},
        }

    def _uncertain(self, pose_id: str, compound_id: str, *, reason: str, seeds: int = 0) -> dict[str, Any]:
        result = self._not_evaluated(pose_id, compound_id, reason=reason)
        result.update(
            {
                "status": "uncertain",
                "aqueous_accessibility": UNCERTAIN,
                "membrane_lateral_accessibility": UNCERTAIN,
                "evidence_level": "geometric_proxy",
                "seed_cell_count": seeds,
            }
        )
        return result

    # ------------------------------------------------------------------ #
    # outputs
    # ------------------------------------------------------------------ #
    def write_outputs(self, directory: str, analyses: list[dict[str, Any]]) -> dict[str, str]:
        summaries = []
        path_entries = []
        for item in analyses:
            summary = {key: value for key, value in item.items() if key != "paths"}
            summaries.append(summary)
            path_entry = {
                "pose_id": item.get("pose_id"),
                "compound_id": item.get("compound_id"),
                "paths": item.get("paths") or {},
            }
            path_entries.append(path_entry)
        analysis_path = os.path.join(directory, "accessibility_analysis.json")
        atomic_write_json(
            analysis_path,
            {
                "schema": "ymavs-accessibility-collection/1",
                "algorithm_version": ALGORITHM_VERSION,
                "generated_at": now_iso(),
                "method": "grid_connectivity",
                "evidence_level": "geometric_proxy",
                "grid": self._prepare_summary(),
                "count": len(summaries),
                "analyses": summaries,
            },
        )
        paths_path = os.path.join(directory, "accessibility_paths.json")
        atomic_write_json(
            paths_path,
            {
                "schema": "ymavs-accessibility-paths/1",
                "algorithm_version": ALGORITHM_VERSION,
                "generated_at": now_iso(),
                "coordinate_frame": self.geometry.coordinate_frame,
                "note": "Coordinates are in the docking receptor frame. They are visualisation geometry, not a permeation mechanism.",
                "paths": path_entries,
            },
        )
        csv_path = os.path.join(directory, "accessibility_analysis.csv")
        fields = [
            "pose_id", "compound_id", "status", "aqueous_accessibility", "membrane_lateral_accessibility",
            "geometric_path_length_A", "aqueous_path_length_A", "membrane_lateral_path_length_A",
            "probe_radius_A", "grid_spacing_A", "reason",
        ]
        _write_csv(csv_path, fields, summaries)
        return {"accessibility_analysis": analysis_path, "accessibility_paths": paths_path, "accessibility_csv": csv_path}


def _write_csv(path: str, fields: list[str], rows: list[dict[str, Any]]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
