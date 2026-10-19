"""Membrane compatibility descriptors for Y-MAVS.

This module describes where a docked ligand sits relative to the membrane and
which local residue classes it contacts.  It also reports a van der Waals
overlap proxy and, only when explicit lipid coordinates are supplied in the same
verified frame, a lipid collision check.

It deliberately does not convert LogP, TPSA or membrane depth into an activity
score.  Those descriptors can explain a hypothesis; they are not calibrated
binding predictions and are not used as a favorable ranking weight.
"""

from __future__ import annotations

import csv
import math
import os
from typing import Any, Iterable

import numpy as np

from ..util import atomic_write_json, now_iso
from .geom import (
    as_coords,
    element_radius,
    nearest_distances,
    receptor_contacts_by_residue,
    residue_class,
)
from .regions import MembraneGeometry, MembraneRegionAnalyzer

ALGORITHM_VERSION = "ymavs-compatibility/1.0.0"


class MembraneCompatibilityAnalyzer:
    def __init__(
        self,
        geometry: MembraneGeometry,
        region_analyzer: MembraneRegionAnalyzer,
        receptor_atoms: Iterable[Any],
        *,
        contact_cutoff_A: float = 4.5,
        explicit_lipid_xyz: np.ndarray | None = None,
        explicit_lipid_source: str = "",
        explicit_lipid_frame: str = "",
        lipid_clash_cutoff_A: float = 2.0,
    ):
        self.geometry = geometry
        self.regions = region_analyzer
        self.receptor_atoms = list(receptor_atoms)
        self.contact_cutoff_A = float(contact_cutoff_A)
        self.explicit_lipid_xyz = (
            None if explicit_lipid_xyz is None else np.asarray(explicit_lipid_xyz, dtype=float).reshape((-1, 3))
        )
        self.explicit_lipid_source = explicit_lipid_source
        self.explicit_lipid_frame = explicit_lipid_frame
        self.lipid_clash_cutoff_A = float(lipid_clash_cutoff_A)

    def analyze(
        self,
        ligand_xyz: np.ndarray,
        *,
        pose_id: str = "",
        compound_id: str = "",
        mol: Any = None,
        ligand_elements: Iterable[str] | None = None,
        existing_metrics: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        ligand = np.asarray(ligand_xyz, dtype=float).reshape((-1, 3))
        if len(ligand) == 0:
            return self._empty(pose_id, compound_id, reason="ligand has no coordinates")
        elements = [str(value) for value in (ligand_elements or [])]
        depth = self.regions.analyze_points(ligand)
        contacts = receptor_contacts_by_residue(ligand, self.receptor_atoms, self.contact_cutoff_A)
        local = self._local_contacts(contacts, existing_metrics)
        steric = self._steric(ligand, elements)
        descriptors = self._descriptors(mol)
        lipid = self._lipid(ligand)
        warnings = []
        if depth["core_fraction"] is not None and depth["core_fraction"] >= 0.5:
            warnings.append("ligand_majority_in_membrane_core_depth_range")
        if (
            depth["core_fraction"] is not None
            and depth["outside_core_fraction"] is not None
            and depth["core_fraction"] >= 0.2
            and depth["outside_core_fraction"] >= 0.2
        ):
            warnings.append("ligand_spans_core_and_outside_depth_ranges")
        if descriptors.get("status") != "evaluated":
            warnings.append("molecular_descriptors_not_available")
        if lipid.get("lipid_clash_status") == "not_evaluated":
            warnings.append("explicit_lipid_atoms_not_evaluated")
        warnings.append("membrane_depth_and_logp_are_not_binding_affinity")
        status = "evaluated" if self.geometry.usable else "needs_input"
        result = {
            "schema": "ymavs-membrane-compatibility/1",
            "algorithm_version": ALGORITHM_VERSION,
            "generated_at": now_iso(),
            "status": status,
            "pose_id": pose_id,
            "compound_id": compound_id,
            "ligand_center_depth_A": depth["center_depth_A"],
            "ligand_depth_min_A": depth["depth_min_A"],
            "ligand_depth_max_A": depth["depth_max_A"],
            "ligand_depth_median_A": depth["depth_median_A"],
            "ligand_core_fraction": depth["core_fraction"],
            "ligand_interface_fraction": depth["interface_fraction"],
            "ligand_outside_core_fraction": depth["outside_core_fraction"],
            "ligand_region": depth["region"],
            "local_hydrophobic_contact_count": local["hydrophobic"],
            "local_polar_contact_count": local["polar"],
            "local_charged_contact_count": local["charged"],
            "local_other_contact_count": local["other"],
            "local_contact_total": local["total"],
            "local_contact_source": local["source"],
            "local_existing_mechanism": local.get("existing_mechanism_counts"),
            "contact_cutoff_A": self.contact_cutoff_A,
            "steric_compatibility_status": steric["status"],
            "vdw_overlap_pairs": steric["overlap_pairs"],
            "vdw_severe_overlap_pairs": steric["severe_overlap_pairs"],
            "ligand_receptor_min_distance_A": steric["min_distance_A"],
            "steric_method": "element_vdw_radii_geometric_proxy",
            "molecular_descriptors": descriptors,
            "lipid_clash_status": lipid["lipid_clash_status"],
            "lipid_clash_count": lipid["lipid_clash_count"],
            "lipid_minimum_distance_A": lipid["lipid_minimum_distance_A"],
            "lipid_source": lipid["lipid_source"],
            "lipid_coordinate_frame": lipid["lipid_coordinate_frame"],
            "orientation_status": self.geometry.orientation_status,
            "coordinate_frame": self.geometry.coordinate_frame,
            "evidence_level": "sourced_geometry" if self.geometry.verified else "geometric_proxy",
            "warnings": warnings,
            "limitations": [
                "membrane depth describes location, not affinity or permeability",
                "local residue composition is not a solvent-exposure measurement",
                "van der Waals overlap is a geometry proxy and does not replace PoseBusters",
                "an explicit lipid collision is only evaluated when same-frame lipid coordinates are supplied",
                "LogP/TPSA and related descriptors are reported for interpretation, not as binding scores",
            ],
        }
        return result

    def _empty(self, pose_id: str, compound_id: str, *, reason: str) -> dict[str, Any]:
        return {
            "schema": "ymavs-membrane-compatibility/1",
            "algorithm_version": ALGORITHM_VERSION,
            "generated_at": now_iso(),
            "status": "not_evaluated",
            "pose_id": pose_id,
            "compound_id": compound_id,
            "reason": reason,
            "orientation_status": self.geometry.orientation_status,
            "coordinate_frame": self.geometry.coordinate_frame,
            "molecular_descriptors": {"status": "not_evaluated"},
            "lipid_clash_status": "not_evaluated",
            "steric_compatibility_status": "not_evaluated",
        }

    def _local_contacts(self, contacts: list[dict[str, Any]], existing: dict[str, Any] | None) -> dict[str, Any]:
        counts = {"hydrophobic": 0, "polar": 0, "charged": 0, "other": 0}
        for row in contacts:
            category = residue_class(row.get("resname"))
            counts[category] = counts.get(category, 0) + 1
        result: dict[str, Any] = {
            **counts,
            "total": len(contacts),
            "source": "receptor_contacts_by_residue at the compatibility contact cutoff",
        }
        if existing and existing.get("hydrophobic_contacts") is not None:
            result["existing_mechanism_counts"] = {
                "hydrophobic": existing.get("hydrophobic_contacts"),
                "polar": existing.get("polar_contacts"),
                "total": existing.get("contact_count"),
                "note": "retained for comparison; not mixed into the Y-MAVS counts above",
            }
        return result

    def _steric(self, ligand: np.ndarray, elements: list[str]) -> dict[str, Any]:
        atoms = [a for a in self.receptor_atoms if not (hasattr(a, "is_hydrogen") and a.is_hydrogen())]
        if not atoms:
            return {
                "status": "not_evaluated",
                "overlap_pairs": None,
                "severe_overlap_pairs": None,
                "min_distance_A": None,
            }
        rec_xyz = as_coords(atoms)
        rec_radii = np.array([element_radius(getattr(atom, "element", "")) for atom in atoms], dtype=float)
        ligand_radii = np.array(
            [element_radius(elements[index] if index < len(elements) else "") for index in range(len(ligand))],
            dtype=float,
        )
        distances = nearest_distances(ligand, rec_xyz)
        minimum = float(distances.min()) if len(distances) else None
        try:
            from scipy.spatial import cKDTree

            tree = cKDTree(rec_xyz)
            overlap_pairs = 0
            severe_pairs = 0
            for index, point in enumerate(ligand):
                # Query a generous shell; the pair-specific radius test is done
                # below with the actual element radii.
                for other in tree.query_ball_point(point, float(ligand_radii[index] + rec_radii.max())):
                    distance = float(np.linalg.norm(point - rec_xyz[other]))
                    radius_sum = float(ligand_radii[index] + rec_radii[other])
                    if radius_sum <= 0:
                        continue
                    if distance < 0.85 * radius_sum:
                        overlap_pairs += 1
                    if distance < 0.65 * radius_sum:
                        severe_pairs += 1
        except ImportError:  # pragma: no cover
            overlap_pairs = None
            severe_pairs = None
        if severe_pairs is None or overlap_pairs is None:
            status = "not_evaluated"
        elif severe_pairs > 0:
            status = "clash"
        elif overlap_pairs > 0:
            status = "minor_overlap"
        else:
            status = "compatible"
        return {
            "status": status,
            "overlap_pairs": overlap_pairs,
            "severe_overlap_pairs": severe_pairs,
            "min_distance_A": None if minimum is None or not math.isfinite(minimum) else round(minimum, 3),
        }

    def _descriptors(self, mol: Any) -> dict[str, Any]:
        if mol is None:
            return {"status": "not_evaluated", "reason": "no RDKit molecule supplied"}
        try:
            from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

            values = {
                "MolLogP": float(Crippen.MolLogP(mol)),
                "TPSA": float(rdMolDescriptors.CalcTPSA(mol)),
                "FormalCharge": int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms())),
                "HBD": int(Lipinski.NumHDonors(mol)),
                "HBA": int(Lipinski.NumHAcceptors(mol)),
                "MolecularWeight": float(Descriptors.MolWt(mol)),
                "RotatableBonds": int(Lipinski.NumRotatableBonds(mol)),
            }
        except Exception as exc:  # noqa: BLE001 - descriptor availability is not a run blocker
            return {"status": "unsupported", "reason": f"{type(exc).__name__}: {exc}"}
        return {
            "status": "evaluated",
            "method": "RDKit descriptors",
            "values": values,
            "units": {
                "MolLogP": "dimensionless (Crippen logP)",
                "TPSA": "angstrom^2",
                "FormalCharge": "elementary charge",
                "MolecularWeight": "g/mol",
                "RotatableBonds": "count",
            },
            "usage": "descriptive only; not used as an affinity score in this version",
        }

    def _lipid(self, ligand: np.ndarray) -> dict[str, Any]:
        if self.explicit_lipid_xyz is None or len(self.explicit_lipid_xyz) == 0:
            return {
                "lipid_clash_status": "not_evaluated",
                "lipid_clash_count": None,
                "lipid_minimum_distance_A": None,
                "lipid_source": None,
                "lipid_coordinate_frame": None,
            }
        if not self.explicit_lipid_frame or self.explicit_lipid_frame != self.geometry.coordinate_frame:
            return {
                "lipid_clash_status": "needs_input",
                "lipid_clash_count": None,
                "lipid_minimum_distance_A": None,
                "lipid_source": self.explicit_lipid_source or None,
                "lipid_coordinate_frame": self.explicit_lipid_frame or None,
                "reason": "explicit lipid coordinates are not in the verified docking-receptor frame",
            }
        try:
            from scipy.spatial import cKDTree

            distances, _ = cKDTree(self.explicit_lipid_xyz).query(ligand, k=1, workers=-1)
            minimum = float(distances.min())
            clashes = int((distances < self.lipid_clash_cutoff_A).sum())
            status = "evaluated"
        except ImportError:  # pragma: no cover
            minimum = float("nan")
            clashes = None
            status = "unsupported"
        return {
            "lipid_clash_status": status,
            "lipid_clash_count": clashes,
            "lipid_minimum_distance_A": None if not math.isfinite(minimum) else round(minimum, 3),
            "lipid_source": self.explicit_lipid_source or None,
            "lipid_coordinate_frame": self.explicit_lipid_frame or None,
            "lipid_clash_cutoff_A": self.lipid_clash_cutoff_A,
        }

    def write_outputs(self, directory: str, analyses: list[dict[str, Any]]) -> dict[str, str]:
        path = os.path.join(directory, "membrane_compatibility.json")
        atomic_write_json(
            path,
            {
                "schema": "ymavs-membrane-compatibility-collection/1",
                "algorithm_version": ALGORITHM_VERSION,
                "generated_at": now_iso(),
                "count": len(analyses),
                "analyses": analyses,
            },
        )
        fields = [
            "pose_id", "compound_id", "status", "ligand_center_depth_A", "ligand_depth_min_A",
            "ligand_depth_max_A", "ligand_core_fraction", "ligand_interface_fraction",
            "ligand_outside_core_fraction", "ligand_region", "local_hydrophobic_contact_count",
            "local_polar_contact_count", "local_charged_contact_count", "local_other_contact_count",
            "steric_compatibility_status", "vdw_overlap_pairs", "vdw_severe_overlap_pairs",
            "ligand_receptor_min_distance_A", "lipid_clash_status", "lipid_clash_count",
            "lipid_minimum_distance_A", "orientation_status", "coordinate_frame",
        ]
        csv_path = os.path.join(directory, "membrane_compatibility.csv")
        _write_csv(csv_path, fields, analyses)
        return {"membrane_compatibility": path, "membrane_compatibility_csv": csv_path}


def _write_csv(path: str, fields: list[str], rows: list[dict[str, Any]]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
