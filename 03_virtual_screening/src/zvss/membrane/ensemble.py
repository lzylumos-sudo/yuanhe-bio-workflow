"""Multi-conformation evidence for Y-MAVS.

The ensemble layer reads real per-state docking evidence and keeps each state's
scores, pose identity, QA and membrane/accessibility evidence separate.  It never
averages GNINA scores across states and never fills a missing state with zero.

Two perspectives are intentionally reported side by side:

* ``state_selective_ranking`` -- ordering within one conformational state;
* ``cross_state_support`` -- how many states actually contain measured evidence
  for a compound.  This is support breadth, not a potency estimate, and is also
  exposed under the requested key ``cross_state_robustness`` with that caveat.

Structural alignment is a C-alpha Kabsch fit with residue identity, chain and
missing-residue checks.  An alignment is a coordinate aid; it does not make two
different conformational states experimentally comparable.
"""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from ..ids import sha256_file
from ..util import atomic_write_json, now_iso, parse_pdb_text, read_text
from .regions import MembraneGeometry
from .settings import ensemble_rules

ALGORITHM_VERSION = "ymavs-ensemble/1.0.0"
EVIDENCE_SCHEMA = "ymavs-ensemble-state-evidence/1"


class MembraneEnsembleAnalyzer:
    def __init__(
        self,
        ctx,
        *,
        settings: dict[str, Any],
        primary_geometry: MembraneGeometry | None = None,
        primary_receptor_sha256: str = "",
        primary_structure_path: str = "",
    ):
        self.ctx = ctx
        self.settings = settings
        self.rules = ensemble_rules(settings)
        self.primary_geometry = primary_geometry
        self.primary_receptor_sha256 = primary_receptor_sha256
        self.primary_structure_path = primary_structure_path

    def run(self) -> dict[str, Any]:
        if not self.rules.get("enabled"):
            return {
                "schema": "ymavs-ensemble-analysis/1",
                "algorithm_version": ALGORITHM_VERSION,
                "generated_at": now_iso(),
                "enabled": False,
                "ensemble_status": "disabled",
                "states": [],
                "compounds": {},
                "notes": ["ensemble screening is disabled for this target; no state evidence was combined"],
            }
        receptors = list(self.rules.get("receptors") or [])
        if not receptors:
            return {
                "schema": "ymavs-ensemble-analysis/1",
                "algorithm_version": ALGORITHM_VERSION,
                "generated_at": now_iso(),
                "enabled": True,
                "ensemble_status": "needs_input",
                "states": [],
                "compounds": {},
                "reason": "ensemble is enabled but no receptors/states are configured",
            }
        primary_pdb = str(getattr(self.ctx.config.target.primary, "pdb_id", "") or "").upper()
        primary_chain = str(getattr(self.ctx.config.target.primary, "chain", "") or "")
        states: list[dict[str, Any]] = []
        for index, spec in enumerate(receptors):
            state = self._analyze_state(
                spec,
                index=index,
                primary_pdb=primary_pdb,
                primary_chain=primary_chain,
            )
            states.append(state)
        compounds = self._integrate(states)
        evaluated = sum(1 for state in states if state.get("evidence", {}).get("status") == "evaluated")
        status = "complete" if len(states) and evaluated == len(states) else ("partial" if evaluated else "needs_input")
        return {
            "schema": "ymavs-ensemble-analysis/1",
            "algorithm_version": ALGORITHM_VERSION,
            "generated_at": now_iso(),
            "enabled": True,
            "ensemble_status": status,
            "state_count": len(states),
            "states_with_evidence": evaluated,
            "states": states,
            "compounds": compounds,
            "state_selective_ranking": self._state_selective_ranking(states, compounds),
            "cross_state_robustness": {
                "definition": "number/fraction of configured states with measured evidence; not a potency or efficacy prediction",
                "per_compound": {
                    compound_id: {
                        "evaluated_states": item["number_of_evaluated_states"],
                        "missing_states": item["number_of_missing_states"],
                        "coverage_fraction": item["state_coverage"],
                    }
                    for compound_id, item in compounds.items()
                },
            },
            "limitations": [
                "per-state GNINA scores are not averaged; different constructs and states are not directly interchangeable",
                "missing state evidence is reported as missing, never as a zero score",
                "C-alpha alignment does not establish that conformational states are biologically equivalent",
                "single-state binding can be meaningful for a state-selective inhibitor and is reported separately",
            ],
        }

    # ------------------------------------------------------------------ #
    # state loading and structure mapping
    # ------------------------------------------------------------------ #
    def _analyze_state(
        self,
        spec: dict[str, Any],
        *,
        index: int,
        primary_pdb: str,
        primary_chain: str,
    ) -> dict[str, Any]:
        structure_id = str(spec.get("structure_id") or spec.get("pdb_id") or f"state_{index + 1}")
        conformational_state = str(spec.get("conformational_state") or "")
        state = {
            "structure_id": structure_id,
            "pdb_id": str(spec.get("pdb_id") or structure_id).upper(),
            "chain": str(spec.get("chain") or primary_chain),
            "conformational_state": conformational_state,
            "role": str(spec.get("role") or ""),
            "is_primary": bool(spec.get("is_primary", False)) or structure_id.upper() == primary_pdb,
            "structure_path": None,
            "structure_sha256": None,
            "expected_sha256": spec.get("sha256") or None,
            "mapping": {"status": "not_evaluated"},
            "membrane": {"status": "not_evaluated"},
            "evidence": {"status": "missing", "compounds": {}},
        }
        if state["is_primary"] and self.primary_structure_path and os.path.exists(self.primary_structure_path):
            state["structure_path"] = self.primary_structure_path
            state["structure_sha256"] = self.primary_receptor_sha256 or sha256_file(self.primary_structure_path)
            if self.primary_geometry is not None and self.primary_geometry.usable:
                state["membrane"] = {
                    "status": self.primary_geometry.orientation_status,
                    "source": self.primary_geometry.source,
                    "coordinate_frame": self.primary_geometry.coordinate_frame,
                    "verified": self.primary_geometry.verified,
                }
        else:
            path = self._resolve_structure_path(spec, structure_id)
            if path:
                state["structure_path"] = path
                state["structure_sha256"] = sha256_file(path)
        if state["expected_sha256"] and state["structure_sha256"] and state["expected_sha256"] != state["structure_sha256"]:
            state["structure_status"] = "sha256_mismatch"
        elif state["structure_path"]:
            state["structure_status"] = "loaded"
        else:
            state["structure_status"] = "needs_input"
        if state["structure_path"] and state["structure_path"].lower().endswith((".cif", ".mmcif")):
            state["structure_status"] = "unsupported_format"
            state["structure_path"] = None
        if state["structure_path"]:
            mapping = self._map_to_primary(
                state["structure_path"],
                chain=state["chain"],
                primary_pdb=primary_pdb,
                primary_chain=primary_chain,
                allow_mutations=bool(spec.get("allow_mutations", False)),
            )
            state["mapping"] = mapping
            if mapping.get("status") == "valid":
                primary_membrane_already = state["is_primary"] and state.get("membrane", {}).get("status") in {
                    "verified_computational", "proxy_only"
                }
                has_state_specific = bool(spec.get("membrane_orientation") or spec.get("evidence_dir"))
                if not primary_membrane_already or has_state_specific:
                    state["membrane"] = self._state_membrane(spec, state["structure_path"], mapping)
        state["evidence"] = self._load_state_evidence(spec, state)
        return state

    def _resolve_structure_path(self, spec: dict[str, Any], structure_id: str) -> str | None:
        if spec.get("structure_path"):
            candidate = Path(str(spec["structure_path"]))
            if not candidate.is_absolute():
                candidate = Path(self.ctx.config.target.source_path).parent / candidate
            return str(candidate.resolve()) if candidate.exists() else None
        pdb_id = str(spec.get("pdb_id") or structure_id).upper()
        directories = []
        target = self.ctx.config.target
        if getattr(target, "structures_dir", ""):
            directories.append(Path(target.structures_dir))
        directories.append(Path(self.ctx.ref_dir()))
        directories.append(Path(self.ctx.root) / "runtime" / "reference_cache")
        for directory in directories:
            for suffix in (".pdb", ".ent", ".cif", ".mmcif"):
                for name in (pdb_id, pdb_id.lower()):
                    candidate = directory / f"{name}{suffix}"
                    if candidate.exists():
                        return str(candidate)
            if directory.is_dir():
                for candidate in sorted(directory.iterdir()):
                    if candidate.is_file() and candidate.name.upper().startswith(pdb_id) and candidate.suffix.lower() in {".pdb", ".ent", ".cif", ".mmcif"}:
                        return str(candidate)
        return None

    def _map_to_primary(
        self,
        structure_path: str,
        *,
        chain: str,
        primary_pdb: str,
        primary_chain: str,
        allow_mutations: bool,
    ) -> dict[str, Any]:
        reference_path = self.primary_structure_path or str(Path(self.ctx.ref_dir()) / f"{primary_pdb}.pdb")
        if not os.path.exists(reference_path):
            return {"status": "needs_input", "reason": "primary receptor PDB is not available for alignment"}
        try:
            reference_ca = self._ca_map(reference_path, primary_chain)
            member_ca = self._ca_map(structure_path, chain)
        except (ValueError, OSError) as exc:
            return {"status": "invalid", "reason": f"{type(exc).__name__}: {exc}"}
        shared = sorted(set(reference_ca) & set(member_ca))
        mismatches = [
            {"key": list(key), "primary": reference_ca[key][0], "state": member_ca[key][0]}
            for key in shared
            if reference_ca[key][0] != member_ca[key][0]
        ]
        usable = [key for key in shared if allow_mutations or reference_ca[key][0] == member_ca[key][0]]
        missing_in_state = sorted(set(reference_ca) - set(member_ca))
        extra_in_state = sorted(set(member_ca) - set(reference_ca))
        min_ca = int(self.rules.get("alignment_min_ca", 12))
        if len(usable) < max(12, min_ca):
            return {
                "status": "invalid",
                "reason": "insufficient matched C-alpha atoms for a rigid alignment",
                "matched_ca": len(usable),
                "minimum_ca": max(12, min_ca),
                "residue_mismatches": mismatches[:50],
                "missing_in_state": [list(key) for key in missing_in_state[:50]],
                "extra_in_state": [list(key) for key in extra_in_state[:50]],
            }
        primary_coords = np.asarray([reference_ca[key][1] for key in usable], dtype=float)
        state_coords = np.asarray([member_ca[key][1] for key in usable], dtype=float)
        from ..util import kabsch_transform

        rotation, translation, rmsd = kabsch_transform(state_coords, primary_coords)
        max_rmsd = float(self.rules.get("alignment_max_rmsd_A", 2.5))
        deviations = np.sqrt(((state_coords @ rotation.T + translation - primary_coords) ** 2).sum(axis=1))
        status = "valid" if rmsd <= max_rmsd else "invalid"
        return {
            "status": status,
            "method": "Kabsch_Calpha_alignment",
            "matched_ca": len(usable),
            "residue_identity_fraction": round(
                sum(1 for key in usable if reference_ca[key][0] == member_ca[key][0]) / max(len(usable), 1), 6
            ),
            "rmsd_A": round(float(rmsd), 4),
            "max_rmsd_A": max_rmsd,
            "max_residue_deviation_A": round(float(deviations.max()), 4),
            "residue_mismatches": mismatches[:50],
            "missing_in_state": [list(key) for key in missing_in_state[:50]],
            "extra_in_state": [list(key) for key in extra_in_state[:50]],
            "state_to_primary_rotation": rotation.tolist(),
            "state_to_primary_translation": translation.tolist(),
            "state_to_primary_note": "the state structure may be mapped into the primary frame with x @ rotation.T + translation; docking scores remain in their own state frame",
            "membrane_copied_from_primary": False,
            "reason": None if status == "valid" else f"alignment RMSD {rmsd:.3f} A exceeds {max_rmsd:.3f} A",
        }

    def _ca_map(self, path: str, chain: str) -> dict[tuple[str, int, str], tuple[str, tuple[float, float, float]]]:
        structure = parse_pdb_text(read_text(path), keep_hydrogens=False)
        result: dict[tuple[str, int, str], tuple[str, tuple[float, float, float]]] = {}
        for atom in structure.atoms:
            if atom.record.strip() != "ATOM" or atom.name.strip() != "CA":
                continue
            if chain and atom.chain != chain:
                continue
            key = (atom.chain, atom.resseq, atom.icode)
            if key in result:
                raise ValueError(f"ambiguous C-alpha atom for residue {key}")
            result[key] = (atom.resname.upper(), (float(atom.x), float(atom.y), float(atom.z)))
        if not result:
            raise ValueError("no C-alpha atoms found for the requested chain")
        return result

    def _state_membrane(self, spec: dict[str, Any], structure_path: str, mapping: dict[str, Any]) -> dict[str, Any]:
        block = spec.get("membrane_orientation") or {}
        evidence_dir = block.get("evidence_dir") or spec.get("evidence_dir")
        if evidence_dir:
            folder = Path(str(evidence_dir))
            if not folder.is_absolute():
                folder = Path(self.ctx.config.target.source_path).parent / folder
            try:
                from .orientation import transfer

                plane = transfer(folder, structure_path)
                return {
                    "status": "verified_computational",
                    "verified": True,
                    "source": plane.get("source"),
                    "method": plane.get("method"),
                    "coordinate_frame": plane.get("coordinate_frame"),
                    "receptor_sha256": plane.get("receptor_sha256"),
                    "plane": plane,
                    "copied_from_primary": False,
                }
            except (ValueError, OSError, KeyError, TypeError) as exc:
                return {"status": "needs_input", "verified": False, "reason": str(exc), "copied_from_primary": False}
        if block.get("center") and block.get("normal") and block.get("thickness_A"):
            try:
                geometry = MembraneGeometry.from_plane(
                    {
                        "center": block["center"],
                        "normal": block["normal"],
                        "thickness_A": block["thickness_A"],
                        "source": block.get("source", ""),
                        "method": block.get("method", ""),
                        "coordinate_frame": block.get("coordinate_frame", "state_structure"),
                        "receptor_sha256": sha256_file(structure_path),
                    },
                    orientation_status="verified_computational" if block.get("verified") is True else "proxy_only",
                )
            except (ValueError, TypeError) as exc:
                return {"status": "needs_input", "verified": False, "reason": str(exc), "copied_from_primary": False}
            return {
                "status": geometry.orientation_status,
                "verified": geometry.verified,
                "source": geometry.source,
                "method": geometry.method,
                "coordinate_frame": geometry.coordinate_frame,
                "receptor_sha256": geometry.receptor_sha256,
                "geometry": geometry.to_dict(buffer_A=float(block.get("interface_buffer_A", 5.0))),
                "copied_from_primary": False,
            }
        return {
            "status": "needs_input",
            "verified": False,
            "reason": "no state-specific membrane orientation evidence; primary orientation is not copied",
            "copied_from_primary": False,
        }

    # ------------------------------------------------------------------ #
    # evidence loading
    # ------------------------------------------------------------------ #
    def _load_state_evidence(self, spec: dict[str, Any], state: dict[str, Any]) -> dict[str, Any]:
        if state.get("structure_status") != "loaded" or not state.get("structure_path"):
            return {
                "status": "needs_input",
                "reason": "state receptor is missing, malformed or has an unverified SHA; evidence is not used",
                "compounds": {},
            }
        path_value = spec.get("docking_results")
        if state.get("is_primary") and not path_value:
            return self._collect_primary_evidence(state)
        if not path_value:
            return {
                "status": "missing",
                "reason": "no docking_results manifest configured for this state",
                "compounds": {},
                "source": None,
            }
        path = Path(str(path_value))
        if not path.is_absolute():
            path = Path(self.ctx.config.target.source_path).parent / path
        if not path.exists():
            return {"status": "missing", "reason": f"docking evidence manifest not found: {path}", "compounds": {}}
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            return {"status": "invalid", "reason": f"cannot read evidence manifest: {type(exc).__name__}: {exc}", "compounds": {}}
        if payload.get("schema") not in (EVIDENCE_SCHEMA, "ymavs-ensemble-state-evidence/1"):
            return {"status": "invalid", "reason": "unsupported state evidence schema", "compounds": {}}
        if str(payload.get("structure_id")) != str(state["structure_id"]):
            return {"status": "invalid", "reason": "evidence structure_id does not match configured state", "compounds": {}}
        expected = state.get("expected_sha256") or state.get("structure_sha256")
        if payload.get("receptor_sha256") and expected and payload["receptor_sha256"] != expected:
            return {"status": "invalid", "reason": "evidence receptor_sha256 does not match the receptor structure", "compounds": {}}
        compounds = payload.get("compounds") or {}
        if not isinstance(compounds, dict):
            return {"status": "invalid", "reason": "compounds must be an object keyed by compound_id", "compounds": {}}
        normalized: dict[str, dict[str, Any]] = {}
        for compound_id, record in compounds.items():
            if not isinstance(record, dict):
                return {"status": "invalid", "reason": f"compound {compound_id} record is not an object", "compounds": {}}
            normalized[str(compound_id)] = dict(record)
        return {
            "status": "evaluated",
            "source": str(path),
            "sha256": sha256_file(str(path)),
            "structure_id": payload.get("structure_id"),
            "coordinate_frame": payload.get("coordinate_frame"),
            "compounds": normalized,
        }

    def _collect_primary_evidence(self, state: dict[str, Any]) -> dict[str, Any]:
        """Read the current run's real pose/assessment rows for the primary state."""
        from .. import data
        from ..scoring.aggregate import representative_pose

        try:
            poses = data.pose_rows(self.ctx)
        except Exception as exc:  # noqa: BLE001
            return {"status": "missing", "reason": f"{type(exc).__name__}: {exc}", "compounds": {}}
        by_compound: dict[str, list[dict[str, Any]]] = {}
        for pose in poses:
            by_compound.setdefault(pose["compound_id"], []).append(pose)
        compounds: dict[str, dict[str, Any]] = {}
        for compound_id, group in sorted(by_compound.items()):
            chosen = representative_pose(group)
            scores = chosen.get("scores") or {}
            metrics = self.ctx.db.metric_map("pose", chosen["pose_id"], "membrane")
            compounds[compound_id] = {
                "state_id": chosen.get("state_id"),
                "pose_id": chosen["pose_id"],
                "engine": chosen.get("engine"),
                "vina_affinity_kcal_mol":scores.get("vina_affinity_kcal_mol"),
                "gnina_affinity":scores.get("minimizedAffinity") if chosen.get("engine")!="vina" else None,
                "cnn_score": scores.get("CNNscore"),
                "cnn_affinity": scores.get("CNNaffinity"),
                "pose_qc_status": chosen.get("assess_status"),
                "membrane_evidence_status": metrics.get("membrane_evidence_status", {}).get("value_text"),
                "state_specific_contacts": metrics.get("pocket_contact_count", {}).get("value"),
                "state_specific_accessibility": {
                    "aqueous": metrics.get("aqueous_accessibility", {}).get("value_text"),
                    "membrane_lateral": metrics.get("membrane_lateral_accessibility", {}).get("value_text"),
                },
                "source": "current run database",
            }
        return {
            "status": "evaluated" if compounds else "missing",
            "source": "current run database",
            "compounds": compounds,
            "note": "primary-state evidence is read from the same run; no score is synthesized",
        }

    # ------------------------------------------------------------------ #
    # integration
    # ------------------------------------------------------------------ #
    def _integrate(self, states: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
        all_compounds = sorted({compound for state in states for compound in (state.get("evidence", {}).get("compounds") or {})})
        total_states = len(states)
        result: dict[str, dict[str, Any]] = {}
        for compound_id in all_compounds:
            per_state = {}
            for state in states:
                evidence = state.get("evidence", {})
                record = (evidence.get("compounds") or {}).get(compound_id)
                per_state[state["structure_id"]] = {
                    "conformational_state": state.get("conformational_state"),
                    "evidence_status": "evaluated" if record else "missing",
                    "record": record,
                }
            evaluated = [key for key, item in per_state.items() if item["evidence_status"] == "evaluated"]
            missing = [key for key, item in per_state.items() if item["evidence_status"] != "evaluated"]
            best_state = None
            if evaluated:
                best_state = sorted(
                    evaluated,
                    key=lambda key: (
                        0 if (per_state[key]["record"] or {}).get("formal_qualified") is True else 1,
                        -(int((per_state[key]["record"] or {}).get("membrane_evidence_rank", 0) or 0)),
                        key,
                    ),
                )[0]
            result[compound_id] = {
                "compound_id": compound_id,
                "state_specific_evidence": per_state,
                "number_of_evaluated_states": len(evaluated),
                "number_of_missing_states": len(missing),
                "state_coverage": round(len(evaluated) / total_states, 6) if total_states else 0.0,
                "best_supported_state": best_state,
                "missing_states": missing,
                "state_specific_contacts": {
                    key: (item["record"] or {}).get("state_specific_contacts")
                    for key, item in per_state.items()
                    if item["record"] is not None
                },
                "state_specific_accessibility": {
                    key: (item["record"] or {}).get("state_specific_accessibility")
                    for key, item in per_state.items()
                    if item["record"] is not None
                },
                "ensemble_status": "complete" if not missing else ("partial" if evaluated else "missing"),
                "note": "missing states are not zero scores; state selectivity is reported separately from coverage",
            }
        return result

    def _state_selective_ranking(self, states: list[dict[str, Any]], compounds: dict[str, dict[str, Any]]) -> dict[str, Any]:
        rankings: dict[str, Any] = {}
        for state in states:
            records = []
            for compound_id, item in compounds.items():
                entry = item["state_specific_evidence"].get(state["structure_id"], {})
                record = entry.get("record") if entry else None
                if not record:
                    continue
                affinity = record.get("vina_affinity_kcal_mol") if record.get("engine")=="vina" else record.get("gnina_affinity")
                try:
                    affinity_value = float(affinity) if affinity is not None else None
                except (TypeError, ValueError):
                    affinity_value = None
                records.append(
                    {
                        "compound_id": compound_id,
                        "pose_id": record.get("pose_id"),
                        "engine":record.get("engine","gnina"), "docking_affinity_kcal_mol":affinity_value,
                        "gnina_affinity":affinity_value if record.get("engine")!="vina" else None,
                        "pose_qc_status": record.get("pose_qc_status"),
                    }
                )
            if len({r["engine"] for r in records})>1:raise ValueError("Mixed-engine ensemble rankings require separate runs")
            records.sort(key=lambda row: (row["docking_affinity_kcal_mol"] is None, row["docking_affinity_kcal_mol"] if row["docking_affinity_kcal_mol"] is not None else math.inf, row["compound_id"]))
            rankings[state["structure_id"]] = {
                "conformational_state": state.get("conformational_state"),
                "state_coverage": f"{len(records)}/{len(compounds)}" if compounds else "0/0",
                "ranking_basis": "within-state engine-labelled affinity only; missing compounds are not scored",
                "ranked": records,
            }
        return {
            "definition": "Within each conformational state independently; no cross-state score averaging.",
            "states": rankings,
        }

    def write_outputs(self, directory: str, payload: dict[str, Any] | None = None) -> str:
        payload = payload or self.run()
        path = os.path.join(directory, "ensemble_analysis.json")
        atomic_write_json(path, payload)
        csv_path = os.path.join(directory, "ensemble_state_coverage.csv")
        fields = [
            "compound_id", "evaluated_states", "missing_states", "state_coverage",
            "best_supported_state", "ensemble_status",
        ]
        rows = [
            {
                "compound_id": compound_id,
                "evaluated_states": item.get("number_of_evaluated_states"),
                "missing_states": item.get("number_of_missing_states"),
                "state_coverage": item.get("state_coverage"),
                "best_supported_state": item.get("best_supported_state"),
                "ensemble_status": item.get("ensemble_status"),
            }
            for compound_id, item in (payload.get("compounds") or {}).items()
        ]
        write_csv(csv_path, fields, rows)
        return path


def write_csv(path: str, fields: list[str], rows: Iterable[dict[str, Any]]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
