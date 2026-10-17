"""Ligand strain, property filters and diversity grouping.

Strain is a force-field difference between the docked conformer and a locally
minimised conformer of the same molecule (MMFF94 where parameterised, UFF
otherwise, with the force field named in the record).  It is a heuristic, not a
binding contribution, and is labelled as such.
"""

from __future__ import annotations

from typing import Any

from ..chem import rdio
from ..chem.props import alert_burden, compute_properties
from ..constants import AssessmentStatus, Reason

STRAIN_FLAG = 25.0
STRAIN_FAIL = 60.0


def run(ctx) -> dict[str, Any]:
    from .. import data

    poses = data.pose_rows(ctx)
    budget = ctx.config.budget
    from ..scoring.aggregate import representative_pose
    groups: dict[str, list[dict[str, Any]]] = {}
    for pose in poses:
        groups.setdefault(pose["compound_id"], []).append(pose)
    # The strain must belong to exactly the pose used by the ranking.
    by_compound = {key: representative_pose(group) for key, group in groups.items()}

    summary: dict[str, Any] = {"compounds": len(by_compound), "flagged": 0, "failed": 0, "scaffolds": {}}
    for compound_id in sorted(by_compound):
        pose = by_compound[compound_id]
        mol = data.pose_mol(ctx, pose)
        if mol is None:
            continue
        properties = compute_properties(mol)
        for key, value in properties.values.items():
            if isinstance(value, (int, float)) and value is not None:
                ctx.db.record_assessment("compound", compound_id, "properties", key, float(value), AssessmentStatus.PASS)
        ctx.db.record_assessment("compound", compound_id, "properties", "alert_burden", alert_burden(properties), AssessmentStatus.PASS)

        strain = _strain(mol, budget)
        strain["pose_id"] = pose["pose_id"]
        strain_status = str(strain.get("status") or "not_evaluated")
        strain_code = {
            "fail": Reason.EXCESSIVE_STRAIN,
            "not_evaluated": Reason.STRAIN_NOT_CONVERGED,
        }.get(strain_status)
        for key, value in strain.items():
            if key == "status":
                ctx.db.record_assessment("compound", compound_id, "strain", "status", None,
                                         AssessmentStatus(strain_status), strain_code, value_text=strain_status)
            elif isinstance(value, (int, float)) and value is not None:
                # A measured number is only a pass when the *strain* verdict is a
                # pass; a number from a run the code refused to trust is stored
                # as not_evaluated so it can never be read as a clean result.
                ctx.db.record_assessment(
                    "compound", compound_id, "strain", key, float(value),
                    AssessmentStatus.PASS if strain_status == "pass"
                    else AssessmentStatus(strain_status if strain_status in {"flag", "fail"} else "not_evaluated"),
                    strain_code,
                )
            else:
                ctx.db.record_assessment("compound", compound_id, "strain", key, None,
                                         AssessmentStatus(strain_status) if value is not None else AssessmentStatus.NOT_EVALUATED,
                                         strain_code if value is not None else Reason.STRAIN_NOT_CONVERGED,
                                         value_text=str(value) if value is not None else None)
        if strain["status"] == "fail":
            summary["failed"] += 1
            ctx.db.record_assessment("compound", compound_id, "strain", "excessive_strain", 1.0, AssessmentStatus.FAIL, Reason.EXCESSIVE_STRAIN)
        elif strain["status"] == "flag":
            summary["flagged"] += 1
            ctx.db.record_assessment("compound", compound_id, "strain", "excessive_strain", 0.5, AssessmentStatus.FLAG, Reason.EXCESSIVE_STRAIN)
        elif strain["status"] == "pass":
            ctx.db.record_assessment("compound", compound_id, "strain", "excessive_strain", 0.0, AssessmentStatus.PASS)
        else:
            # No usable strain number: "excessive_strain = 0" would read as a
            # clean bill of health, so it is recorded as not evaluated instead.
            summary["not_evaluated"] = int(summary.get("not_evaluated", 0)) + 1
            ctx.db.record_assessment(
                "compound", compound_id, "strain", "excessive_strain", None,
                AssessmentStatus.NOT_EVALUATED, Reason.STRAIN_NOT_CONVERGED,
                note=str(strain.get("reason") or "strain was not measured"),
            )

        violations = _property_violations(properties, budget)
        ctx.db.record_assessment("compound", compound_id, "properties", "violations", violations,
                                 AssessmentStatus.PASS if not violations else AssessmentStatus.FLAG,
                                 Reason.PROPERTY_VIOLATION if violations else None)
        scaffold = rdio.murcko_scaffold(mol)
        summary["scaffolds"].setdefault(scaffold or "acyclic", []).append(compound_id)
    summary["scaffold_count"] = len(summary["scaffolds"])
    return summary


def _strain(mol, budget=None) -> dict[str, Any]:
    """Energy difference between the docked conformer and its local minimum."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    flag = float(getattr(budget, "strain_max_kcal", None) or STRAIN_FLAG)
    fail = flag * 2.0
    import math
    # Preserve every deposited heavy-atom coordinate; regenerate and relax only
    # hydrogens before comparing against an unconstrained local minimum.
    work = Chem.AddHs(Chem.RemoveHs(Chem.Mol(mol)), addCoords=True)
    if work.GetNumConformers() != 1:
        return {"status": "not_evaluated", "strain_kcal_mol": None, "reason": "one 3D conformer required"}
    for name in ("MMFF94", "UFF"):
        def force_field(m):
            if name == "MMFF94":
                if not AllChem.MMFFHasAllMoleculeParams(m): return None
                return AllChem.MMFFGetMoleculeForceField(m, AllChem.MMFFGetMoleculeProperties(m))
            if not AllChem.UFFHasAllMoleculeParams(m): return None
            return AllChem.UFFGetMoleculeForceField(m)
        field = force_field(work)
        if field is None: continue
        raw_energy = float(field.CalcEnergy())
        for atom in work.GetAtoms():
            if atom.GetAtomicNum() > 1: field.AddFixedPoint(atom.GetIdx())
        field.Initialize()
        h_status = field.Minimize(maxIts=2000)
        energy = float(field.CalcEnergy())
        reference = Chem.Mol(work)
        reference_field = force_field(reference)
        if reference_field is None: continue
        minimum_status = reference_field.Minimize(maxIts=2000)
        minimum = float(reference_field.CalcEnergy())
        result = {"force_field": name, "pose_energy": energy, "raw_pose_energy": raw_energy,
                  "minimised_energy": minimum, "hydrogen_minimization_converged": h_status == 0,
                  "reference_minimization_converged": minimum_status == 0,
                  "protocol": "heavy_atoms_fixed_H_relaxed_minus_local_minimum"}
        if h_status != 0 or minimum_status != 0 or not all(math.isfinite(x) and abs(x)<1e5 for x in (energy,minimum)):
            return {**result, "status": "not_evaluated", "strain_kcal_mol": None,
                    "reason": "minimization not converged or energy nonphysical"}
        delta = energy-minimum
        if delta < -0.01:
            return {**result, "status": "not_evaluated", "strain_kcal_mol": None,
                    "reason": "local minimum is higher than fixed-heavy pose"}
        return {**result, "strain_kcal_mol": round(delta,4),
                "status": "fail" if delta>fail else "flag" if delta>flag else "pass"}
    return {"status": "not_evaluated", "strain_kcal_mol": None, "force_field": None,
            "reason": "no force field parameterises the molecule"}


def _property_violations(properties, budget) -> list[str]:
    violations: list[str] = []
    burden = alert_burden(properties)
    limit = float(getattr(budget, "membrane_alert_burden_max", 3) or 3)
    if burden is not None and burden > limit:
        violations.append(f"alert_burden>{limit:g}")
    return violations
