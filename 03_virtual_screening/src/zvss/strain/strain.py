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
    by_compound: dict[str, dict[str, Any]] = {}
    for pose in poses:
        by_compound.setdefault(pose["compound_id"], pose)

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
                                         AssessmentStatus.NOT_EVALUATED, Reason.STRAIN_NOT_CONVERGED,
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
    work = Chem.Mol(mol)
    if work.GetNumAtoms() == work.GetNumHeavyAtoms():
        work = Chem.AddHs(work)
    energy = None
    minimum = None
    used = None
    for name in ("MMFF94", "UFF"):
        # The force-field object and its parameters must come from the same
        # molecule instance; mixing two copies crashes RDKit natively.
        reference = Chem.Mol(work)
        if name == "MMFF94":
            if not AllChem.MMFFHasAllMoleculeParams(work):
                continue
            props = AllChem.MMFFGetMoleculeProperties(work)
            ref_props = AllChem.MMFFGetMoleculeProperties(reference)
            if props is None or ref_props is None:
                continue
            field = AllChem.MMFFGetMoleculeForceField(work, props)
            reference_field = AllChem.MMFFGetMoleculeForceField(reference, ref_props)
        else:
            if not AllChem.UFFHasAllMoleculeParams(work):
                continue
            field = AllChem.UFFGetMoleculeForceField(work)
            reference_field = AllChem.UFFGetMoleculeForceField(reference)
        if field is None or reference_field is None:
            continue
        energy = float(field.CalcEnergy())
        reference_field.Minimize(maxIts=2000)
        minimum = float(reference_field.CalcEnergy())
        used = name
        break
    if energy is None or minimum is None:
        return {
            "strain_kcal_mol": None,
            "force_field": None,
            "status": "not_evaluated",
            "reason": "no force field in this build parameterises the molecule",
        }
    if abs(energy) > 1.0e5 or abs(minimum) > 1.0e5:
        return {
            "strain_kcal_mol": None,
            "force_field": used,
            "status": "not_evaluated",
            "reason": "single-point energy is unphysical; coordinates are probably degenerate",
        }
    delta = energy - minimum
    status = "pass"
    if delta > fail:
        status = "fail"
    elif delta > flag:
        status = "flag"
    return {
        "strain_kcal_mol": round(float(delta), 4),
        "pose_energy": round(float(energy), 4),
        "minimised_energy": round(float(minimum), 4),
        "force_field": used,
        "status": status,
    }


def _property_violations(properties, budget) -> list[str]:
    violations: list[str] = []
    burden = alert_burden(properties)
    limit = float(getattr(budget, "membrane_alert_burden_max", 3) or 3)
    if burden is not None and burden > limit:
        violations.append(f"alert_burden>{limit:g}")
    return violations
