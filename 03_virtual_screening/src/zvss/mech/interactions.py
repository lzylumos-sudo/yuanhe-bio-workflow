"""Local mechanism evidence: which receptor residues does a pose actually touch.

This stage answers a narrow, checkable question -- "does the pose sit in the
configured local region, and with what geometry" -- and deliberately does not
claim a binding mechanism.  Reference-ligand comparison is a contact-set
overlap, not an energy.

Numbering: the target profile speaks native/UniProt numbering; the receptor file
speaks author numbering.  The conversion comes from the verified residue map, and
a region residue that is absent from the receptor is reported as such.
"""

from __future__ import annotations

import os
from typing import Any, Iterable

import numpy as np

from ..constants import AssessmentStatus, Reason
from ..util import read_json

CONTACT_CUTOFF = 4.5
HBOND_CUTOFF = 3.5


def _mapping_rows(ctx) -> list[dict[str, Any]]:
    directory = ctx.ref_dir()
    if not os.path.isdir(directory):
        return []
    for name in sorted(os.listdir(directory)):
        if name.startswith("mapping_") and name.endswith(".json"):
            payload = read_json(os.path.join(directory, name), {}) or {}
            if payload.get("rows"):
                return payload["rows"]
    return []


def native_to_auth(ctx, native_residues: Iterable[int]) -> dict[int, dict[str, Any]]:
    rows = _mapping_rows(ctx)
    by_native: dict[int, dict[str, Any]] = {}
    for row in rows:
        native = row.get("uniprot_resseq")
        if native is None:
            continue
        by_native[int(native)] = {
            "auth_seq_id": row.get("auth_seq_id"),
            "resname": row.get("resname"),
            "construct_match": row.get("construct_match"),
            "label_seq_id": row.get("label_seq_id"),
        }
    return {int(n): by_native.get(int(n), {}) for n in native_residues}


def region_atoms(ctx, receptor_atoms, native_residues: Iterable[int]) -> dict[str, Any]:
    mapping = native_to_auth(ctx, native_residues)
    present: dict[int, list] = {}
    missing: list[int] = []
    mismatched: list[dict[str, Any]] = []
    for native, info in mapping.items():
        auth = info.get("auth_seq_id")
        if auth is None:
            missing.append(native)
            continue
        atoms = [a for a in receptor_atoms if a.resseq == auth and not a.is_hydrogen()]
        if not atoms:
            missing.append(native)
            continue
        if info.get("resname") and atoms[0].resname != info["resname"]:
            mismatched.append({"native": native, "auth": auth, "expected": info.get("resname"), "found": atoms[0].resname})
        present[native] = atoms
    return {
        "residues": present,
        "missing": missing,
        "mismatched": mismatched,
        "mapping_source": "verified_residue_map" if mapping else "none",
        "resolved_auth": {n: (v[0].resseq if v else None) for n, v in present.items()},
    }


def contacts(ligand_xyz: np.ndarray, receptor_atoms, *, cutoff: float = CONTACT_CUTOFF) -> dict[str, Any]:
    """Per-residue contacts and a hydrophobic/polar split."""
    per_residue: dict[tuple[str, int, str], dict[str, Any]] = {}
    for atom in receptor_atoms:
        if atom.is_hydrogen():
            continue
        key = (atom.chain, atom.resseq, atom.resname)
        entry = per_residue.setdefault(
            key, {"chain": atom.chain, "resseq": atom.resseq, "resname": atom.resname, "min_distance": 1e9, "atoms": 0, "polar": False}
        )
        delta = ligand_xyz - np.array([atom.x, atom.y, atom.z])
        distance = float(np.sqrt((delta ** 2).sum(axis=1)).min())
        entry["min_distance"] = min(entry["min_distance"], distance)
        if distance <= cutoff:
            entry["atoms"] += 1
            if atom.element.upper() in {"N", "O", "S"}:
                entry["polar"] = True
    touching = [entry for entry in per_residue.values() if entry["min_distance"] <= cutoff]
    touching.sort(key=lambda e: (e["min_distance"], e["chain"], e["resseq"]))
    for entry in touching:
        entry["min_distance"] = round(entry["min_distance"], 3)
    return {
        "contacts": touching,
        "contact_count": len(touching),
        "polar_contacts": sum(1 for e in touching if e["polar"]),
        "hydrophobic_contacts": sum(1 for e in touching if not e["polar"]),
        "hbond_proxy": sum(
            1 for e in touching if e["polar"] and e["min_distance"] <= HBOND_CUTOFF
        ),
        "fingerprint": sorted(f"{e['resname']}{e['resseq']}" for e in touching),
    }


def run(ctx) -> dict[str, Any]:
    from .. import data

    poses = data.pose_rows(ctx)
    profile = ctx.config.target
    region_cfg = getattr(profile, "region", None)
    sides = list(getattr(region_cfg, "sides", []) or [])
    receptor_atoms, receptor_path = data.receptor_atoms(ctx)
    summary: dict[str, Any] = {
        "poses": len(poses),
        "receptor": receptor_path,
        "cutoff": CONTACT_CUTOFF,
        "regions": {},
        "evaluated": 0,
        "not_evaluated": 0,
        "reference": {},
    }
    if not receptor_atoms:
        for pose in poses:
            ctx.db.record_assessment("pose", pose["pose_id"], "mechanism", "region_contacts", None,
                                     AssessmentStatus.NEEDS_INPUT, Reason.RECEPTOR_UNBOUND)
        summary["status"] = "needs_input"
        summary["reason"] = str(Reason.RECEPTOR_UNBOUND)
        return summary

    region_residues: dict[str, list[int]] = {}
    for side in sides:
        region_residues[side.name] = list(side.native_residues or [])
    resolved = {name: region_atoms(ctx, receptor_atoms, residues) for name, residues in region_residues.items()}
    summary["regions"] = {
        name: {
            "native": region_residues[name],
            "resolved_auth": resolved[name]["resolved_auth"],
            "missing": resolved[name]["missing"],
            "mismatched": resolved[name]["mismatched"],
            "mapping_source": resolved[name]["mapping_source"],
        }
        for name in region_residues
    }

    region_xyz: dict[str, np.ndarray] = {}
    for name, payload in resolved.items():
        atoms = [a for group in payload["residues"].values() for a in group]
        region_xyz[name] = data.coords_of(atoms) if atoms else np.zeros((0, 3))

    reference = _reference_contacts(ctx, receptor_atoms)
    summary["reference"] = reference

    for pose in poses:
        mol = data.pose_mol(ctx, pose)
        if mol is None or mol.GetNumConformers() == 0:
            ctx.db.record_assessment("pose", pose["pose_id"], "mechanism", "contacts", None,
                                     AssessmentStatus.FAIL, Reason.MISSING_COORDINATES)
            summary["not_evaluated"] += 1
            continue
        xyz = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
        contact_payload = contacts(xyz, receptor_atoms)
        ctx.db.record_assessment("pose", pose["pose_id"], "mechanism", "contact_count", contact_payload["contact_count"], AssessmentStatus.PASS)
        ctx.db.record_assessment("pose", pose["pose_id"], "mechanism", "polar_contacts", contact_payload["polar_contacts"], AssessmentStatus.PASS)
        ctx.db.record_assessment("pose", pose["pose_id"], "mechanism", "hbond_proxy", contact_payload["hbond_proxy"], AssessmentStatus.PASS)
        ctx.db.record_assessment("pose", pose["pose_id"], "mechanism", "contact_fingerprint", contact_payload["fingerprint"], AssessmentStatus.PASS)

        hits: dict[str, Any] = {}
        total_region_contacts = 0
        for name, arr in region_xyz.items():
            if arr.size == 0:
                ctx.db.record_assessment("pose", pose["pose_id"], "mechanism", f"region_{name}_contacts", None,
                                         AssessmentStatus.NEEDS_INPUT, Reason.RESIDUE_MAPPING_MISSING)
                hits[name] = {"status": "needs_input"}
                continue
            distances = np.sqrt(((xyz[:, None, :] - arr[None, :, :]) ** 2).sum(axis=2))
            count = int((distances.min(axis=0) <= CONTACT_CUTOFF).sum())
            total_region_contacts += count
            hits[name] = {"contacts": count, "min_distance": round(float(distances.min()), 3)}
            ctx.db.record_assessment("pose", pose["pose_id"], "mechanism", f"region_{name}_contacts", count, AssessmentStatus.PASS)
            ctx.db.record_assessment("pose", pose["pose_id"], "mechanism", f"region_{name}_min_distance", round(float(distances.min()), 3), AssessmentStatus.PASS)
        ctx.db.record_assessment("pose", pose["pose_id"], "mechanism", "region_contacts_total", total_region_contacts, AssessmentStatus.PASS)

        if reference.get("fingerprint"):
            overlap = _jaccard(contact_payload["fingerprint"], reference["fingerprint"])
            ctx.db.record_assessment("pose", pose["pose_id"], "mechanism", "reference_contact_jaccard", round(overlap, 4),
                                     AssessmentStatus.PASS if overlap > 0 else AssessmentStatus.FLAG)
        else:
            ctx.db.record_assessment("pose", pose["pose_id"], "mechanism", "reference_contact_jaccard", None,
                                     AssessmentStatus.NOT_EVALUATED, Reason.REFERENCE_POSE_MISSING)
        summary["evaluated"] += 1
        summary.setdefault("per_pose", {})[pose["pose_id"]] = {
            "contacts": contact_payload["contact_count"],
            "region_contacts": total_region_contacts,
            "regions": hits,
        }
    return summary


def _reference_contacts(ctx, receptor_atoms) -> dict[str, Any]:
    """Contacts of the deposited reference ligands, as an overlap baseline."""
    from .. import data

    ref_atoms, record = data.reference_structure(ctx)
    ligand = record.get("reference_ligand") or {}
    if not ligand.get("atoms"):
        return {
            "status": "not_evaluated",
            "reason": str(Reason.REFERENCE_POSE_MISSING),
            "note": "no reference ligand pose is configured for this target",
        }
    xyz = np.asarray(ligand["atoms"], dtype=float)
    payload = contacts(xyz, receptor_atoms)
    return {
        "status": "ok",
        "label": ligand.get("label"),
        "source": ligand.get("source"),
        "fingerprint": payload["fingerprint"],
        "contact_count": payload["contact_count"],
    }


def _jaccard(a: list[str], b: list[str]) -> float:
    sa, sb = set(a), set(b)
    if not sa or not sb:
        return 0.0
    return len(sa & sb) / len(sa | sb)
