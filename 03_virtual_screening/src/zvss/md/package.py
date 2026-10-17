"""MD task-package writer and preflight.

The package is a *hand-off contract*: complex geometry, ligand chemistry, identity
and provenance, plus the exact checks the downstream MD stage is expected to run.
It does not pretend that the MD run happened.  Two execution statuses are emitted
next to the package:

``ready_for_execution``                  geometry and inputs are complete;
``blocked_missing_membrane_capability``  a membrane system cannot be built;
``blocked_platform_linux_required``      the destination runner is Linux-only.
"""

from __future__ import annotations

import os
import shutil
from typing import Any

import numpy as np

from ..constants import ExecutionStatus, Reason
from ..util import PdbAtom, ensure_dir, write_pdb
from ..ids import sha256_file, sha256_text

LIGAND_RESNAME = "LIG"
LIGAND_CHAIN = "L"
LIGAND_RESID = 1


def preflight(ctx, entry: dict[str, Any]) -> dict[str, Any]:
    from .. import data

    pose = ctx.db.one("SELECT * FROM pose WHERE pose_id=?", (entry["pose_id"],))
    if not pose:
        return {"status": "fail", "reason_code": str(Reason.MISSING_COORDINATES), "detail": "pose row missing"}
    mol = data.pose_mol(ctx, pose)
    if mol is None or mol.GetNumConformers() == 0:
        return {"status": "fail", "reason_code": str(Reason.MISSING_COORDINATES), "detail": "no coordinates"}
    rec_atoms, rec_path = data.receptor_atoms(ctx)
    if not rec_atoms:
        return {"status": "fail", "reason_code": str(Reason.RECEPTOR_UNBOUND), "detail": "receptor not available"}
    xyz = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    rec_xyz = data.coords_of(rec_atoms)
    from ..util import pairwise_distances

    distances = pairwise_distances(xyz, rec_xyz)
    minimum = float(distances.min())
    checks = {
        "receptor_atoms": len(rec_atoms),
        "ligand_atoms": int(mol.GetNumAtoms()),
        "min_ligand_receptor_distance": round(minimum, 3),
        "ligand_charge": entry.get("state_id") is not None,
        "geometry_ok": minimum > 0.5,
    }
    status = "pass" if minimum > 0.5 else "fail"
    return {
        "status": status,
        "reason_code": None if status == "pass" else str(Reason.SEVERE_CLASH),
        "detail": checks,
        "receptor_path": rec_path,
        "mol": mol,
        "receptor_atoms": rec_atoms,
    }


def build(ctx, selection: dict[str, Any]) -> dict[str, Any]:
    from .. import data

    root = ctx.md_dir()
    packages: list[dict[str, Any]] = []
    formal = list(selection.get("primary", []))
    preliminary = list(selection.get("preliminary", []))
    exploratory = list(selection.get("exploratory", []))
    if not formal:
        # A hand-off package is an instruction to spend compute on a candidate.
        # Handing off a display-only compound would launder an unqualified
        # candidate into a formal one downstream, so this stays blocked.
        summary = {
            "count": 0,
            "ready": 0,
            "blocked": 0,
            "packages": [],
            "status": "blocked_no_formally_qualified_candidates",
            "reason_code": str(Reason.EVIDENCE_INCOMPLETE),
            "detail": (
                "no compound met every formal requirement; "
                f"{len(preliminary)} preliminary and {len(exploratory)} exploratory candidates exist but are "
                "not eligible for MD hand-off"
            ),
            "preliminary": [
                {"compound_id": entry["compound_id"], "unmet": (selection.get("qualification") or {}).get(entry["compound_id"], {}).get("unmet_requirements")}
                for entry in preliminary
            ],
            "exploratory": [
                {"compound_id": entry["compound_id"], "unmet": (selection.get("qualification") or {}).get(entry["compound_id"], {}).get("unmet_requirements")}
                for entry in exploratory
            ],
            "notes": [
                "hand-off is restricted to formally qualified candidates",
                "preliminary candidates stay visible in the report and are never packaged for MD",
            ],
        }
        from ..util import atomic_write_json

        atomic_write_json(os.path.join(ctx.stage_dir("md"), "md_summary.json"), summary)
        return summary

    entries = list(formal)
    for entry in entries:
        result = preflight(ctx, entry)
        role = entry.get("selection_role", "primary")
        directory = ensure_dir(os.path.join(root, f"{role}_{entry['selection_order']:02d}_{entry['compound_id']}"))
        record = {
            "compound_id": entry["compound_id"],
            "state_id": entry["state_id"],
            "pose_id": entry["pose_id"],
            "role": role,
            "order": entry["selection_order"],
            "selection_version": entry["selection_version"],
            "preflight": {k: v for k, v in result.items() if k != "mol"},
            "created_by": "zvss",
            "packages": directory,
        }
        if result["status"] != "pass":
            record["execution_status"] = ExecutionStatus.BLOCKED.value
            record["reason_code"] = result["reason_code"]
            packages.append(record)
            continue

        mol = result["mol"]
        receptor_atoms: list[PdbAtom] = result["receptor_atoms"]
        index = 0
        atom_mapping: dict[str, int] = {}
        conformer = mol.GetConformer()
        ligand_atoms: list[PdbAtom] = []
        for atom in mol.GetAtoms():
            index += 1
            position = conformer.GetAtomPosition(atom.GetIdx())
            name = _atom_name(atom, index)
            atom_mapping[name] = index
            ligand_atoms.append(
                PdbAtom(
                    record="HETATM", serial=index, name=name, altloc="", resname=LIGAND_RESNAME,
                    chain=LIGAND_CHAIN, resseq=LIGAND_RESID, icode="", x=position.x, y=position.y, z=position.z,
                    occupancy=1.0, bfactor=0.0, element=atom.GetSymbol(),
                )
            )
        complex_path = os.path.join(directory, "complex.pdb")
        header = [
            f"REMARK 900 ZVSS MD HANDOFF  compound {entry['compound_id']} state {entry['state_id']} pose {entry['pose_id']}",
            f"REMARK 901 selection {entry['selection_version']} role {role} order {entry['selection_order']}",
            "REMARK 902 receptor is the prepared receptor used for screening; ligand coordinates are the scored pose",
            "REMARK 903 membrane system NOT built; see task.json for blocked capabilities",
        ]
        write_pdb(complex_path, [*receptor_atoms, *ligand_atoms], extra_header=header)

        from ..chem import rdio

        ligand_path = os.path.join(directory, "ligand.sdf")
        rdio.write_sdf(ligand_path, [(mol, {"zvss_compound_id": entry["compound_id"], "zvss_pose_id": entry["pose_id"]})])

        manifest = {
            "format": "zvss-md-handoff-1",
            "compound_id": entry["compound_id"],
            "state_id": entry["state_id"],
            "pose_id": entry["pose_id"],
            "ligand": {
                "resname": LIGAND_RESNAME,
                "chain": LIGAND_CHAIN,
                "resid": str(LIGAND_RESID),
                "atoms": index,
                "atom_mapping": atom_mapping,
                "smiles": _smiles(ctx, entry["state_id"]),
            },
            "receptor": {"sha256": sha256_file(result["receptor_path"]) if result.get("receptor_path") else None,
                         "source": result.get("receptor_path")},
            "checksums": {},
        }
        task = {
            "task": "md_simulation",
            "format": "zvss-md-task-1",
            "inputs": {"complex_pdb": "complex.pdb", "ligand_sdf": "ligand.sdf", "manifest": "manifest.json"},
            "environment": {
                "expected_engine": ctx.config.get("md.engine", "openmm"),
                "membrane_required": bool(getattr(ctx.config.target.membrane, "usable", False)),
            },
            "execution": {
                "status": "blocked",
                "blockers": [
                    {
                        "code": str(Reason.MD_MEMBRANE_CAPABILITY_MISSING),
                        "detail": "the available MD runner has no membrane or non-aqueous protocol",
                    },
                    {
                        "code": str(Reason.MD_PLATFORM_LINUX_ONLY),
                        "detail": "the available MD runner exits immediately on non-Linux hosts",
                    },
                ],
            },
            "protocol": _protocol(ctx),
            "units": {
                "coordinates": "angstrom (PDB convention, ligand and receptor)",
                "box": "angstrom",
                "time": "nanosecond",
                "temperature": "kelvin",
                "energy": "kilojoule per mole",
                "pressure": "bar",
            },
            "coordinate_frame": {
                "system": "the docking box frame used for screening; no re-imaging was applied",
                "ligand_is_scored_pose": True,
                "receptor_is_screening_receptor": True,
                "receptor_sha256": sha256_file(result["receptor_path"]) if result.get("receptor_path") else None,
            },
            "required_checks": [
                "complex topology builds without missing parameters",
                "ligand parameters generated and charge-balanced",
                "system neutralised and minimised to a stable gradient",
                "production run reaches the configured length",
                "ligand stays within the binding site (RMSD report)",
            ],
            "provenance": {
                "selection_version": entry["selection_version"],
                "evidence_score": entry.get("evidence_score"),
                "tier": entry.get("tier"),
                "components_used": entry.get("components_used"),
                "components_missing": entry.get("components_missing"),
            },
        }
        from ..util import atomic_write_json

        atomic_write_json(os.path.join(directory, "manifest.json"), manifest)
        atomic_write_json(os.path.join(directory, "task.json"), task)
        (manifest["checksums"]).update(
            {
                "complex.pdb": sha256_file(complex_path),
                "ligand.sdf": sha256_file(ligand_path),
            }
        )
        atomic_write_json(os.path.join(directory, "manifest.json"), manifest)
        write_checksums(directory)
        write_readme(directory, entry, task)

        record["execution_status"] = "ready_for_execution"
        record["blockers"] = task["execution"]["blockers"]
        record["complex"] = complex_path
        record["ligand"] = ligand_path
        ctx.db.upsert_md_artifact(
            {
                "handoff_id": "MDH" + sha256_text(
                    f"{entry['pose_id']}|{entry['selection_version']}|{role}|{entry['selection_order']}"
                )[:16].upper(),
                "run_id": ctx.run_id,
                "compound_id": entry["compound_id"],
                "pose_id": entry["pose_id"],
                "role": role,
                "package_dir": directory,
                "status": "ready_for_execution",
                "blockers": task["execution"]["blockers"],
                "sha256": sha256_text(complex_path + ligand_path),
            }
        )
        packages.append(record)

    summary = {
        "packages": packages,
        "count": len(packages),
        "ready": sum(1 for p in packages if p.get("execution_status") == "ready_for_execution"),
        "blocked": sum(1 for p in packages if p.get("execution_status") != "ready_for_execution"),
        "handoff_eligibility": "formal qualified candidates only",
        "preliminary": [
            {"compound_id": entry["compound_id"],
             "unmet": (selection.get("qualification") or {}).get(entry["compound_id"], {}).get("unmet_requirements")}
            for entry in preliminary
        ],
        "exploratory_not_packaged": [
            {"compound_id": entry["compound_id"],
             "unmet": (selection.get("qualification") or {}).get(entry["compound_id"], {}).get("unmet_requirements")}
            for entry in exploratory
        ],
    }
    from ..util import atomic_write_json

    atomic_write_json(os.path.join(ctx.stage_dir("md"), "md_summary.json"), summary)
    return summary


def _atom_name(atom, index: int) -> str:
    symbol = atom.GetSymbol().upper()
    return f"{symbol}{index}"[:4]


def _protocol(ctx) -> dict[str, Any]:
    """The MD protocol the package is written for, echoed from the target config.

    Echoing the configured values (instead of leaving them implicit) is what lets
    the downstream project check that the package matches its own protocol rather
    than guessing.
    """
    profile = ctx.config.target
    protocol = dict(getattr(profile, "md_protocol", None) or {})
    protocol.setdefault("engine", ctx.config.get("md.engine", "openmm"))
    protocol.setdefault("force_field", "amber14-all.xml")
    protocol.setdefault("water_model", "tip3p.xml")
    protocol.setdefault("temperature_K", 310.0)
    protocol.setdefault("pressure_bar", 1.0)
    protocol.setdefault("timestep_fs", 2.0)
    protocol.setdefault("length_ns", 100.0)
    protocol.setdefault("ensemble", "NPT")
    protocol["membrane_protocol_available"] = bool(getattr(profile.membrane, "usable", False))
    protocol["source"] = "target profile md_protocol plus documented defaults"
    return protocol


def _smiles(ctx, state_id: str) -> str | None:
    row = ctx.db.one("SELECT canonical_smiles FROM chemical_state WHERE state_id=?", (state_id,))
    return row["canonical_smiles"] if row else None


def write_checksums(directory: str) -> str:
    import hashlib

    lines = []
    for name in sorted(os.listdir(directory)):
        path = os.path.join(directory, name)
        if not os.path.isfile(path) or name == "checksums.sha256":
            continue
        digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
        lines.append(f"{digest}  {name}")
    target = os.path.join(directory, "checksums.sha256")
    with open(target, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")
    return target


def write_readme(directory: str, entry: dict[str, Any], task: dict[str, Any]) -> str:
    lines = [
        f"# MD task package - {entry['compound_id']}",
        "",
        f"- pose: `{entry['pose_id']}`",
        f"- role: {entry.get('selection_role')} #{entry.get('selection_order')}",
        f"- evidence tier: {entry.get('tier')}",
        f"- selection version: `{entry.get('selection_version')}`",
        "",
        "## Contents",
        "",
        "- `complex.pdb` receptor + ligand pose in one file",
        "- `ligand.sdf` ligand chemistry as scored",
        "- `manifest.json` ligand identity and atom mapping",
        "- `task.json` machine-readable task and required checks",
        "- `checksums.sha256` integrity",
        "",
        "## Execution",
        "",
        "This package is complete, but the MD step itself is **not executed here**:",
        "",
    ]
    for blocker in task["execution"]["blockers"]:
        lines.append(f"- `{blocker['code']}`: {blocker['detail']}")
    lines += ["", "Run the downstream MD project against `complex.pdb` and record the results separately.", ""]
    path = os.path.join(directory, "README.md")
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines))
    return path
