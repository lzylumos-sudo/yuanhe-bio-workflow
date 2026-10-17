"""Adapters for the two downstream projects.

``adapter04`` maps a ZVSS MD hand-off package onto the existing OpenMM MD project
(`05_molecular_dynamics`), which accepts a ``complex.pdb`` plus an adjacent
``ligand.sdf``/``manifest.json``.  ``adapter05`` maps the same reference onto the
ligand-interaction project (`04_interaction_analysis`, ZLIS) command line.

Both adapters *plan* commands and validate inputs.  They never claim a run
happened: the returned record carries ``executed: false`` and the concrete
blockers (Linux-only runner, missing membrane protocol).
"""

from __future__ import annotations

import os
import shlex
from typing import Any

from ..constants import Reason
from ..util import atomic_write_json, read_json


def adapter04(ctx, package_dir: str, *, project_root: str | None = None) -> dict[str, Any]:
    root = project_root or (ctx.config.get("md.adapter04_root") or "")
    complex_pdb = os.path.join(package_dir, "complex.pdb")
    ligand_sdf = os.path.join(package_dir, "ligand.sdf")
    manifest = os.path.join(package_dir, "manifest.json")
    checks = {
        "complex_pdb": os.path.exists(complex_pdb),
        "ligand_sdf": os.path.exists(ligand_sdf),
        "manifest": os.path.exists(manifest),
        "project_root_present": bool(root) and os.path.isdir(root),
    }
    entry_script = os.path.join(root, "run_md.sh") if root else ""
    command = ["bash", "run_md.sh", "--input", "complex.pdb"] if root else []
    record = {
        "adapter": "04_universal_md",
        "package": package_dir,
        "checks": checks,
        "command_template": command,
        "command_readable": " ".join(shlex.quote(part) for part in command),
        "contract": {
            "required_input": "complex.pdb",
            "optional": ["ligand.sdf", "manifest.json"],
            "manifest_fields_used": ["ligand.resname", "ligand.chain", "ligand.resid", "ligand.atom_mapping"],
            "embedded_remark_format": "REMARK 950",
            "note": "the existing runner also accepts an embedded manifest; ZVSS writes a sidecar manifest instead",
        },
        "executed": False,
        "blockers": [
            {
                "code": str(Reason.MD_PLATFORM_LINUX_ONLY),
                "detail": "run_md.sh exits with code 2 on non-Linux hosts; this environment is Windows",
            },
            {
                "code": str(Reason.MD_MEMBRANE_CAPABILITY_MISSING),
                "detail": "the existing MD project declares no membrane or non-aqueous protocol",
            },
        ],
        "verified_offline": checks["complex_pdb"] and checks["ligand_sdf"] and checks["manifest"],
    }
    atomic_write_json(os.path.join(package_dir, "adapter04.json"), record)
    return record


def adapter05(ctx, package_dir: str, *, project_root: str | None = None) -> dict[str, Any]:
    root = project_root or (ctx.config.get("analysis.adapter05_root") or "")
    manifest = read_json(os.path.join(package_dir, "manifest.json"), {}) or {}
    ligand = manifest.get("ligand") or {}
    selector = f"{ligand.get('chain', 'B')}:{ligand.get('resname', 'LIG')}{ligand.get('resid', '1')}"
    command = [
        "bash", "run.sh", "analyze", "complex.pdb",
        "--ligand", selector,
        "--template", "ligand.sdf",
    ]
    record = {
        "adapter": "05_ligand_interaction_studio",
        "package": package_dir,
        "command_template": command,
        "command_readable": " ".join(shlex.quote(part) for part in command),
        "selector": selector,
        "contract": {
            "entry": "bash run.sh analyze <pdb> [--ligand A:BEN1] [--template lig.sdf]",
            "engine": "ProLIF",
            "outputs": ["interactions.csv", "interactions.json", "2D SVG", "scene.pse", "ligand.sdf"],
        },
        "project_root_present": bool(root) and os.path.isdir(root),
        "executed": False,
        "blockers": [
            {
                "code": str(Reason.INTERACTION_ENGINE_FAILED),
                "detail": "the interaction project requires its own Python environment; ZVSS plans the command but does not run it",
            }
        ] if not root else [],
    }
    atomic_write_json(os.path.join(package_dir, "adapter05.json"), record)
    return record


def build_adapters(ctx, md_summary: dict[str, Any]) -> dict[str, Any]:
    results = []
    for package in md_summary.get("packages") or []:
        directory = package.get("packages")
        if not directory or not os.path.isdir(directory):
            continue
        if package.get("execution_status") != "ready_for_execution":
            results.append({"package": directory, "skipped": package.get("reason_code")})
            continue
        results.append(
            {
                "package": directory,
                "adapter04": adapter04(ctx, directory),
                "adapter05": adapter05(ctx, directory),
            }
        )
    summary = {"packages": len(results), "records": results}
    atomic_write_json(os.path.join(ctx.stage_dir("md"), "adapters.json"), summary)
    return summary
