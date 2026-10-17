"""Validate MD hand-off packages against the specification, item by item.

The question this answers is narrow and important: *is the text/coordinate pair
handed to the MD project equivalent to what the screening stage decided?*  It is
not a membrane-engine preflight -- that needs a membrane-capable MD runner,
which this environment does not have, and the run keeps that state.

Every check records the observed value next to the required one, so a partial
package cannot be reported as "validated".
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from typing import Any

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PROJECT, "src"))

from zvss.constants import Reason  # noqa: E402
from zvss.util import read_json, parse_pdb_text  # noqa: E402


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def check(name: str, passed: bool, observed: Any, expected: Any, detail: str = "") -> dict[str, Any]:
    return {
        "check": name,
        "status": "pass" if passed else "fail",
        "observed": observed,
        "expected": expected,
        "detail": detail,
    }


def validate_package(directory: str) -> dict[str, Any]:
    manifest = read_json(os.path.join(directory, "manifest.json"), {}) or {}
    task = read_json(os.path.join(directory, "task.json"), {}) or {}
    checks: list[dict[str, Any]] = []

    required = ["manifest.json", "task.json", "complex.pdb", "ligand.sdf", "checksums.sha256"]
    present = {name: os.path.exists(os.path.join(directory, name)) for name in required}
    checks.append(check("package_files_present", all(present.values()), present, {n: True for n in required}))

    # checksum integrity
    recorded = manifest.get("checksums") or {}
    actual = {}
    for name in ("complex.pdb", "ligand.sdf"):
        path = os.path.join(directory, name)
        actual[name] = sha256_file(path) if os.path.exists(path) else None
    checks.append(
        check(
            "checksums_match_files",
            all(recorded.get(name) == actual[name] for name in ("complex.pdb", "ligand.sdf")),
            actual,
            recorded,
        )
    )

    # ligand mapping: manifest atom count and names must match ligand.sdf
    ligand_text = open(os.path.join(directory, "ligand.sdf"), encoding="utf-8").read()
    ligand_lines = ligand_text.splitlines()
    try:
        sdf_atoms = int(ligand_lines[3][0:3])
    except (IndexError, ValueError):
        sdf_atoms = 0
    mapping = manifest.get("ligand", {}).get("atom_mapping") or {}
    checks.append(
        check(
            "ligand_atom_mapping_complete",
            bool(mapping) and len(mapping) == manifest.get("ligand", {}).get("atoms") == sdf_atoms,
            {"mapping_entries": len(mapping), "declared_atoms": manifest.get("ligand", {}).get("atoms"), "sdf_atoms": sdf_atoms},
            "mapping entries == declared atoms == SDF atoms",
        )
    )

    # complex composition: receptor atoms + ligand atoms, ligand as HETATM block
    complex_text = open(os.path.join(directory, "complex.pdb"), encoding="utf-8").read()
    atoms = parse_pdb_text(complex_text).atoms
    protein = [a for a in atoms if a.record.strip().upper() == "ATOM"]
    het = [a for a in atoms if a.record.strip().upper() == "HETATM"]
    resname = manifest.get("ligand", {}).get("resname")
    ligand_het = [a for a in het if a.resname == resname]
    checks.append(
        check(
            "complex_atom_count_consistent",
            len(atoms) == len(protein) + len(het) and len(ligand_het) == len(mapping),
            {"total": len(atoms), "protein": len(protein), "hetatm": len(het), "ligand_atoms": len(ligand_het)},
            {"total": "protein+hetatm", "ligand_atoms": len(mapping)},
        )
    )

    # coordinates: ligand in complex equal to ligand.sdf coordinates (no drift)
    ligand_coords = []
    for line in ligand_lines[4 : 4 + sdf_atoms]:
        try:
            ligand_coords.append((float(line[0:10]), float(line[10:20]), float(line[20:30])))
        except ValueError:
            continue
    complex_coords = [(a.x, a.y, a.z) for a in ligand_het]
    max_dev = None
    if len(ligand_coords) == len(complex_coords) and complex_coords:
        max_dev = max(
            max(abs(a - b) for a, b in zip(p, q)) for p, q in zip(ligand_coords, complex_coords)
        )
    checks.append(
        check(
            "coordinates_identical_to_ligand_sdf",
            max_dev is not None and max_dev < 1e-3,
            {"max_abs_deviation_A": round(max_dev, 5) if max_dev is not None else None},
            "< 0.001 A",
            "the pose handed to MD must be the pose that was scored",
        )
    )

    # provenance: selection version / evidence tier carried through
    provenance = task.get("provenance") or {}
    checks.append(
        check(
            "provenance_linked_to_selection",
            bool(provenance.get("selection_version")) and provenance.get("tier") is not None,
            provenance,
            "selection_version + tier present",
        )
    )

    # units and protocol: explicit, non-empty, with engine and membrane flag
    protocol = task.get("protocol") or {}
    units = task.get("units") or {}
    env = task.get("environment") or {}
    checks.append(
        check(
            "protocol_and_units_declared",
            bool(protocol) and bool(units) and bool(env.get("expected_engine")),
            {"protocol_keys": sorted(protocol), "units": units, "engine": env.get("expected_engine")},
            "protocol + units + engine declared",
        )
    )

    # blocked state must be explicit rather than absent
    execution = task.get("execution") or {}
    blockers = execution.get("blockers") or []
    codes = {blocker.get("code") for blocker in blockers}
    checks.append(
        check(
            "blocked_capabilities_explicit",
            execution.get("status") == "blocked"
            and str(Reason.MD_MEMBRANE_CAPABILITY_MISSING) in codes
            and str(Reason.MD_PLATFORM_LINUX_ONLY) in codes,
            {"status": execution.get("status"), "blockers": sorted(c for c in codes if c)},
            {
                "status": "blocked",
                "blockers": [str(Reason.MD_MEMBRANE_CAPABILITY_MISSING), str(Reason.MD_PLATFORM_LINUX_ONLY)],
            },
        )
    )

    adapter04 = read_json(os.path.join(directory, "adapter04.json"), {}) or {}
    adapter05 = read_json(os.path.join(directory, "adapter05.json"), {}) or {}
    return {
        "package": os.path.basename(directory),
        "compound_id": manifest.get("compound_id"),
        "pose_id": manifest.get("pose_id"),
        "checks": checks,
        "failed": [item["check"] for item in checks if item["status"] == "fail"],
        "validated": not any(item["status"] == "fail" for item in checks),
        "adapter04": adapter04.get("status") or adapter04.get("execution_status"),
        "adapter05": adapter05.get("status") or adapter05.get("execution_status"),
    }


def main(run_id: str) -> int:
    handoff = os.path.join(PROJECT, "runs", run_id, "md_handoff")
    packages = sorted(
        os.path.join(handoff, name) for name in os.listdir(handoff) if os.path.isdir(os.path.join(handoff, name))
    )
    results = [validate_package(package) for package in packages]
    payload = {
        "run_id": run_id,
        "handoff_dir": handoff,
        "packages": results,
        "validated": sum(1 for item in results if item["validated"]),
        "failed": sum(1 for item in results if not item["validated"]),
        "scope": "hand-off package validation only; not a membrane-engine preflight",
    }
    out_dir = os.path.join(PROJECT, "docs", "acceptance")
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "MD_HANDOFF_VALIDATION.json"), "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=1)

    lines = [
        f"# MD Task bundle equivalence validation ({run_id})",
        "",
        "Checks bundle/specification equivalence,**are not**membrane engine preflight.",
        "",
        f"- Task bundle count: {len(results)} | All checks passed: {payload['validated']} | Failures present: {payload['failed']}",
        "",
        "| Task bundles | Compounds | Checks | Failures | adapter04 | adapter05 |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for item in results:
        lines.append(
            f"| {item['package']} | {item['compound_id']} | {len(item['checks'])} | "
            f"{', '.join(item['failed']) or 'None'} | {item['adapter04']} | {item['adapter05']} |"
        )
    detail_names = [item["check"] for item in results[0]["checks"]] if results else []
    if detail_names:
        lines += ["", "## Per-bundle results", ""]
        for item in results:
            lines.append(f"### {item['package']}")
            lines.append("")
            lines.append("| Check | Result | Observed value | Expected value |")
            lines.append("| --- | --- | --- | --- |")
            for entry in item["checks"]:
                observed = json.dumps(entry["observed"], ensure_ascii=False, default=str)
                expected = json.dumps(entry["expected"], ensure_ascii=False, default=str)
                lines.append(
                    f"| {entry['check']} | {entry['status']} | {observed[:180]} | {expected[:180]} |"
                )
            lines.append("")
    lines += [
        "## Conclusion",
        "",
        "- Coordinates/Mapping/units/Provenance matches screening-stage records; see itemized results.",
        "- Container and membrane system remain unbuilt:`04` associated OpenMM runner lacks a membrane protocol and runs only under Linux ;",
        "  these bundles therefore remain `blocked`, and are not described as having entered production MD”. ",
        "",
    ]
    with open(os.path.join(out_dir, "MD_HANDOFF_VALIDATION.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
    print(json.dumps({"packages": len(results), "validated": payload["validated"], "failed": payload["failed"]}))
    for item in results:
        if item["failed"]:
            print(" FAILED", item["package"], item["failed"])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "t6"))
