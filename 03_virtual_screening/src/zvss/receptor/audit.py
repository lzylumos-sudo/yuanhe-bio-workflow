"""Receptor audit: prove which file is actually used, and what was removed.

The audit answers four questions that a bare "receptor: 4598 atoms" cannot:

1. what the deposited entry contains per chain, split into fusion / native
   domain / C-terminal tag using the verified UniProt mapping;
2. what the file handed to later stages actually contains, verified by
   re-reading it and comparing residue identity with the source;
3. which residues and atoms were deleted, and why;
4. whether the deposited mutations were left as crystallised, restored to the
   wild type, or not resolved at all.

Nothing here assumes a specific atom count: the counts are measured, and a
disagreement with the deposit is reported instead of being normalised away.
"""

from __future__ import annotations

import collections
import csv
import os
from typing import Any, Iterable

from ..constants import AssessmentStatus
from ..ids import sha256_file
from ..util import atomic_write_json, atomic_write_text
from .pdbio import load_structure, protein_atoms

ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
    "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
    "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


def _residue_map(structure, chain: str) -> dict[int, dict[str, Any]]:
    out: dict[int, dict[str, Any]] = {}
    for atom in protein_atoms(structure):
        if atom.chain != chain:
            continue
        slot = out.setdefault(int(atom.resseq), {"resseq": int(atom.resseq), "resname": atom.resname, "atoms": 0})
        slot["atoms"] += 1
        # heavy-atom-only name inventory is enough to spot a rebuilt side chain
        slot.setdefault("atom_names", set()).add(atom.name.strip())
    return out


def _chain_inventory(structure) -> list[dict[str, Any]]:
    groups: dict[str, dict[str, Any]] = collections.OrderedDict()
    for atom in structure.atoms:
        chain = atom.chain or " "
        slot = groups.setdefault(
            chain,
            {"chain": chain, "atoms": 0, "protein_atoms": 0, "hetero_atoms": 0,
             "residues": set(), "hetero": collections.Counter()},
        )
        slot["atoms"] += 1
        if atom.record.strip().upper() == "ATOM" and not atom.is_hydrogen():
            slot["protein_atoms"] += 1
            slot["residues"].add(int(atom.resseq))
        elif atom.record.strip().upper() == "HETATM":
            slot["hetero_atoms"] += 1
            slot["hetero"][atom.resname] += 1
    out = []
    for slot in groups.values():
        nums = sorted(slot["residues"])
        out.append(
            {
                "chain": slot["chain"],
                "atoms": slot["atoms"],
                "protein_atoms": slot["protein_atoms"],
                "hetero_atoms": slot["hetero_atoms"],
                "protein_residues": len(nums),
                "auth_first": nums[0] if nums else None,
                "auth_last": nums[-1] if nums else None,
                "hetero": dict(slot["hetero"]),
            }
        )
    return out


def _segment_of(auth_seq_id: int, native_range: tuple[int, int] | None) -> str:
    if native_range is None:
        return "unknown"
    low, high = native_range
    if auth_seq_id < low:
        return "n_terminal_fusion"
    if auth_seq_id > high:
        return "c_terminal_tag"
    return "native_domain"


def audit_receptor(
    *,
    original_path: str,
    receptor_path: str,
    reference_complex_path: str,
    chain: str,
    native_range: tuple[int, int] | None,
    mapping_rows: Iterable[dict[str, Any]] = (),
    mutations: Iterable[Any] = (),
    c_terminal_tag: dict[str, Any] | None = None,
    strip_log: dict[str, Any] | None = None,
) -> dict[str, Any]:
    before = load_structure(original_path)
    after = load_structure(receptor_path)
    before_res = _residue_map(before, chain)
    after_res = _residue_map(after, chain)

    native_of: dict[int, int | None] = {}
    native_resname: dict[int, str] = {}
    for row in mapping_rows:
        auth = row.get("auth_seq_id")
        if auth is None:
            continue
        native_of[int(auth)] = row.get("uniprot_resseq")
        native_resname[int(auth)] = str(row.get("uniprot_resname") or "")

    segments: dict[str, dict[str, Any]] = {}
    for resseq, info in before_res.items():
        name = _segment_of(resseq, native_range)
        slot = segments.setdefault(name, {"residues": 0, "atoms": 0, "first": resseq, "last": resseq})
        slot["residues"] += 1
        slot["atoms"] += info["atoms"]
        slot["first"] = min(slot["first"], resseq)
        slot["last"] = max(slot["last"], resseq)

    # Numbering must round-trip: every residue that survives into the receptor
    # keeps its deposited author number and residue name, and no new number
    # appears. Deliberate removals are accounted for separately in ``deletions``;
    # counting them here would hide a real renumbering behind expected stripping.
    dropped = sorted(set(before_res) - set(after_res))
    added = sorted(set(after_res) - set(before_res))
    renamed = [
        {"resseq": resseq, "deposit": before_res[resseq]["resname"], "receptor": after_res[resseq]["resname"]}
        for resseq in sorted(set(before_res) & set(after_res))
        if before_res[resseq]["resname"] != after_res[resseq]["resname"]
    ]
    numbering_ok = not added and not renamed
    expected_native = {
        resseq for resseq in before_res if native_range and native_range[0] <= resseq <= native_range[1]
    }
    foreign_in_receptor = sorted(set(after_res) - expected_native) if native_range else []
    missing_from_receptor = sorted(expected_native - set(after_res))

    mutation_records = []
    for mutation in mutations:
        auth = int(getattr(mutation, "auth_resseq"))
        declared_native = str(getattr(mutation, "native_resname", "") or "").upper()
        declared_construct = str(getattr(mutation, "construct_resname", "") or "").upper()
        observed = before_res.get(auth, {}).get("resname")
        observed_after = after_res.get(auth, {}).get("resname")
        mapped = native_resname.get(auth) or ONE_TO_THREE.get(str(native_of.get(auth) and ""), "")
        mutation_records.append(
            {
                "auth_seq_id": auth,
                "native_resseq": native_of.get(auth),
                "declared_native_resname": declared_native,
                "declared_construct_resname": declared_construct,
                "mapping_native_resname": mapped,
                "observed_resname_before": observed,
                "observed_resname_after": observed_after,
                "wt_restored": False,
                "status": "as_deposited_mutant",
                "reason": (
                    "no rotamer rebuilding is performed; the receptor keeps the crystallised "
                    "side chain so every downstream number refers to this exact structure"
                ),
                "consistent_with_config": observed == declared_construct if observed else None,
            }
        )

    checks = [
        {
            "check": "receptor_numbering_roundtrip",
            "status": AssessmentStatus.PASS.value if numbering_ok else AssessmentStatus.FAIL.value,
            "detail": (
                f"all {len(set(before_res) & set(after_res))} retained residues keep their deposited number and name; "
                f"{len(dropped)} residues were removed by design"
                if numbering_ok
                else f"{len(added)} unexpected residues and {len(renamed)} renamed residues in the receptor file"
            ),
            "evidence": {"unexpected": added[:40], "renamed": renamed[:40], "removed_by_design": len(dropped)},
        },
        {
            "check": "native_domain_only",
            "status": AssessmentStatus.PASS.value if not foreign_in_receptor else AssessmentStatus.FAIL.value,
            "detail": (
                "no atom outside the declared native range reaches later stages"
                if not foreign_in_receptor
                else "receptor still contains residues outside the native range"
            ),
            "evidence": foreign_in_receptor[:40],
        },
        {
            "check": "complete_native_coverage",
            "status": AssessmentStatus.PASS.value if not missing_from_receptor else AssessmentStatus.FLAG.value,
            "detail": (
                "all modelled native residues are present in the receptor"
                if not missing_from_receptor
                else f"{len(missing_from_receptor)} modelled native residues are absent from the receptor"
            ),
            "evidence": missing_from_receptor[:40],
        },
        {
            "check": "reference_complex_archived",
            "status": AssessmentStatus.PASS.value if os.path.exists(reference_complex_path) else AssessmentStatus.FAIL.value,
            "detail": "the untouched deposited coordinates stay available for mechanism reference",
            "evidence": os.path.basename(reference_complex_path),
        },
        {
            "check": "wild_type_sidechain_restoration",
            "status": AssessmentStatus.NOT_EVALUATED.value,
            "detail": "not performed: no rotamer rebuild is attempted, deposited mutant side chains are kept",
            "evidence": [record["auth_seq_id"] for record in mutation_records],
        },
    ]

    report = {
        "chain": chain,
        "native_range": list(native_range) if native_range else None,
        "original": {
            "path": original_path,
            "sha256": sha256_file(original_path),
            "atoms": len(before.atoms),
            "chains": _chain_inventory(before),
            "segments": segments,
        },
        "receptor": {
            "path": receptor_path,
            "sha256": sha256_file(receptor_path),
            "atoms": len(after.atoms),
            "residues": len(after_res),
            "auth_first": min(after_res) if after_res else None,
            "auth_last": max(after_res) if after_res else None,
            "chains": _chain_inventory(after),
        },
        "reference_complex": {
            "path": reference_complex_path,
            "sha256": sha256_file(reference_complex_path) if os.path.exists(reference_complex_path) else "",
            "chains": _chain_inventory(before),
            "role": "mechanism reference only; never used as the docking receptor",
        },
        "deletions": dict(strip_log or {}),
        "removed_residues": {
            name: sorted({r for r in before_res if _segment_of(r, native_range) == name})
            for name in segments
            if name != "native_domain"
        },
        "mutations": mutation_records,
        "c_terminal_tag": c_terminal_tag or {},
        "checks": checks,
        "notes": [
            "atom counts are measured, not asserted; the deposit is the reference, not a target number",
            "chain B of 9NU4 is a 29-residue peptide (Sgl) and is kept only in the reference complex",
        ],
    }
    report["status"] = (
        AssessmentStatus.PASS.value
        if all(check["status"] == AssessmentStatus.PASS.value for check in checks[:4])
        else AssessmentStatus.FAIL.value
    )
    return report


def write_audit(report: dict[str, Any], directory: str) -> dict[str, str]:
    os.makedirs(directory, exist_ok=True)
    paths: dict[str, str] = {}
    paths["json"] = atomic_write_json(os.path.join(directory, "receptor_audit.json"), report)

    chain = report["chain"]
    before = {row["resseq"]: row for row in _residue_rows(report["original"]["path"], chain)}
    after = {row["resseq"]: row for row in _residue_rows(report["receptor"]["path"], chain)}
    paths["residues_before"] = _write_csv(
        os.path.join(directory, "residues_before.csv"),
        ["resseq", "resname", "atoms", "segment"],
        [
            {**row, "segment": _segment_of(row["resseq"], tuple(report["native_range"])) if report["native_range"] else "unknown"}
            for row in before.values()
        ],
    )
    paths["residues_after"] = _write_csv(
        os.path.join(directory, "residues_after.csv"),
        ["resseq", "resname", "atoms", "present_in_deposit"],
        [{**row, "present_in_deposit": row["resseq"] in before} for row in after.values()],
    )
    removed: list[dict[str, Any]] = []
    for name, residues in (report.get("removed_residues") or {}).items():
        for resseq in residues:
            info = before.get(resseq, {})
            removed.append({"resseq": resseq, "resname": info.get("resname", ""), "segment": name})
    paths["deletions"] = _write_csv(os.path.join(directory, "deletions.csv"), ["resseq", "resname", "segment"], removed)
    paths["mutations"] = _write_csv(
        os.path.join(directory, "mutations.csv"),
        ["auth_seq_id", "declared_native_resname", "declared_construct_resname",
         "mapping_native_resname", "observed_resname_before", "observed_resname_after", "wt_restored", "status"],
        report.get("mutations") or [],
    )
    paths["markdown"] = atomic_write_text(
        os.path.join(directory, "receptor_audit.md"), render_markdown(report)
    )
    return paths


def _residue_rows(path: str, chain: str) -> list[dict[str, Any]]:
    structure = load_structure(path)
    rows = []
    for resseq, info in sorted(_residue_map(structure, chain).items()):
        rows.append({"resseq": resseq, "resname": info["resname"], "atoms": info["atoms"]})
    return rows


def _write_csv(path: str, fields: list[str], rows: Iterable[dict[str, Any]]) -> str:
    rows = list(rows)
    lines = [",".join(fields)]
    for row in rows:
        values = []
        for field in fields:
            value = row.get(field, "")
            text = "" if value is None else str(value)
            if "," in text:
                text = '"' + text.replace('"', '""') + '"'
            values.append(text)
        lines.append(",".join(values))
    return atomic_write_text(path, "\n".join(lines) + "\n")


def render_markdown(report: dict[str, Any]) -> str:
    lines = ["# Receptor audit", ""]
    original = report["original"]
    receptor = report["receptor"]
    lines += [
        f"- Chain: `{report['chain']}` · Native interval: `{report['native_range']}`",
        f"- Input before processing: `{os.path.basename(original['path'])}` sha256 `{original['sha256'][:16]}` · {original['atoms']} Atom",
        f"- Actual downstream receptor: `{os.path.basename(receptor['path'])}` sha256 `{receptor['sha256'][:16]}` · "
        f"{receptor['atoms']} Atom / {receptor['residues']} Residue / auth {receptor['auth_first']}..{receptor['auth_last']}",
        f"- Mechanistic reference complex: `{os.path.basename(report['reference_complex']['path'])}` "
        f"sha256 `{(report['reference_complex']['sha256'] or '')[:16]}`",
        "",
        "## Composition before processing",
        "",
        "| Chain | Protein atoms | Protein residues | auth Range | Nonprotein atoms |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in original["chains"]:
        lines.append(
            f"| {row['chain']} | {row['protein_atoms']} | {row['protein_residues']} | "
            f"{row['auth_first']}..{row['auth_last']} | {row['hetero_atoms']} |"
        )
    lines += ["", "## Backbone segments", "", "| Segment | Residue count | Atom count | auth Range |", "| --- | --- | --- | --- |"]
    for name, slot in sorted(report["original"]["segments"].items()):
        lines.append(f"| {name} | {slot['residues']} | {slot['atoms']} | {slot['first']}..{slot['last']} |")
    lines += [
        "",
        "## Actual receptor composition",
        "",
        "| Chain | Protein atoms | Protein residues | auth Range |",
        "| --- | --- | --- | --- |",
    ]
    for row in receptor["chains"]:
        lines.append(
            f"| {row['chain']} | {row['protein_atoms']} | {row['protein_residues']} | {row['auth_first']}..{row['auth_last']} |"
        )
    lines += ["", "## Check", "", "| Check | Status | Description |", "| --- | --- | --- |"]
    for check in report["checks"]:
        lines.append(f"| {check['check']} | {check['status']} | {check['detail']} |")
    if report.get("mutations"):
        lines += [
            "",
            "## Mutation sites",
            "",
            "| auth | Declared construct | Repository residues | After processing | WT Resume | Status |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
        for row in report["mutations"]:
            lines.append(
                f"| {row['auth_seq_id']} | {row['declared_construct_resname']} | {row['observed_resname_before']} | "
                f"{row['observed_resname_after']} | {row['wt_restored']} | {row['status']} |"
            )
    lines += ["", "## Description", ""] + [f"- {note}" for note in report.get("notes", [])]
    return "\n".join(lines) + "\n"
