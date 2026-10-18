"""Import stage: turn discovered files into identity + pose rows.

Four identity layers are created here and nothing downstream is allowed to merge
them:

* ``compound_id``      - parent compound, stereochemistry preserved;
* ``state_id``         - concrete protonation/charge state;
* ``job_id``           - receptor version + box + protocol + seed + state;
* ``pose_id``          - one concrete pose, bound to its own score record.

A molecule whose chemical graph cannot be recovered keeps its coordinates and
scores but is flagged ``needs_input`` and can never reach the MD export list.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any

from rdkit import Chem

from ..chem import rdio
from ..chem.recover import RecoveredLigand, recover_from_template
from ..constants import AssessmentStatus, ExecutionStatus, Reason, SelectionStatus
from ..ids import chemical_state_id, compound_id, docking_job_id, pose_id, receptor_version, sha256_text
from ..util import atomic_write_json, now_iso, read_json, relpath, safe_int
from . import poses as pose_io
from .discover import InputBundle, box_from_log

# --------------------------------------------------------------------------- #
# records
# --------------------------------------------------------------------------- #


@dataclass
class ImportedPose:
    pose_row: dict[str, Any]
    recovered: RecoveredLigand | None
    compound_row: dict[str, Any] | None
    state_row: dict[str, Any] | None
    notes: list[str] = field(default_factory=list)


@dataclass
class ImportReport:
    run_id: str
    poses_seen: int = 0
    poses_imported: int = 0
    poses_failed: int = 0
    parents: int = 0
    states: int = 0
    jobs: int = 0
    needs_input: list[dict[str, Any]] = field(default_factory=list)
    conflicts: list[dict[str, Any]] = field(default_factory=list)
    files: list[dict[str, Any]] = field(default_factory=list)
    table_binding: dict[str, Any] | None = None
    log_binding: list[dict[str, Any]] = field(default_factory=list)
    receptor: dict[str, Any] = field(default_factory=dict)
    box: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "run_id": self.run_id,
            "poses_seen": self.poses_seen,
            "poses_imported": self.poses_imported,
            "poses_failed": self.poses_failed,
            "parents": self.parents,
            "states": self.states,
            "jobs": self.jobs,
            "needs_input": self.needs_input,
            "conflicts": self.conflicts,
            "files": self.files,
            "table_binding": self.table_binding,
            "log_binding": self.log_binding,
            "receptor": self.receptor,
            "box": self.box,
            "warnings": self.warnings,
        }


# --------------------------------------------------------------------------- #
# receptor + box resolution
# --------------------------------------------------------------------------- #


def resolve_receptor(ctx, bundle: InputBundle, hint: dict[str, Any] | None = None) -> dict[str, Any]:
    hint = hint or {}
    candidate = hint.get("receptor_file") or bundle.receptor_file
    if candidate and not os.path.isabs(candidate):
        candidate = os.path.join(bundle.root, candidate)
    if not candidate or not os.path.exists(candidate):
        return {
            "status": "unbound",
            "reason_code": str(Reason.RECEPTOR_UNBOUND),
            "note": "no receptor file was found or declared; pose scores keep their original receptor hash only",
        }
    from ..ids import sha256_file

    digest = sha256_file(candidate)
    declared = hint.get("receptor_sha256")
    if declared and declared != digest:
        raise ValueError("declared docking receptor sha256 does not match the file")
    recipe = str(hint.get("preparation") or "as_provided")
    frame = str(hint.get("coordinate_frame") or "receptor_input_frame")
    version = receptor_version(digest, recipe, frame, str(hint.get("mutant_recipe") or ""))
    archived = ctx.archive_input(candidate)
    return {
        "status": "bound",
        "path": candidate,
        "archived": archived,
        "sha256": digest,
        "version": version,
        "preparation": recipe,
        "coordinate_frame": frame,
    }


def resolve_box(bundle: InputBundle, hint: dict[str, Any] | None = None, log_texts: list[str] | None = None) -> dict[str, Any]:
    from .discover import detect_box

    if hint and hint.get("box"):
        box = detect_box(hint["box"])
        if box:
            box["source"] = box.get("source") or "manifest"
            return box
    if bundle.box_file and os.path.exists(bundle.box_file):
        box = detect_box(read_json(bundle.box_file, {}))
        if box:
            box["source"] = bundle.box_file
            return box
    for text in log_texts or []:
        box = box_from_log(text)
        if box:
            return box
    return {}


# --------------------------------------------------------------------------- #
# chemical template resolution
# --------------------------------------------------------------------------- #


def _template_candidates(bundle: InputBundle, hints: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Map parent label / resname -> template source declared in the manifest."""
    out: dict[str, dict[str, Any]] = {}
    templates = hints.get("templates") or {}
    if isinstance(templates, dict):
        for key, value in templates.items():
            if isinstance(value, str):
                out[str(key)] = {"smiles": value} if not value.lower().endswith((".sdf", ".mol")) else {"sdf": value}
            elif isinstance(value, dict):
                out[str(key)] = dict(value)
    return out


def _load_template(spec: dict[str, Any], bundle: InputBundle) -> tuple[Chem.Mol | None, str]:
    if not spec:
        return None, ""
    if spec.get("sdf"):
        path = str(spec["sdf"])
        if not os.path.isabs(path):
            path = os.path.join(bundle.root, path)
        if os.path.exists(path):
            for mol, _props in rdio.iter_sdf(path):
                if mol is not None:
                    return mol, f"sdf:{path}"
        return None, f"sdf missing: {path}"
    if spec.get("smiles"):
        mol = rdio.mol_from_smiles(str(spec["smiles"]))
        return (mol, f"smiles:{spec['smiles']}") if mol is not None else (None, "invalid smiles")
    return None, "template spec has neither sdf nor smiles"


def _adjacent_ligand_template(pose_path: str) -> tuple[Chem.Mol | None, str]:
    directory = os.path.dirname(pose_path)
    for name in ("ligand.sdf", "ligand_template.sdf", "input.sdf", "template.sdf"):
        candidate = os.path.join(directory, name)
        if os.path.exists(candidate):
            for mol, _props in rdio.iter_sdf(candidate):
                if mol is not None:
                    return mol, candidate
    return None, ""


# --------------------------------------------------------------------------- #
# identity helpers
# --------------------------------------------------------------------------- #


def _identity(mol: Chem.Mol, *, note: str = "") -> tuple[dict[str, Any], dict[str, Any]]:
    smiles = rdio.canonical_isomeric_smiles(mol)
    if not smiles:
        raise ValueError("molecule has no canonical SMILES")
    charge = rdio.net_charge(mol)
    cid = compound_id(smiles)
    sid = chemical_state_id(smiles, charge, note)
    compound_row = {
        "compound_id": cid,
        "canonical_smiles": smiles,
        "inchikey": _inchikey(mol),
        "formula": rdio.mol_formula(mol),
        "heavy_atoms": rdio.heavy_atom_count(mol),
        "net_charge": charge,
        "parent_label": note or None,
        "meta": {"stereo_note": "isomeric SMILES retained"},
    }
    state_row = {
        "state_id": sid,
        "compound_id": cid,
        "canonical_smiles": smiles,
        "formal_charge": charge,
        "protonation_note": note or "as_provided",
        "tautomer_note": "",
        "stereo_note": "from input chemistry",
        "source": "import",
    }
    return compound_row, state_row


def _inchikey(mol: Chem.Mol) -> str | None:
    try:
        return Chem.MolToInchiKey(mol)
    except Exception:
        return None


# --------------------------------------------------------------------------- #
# main import
# --------------------------------------------------------------------------- #


def import_inputs(
    ctx,
    bundle: InputBundle,
    *,
    hints: dict[str, Any] | None = None,
    selector: dict[str, str] | None = None,
    resume: bool = False,
) -> ImportReport:
    hints = {**bundle.hints, **(hints or {})}
    report = ImportReport(run_id=ctx.run_id)
    report.warnings.extend(bundle.warnings)
    report.files = [item.to_dict() for item in bundle.files]

    tables = [item for item in bundle.by_kind("table")]
    logs = [item for item in bundle.by_kind("log") if item.note == "gnina_log"]

    table = None
    table_report: dict[str, Any] | None = None
    if tables:
        table = pose_io.read_table(tables[0].path)
        table_report = {"file": tables[0].path, "columns": table.columns, "rows": len(table.rows)}

    # GNINA logs: parse early so the box and seed can be recovered from them.
    log_texts: list[str] = []
    parsed_logs = []
    from ..util import read_text

    for item in logs:
        text = read_text(item.path)
        log_texts.append(text)
        parsed = pose_io.parse_gnina_log(text, item.path)
        parsed_logs.append(parsed)
        report.log_binding.append(
            {
                "file": item.path,
                "version": parsed.version,
                "seed": parsed.seed,
                "model": parsed.model,
                "modes": len(parsed.modes),
                "warnings": parsed.warnings[:20],
            }
        )

    # Protocol values recovered from a real GNINA log take precedence over
    # generic defaults; the log is evidence, the default is only a placeholder.
    if parsed_logs:
        first = parsed_logs[0]
        hints.setdefault("engine", "gnina")
        if first.scoring:
            hints.setdefault("scoring_function", first.scoring)
        if first.model:
            hints.setdefault("cnn_model", first.model)
        if first.seed is not None:
            hints.setdefault("seed", first.seed)
        if first.version:
            hints.setdefault("gnina_version", first.version)

    report.receptor = resolve_receptor(ctx, bundle, hints)
    report.box = resolve_box(bundle, hints, log_texts)

    pose_files = [item for item in bundle.files if item.kind in {"pose_sdf", "pose_pdb"}]
    collected = pose_io.collect_poses(pose_files, table=table, selector=selector,
                                     titles={str(k): str(v) for k, v in (hints.get("pose_titles") or {}).items()})
    raw_poses: list[pose_io.RawPose] = collected["poses"]
    report.poses_seen = len(raw_poses)
    report.table_binding = table_report
    for entry in collected["reports"]:
        if isinstance(entry, dict) and entry.get("table_binding"):
            report.table_binding = {**(table_report or {}), **entry["table_binding"]}

    # Cross-check GNINA log mode tables against parsed SDF records (per file).
    for parsed in parsed_logs:
        matching = [p for p in raw_poses if os.path.dirname(p.source_file) == os.path.dirname(parsed.path)]
        if not matching:
            matching = raw_poses
        records = [{"scores": p.scores} for p in matching]
        comparison = pose_io.match_log_modes_to_records(parsed, records)
        if comparison["log_modes"]:
            report.log_binding.append(
                {
                    "file": parsed.path,
                    "cross_check": comparison,
                    "note": "log mode rows paired to SDF records by score values, not by order",
                }
            )

    templates = _template_candidates(bundle, hints)
    label_to_smiles: dict[str, set[str]] = {}
    recovered_by_pose: dict[str, RecoveredLigand | None] = {}

    for raw in raw_poses:
        label = str(raw.extras.get("table_parent_id") or raw.parent_label or "").strip()
        if label:
            raw.parent_label = label
        try:
            (pose_row, recovered, compound_row, state_row, notes) = _import_single(
                ctx, raw, bundle, hints, templates, report
            )
        except Exception as exc:  # isolate molecule-level failure
            report.poses_failed += 1
            report.needs_input.append(
                {
                    "source_file": raw.source_file,
                    "record_index": raw.record_index,
                    "reason": f"{type(exc).__name__}: {exc}",
                    "reason_code": str(Reason.CHEM_TEMPLATE_MISSING),
                }
            )
            ctx.log(
                f"pose import failed: {exc}",
                stage="import",
                level="warning",
                code=str(Reason.CHEM_TEMPLATE_MISSING),
                payload={"source": raw.source_file, "record": raw.record_index},
            )
            continue

        if pose_row is None:
            report.poses_failed += 1
            continue

        if compound_row and state_row:
            ctx.db.upsert_compound(compound_row)
            ctx.db.upsert_state(state_row)
            if label:
                label_to_smiles.setdefault(label, set()).add(state_row["canonical_smiles"])
        stored_id = ctx.db.upsert_pose(pose_row)
        pose_row["pose_id"] = stored_id
        recovered_by_pose[stored_id] = recovered
        report.poses_imported += 1
        if notes:
            for note in notes:
                report.warnings.append(f"{os.path.basename(raw.source_file)}#{raw.record_index}: {note}")

        # Persist the pose coordinates that were actually scored.
        if raw.mol is not None:
            coords_path = os.path.join(ctx.poses_dir(), f"{stored_id}.sdf")
            rdio.write_sdf(
                coords_path,
                [(raw.mol, {
                    "zvss_pose_id": stored_id,
                    "zvss_source_file": os.path.basename(raw.source_file),
                    "zvss_record_index": raw.record_index,
                    "zvss_parent_label": raw.parent_label,
                })],
            )
            ctx.db.execute("UPDATE pose SET coords_path=? WHERE pose_id=?", (coords_path, stored_id))

    for label, smiles_set in label_to_smiles.items():
        if len(smiles_set) > 1:
            report.conflicts.append(
                {
                    "parent_label": label,
                    "distinct_structures": sorted(smiles_set),
                    "note": "one parent label maps to several chemical structures; ids kept separate, not overwritten",
                }
            )

    counts = ctx.db.query(
        "SELECT (SELECT COUNT(*) FROM compound) AS compounds,"
        " (SELECT COUNT(*) FROM chemical_state) AS states,"
        " (SELECT COUNT(*) FROM docking_job) AS jobs"
    )[0]
    report.parents = counts["compounds"]
    report.states = counts["states"]
    report.jobs = counts["jobs"]

    # Cache the recovered molecules for downstream stages in one validated file.
    payload = {
        "run_id": ctx.run_id,
        "created": now_iso(),
        "poses": {
            pid: {
                "status": rec.status,
                "reason_code": rec.reason_code,
                "notes": rec.notes,
                "evidence": rec.evidence,
                "atom_map": rec.atom_map,
                "smiles": rec.smiles,
            }
            for pid, rec in recovered_by_pose.items()
            if rec is not None
        },
    }
    atomic_write_json(os.path.join(ctx.stage_dir("import"), "recovery_index.json"), payload)
    atomic_write_json(os.path.join(ctx.stage_dir("import"), "import_report.json"), report.to_dict())
    # The receptor the poses were docked against is stage-import state; the
    # re-dock stage (data.receptor_record) reads exactly this file.
    atomic_write_json(os.path.join(ctx.stage_dir("import"), "receptor.json"), report.receptor)
    ctx.db.upsert_run(ctx.run_id, status="imported")
    return report


def _import_single(ctx, raw, bundle, hints, templates, report):
    notes: list[str] = list(raw.issues)
    recovered: RecoveredLigand | None = None
    compound_row = None
    state_row = None

    template_spec = templates.get(raw.parent_label) or templates.get(raw.resname) or {}
    template_mol, template_source = _load_template(template_spec, bundle) if template_spec else (None, "")

    mol: Chem.Mol | None = raw.mol
    if raw.needs_template or mol is None:
        if template_mol is None:
            template_mol, template_source = _adjacent_ligand_template(raw.source_file)
        if template_mol is None:
            report.needs_input.append(
                {
                    "source_file": raw.source_file,
                    "record_index": raw.record_index,
                    "model_index": raw.model_index,
                    "reason": "coordinates present but no chemical template; provide the original SDF/SMILES",
                    "reason_code": str(Reason.TEMPLATE_NOT_FOUND),
                    "candidates": raw.extras.get("candidates"),
                }
            )
            notes.append("chemical graph not recovered")
            pose_row = _base_pose_row(
                ctx, raw, None, None, None, ExecutionStatus.COMPLETED, AssessmentStatus.NEEDS_INPUT,
                SelectionStatus.EXCLUDED, Reason.TEMPLATE_NOT_FOUND, notes,
            )
            return pose_row, None, None, None, notes
        recovered = recover_from_template(
            template_mol, raw.ligand_block, atom_names=raw.atom_names, resname=raw.resname
        )
        notes.extend(recovered.notes)
        if recovered.status in {"fail", "needs_input"}:
            pose_row = _base_pose_row(
                ctx, raw, recovered.mol, None, None, ExecutionStatus.COMPLETED,
                AssessmentStatus.FAIL if recovered.status == "fail" else AssessmentStatus.NEEDS_INPUT,
                SelectionStatus.EXCLUDED, recovered.reason_code, notes,
            )
            return pose_row, recovered, None, None, notes
        mol = recovered.mol
        raw.mol = mol
        notes.append(f"chemistry recovered from {template_source}")
    else:
        # SDF chemistry is already authoritative; validate it rather than trusting it.
        try:
            Chem.SanitizeMol(Chem.Mol(mol))
        except Exception as exc:
            report.needs_input.append(
                {
                    "source_file": raw.source_file,
                    "record_index": raw.record_index,
                    "reason": f"SDF chemistry failed sanitisation: {exc}",
                    "reason_code": str(Reason.CHEM_BOND_ORDER_UNRECOVERED),
                }
            )
            pose_row = _base_pose_row(
                ctx, raw, mol, None, None, ExecutionStatus.COMPLETED, AssessmentStatus.FAIL,
                SelectionStatus.EXCLUDED, Reason.CHEM_BOND_ORDER_UNRECOVERED, notes,
            )
            return pose_row, None, None, None, notes
        if template_mol is not None:
            from ..chem.recover import validate_bound_against_template

            check = validate_bound_against_template(mol, template_mol)
            if not check["graph_match"]:
                report.conflicts.append(
                    {
                        "source_file": raw.source_file,
                        "record_index": raw.record_index,
                        "template_source": template_source,
                        "detail": check,
                        "note": "pose chemistry disagrees with the declared template; pose kept but flagged",
                    }
                )
                notes.append("pose chemistry disagrees with declared template")

    compound_row, state_row = _identity(mol, note=raw.parent_label)
    # Foreign keys are enforced: compound -> state -> job -> pose, in that order.
    ctx.db.upsert_compound(compound_row)
    ctx.db.upsert_state(state_row)

    receptor_hash = str((hints.get("receptor_sha256") or report.receptor.get("sha256") or ""))
    box = report.box or {}
    box_hash = sha256_text(str(box.get("center")) + str(box.get("size")))
    protocol = {
        "engine": (hints.get("engine") or "gnina"),
        "scoring_function": (hints.get("scoring_function") or "vina"),
        "cnn_model": hints.get("cnn_model") or "crossdock_default2018_1_3",
        "exhaustiveness": hints.get("exhaustiveness"),
        "num_modes": hints.get("num_modes"),
        "seed": hints.get("seed"),
    }
    seed = safe_int(raw.extras.get("seed")) if raw.extras.get("seed") is not None else safe_int(hints.get("seed"))
    job_key = docking_job_id(
        state_row["state_id"],
        str(report.receptor.get("version") or "unknown"),
        box_hash,
        sha256_text(str(sorted(protocol.items()))),
        seed,
        str(protocol["engine"]),
        str(protocol["cnn_model"]),
    )

    ctx.db.upsert_job(
        {
            "job_id": job_key,
            "state_id": state_row["state_id"],
            "receptor_version": report.receptor.get("version"),
            "receptor_sha256": receptor_hash or None,
            "box": box or None,
            "box_hash": box_hash,
            "protocol": protocol,
            "protocol_hash": sha256_text(str(sorted(protocol.items()))),
            "seed": seed,
            "engine": protocol["engine"],
            "engine_version": hints.get("engine_version") or _engine_version(hints),
            "model_id": protocol["cnn_model"],
            "scoring_function": protocol["scoring_function"],
            "source": os.path.basename(raw.source_file),
        }
    )

    pose_row = _base_pose_row(
        ctx, raw, mol, compound_row, state_row, ExecutionStatus.COMPLETED,
        AssessmentStatus.NOT_EVALUATED, SelectionStatus.ELIGIBLE, None, notes,
        job_id=job_key,
    )
    return pose_row, recovered, compound_row, state_row, notes


def parsed_default(key: str) -> Any:  # pragma: no cover - placeholder kept for clarity
    return None


def _engine_version(hints: dict[str, Any]) -> str | None:
    version = hints.get("gnina_version") or hints.get("engine_version")
    return str(version) if version else None


def _base_pose_row(
    ctx,
    raw: pose_io.RawPose,
    mol: Chem.Mol | None,
    compound_row: dict | None,
    state_row: dict | None,
    exec_status: ExecutionStatus,
    assess_status: AssessmentStatus,
    select_status: SelectionStatus,
    reason_code: Reason | None,
    notes: list[str],
    *,
    job_id: str | None = None,
) -> dict[str, Any]:
    coords_hash = raw.coords_digest or sha256_text(f"{raw.source_file}:{raw.record_index}")
    job = job_id or "JOB_UNRESOLVED"
    pid = pose_id(job, raw.record_index, coords_hash, None)
    return {
        "pose_id": pid,
        "job_id": job,
        "state_id": state_row["state_id"] if state_row else "STA_UNRESOLVED",
        "compound_id": compound_row["compound_id"] if compound_row else "CMP_UNRESOLVED",
        "rank_in_job": raw.extras.get("rank") or raw.record_index + 1,
        "record_index": raw.record_index,
        "model_index": raw.model_index,
        "source_file": raw.source_file,
        "source_hash": raw.source_hash,
        "coords_hash": coords_hash,
        "scores": raw.scores,
        "exec_status": str(exec_status),
        "assess_status": str(assess_status),
        "select_status": str(select_status),
        "reason_code": str(reason_code) if reason_code else None,
        "reason": "; ".join(notes) if notes else None,
        "meta": {
            "score_provenance": raw.score_provenance,
            "parent_label": raw.parent_label,
            "resname": raw.resname,
            "ligand_chain": raw.chain,
            "ligand_resid": raw.resid,
            "notes": notes,
            "atom_count": mol.GetNumAtoms() if mol is not None else 0,
        },
    }


def recovery_index(ctx) -> dict[str, Any]:
    return read_json(os.path.join(ctx.stage_dir("import"), "recovery_index.json"), {}) or {}


def relative(path: str, ctx) -> str:
    return relpath(path, ctx.run_dir)
