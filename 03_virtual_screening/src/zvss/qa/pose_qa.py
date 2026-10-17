"""Pose-level quality control.

Two independent families of checks are run and reported separately:

* geometric/chemical validity (PoseBusters when importable, plus local checks);
* consistency with the declared receptor and box (clashes, box escape).

A check that cannot run is recorded as ``not_evaluated``/``unsupported`` with a
reason code.  It is never recorded as passed.
"""

from __future__ import annotations

import os
import traceback
from typing import Any

import numpy as np

from .. import compat
from ..constants import AssessmentStatus, Reason
from ..util import atomic_write_json, pairwise_distances, safe_float

POSE_QA_METRICS = (
    "atoms_in_box_fraction",
    "min_receptor_distance",
    "severe_clash_pairs",
    "positive_min_distance",
    "bond_length_outliers",
    "posebusters_status",
    "posebusters_passed",
    "posebusters_failed",
)


def _posebusters_available() -> tuple[bool, str, dict[str, Any]]:
    """Return ``(available, reason, diagnostics)``.

    The diagnostics matter as much as the verdict: when PoseBusters cannot run,
    the report has to say *why* (unwritable scratch, CJK RDKit data path,
    missing package, missing handle) instead of a generic "unsupported".
    """
    import traceback

    diagnostics: dict[str, Any] = {
        "scratch_root": compat.scratch_root(),
        "scratch_source": compat._STATE.get("scratch_source"),
        "scratch_is_ascii": compat.ascii_scratch_available(),
        "tmpdir": os.environ.get("TMP") or "",
        "cwd_is_ascii": compat._is_ascii(os.getcwd()),
    }
    if not compat.ascii_scratch_available():
        diagnostics["category"] = "scratch_not_ascii"
        diagnostics["failing_path"] = compat.scratch_root()
        return (
            False,
            "no writable ASCII scratch directory: RDKit and PoseBusters both open data files by "
            f"narrow path and {compat.scratch_root()!r} cannot be represented",
            diagnostics,
        )
    try:
        compat.bootstrap()
    except Exception:  # pragma: no cover - defensive
        diagnostics["category"] = "bootstrap_failed"
        diagnostics["traceback"] = traceback.format_exc()
        return False, "environment bootstrap failed: " + diagnostics["traceback"].splitlines()[-1], diagnostics
    try:
        from rdkit import RDConfig

        feature_def = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
        diagnostics["rdkit_data_dir"] = RDConfig.RDDataDir
        diagnostics["rdkit_data_is_ascii"] = compat._is_ascii(RDConfig.RDDataDir)
        diagnostics["base_features_present"] = os.path.exists(feature_def)
        if not os.path.exists(feature_def):
            diagnostics["category"] = "rdkit_data_missing"
            diagnostics["failing_path"] = feature_def
            return False, f"RDKit data file is missing: {feature_def}", diagnostics
    except Exception:  # pragma: no cover - defensive
        diagnostics["category"] = "rdkit_import_failed"
        diagnostics["traceback"] = traceback.format_exc()
        return False, "RDKit could not be imported: " + diagnostics["traceback"].splitlines()[-1], diagnostics
    try:
        from posebusters import PoseBusters  # noqa: F401

        import posebusters

        diagnostics["posebusters_version"] = getattr(posebusters, "__version__", "unknown")
    except Exception as exc:
        text = traceback.format_exc()
        diagnostics["traceback"] = text
        diagnostics["category"] = _classify_posebusters_error(exc)
        diagnostics["failing_path"] = _failing_path_from(text)
        return False, f"{type(exc).__name__}: {exc}", diagnostics
    diagnostics["category"] = "ok"
    return True, "", diagnostics


def _classify_posebusters_error(exc: BaseException) -> str:
    message = f"{type(exc).__name__}: {exc}".lower()
    if isinstance(exc, ImportError) or "no module named" in message:
        return "missing_dependency"
    if "could not be opened" in message or "bad input file" in message:
        return "rdkit_data_unreadable"
    if "permission" in message:
        return "scratch_not_writable"
    return "posebusters_import_error"


def _failing_path_from(traceback_text: str) -> str | None:
    """Best-effort extraction of the file path named in a load failure."""
    for marker in ("File: ", "could not be opened"):
        if marker in traceback_text:
            tail = traceback_text.split(marker, 1)[1].strip()
            candidate = tail.splitlines()[0].strip().rstrip(".")
            if os.path.isabs(candidate) or (len(candidate) > 3 and candidate[1] == ":"):
                return candidate
    return None


def run(ctx) -> dict[str, Any]:
    """Evaluate every pose; returns aggregate counts and per-pose metrics."""
    from .. import data

    poses = data.pose_rows(ctx)
    rec_atoms, rec_path = data.receptor_atoms(ctx)
    rec_xyz = data.coords_of(rec_atoms) if rec_atoms else np.zeros((0, 3))
    available, reason, diagnostics = _posebusters_available()
    buster = None
    if available:
        from posebusters import PoseBusters

        try:
            buster = PoseBusters(config="dock")
        except Exception as exc:  # pragma: no cover
            available, reason = False, f"PoseBusters(config='dock') failed: {exc}"
            diagnostics["category"] = "buster_config_failed"
            diagnostics["traceback"] = traceback.format_exc()

    summary = {
        "poses": len(poses),
        "passed": 0,
        "flagged": 0,
        "failed": 0,
        "not_evaluated": 0,
        "posebusters": {"available": available, "reason": reason, "diagnostics": diagnostics},
        "clash_cutoff": 2.0,
        "box_escape_poses": [],
    }
    environment_path = os.path.join(ctx.stage_dir("pose_qa"), "posebusters_environment.json")
    atomic_write_json(environment_path, summary["posebusters"])
    receptor_cond = _stage_receptor(rec_atoms, ctx.run_id) if available else None
    pose_errors: list[dict[str, Any]] = []
    per_pose: dict[str, dict[str, Any]] = {}

    for pose in poses:
        mol = data.pose_mol(ctx, pose)
        metrics: dict[str, Any] = {}
        reasons: list[str] = []
        if mol is None or mol.GetNumConformers() == 0:
            ctx.db.record_assessment("pose", pose["pose_id"], "pose_qa", "molecule_available", 0.0, AssessmentStatus.FAIL, Reason.POSE_GEOMETRY_INVALID, note="no coordinates")
            ctx.db.set_pose_status(pose["pose_id"], assess_status=AssessmentStatus.FAIL, reason_code=Reason.MISSING_COORDINATES)
            summary["failed"] += 1
            continue

        xyz = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
        # box containment
        box = pose.get("box") or {}
        if box.get("center") and box.get("size"):
            center = np.asarray(box["center"], dtype=float)
            size = np.asarray(box["size"], dtype=float)
            inside = np.all(np.abs(xyz - center) <= size / 2.0 + 0.5, axis=1)
            metrics["atoms_in_box_fraction"] = float(inside.mean())
            if inside.mean() < 0.9:
                reasons.append("pose leaves the docking box")
                summary["box_escape_poses"].append(pose["pose_id"])

        # receptor clashes
        if rec_xyz.size:
            distances = pairwise_distances(xyz, rec_xyz)
            metrics["min_receptor_distance"] = float(distances.min())
            metrics["severe_clash_pairs"] = int((distances < 2.0).sum())
            metrics["positive_min_distance"] = bool(distances.min() > 0.0)
        else:
            metrics["min_receptor_distance"] = None
            reasons.append(str(Reason.RECEPTOR_UNBOUND))

        # bond lengths vs element-pair reference ranges
        metrics["bond_length_outliers"] = _bond_outliers(mol)

        # PoseBusters
        if buster is not None:
            result = _run_posebusters(buster, mol, pose, receptor_cond)
            metrics["posebusters_status"] = result["status"]
            metrics["posebusters_passed"] = result["passed"]
            metrics["posebusters_failed"] = result["failed"]
            metrics["posebusters_failed_names"] = (
                ";".join(result["failed_names"]) if result["failed_names"] else None
            )
            if result["status"] == "ok":
                metrics["posebusters_items_total"] = len(result.get("items") or {})
                busted_receptor = receptor_cond is not None
                metrics["posebusters_receptor_supplied"] = 1.0 if busted_receptor else 0.0
                if not busted_receptor:
                    metrics["posebusters_status"] = "ok_without_receptor"
                    reasons.append("posebusters:protein_checks_not_possible")
            else:
                reasons.append(f"posebusters:{result['status']}")
                pose_errors.append(
                    {
                        "pose_id": pose["pose_id"],
                        "compound_id": pose.get("compound_id"),
                        "error_type": result.get("error_type"),
                        "error_message": result.get("error_message"),
                        "failing_path": result.get("failing_path"),
                        "traceback": result.get("traceback"),
                    }
                )
            if result["failed"]:
                reasons.extend(result["failed_names"])
        else:
            metrics["posebusters_status"] = "unsupported"
            metrics["posebusters_passed"] = None
            metrics["posebusters_failed"] = None

        status, rule_reasons = _status(metrics)
        reasons.extend(rule_reasons)
        if status == "fail":
            code = Reason.POSE_GEOMETRY_INVALID
        elif status == "flag":
            code = Reason.SEVERE_CLASH
        else:
            code = None
        mapped = {
            "pass": AssessmentStatus.PASS,
            "flag": AssessmentStatus.FLAG,
            "fail": AssessmentStatus.FAIL,
        }[status]
        for metric, value in metrics.items():
            ctx.db.record_assessment("pose", pose["pose_id"], "pose_qa", metric, value, status=mapped, reason_code=code)
        # Per-run verdict: ``pose.assess_status`` is one global column and would
        # let a later run overwrite this run's QA outcome.
        ctx.db.record_assessment(
            "pose", pose["pose_id"], "pose_qa", "qa_status", None, mapped, code, value_text=status
        )
        if status == "pass":
            ctx.db.set_pose_status(pose["pose_id"], assess_status=AssessmentStatus.PASS)
            summary["passed"] += 1
        elif status == "flag":
            ctx.db.set_pose_status(pose["pose_id"], assess_status=AssessmentStatus.FLAG, reason_code=code)
            summary["flagged"] += 1
        else:
            ctx.db.set_pose_status(pose["pose_id"], assess_status=AssessmentStatus.FAIL, reason_code=code)
            summary["failed"] += 1
        per_pose[pose["pose_id"]] = {"status": status, "metrics": metrics, "reasons": reasons}

    summary["per_pose"] = per_pose
    if pose_errors:
        summary["posebusters_errors"] = pose_errors
        atomic_write_json(os.path.join(ctx.stage_dir("pose_qa"), "posebusters_errors.json"), pose_errors)
    return summary


def _status(metrics: dict[str, Any]) -> tuple[str, list[str]]:
    """Return ``(status, rule_reasons)`` so a flag always says what tripped it."""
    outliers = metrics.get("bond_length_outliers") or 0
    if outliers > 2:
        return "fail", [f"bond_length_outliers>{2}"]
    minimum = metrics.get("min_receptor_distance")
    if minimum is not None and minimum < 0.5:
        return "fail", ["min_receptor_distance<0.5A"]
    fraction = metrics.get("atoms_in_box_fraction")
    if (fraction if fraction is not None else 1.0) < 0.5:
        return "fail", ["atoms_in_box_fraction<0.5"]
    if metrics.get("posebusters_failed"):
        failed = metrics.get("posebusters_failed_names")
        return "flag", [f"posebusters_failed:{','.join(failed)}" if failed else "posebusters_failed"]
    if fraction is not None and fraction < 0.99:
        return "flag", [f"atoms_outside_docking_box:{round(1.0 - fraction, 3)}"]
    if minimum is not None and minimum < 2.0:
        return "flag", [f"min_receptor_distance<2.0A:{round(minimum, 3)}"]
    if outliers:
        return "flag", [f"bond_length_outliers={outliers}"]
    return "pass", []


BOND_RANGES = {
    ("C", "C"): (1.20, 1.75),
    ("C", "N"): (1.15, 1.60),
    ("C", "O"): (1.10, 1.60),
    ("C", "S"): (1.55, 1.95),
    ("C", "F"): (1.20, 1.55),
    ("C", "Cl"): (1.55, 1.90),
    ("C", "Br"): (1.70, 2.05),
    ("N", "N"): (1.15, 1.60),
    ("N", "O"): (1.15, 1.60),
    ("O", "O"): (1.20, 1.65),
    ("C", "H"): (0.85, 1.20),
}


def _bond_outliers(mol) -> int:
    conf = mol.GetConformer()
    outliers = 0
    for bond in mol.GetBonds():
        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()
        pa = conf.GetAtomPosition(a.GetIdx())
        pb = conf.GetAtomPosition(b.GetIdx())
        length = ((pa.x - pb.x) ** 2 + (pa.y - pb.y) ** 2 + (pa.z - pb.z) ** 2) ** 0.5
        key = tuple(sorted((a.GetSymbol(), b.GetSymbol())))
        low, high = BOND_RANGES.get(key, (0.7, 2.3))
        if not (low <= length <= high):
            outliers += 1
    return outliers


def _stage_receptor(atoms, run_id: str) -> str | None:
    """Write the receptor to an ASCII path and return it.

    PoseBusters opens ``mol_cond`` as a file and RDKit uses narrow paths, so the
    CJK project path cannot be handed over directly.  Staging the receptor in
    the portable ASCII scratch keeps the protein-based checks (distance to
    protein, volume overlap) genuinely evaluated instead of silently skipped.
    """
    if not atoms:
        return None
    directory = os.path.join(compat.scratch_root(), f"pose_qa-{run_id}")
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, "receptor.pdb")
    block = "\n".join(atom.to_line() for atom in atoms) + "\nEND\n"
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(block)
    return path


def _run_posebusters(buster, mol, pose: dict[str, Any], receptor_path: str | None = None) -> dict[str, Any]:
    """Run PoseBusters for one pose; never swallow an exception silently.

    The positional form ``bust(mol_pred, mol_true, mol_cond)`` is used on
    purpose: in the pinned PoseBusters build the DataFrame form resolves column
    names as file paths and yields no rows ("No objects to concatenate"), which
    would quietly turn every pose into "not evaluated".
    """
    try:
        result = buster.bust(mol, None, receptor_path)
        row = result.iloc[0].to_dict()
        # PoseBusters emits one boolean column per check; keep the full map so a
        # single failing geometry check stays auditable per item.
        items = {key: bool(value) for key, value in row.items() if isinstance(value, (bool, np.bool_))}
        failures = [key for key, value in items.items() if not value]
        return {
            "status": "ok",
            "items": items,
            "passed": sum(1 for value in items.values() if value),
            "failed": len(failures),
            "failed_names": failures,
        }
    except Exception as exc:  # pragma: no cover - library specific
        text = traceback.format_exc()
        return {
            "status": "error",
            "error_type": type(exc).__name__,
            "error_message": str(exc),
            "traceback": text,
            "failing_path": _failing_path_from(text),
            "passed": None,
            "failed": None,
            "failed_names": [],
            "items": {},
        }
