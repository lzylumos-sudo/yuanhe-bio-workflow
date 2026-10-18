"""GNINA capability detection and multi-seed re-docking with pose clustering.

Nothing here pretends a tool exists.  ``detect`` runs the configured command
with a short timeout and records exactly what happened.  When GNINA is
unavailable the stage returns ``blocked_external_tool`` and the selected
candidates keep their evidence tier without a re-dock component, rather than
receiving an invented score.

When GNINA *is* available the stage runs one docking per configured seed and
answers the only question a re-dock can answer: do independent stochastic runs
land in the same pose cluster as the imported pose?  The counts that feed formal
candidate qualification are ``seeds_planned``, ``seeds_attempted``,
``seeds_succeeded`` and ``supporting_seeds_same_cluster``; they are recorded
separately and never collapsed into a single "reproducible: yes".
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from typing import Any

from ..constants import AssessmentStatus, ExecutionStatus, Reason
from ..ids import sha256_text
from ..util import atomic_write_json, ensure_dir, read_text, safe_float

VERSION_RE = re.compile(r"gnina\s*v?([0-9]+\.[0-9]+\.[0-9]+)", re.IGNORECASE)


@dataclass
class Capability:
    tool: str
    available: bool
    reason: str = ""
    version: str | None = None
    command: list[str] = field(default_factory=list)
    detail: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "tool": self.tool,
            "available": self.available,
            "reason": self.reason,
            "version": self.version,
            "command": self.command,
            "detail": self.detail,
        }


def _engine_environment(ctx):
 raw=getattr(ctx.config.host,'settings',None) or getattr(ctx.config.host,'raw',None) or {}
 extra=((raw.get('tools') or {}).get('gnina_settings') or {}).get('environment') or {}
 return {**os.environ,**{k:str(v) for k,v in extra.items() if k in {'LD_LIBRARY_PATH'}}}


def detect_gnina(ctx) -> Capability:
    host = ctx.config.host
    settings = getattr(host, "settings", None) or getattr(host, "raw", None) or {}
    tools = settings.get("tools") or {} if isinstance(settings, dict) else {}
    configured = tools.get("gnina") if isinstance(tools, dict) else None
    if configured is None and isinstance(settings, dict):
        configured = settings.get("gnina")
    candidates: list[list[str]] = []
    if isinstance(configured, list) and configured:
        candidates.append([str(part) for part in configured])
    elif isinstance(configured, str) and configured.strip():
        import shlex

        candidates.append(shlex.split(configured))
    executable = shutil.which("gnina")
    if executable:
        candidates.append([executable])
    if not candidates:
        return Capability("gnina", False, "no gnina command configured and none on PATH")
    last_detail: dict[str, Any] = {}
    for command in candidates:
        for flag in ("--version", "--help"):
            try:
                started = time.time()
                probe = subprocess.run(
                    [*command, flag], capture_output=True, text=True, timeout=120, check=False, env=_engine_environment(ctx)
                )
            except FileNotFoundError as exc:
                last_detail = {"command": command, "error": f"FileNotFoundError: {exc}"}
                break
            except subprocess.TimeoutExpired:
                last_detail = {"command": command, "error": "timeout after 120 s"}
                continue
            except Exception as exc:  # noqa: BLE001
                last_detail = {"command": command, "error": f"{type(exc).__name__}: {exc}"}
                break
            output = f"{probe.stdout}\n{probe.stderr}"
            match = VERSION_RE.search(output)
            ok = probe.returncode == 0 or bool(match)
            last_detail = {
                "command": command,
                "flag": flag,
                "returncode": probe.returncode,
                "seconds": round(time.time() - started, 3),
                "stdout_head": probe.stdout[:2000],
                "stderr_head": probe.stderr[:2000],
            }
            if ok:
                return Capability(
                    "gnina", True, version=match.group(1) if match else None, command=command, detail=last_detail
                )
    return Capability("gnina", False, "every configured gnina command failed to run", detail=last_detail)


def to_linux_path(path: str) -> str:
    """Translate a Windows path for a WSL-hosted binary (``C:\\x`` -> ``/mnt/c/x``)."""
    if not path:
        return path
    if path.startswith("/"):
        return path
    drive, _, rest = path.partition(":")
    if len(drive) == 1 and rest:
        return "/mnt/" + drive.lower() + rest.replace("\\", "/")
    return path.replace("\\", "/")


def _is_wsl_command(command: list[str]) -> bool:
    return bool(command) and os.path.basename(command[0]).lower() in {"wsl", "wsl.exe"}


def redock_argv(
    *,
    command: list[str],
    receptor: str,
    ligand: str,
    center: list[float],
    size: list[float],
    seed: int,
    out_prefix: str,
    exhaustiveness: int = 8,
    num_modes: int = 3,
    scoring: str = "vina",
    cnn: str = "rescore",
    cpu: int = 1,
    cnn_model: str | None = None,
    device: str | None = None,
) -> list[str]:
    translate = _is_wsl_command(command)
    receptor_arg = to_linux_path(receptor) if translate else receptor
    ligand_arg = to_linux_path(ligand) if translate else ligand
    out_arg = to_linux_path(out_prefix) if translate else out_prefix
    argv = [
        *command,
        "--receptor", receptor_arg,
        "--ligand", ligand_arg,
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", f"{size[0]:.3f}",
        "--size_y", f"{size[1]:.3f}",
        "--size_z", f"{size[2]:.3f}",
        "--scoring", scoring,
        "--cnn_scoring", cnn,
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--seed", str(seed),
        "--cpu", str(cpu),
        "--out", out_arg,
    ]
    if cnn_model: argv.extend(["--cnn", cnn_model])
    if device == "gpu": argv.extend(["--device", "0"])
    elif device == "cpu": argv.append("--no_gpu")
    elif device is not None: raise ValueError("Invalid GNINA device")
    return argv


def run_redock(ctx, plan: list[dict[str, Any]]) -> dict[str, Any]:
    """Run re-docking for the supplied plan and cluster the resulting poses."""
    capability = detect_gnina(ctx)
    summary: dict[str, Any] = {
        "capability": capability.to_dict(),
        "requested": len(plan),
        "completed": 0,
        "failed": 0,
        "skipped": 0,
        "results": [],
        "compounds": {},
    }
    atomic_write_json(os.path.join(ctx.stage_dir("redock"), "capability.json"), capability.to_dict())
    if not capability.available:
        summary["status"] = "blocked_external_tool"
        summary["reason_code"] = str(Reason.REDOCK_TOOL_UNAVAILABLE)
        summary["note"] = (
            "GNINA could not be executed from this environment; re-docking evidence is absent and "
            "the ranking keeps using the imported single-seed scores"
        )
        for item in plan:
            for pose in item.get("pose_ids", []):
                ctx.db.record_assessment("pose", pose, "redock", "redock_completed", None,
                                         AssessmentStatus.NEEDS_INPUT, Reason.REDOCK_TOOL_UNAVAILABLE)
        return summary

    receptor = _receptor_path(ctx)
    if not receptor:
        summary["status"] = "blocked_missing_input"
        summary["reason_code"] = str(Reason.RECEPTOR_UNBOUND)
        return summary

    host_settings = getattr(ctx.config.host, "settings", None) or getattr(ctx.config.host, "raw", None) or {}
    engine = ((host_settings.get("tools") or {}).get("gnina_settings", {})) if isinstance(host_settings, dict) else {}
    cluster_cutoff = float(getattr(ctx.config.budget, "cluster_rmsd_A", 2.0) or 2.0)
    workdir = ensure_dir(os.path.join(ctx.stage_dir("redock"), "runs"))
    support_note = (
        "a seed supports the imported pose when at least one of its modes is within "
        f"{cluster_cutoff} A symmetry-corrected RMSD in the unchanged receptor frame of the imported pose"
    )

    for item in plan:
        state_id = item["state_id"]
        compound_id = item.get("compound_id")
        ligand_path = item.get("ligand_path")
        box = item.get("box") or {}
        seeds = [int(seed) for seed in item.get("seeds", [])]
        record: dict[str, Any] = {
            "compound_id": compound_id,
            "reference_pose_ids": item.get("pose_ids", []),
            "state_id": state_id,
            "seeds_planned": len(seeds),
            "seeds_attempted": 0,
            "seeds_succeeded": 0,
            "supporting_seeds_same_cluster": 0,
            "cluster_rmsd_A": cluster_cutoff,
            "support_note": support_note,
            "per_seed": [],
        }
        summary["compounds"][compound_id or state_id] = record
        if not ligand_path or not os.path.exists(ligand_path) or not box.get("center"):
            record["status"] = "skipped_missing_input"
            record["detail"] = "ligand coordinates or docking box are missing"
            summary["skipped"] += 1
            item["status"] = "skipped_missing_input"
            continue
        reference = _reference_mol(ligand_path)
        for seed in seeds:
            record["seeds_attempted"] += 1
            prefix = os.path.join(workdir, f"{state_id}_seed{seed}")
            out_path = prefix + ".sdf"
            log_path = prefix + ".log"
            job_path = prefix + ".json"
            if os.path.exists(out_path) and os.path.exists(job_path):
                job = _read_json(job_path)
                seed_record = _seed_record(
                    out_path, reference, seed, cluster_cutoff, cached=True, job=job, log_path=log_path
                )
                record["per_seed"].append(seed_record)
                if seed_record["status"] == "ok" and seed_record.get("returncode") == 0:
                    record["supporting_seeds_same_cluster"] += int(seed_record.get("supports_imported_pose", False))
                    record["seeds_succeeded"] += 1
                    summary["completed"] += 1
                else:
                    summary["failed"] += 1
                continue
            argv = redock_argv(
                command=capability.command,
                receptor=receptor,
                ligand=ligand_path,
                center=box["center"],
                size=box["size"],
                seed=seed,
                out_prefix=out_path,
                exhaustiveness=int(engine.get("exhaustiveness", 8)),
                num_modes=int(engine.get("num_modes", 3)),
                scoring=str(engine.get("scoring", "vina")),
                cnn=str(engine.get("cnn_scoring", "rescore")),
                cpu=int(engine.get("cpu", 1)),
                cnn_model=engine.get("cnn_model"),
                device=engine.get("device"),
            )
            started = time.time()
            try:
                probe = subprocess.run(
                    argv, capture_output=True, text=True, timeout=int(engine.get("timeout_seconds", 3600)), env=_engine_environment(ctx)
                )
                payload = {
                    "argv": argv,
                    "returncode": probe.returncode,
                    "seconds": round(time.time() - started, 3),
                    "executed_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
                    "stdout_tail": probe.stdout[-8000:],
                    "stderr_tail": probe.stderr[-4000:],
                }
                atomic_write_json(job_path, payload)
                with open(log_path, "w", encoding="utf-8", newline="\n") as handle:
                    handle.write(f"$ {' '.join(argv)}\n")
                    handle.write(f"# exit={probe.returncode} seconds={payload['seconds']}\n")
                    handle.write(probe.stdout)
                    handle.write("\n[stderr]\n")
                    handle.write(probe.stderr)
            except Exception as exc:  # noqa: BLE001
                summary["failed"] += 1
                record["per_seed"].append(
                    {"seed": seed, "status": "error", "error": f"{type(exc).__name__}: {exc}"}
                )
                continue
            seed_record = _seed_record(
                out_path, reference, seed, cluster_cutoff, cached=False, job=_read_json(job_path), log_path=log_path
            )
            seed_record["returncode"] = int(_read_json(job_path).get("returncode", -1))
            record["per_seed"].append(seed_record)
            if seed_record["status"] == "ok" and seed_record.get("returncode") == 0:
                record["seeds_succeeded"] += 1
                summary["completed"] += 1
                if seed_record.get("supports_imported_pose"):
                    record["supporting_seeds_same_cluster"] += 1
            else:
                summary["failed"] += 1
        record["status"] = (
            "ok" if record["seeds_succeeded"] == record["seeds_planned"] > 0 else "partial_or_failed"
        )
        _record_compound_evidence(ctx, record)

    summary["status"] = "completed"
    summary["scored_compounds"] = len(summary["compounds"])
    atomic_write_json(os.path.join(ctx.stage_dir("redock"), "summary.json"), summary)
    return summary


def _record_compound_evidence(ctx, record: dict[str, Any]) -> None:
    compound_id = record.get("compound_id")
    if not compound_id:
        return
    attempts = record["seeds_attempted"]
    successes = record["seeds_succeeded"]
    support = record["supporting_seeds_same_cluster"]
    planned = record["seeds_planned"]
    ctx.db.record_assessment(
        "compound", compound_id, "redock", "seeds_planned", float(planned),
        AssessmentStatus.PASS if planned else AssessmentStatus.NOT_EVALUATED,
    )
    ctx.db.record_assessment(
        "compound", compound_id, "redock", "seeds_attempted", float(attempts),
        AssessmentStatus.PASS if attempts else AssessmentStatus.NOT_EVALUATED,
    )
    ctx.db.record_assessment(
        "compound", compound_id, "redock", "seeds_succeeded", float(successes),
        AssessmentStatus.PASS if successes else AssessmentStatus.FAIL,
        None if successes else Reason.REDOCK_SEED_INCOMPLETE,
    )
    ctx.db.record_assessment(
        "compound", compound_id, "redock", "supporting_seeds_same_cluster", float(support),
        AssessmentStatus.PASS if successes else AssessmentStatus.NOT_EVALUATED,
    )


def _reference_mol(path: str):
    from ..chem import rdio

    for mol, _props in rdio.iter_sdf(path):
        if mol is not None and mol.GetNumConformers():
            return mol
    return None


def _seed_record(
    out_path: str,
    reference,
    seed: int,
    cluster_cutoff: float,
    *,
    cached: bool,
    job: dict[str, Any] | None = None,
    log_path: str | None = None,
) -> dict[str, Any]:
    from ..chem import rdio

    if not os.path.exists(out_path):
        return {
            "seed": seed,
            "status": "no_output",
            "error": "GNINA produced no SDF file",
            "returncode": int((job or {}).get("returncode", -1)),
        }
    best_rmsd = None
    scores: list[dict[str, Any]] = []
    for mol, props in rdio.iter_sdf(out_path):
        if mol is None:
            continue
        scores.append(_normalise(props))
        if reference is not None and mol.GetNumConformers():
            rmsd = _best_rmsd(mol, reference)
            if rmsd is not None and (best_rmsd is None or rmsd < best_rmsd):
                best_rmsd = rmsd
    payload = {
        "seed": seed,
        "status": "ok" if scores and best_rmsd is not None else "comparison_unavailable" if scores else "empty",
        "modes": len(scores),
        "min_rmsd_to_imported_A": round(best_rmsd, 3) if best_rmsd is not None else None,
        "supports_imported_pose": bool(best_rmsd is not None and best_rmsd <= cluster_cutoff),
        "cached": cached,
        "file": out_path,
        "sha256": __import__("hashlib").sha256(open(out_path, "rb").read()).hexdigest(),
        "best": scores[0] if scores else None,
        "returncode": int((job or {}).get("returncode", 0)),
    }
    if job:
        payload["argv"] = job.get("argv")
        payload["seconds"] = job.get("seconds")
    if log_path and os.path.exists(log_path):
        payload["log_head"] = read_text(log_path)[:600]
    return payload


def _best_rmsd(mol, reference) -> float | None:
    """Symmetry-corrected RMSD in the *receptor* frame.

    ``GetBestRMS`` would superpose the ligand onto the reference and report ~0
    for any two poses of the same molecule, so it cannot answer "same binding
    site?".  ``CalcRMS`` keeps the coordinates fixed and only permutes
    symmetry-equivalent atoms, which is the comparison the specification asks
    for.
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    try:
        probe, target = Chem.RemoveHs(mol), Chem.RemoveHs(reference)
        if probe.GetNumAtoms() != target.GetNumAtoms():
            return None
        if Chem.MolToSmiles(probe, isomericSmiles=True) != Chem.MolToSmiles(target, isomericSmiles=True):
            return None
        # Let RDKit enumerate equivalent atom mappings. CalcRMS takes a list
        # of mappings, not the flat atom-pair list formerly passed here.
        return float(rdMolAlign.CalcRMS(probe, target, symmetrizeConjugatedTerminalGroups=True))
    except (ValueError, RuntimeError):
        return None


def _read_json(path: str) -> dict[str, Any]:
    from ..util import read_json

    return read_json(path, {}) or {}


def _receptor_path(ctx) -> str | None:
    from .. import data

    record = data.receptor_record(ctx)
    return record.get("archived") or record.get("path")


def _normalise(props: dict[str, Any]) -> dict[str, Any]:
    from ..ingest.poses import normalise_scores

    return normalise_scores({k: v for k, v in props.items() if not k.startswith("_")})[0]


def summarise_seeds(results: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Aggregate per-seed best scores into a reproducibility summary per state."""
    out: dict[str, dict[str, Any]] = {}
    for item in results:
        if item.get("status") != "ok":
            continue
        state = item["state_id"]
        entry = out.setdefault(state, {"seeds": [], "affinities": [], "cnn": [], "poses": 0})
        entry["seeds"].append(item["seed"])
        entry["poses"] += item.get("poses", 0)
        if item.get("best"):
            affinity = safe_float(item["best"].get("minimizedAffinity"))
            cnn = safe_float(item["best"].get("CNNscore"))
            if affinity is not None:
                entry["affinities"].append(affinity)
            if cnn is not None:
                entry["cnn"].append(cnn)
    for state, entry in out.items():
        affinities = entry["affinities"]
        entry["seed_count"] = len(entry["seeds"])
        entry["median_affinity"] = _median(affinities)
        entry["std_affinity"] = _std(affinities)
        entry["consistency"] = round(1.0 / (1.0 + (entry["std_affinity"] or 0.0)), 4)
    return out


def _median(values: list[float]) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    mid = len(ordered) // 2
    return ordered[mid] if len(ordered) % 2 else (ordered[mid - 1] + ordered[mid]) / 2.0


def _std(values: list[float]) -> float | None:
    if len(values) < 2:
        return None
    mean = sum(values) / len(values)
    return (sum((v - mean) ** 2 for v in values) / (len(values) - 1)) ** 0.5
