#!/usr/bin/env python3
"""Command-line entry point for the isolated MurJ v2 pipeline.

Every stage writes ``stage_success.json`` only after it succeeds. Model stages
fail closed: no legacy molecule or hand-written structure is ever substituted
for a failed Peptide2Mol / REINVENT4 run.
"""
from __future__ import annotations

import argparse
import re
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

# Conda environments contain absolute prefixes and symlinks, so they stay on the
# internal disk (the project may live on an exFAT drive). Override with MURJ_ENV_ROOT.
import os
ENV_ROOT = Path(os.environ.get("MURJ_ENV_ROOT",
                               Path.home() / ".local/share/yuanhe/envs"))

STAGES = ["prepare", "hypotheses", "generate-seeds", "select-cores", "expand",
          "synthesis-check", "geometry-check", "select", "export", "validate"]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _stage_dir(run_id: str, stage: str) -> Path:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,79}", run_id):
        raise ValueError("unsafe run ID")
    return ROOT / "runs" / run_id / stage


def _mark(run_id: str, stage: str, outputs: list[Path], extra: dict | None = None) -> None:
    missing = [str(p) for p in outputs if not p.exists()]
    if missing:
        raise RuntimeError(f"stage {stage} missing outputs: {missing}")
    d = _stage_dir(run_id, stage)
    d.mkdir(parents=True, exist_ok=True)
    payload = {
        "stage": stage, "status": "success",
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "outputs": [{"path": str(p), "sha256": _sha256(p)} for p in outputs if p.is_file()],
        **(extra or {}),
    }
    (d / "stage_success.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def _require(run_id: str, stage: str) -> None:
    marker = _stage_dir(run_id, stage) / "stage_success.json"
    if not marker.is_file():
        raise RuntimeError(f"required stage '{stage}' has not succeeded for run {run_id}")
    record = json.loads(marker.read_text(encoding="utf-8"))
    if record.get("status") != "success" or not record.get("outputs"):
        raise RuntimeError("invalid success marker")
    for item in record["outputs"]:
        path = Path(item["path"])
        if not path.is_file() or _sha256(path) != item["sha256"]:
            raise RuntimeError(f"upstream output changed or missing: {path}")


def cmd_prepare(args, cfg):
    from src.structure import prepare_structure
    ps = cfg["primary_structure"]
    pdb = (ROOT / ps["path"]).resolve()
    if _sha256(pdb) != ps["sha256"]:
        raise RuntimeError("primary PDB SHA-256 differs from configs/project.json")
    out = _stage_dir(args.run_id, "prepare")
    res = prepare_structure(
        pdb, out, structure_id=ps["structure_id"], receptor_chain=ps["receptor_chain"],
        peptide_chain=ps["peptide_chain"], receptor_start=ps["receptor_target_range"][0],
        receptor_end=ps["receptor_target_range"][1], crop_distance_A=float(cfg["pocket"]["radius_A"]),
    )
    paths = [Path(res[k]) for k in ("receptor_full", "receptor_prepared", "peptide_reference", "structures_manifest",
                                    "residue_map", "preparation_report", "checksums")]
    _mark(args.run_id, "prepare", [p if p.is_absolute() else ROOT / p if (ROOT / p).exists() else p for p in paths],
          {"source_sha256": ps["sha256"]})


def cmd_hypotheses(args, cfg):
    from src.hypotheses import generate_hypotheses
    _require(args.run_id, "prepare")
    ps = cfg["primary_structure"]
    out = _stage_dir(args.run_id, "hypotheses") / "hypotheses.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    generate_hypotheses((ROOT / ps["path"]).resolve(), out, structure_id=ps["structure_id"],
                        receptor_chain=ps["receptor_chain"], peptide_chain=ps["peptide_chain"],
                        feature_radius_A=float(cfg["geometry"]["feature_radius_A"]))
    _mark(args.run_id, "hypotheses", [out])


def cmd_generate(args, cfg):
    from src.adapters import peptide2mol
    _require(args.run_id, "hypotheses")
    mc = cfg["models"]["peptide2mol"]
    pocket = _stage_dir(args.run_id, "prepare") / "receptor_prepared.pdb"
    exe = Path(args.executable or ENV_ROOT / "peptide2mol/bin/python")
    manifests = []
    for seed in mc["seeds"]:
        out = _stage_dir(args.run_id, "generate-seeds") / f"seed_{seed}"
        manifests.append(peptide2mol.run(
            pocket_pdb=pocket, checkpoint=ROOT / mc["checkpoint"], checkpoint_sha256=mc["checkpoint_sha256"],
            output_dir=out, budget=args.budget or mc["sample_budget"], seed=seed, executable=exe,
            repo=ROOT / "vendor/Peptide2Mol"))
    _mark(args.run_id, "generate-seeds",
          [_stage_dir(args.run_id, "generate-seeds") / f"seed_{s}" / "manifest.json" for s in mc["seeds"]],
          {"real_model": "Peptide2Mol", "reports": [m["report"] for m in manifests]})


def cmd_expand(args, cfg):
    from src.adapters import reinvent4
    mc = cfg["models"]["reinvent4"]
    exe = Path(args.executable or ENV_ROOT / "reinvent4/bin/reinvent")
    _require(args.run_id, "generate-seeds")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,79}", args.name):
        raise ValueError("unsafe expansion name")
    out = _stage_dir(args.run_id, "expand") / args.name
    m = reinvent4.run(scaffold_smiles=args.scaffold, checkpoint=ROOT / mc["prior"],
                      checkpoint_sha256=mc["prior_sha256"], output_dir=out,
                      budget=args.budget or mc["sample_budget_per_core"], seed=mc["seed"],
                      executable=exe, repo=ROOT / "vendor/REINVENT4", device=mc["device"])
    _mark(args.run_id, "expand", [out / "manifest.json", out / "candidates.json"],
          {"real_model": "REINVENT4/LibInvent", "report": m["report"]})


def cmd_validate(args, cfg):
    from src.export import validate_export
    report = validate_export(Path(args.export_dir))
    if not report.get("valid", report.get("ok")):
        raise RuntimeError(json.dumps(report, ensure_ascii=False, default=str))
    d = _stage_dir(args.run_id, "validate"); d.mkdir(parents=True, exist_ok=True)
    p = d / "validation.json"
    p.write_text(json.dumps(report, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    _mark(args.run_id, "validate", [p])


def cmd_status(args, cfg):
    base = ROOT / "runs" / args.run_id
    st = {s: (base / s / "stage_success.json").exists() for s in STAGES}
    st["platform"] = {"python": sys.version.split()[0], "system": platform.platform(), "machine": platform.machine()}
    print(json.dumps(st, indent=2, ensure_ascii=False))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", default=str(ROOT / "configs/project.json"))
    p.add_argument("--run-id", default="murj_9nu4_pilot")
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("prepare"); sub.add_parser("hypotheses"); sub.add_parser("status")
    g = sub.add_parser("generate-seeds"); g.add_argument("--budget", type=int); g.add_argument("--executable")
    e = sub.add_parser("expand"); e.add_argument("--scaffold", required=True); e.add_argument("--name", required=True)
    e.add_argument("--budget", type=int); e.add_argument("--executable")
    v = sub.add_parser("validate"); v.add_argument("--export-dir", required=True)
    args = p.parse_args()
    _stage_dir(args.run_id, "validate-id")
    cfg = json.loads(Path(args.config).read_text(encoding="utf-8"))
    {"prepare": cmd_prepare, "hypotheses": cmd_hypotheses, "generate-seeds": cmd_generate,
     "expand": cmd_expand, "validate": cmd_validate, "status": cmd_status}[args.command](args, cfg)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # fail closed, non-zero exit
        print(f"FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
        sys.exit(1)
