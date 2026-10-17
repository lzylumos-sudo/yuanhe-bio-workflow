"""Strict adapter for the pinned Peptide2Mol checkout.

This module orchestrates the official preparation script and inference shell script.
It never substitutes a legacy implementation or fabricates molecules.
"""
from __future__ import annotations

import argparse
import hashlib
import re
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

REPO_COMMIT = "2c5af9628aeeab5ed2dff5dc5326ce2463316291"
SUPPORTED_ELEMENTS = {"C", "N", "O", "F", "P", "S", "CL", "BR"}


def hydra_string(value) -> str:
    # Quotes belong to Hydra's grammar, not a shell. Escape backslashes/quotes.
    return json.dumps(str(value), ensure_ascii=False)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def platform_readiness(repo: str | Path, executable: str | Path) -> dict[str, Any]:
    """Return explicit readiness information; no fallback is attempted."""
    repo_path, exe_path = Path(repo).resolve(), Path(executable).absolute()
    prepare = repo_path / "notebooks" / "deal_with_mol_test_from_pdb_pal_5A.py"
    inference = repo_path / "src" / "eval.py"
    issues: list[str] = []
    if platform.system() != "Linux":
        issues.append("Peptide2Mol's released CUDA environment is supported on Linux only")
    if not repo_path.is_dir():
        issues.append(f"vendor repository missing: {repo_path}")
    if not prepare.is_file():
        issues.append(f"official preparation script missing: {prepare}")
    if not inference.is_file():
        issues.append(f"official inference entry point missing: {inference}")
    if not exe_path.is_file() or not os.access(exe_path, os.X_OK):
        issues.append(f"Python executable missing or not executable: {exe_path}")
    return {
        "ready": not issues,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "issues": issues,
        "legacy_fallback": False,
    }


def _pdb_elements(path: Path) -> set[str]:
    elements: set[str] = set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element:
            atom_name = line[12:16].strip().lstrip("0123456789")
            element = (atom_name[:2] if atom_name[:2].upper() == "CL" else atom_name[:1]).upper()
        if element and element != "H":
            elements.add(element)
    return elements


def _run(command: list[str], cwd: Path, env: dict[str, str] | None = None) -> dict[str, Any]:
    import threading
    started = _now()
    environment = {**os.environ, **(env or {}), "PYTHONUNBUFFERED": "1"}
    proc = subprocess.Popen(command, cwd=cwd, env=environment, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1)
    captured = {"stdout": [], "stderr": []}
    def consume(name, stream):
        for line in stream:
            captured[name].append(line)
            print(f"[Peptide2Mol {proc.pid} {name}] {line}", end="", flush=True)
    readers=[threading.Thread(target=consume,args=(name,getattr(proc,name)),daemon=True) for name in captured]
    for t in readers:t.start()
    try:proc.wait(timeout=86400)
    except subprocess.TimeoutExpired:
        proc.kill();proc.wait();raise
    finally:
        for t in readers:t.join()
    return {"command":command,"cwd":str(cwd),"started_at":started,"finished_at":_now(),
            "returncode":proc.returncode,**{k:"".join(v) for k,v in captured.items()}}


def _rdkit():
    try:
        from rdkit import Chem
    except ImportError as exc:  # pragma: no cover - environment-dependent
        raise RuntimeError("RDKit is required to validate Peptide2Mol SDF output") from exc
    return Chem


def _collect_sdfs(search_root: Path, archive: Path, requested: int) -> tuple[list[dict[str, Any]], dict[str, int]]:
    Chem = _rdkit()
    archive.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    reasons: dict[str, int] = {}
    seen: set[str] = set()
    candidates = [p for p in search_root.rglob("*.sdf") if archive not in p.parents]
    candidates = [p for p in candidates if "origin" not in p.stem.lower()]
    for index, source in enumerate(sorted(candidates)):
        target = archive / f"{index:05d}_{source.name}"
        shutil.copy2(source, target)
        mol = next((m for m in Chem.SDMolSupplier(str(source), removeHs=False) if m is not None), None)
        reason = None
        canonical = None
        if mol is None:
            reason = "sdf_parse_failed"
        else:
            try:
                Chem.SanitizeMol(mol)
                if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
                    raise ValueError("disconnected or dummy-containing molecule")
                if not mol.GetNumConformers():
                    raise ValueError("missing coordinates")
                import math
                if not all(math.isfinite(float(x)) for xyz in mol.GetConformer().GetPositions() for x in xyz):
                    raise ValueError("nonfinite coordinates")
                canonical = Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True)
            except Exception:
                reason = "sanitization_failed"
        if canonical and canonical in seen:
            reason = "duplicate"
        elif canonical:
            seen.add(canonical)
        if reason:
            reasons[reason] = reasons.get(reason, 0) + 1
        rows.append({"raw_sdf": str(target), "smiles": canonical, "valid": reason is None, "failure_reason": reason})
        if len(rows) >= requested:
            break
    if len(rows) < requested:
        reasons["model_produced_fewer_than_requested"] = requested - len(rows)
    return rows, reasons


def run(
    *, pocket_pdb: str | Path, checkpoint: str | Path, checkpoint_sha256: str,
    output_dir: str | Path, budget: int, seed: int, executable: str | Path = sys.executable,
    repo: str | Path | None = None, batch_size: int = 2,
) -> dict[str, Any]:
    """Prepare a pocket, invoke official inference, validate and archive raw SDFs."""
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    if budget <= 0:
        raise ValueError("budget must be positive")
    root = Path(__file__).resolve().parents[2]
    repo_path = Path(repo or root / "vendor" / "Peptide2Mol").resolve()
    out = Path(output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise FileExistsError("model output directory must be empty; use a new attempt directory")
    manifest_path = out / "manifest.json"
    pocket, ckpt, exe = Path(pocket_pdb).resolve(), Path(checkpoint).resolve(), Path(executable).absolute()
    manifest: dict[str, Any] = {
        "adapter": "Peptide2Mol", "repo": str(repo_path), "repo_commit": REPO_COMMIT,
        "executable": str(exe), "checkpoint": {"path": str(ckpt), "required_sha256": checkpoint_sha256},
        "inputs": {"pocket_pdb": str(pocket)}, "config": {"batch_size": batch_size},
        "seed": seed, "budget": budget, "started_at": _now(), "commands": [], "status": "running",
    }
    try:
        if len(checkpoint_sha256) != 64:
            raise ValueError("checkpoint SHA256 must be exactly 64 hexadecimal characters")
        int(checkpoint_sha256, 16)
        if not pocket.is_file():
            raise FileNotFoundError(f"pocket PDB missing: {pocket}")
        if not ckpt.is_file():
            raise FileNotFoundError(f"checkpoint missing: {ckpt}")
        actual_hash = _sha256(ckpt)
        manifest["checkpoint"]["actual_sha256"] = actual_hash
        if actual_hash.lower() != checkpoint_sha256.lower():
            raise ValueError(f"checkpoint SHA256 mismatch: expected {checkpoint_sha256}, got {actual_hash}")
        readiness = platform_readiness(repo_path, exe)
        manifest["platform_readiness"] = readiness
        manifest["input_sha256"] = _sha256(pocket)
        from .provenance import source_fingerprint
        manifest["source_provenance"] = source_fingerprint(repo_path, REPO_COMMIT)
        if not readiness["ready"]:
            raise RuntimeError("platform not ready: " + "; ".join(readiness["issues"]))
        unsupported = sorted(_pdb_elements(pocket) - SUPPORTED_ELEMENTS)
        if unsupported:
            raise ValueError("unsupported pocket elements: " + ", ".join(unsupported))
        input_dir = out / "prepared_input"
        input_dir.mkdir(exist_ok=True)
        prepared_pdb = input_dir / (pocket.stem if "poc" in pocket.stem.lower() else pocket.stem + "_poc")
        prepared_pdb = prepared_pdb.with_suffix(".pdb")
        shutil.copy2(pocket, prepared_pdb)
        prepare_command = [str(exe), str(repo_path / "notebooks" / "deal_with_mol_test_from_pdb_pal_5A.py"), str(input_dir)]
        prepare_result = _run(prepare_command, repo_path)
        manifest["commands"].append(prepare_result)
        if prepare_result["returncode"]:
            raise RuntimeError(f"Peptide2Mol preparation failed with exit code {prepare_result['returncode']}")
        pt = input_dir / f"{prepared_pdb.stem}.pt"
        if not pt.is_file():
            raise RuntimeError(f"preparation completed without expected file: {pt}")
        raw_root = out / "model_output"
        inference_command = [
            str(exe), str(repo_path / "src" / "eval.py"),
            "experiment=mol_test", "ckpt_path=" + hydra_string(ckpt), "++paths.data_dir=" + hydra_string(input_dir),
            "+data.lmdb_fn=" + hydra_string(pt.name), "data=mol_test_true", "model=Moldiff_test",
            f"data.infer_batch_size={batch_size}", "trainer.devices=1", "++paths.log_dir=" + hydra_string(out / "logs"),
            "++model.net.sample.log_dir=" + hydra_string(raw_root), "++model.net.sample.pdb_dir=" + hydra_string(input_dir),
            "++model.net.sample.batch_size=1", f"++model.net.sample.num_mols={budget}", f"++seed={seed}",
        ]
        env = os.environ.copy()
        env.update({"PYTHONHASHSEED": str(seed), "PL_GLOBAL_SEED": str(seed)})
        inference_result = _run(inference_command, repo_path, env)
        manifest["commands"].append(inference_result)
        if inference_result["returncode"]:
            raise RuntimeError(f"Peptide2Mol inference failed with exit code {inference_result['returncode']}")
        generated_sdf_dir = Path(str(raw_root) + "_SDF")
        if not generated_sdf_dir.is_dir():
            raise RuntimeError(f"Peptide2Mol completed without expected SDF directory: {generated_sdf_dir}")
        rows, reasons = _collect_sdfs(generated_sdf_dir, out / "raw_sdf", budget)
        valid = [row for row in rows if row["valid"]]
        if not valid:
            raise RuntimeError("Peptide2Mol produced zero valid unique molecules")
        report = {"requested": budget, "completed": len(rows), "valid": len(valid), "unique": len({r["smiles"] for r in valid}), "failure_reasons": reasons}
        _write_json(out / "candidates.json", {"candidates": rows, "report": report})
        manifest.update({"report": report, "status": "completed"})
    except Exception as exc:
        manifest.update({"status": "failed", "error": f"{type(exc).__name__}: {exc}"})
        raise
    finally:
        manifest["finished_at"] = _now()
        _write_json(manifest_path, manifest)
    return manifest


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the real Peptide2Mol adapter")
    parser.add_argument("--pocket-pdb", required=True); parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--checkpoint-sha256", required=True); parser.add_argument("--output-dir", required=True)
    parser.add_argument("--budget", required=True, type=int); parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--executable", default=sys.executable); parser.add_argument("--repo"); parser.add_argument("--batch-size", type=int, default=2)
    args = parser.parse_args(argv)
    try:
        run(**vars(args))
    except Exception as exc:
        print(str(exc), file=sys.stderr); return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
