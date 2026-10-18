"""Strict adapter for sampling with the pinned official REINVENT4 CLI."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

REPO_COMMIT = "3cf206ea5b385ac955dc44bbae6063b55ee31e0a"
CONFIG_VERSION = 4


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_json(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _toml_string(value: str | Path) -> str:
    return json.dumps(str(value))


def _rdkit():
    try:
        from rdkit import Chem
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("RDKit is required to validate REINVENT4 chemistry") from exc
    return Chem


def validate_scaffold(smiles: str) -> dict[str, Any]:
    """Require one uniquely labelled dummy atom in a parseable scaffold."""
    Chem = _rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("scaffold is not valid SMILES")
    dummies = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
    if not dummies:
        raise ValueError("scaffold must contain an attachment dummy atom")
    labels = [atom.GetAtomMapNum() or atom.GetIsotope() for atom in dummies]
    if any(label <= 0 for label in labels):
        raise ValueError("every scaffold dummy must have an isotope or atom-map attachment label")
    if len(labels) != len(set(labels)):
        raise ValueError("scaffold attachment labels must be unique")
    if len(dummies) != 1:
        raise ValueError("LibInvent sampling requires exactly one attachment dummy")
    core = Chem.DeleteSubstructs(mol, Chem.MolFromSmarts("[#0]"))
    if core.GetNumAtoms() == 0:
        raise ValueError("scaffold has no non-dummy core")
    return {"canonical": Chem.MolToSmiles(mol, canonical=True), "labels": labels, "core": Chem.MolToSmiles(core, canonical=True)}


def render_sampling_toml(*, checkpoint: Path, scaffold_file: Path, output_csv: Path, budget: int, device: str) -> str:
    """Generate a REINVENT4 v4 sampling configuration for LibInvent."""
    return "\n".join([
        'run_type = "sampling"', f'device = {_toml_string(device)}', "",
        "[parameters]", f'model_file = {_toml_string(checkpoint)}', f'smiles_file = {_toml_string(scaffold_file)}',
        f'output_file = {_toml_string(output_csv)}', f"num_smiles = {budget}",
        "unique_molecules = false", "randomize_smiles = false", "",])


def _bracket_dummy_smiles(smiles):
    """Rewrite each attachment dummy to the single bracket token "[*]"."""
    chem = _rdkit()
    mol = chem.MolFromSmiles(smiles)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetAtomMapNum(0)
            atom.SetIsotope(0)
    return chem.MolToSmiles(mol, isomericSmiles=True, canonical=False).replace("*", "[*]")


def scaffold_identity(smiles):
    Chem = _rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
        if atom.GetAtomicNum() == 0:
            atom.SetIsotope(0)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _core_retained(scaffold: str, candidate: str) -> tuple[bool, str | None]:
    Chem = _rdkit()
    candidate_mol = Chem.MolFromSmiles(candidate)
    if candidate_mol is None:
        return False, "invalid_smiles"
    if len(Chem.GetMolFrags(candidate_mol)) != 1:
        return False, "disconnected_molecule"
    if any(atom.GetAtomicNum() == 0 for atom in candidate_mol.GetAtoms()):
        return False, "attachment_not_consumed"
    scaffold_mol = Chem.MolFromSmiles(scaffold)
    core = Chem.DeleteSubstructs(scaffold_mol, Chem.MolFromSmarts("[#0]"))
    if not candidate_mol.HasSubstructMatch(core, useChirality=True):
        return False, "core_not_retained"
    return True, None


def run(
    *, scaffold_smiles: str, checkpoint: str | Path, checkpoint_sha256: str,
    output_dir: str | Path, budget: int, seed: int, executable: str | Path,
    repo: str | Path | None = None, device: str = "cpu",
) -> dict[str, Any]:
    if budget <= 0:
        raise ValueError("budget must be positive")
    root = Path(__file__).resolve().parents[2]
    repo_path = Path(repo or root / "vendor" / "REINVENT4").resolve()
    out = Path(output_dir).resolve(); out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise FileExistsError("model output directory must be empty; use a new attempt directory")
    manifest_path = out / "manifest.json"
    ckpt, exe = Path(checkpoint).resolve(), Path(executable).resolve()
    manifest: dict[str, Any] = {
        "adapter": "REINVENT4", "repo": str(repo_path), "repo_commit": REPO_COMMIT,
        "executable": str(exe), "checkpoint": {"path": str(ckpt), "required_sha256": checkpoint_sha256},
        "inputs": {"scaffold_smiles": scaffold_smiles}, "config": {"version": CONFIG_VERSION, "device": device},
        "seed": seed, "budget": budget, "started_at": _now(), "commands": [], "status": "running",
    }
    try:
        if len(checkpoint_sha256) != 64:
            raise ValueError("checkpoint SHA256 must be exactly 64 hexadecimal characters")
        int(checkpoint_sha256, 16)
        if not ckpt.is_file(): raise FileNotFoundError(f"checkpoint missing: {ckpt}")
        actual_hash = _sha256(ckpt); manifest["checkpoint"]["actual_sha256"] = actual_hash
        if actual_hash.lower() != checkpoint_sha256.lower():
            raise ValueError(f"checkpoint SHA256 mismatch: expected {checkpoint_sha256}, got {actual_hash}")
        if not repo_path.is_dir(): raise FileNotFoundError(f"vendor repository missing: {repo_path}")
        if not exe.is_file() or not os.access(exe, os.X_OK):
            raise RuntimeError(f"REINVENT executable missing or not executable: {exe}")
        from .provenance import source_fingerprint
        manifest["source_provenance"] = source_fingerprint(repo_path, REPO_COMMIT)
        scaffold = validate_scaffold(scaffold_smiles); manifest["inputs"]["validated_scaffold"] = scaffold
        # The CLI only knows the bracket token "[*]" (chemistry/tokens.py +
        # libinvent vocabulary), so the dummy is written without map numbers;
        # the original traversal order is kept.
        scaffold_for_model = _bracket_dummy_smiles(scaffold_smiles)
        manifest["inputs"]["model_scaffold"] = scaffold_for_model
        scaffold_file, raw_csv = out / "scaffold.smi", out / "raw_sampling.csv"
        scaffold_file.write_text(scaffold_for_model + "\n", encoding="utf-8")
        config = out / "sampling.toml"
        config.write_text(render_sampling_toml(checkpoint=ckpt, scaffold_file=scaffold_file, output_csv=raw_csv, budget=budget, device=device), encoding="utf-8")
        command = [str(exe), "--seed", str(seed), str(config)]
        started = _now(); completed = subprocess.run(command, cwd=repo_path, text=True, capture_output=True, check=False, timeout=86400)
        manifest["commands"].append({"command": command, "cwd": str(repo_path), "started_at": started, "finished_at": _now(), "returncode": completed.returncode, "stdout": completed.stdout, "stderr": completed.stderr})
        if completed.returncode: raise RuntimeError(f"REINVENT4 failed with exit code {completed.returncode}")
        if not raw_csv.is_file(): raise RuntimeError("REINVENT4 completed without its configured CSV output")
        archive = out / "raw"; archive.mkdir(exist_ok=True); shutil.copy2(raw_csv, archive / raw_csv.name)
        with raw_csv.open(newline="", encoding="utf-8") as handle: rows = list(csv.DictReader(handle))
        if not rows: raise RuntimeError("REINVENT4 CSV contains no samples")
        required = {"SMILES", "NLL"}
        if not required.issubset(rows[0]): raise RuntimeError("REINVENT4 CSV must contain SMILES and NLL columns")
        input_column = "Input_SMILES" if "Input_SMILES" in rows[0] else "Scaffold" if "Scaffold" in rows[0] else None
        if input_column is None: raise RuntimeError("REINVENT4 CSV must preserve Input_SMILES or Scaffold")
        accepted, rejected, seen = [], [], set()
        for raw in rows:
            preserved = {"SMILES": raw.get("SMILES"), "Input_SMILES": raw.get(input_column), "NLL": raw.get("NLL")}
            if scaffold_identity(raw.get(input_column, "")) != scaffold_identity(scaffold_smiles):
                preserved["failure_reason"] = "input_scaffold_mismatch"; rejected.append(preserved); continue
            ok, reason = _core_retained(scaffold_smiles, raw.get("SMILES", ""))
            if not ok: preserved["failure_reason"] = reason; rejected.append(preserved); continue
            canonical = _rdkit().MolToSmiles(_rdkit().MolFromSmiles(raw["SMILES"]), canonical=True)
            if canonical in seen: preserved["failure_reason"] = "duplicate"; rejected.append(preserved); continue
            seen.add(canonical); preserved["canonical_smiles"] = canonical; accepted.append(preserved)
        if not accepted: raise RuntimeError("REINVENT4 produced no valid unique molecules retaining the scaffold core")
        report = {"requested": budget, "completed": len(rows), "valid": len(accepted), "unique": len(accepted), "rejected": len(rejected)}
        _write_json(out / "candidates.json", {"candidates": accepted, "rejected": rejected, "report": report})
        manifest.update({"report": report, "raw_output": str(raw_csv), "status": "completed"})
    except Exception as exc:
        manifest.update({"status": "failed", "error": f"{type(exc).__name__}: {exc}"}); raise
    finally:
        manifest["finished_at"] = _now(); _write_json(manifest_path, manifest)
    return manifest


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the real REINVENT4 adapter")
    parser.add_argument("--scaffold-smiles", required=True); parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--checkpoint-sha256", required=True); parser.add_argument("--output-dir", required=True)
    parser.add_argument("--budget", required=True, type=int); parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--executable", required=True); parser.add_argument("--repo"); parser.add_argument("--device", default="cpu")
    args = parser.parse_args(argv)
    try: run(**vars(args))
    except Exception as exc: print(str(exc), file=sys.stderr); return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
