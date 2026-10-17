"""GNINA job adapter (plan §17-§18). One job = one receptor revision x one ligand conformer x one seed, own directory.

* CNN is never disabled: cnn_scoring=none is refused for a GNINA protocol (use the Vina adapter for Vina).
* device 'gpu' requires a visible GPU; 'cpu' adds --no_gpu. There is no silent fallback: a GPU request that
  cannot be satisfied FAILS and says so (Protocol.allow_cpu_fallback only permits an *explicit, recorded* retry by the caller).
* Scores are read ONLY from the engine's own SDF output properties. Missing/non-finite => status 'failed'.
"""
from __future__ import annotations
import json, math, os, re, subprocess, time
from pathlib import Path

from ..hashing import atomic_write_json, sha256_file, sha256_bytes
from .execute import run_argv, Cancelled
from .protocol import Protocol

SCORE_KEYS = {"minimizedAffinity": "empirical_score_kcal_mol", "CNNscore": "cnn_pose_score",
              "CNNaffinity": "cnn_affinity_pK", "CNNvariance": "cnn_variance"}
REQUIRED = ("minimizedAffinity", "CNNscore", "CNNaffinity")


class GninaError(Exception):
    pass


def gnina_version(exe: str, env_extra=None) -> str:
    p = subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=60, env={**os.environ, **(env_extra or {})})
    txt = (p.stdout + p.stderr).strip()
    if p.returncode != 0 or not txt:
        raise GninaError(f"GNINA_NOT_RUNNABLE rc={p.returncode}: {txt[-300:]}")
    return txt.splitlines()[0]


def build_argv(exe, proto: Protocol, receptor, ligand, box, seed, out_sdf, *, device, cpu, gpu_index=0,
               operation="redock", log=None):
    if proto.engine != "gnina":
        raise GninaError("PROTOCOL_ENGINE_MISMATCH")
    if proto.cnn_scoring in ("none", "0", ""):
        raise GninaError("CNN_DISABLED_NOT_ALLOWED: a GNINA run with CNN off is a different protocol; do not label it GNINA-CNN")
    if device not in ("cpu", "gpu"):
        raise GninaError("DEVICE_MUST_BE_EXPLICIT: 'cpu' or 'gpu'")
    a = [exe, "-r", receptor, "-l", ligand, "-o", out_sdf,
         "--center_x", box.center[0], "--center_y", box.center[1], "--center_z", box.center[2],
         "--size_x", box.size[0], "--size_y", box.size[1], "--size_z", box.size[2],
         "--cnn_scoring", proto.cnn_scoring, "--cnn", proto.cnn_model, "--cpu", cpu, "--seed", seed]
    if log:
        a += ["--log", log]
    if operation == "redock":
        a += ["--exhaustiveness", proto.exhaustiveness, "--num_modes", proto.num_modes]
    elif operation == "score_only":
        a += ["--score_only"]
    elif operation == "minimize":
        a += ["--minimize"]
    else:
        raise GninaError(f"UNKNOWN_OPERATION {operation}")
    a += ["--no_gpu"] if device == "cpu" else ["--device", gpu_index]
    return [str(x) for x in a]


def parse_gnina_sdf(path: Path) -> list[dict]:
    from rdkit import Chem
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise GninaError("EMPTY_OUTPUT: engine produced no pose file content")
    try:
        sup = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
        mols = list(sup)
    except Exception as e:   # noqa
        raise GninaError(f"INVALID_SDF: {e}")
    poses = []
    for i, m in enumerate(mols):
        if m is None:
            raise GninaError(f"INVALID_SDF_RECORD {i}")
        rec = {"pose_index": i + 1}
        for k, name in SCORE_KEYS.items():
            if m.HasProp(k):
                v = float(m.GetProp(k))
                if not math.isfinite(v):
                    raise GninaError(f"NONFINITE_{k}")
                rec[name] = v
            elif k in REQUIRED:
                raise GninaError(f"MISSING_{k}")
            else:
                rec[name] = None     # NA, never 0
        poses.append(rec)
    if not poses:
        raise GninaError("EMPTY_OUTPUT")
    return poses


def run_gnina_job(exe, receptor_pdb: Path, ligand_sdf: Path, box, proto: Protocol, seed: int, job_dir: Path, *,
                  device: str, cpu: int, timeout: float, gpu_index=0, operation="redock", env_extra=None,
                  nvidia_smi=None, cancel=None, cancel_file=None, meta=None) -> dict:
    job_dir = Path(job_dir)
    job_dir.mkdir(parents=True, exist_ok=False)           # a job dir is never reused or overwritten
    rec = {"job": job_dir.name, "engine": "gnina", "operation": operation, "device": device, "seed": seed,
           "protocol_hash": proto.hash(), "box_hash": box.hash(), "status": "running", "error": None, "meta": meta or {},
           "started": time.time()}
    try:
        rec["engine_version"] = gnina_version(exe, env_extra)
        rec["engine_sha256"] = sha256_file(exe)
        rec["receptor_sha256"] = sha256_file(receptor_pdb)
        lig = job_dir / "ligand.sdf"; lig.write_bytes(Path(ligand_sdf).read_bytes())
        rec["ligand_sha256"] = sha256_file(lig)
        out = job_dir / "poses.sdf"
        argv = build_argv(exe, proto, str(receptor_pdb), str(lig), box, seed, str(out), device=device, cpu=cpu,
                          gpu_index=gpu_index, operation=operation, log=str(job_dir / "gnina.log"))
        rec["argv"] = argv; atomic_write_json(job_dir / "job_record.json", rec)
        smi = nvidia_smi if device == "gpu" else None
        rec["audit"] = run_argv(argv, job_dir, timeout=timeout, env_extra=env_extra, cancel=cancel, cancel_file=cancel_file, sample_gpu=smi)
        if rec["audit"]["returncode"] != 0:
            raise GninaError(f"ENGINE_EXIT_{rec['audit']['returncode']}: " + (job_dir / "stderr.log").read_text(errors="replace")[-600:])
        rec["poses"] = parse_gnina_sdf(out)
        rec["output_sha256"] = sha256_file(out)
        if device == "gpu":
            rows = [r for s in json.loads((job_dir / "gpu_samples.json").read_text()) if (job_dir / "gpu_samples.json").exists() for r in s.get("our_process_rows", [])] if (job_dir / "gpu_samples.json").exists() else []
            rec["gpu_process_seen_in_nvidia_smi"] = bool(rows)
        rec["status"] = "ok"
    except Cancelled:
        rec["status"] = "cancelled"
    except TimeoutError as e:
        rec["status"] = "timeout"; rec["error"] = str(e)
    except Exception as e:      # noqa: failure is recorded, never swallowed into a fake score
        rec["status"] = "failed"; rec["error"] = f"{type(e).__name__}: {e}"
    rec["finished"] = time.time()
    atomic_write_json(job_dir / "job_record.json", rec)
    return rec
