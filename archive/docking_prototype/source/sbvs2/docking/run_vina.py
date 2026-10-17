"""Real AutoDock Vina execution through the python bindings, with per-job isolation (plan §17, §25).

Every job: own work dir, raw output file, stdout/err log, exit status. Scores come only from the engine output.
"""
from __future__ import annotations
import json, re, time, traceback
from pathlib import Path

from ..hashing import atomic_write_json, sha256_file
from .box import Box
from .protocol import Protocol


def parse_vina_pdbqt(text: str) -> list[dict]:
    poses = []
    for blk in text.split("ENDMDL"):
        m = re.search(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", blk)
        mo = re.search(r"MODEL\s+(\d+)", blk)
        if m and mo:
            poses.append({"pose_id": int(mo.group(1)), "vina_affinity": float(m.group(1)),
                          "rmsd_lb": float(m.group(2)), "rmsd_ub": float(m.group(3))})
    return poses


def run_vina_job(receptor_pdbqt: Path, ligand_pdbqt_text: str, box: Box, proto: Protocol, seed: int,
                 job_dir: Path, cpu: int = 0) -> dict:
    from vina import Vina
    job_dir.mkdir(parents=True, exist_ok=False)          # a job dir is never reused/overwritten
    lig = job_dir / "ligand.pdbqt"; lig.write_text(ligand_pdbqt_text)
    out = job_dir / "poses.pdbqt"
    rec = {"job": job_dir.name, "seed": seed, "receptor_sha256": sha256_file(receptor_pdbqt),
           "ligand_sha256": sha256_file(lig), "protocol_hash": proto.hash(), "box_hash": box.hash(),
           "status": "running", "engine": "vina-python-1.2.7", "started": time.time()}
    try:
        v = Vina(sf_name="vina", seed=seed, cpu=cpu, verbosity=0)
        v.set_receptor(str(receptor_pdbqt))
        v.set_ligand_from_file(str(lig))
        v.compute_vina_maps(center=list(box.center), box_size=list(box.size))
        v.dock(exhaustiveness=proto.exhaustiveness, n_poses=proto.num_modes)
        v.write_poses(str(out), n_poses=proto.num_modes, overwrite=True)
        poses = parse_vina_pdbqt(out.read_text())
        if not poses:
            raise RuntimeError("EMPTY_OUTPUT: engine produced no poses")       # never counted as success
        rec.update(status="ok", poses=poses, output_sha256=sha256_file(out))
    except Exception as e:                                                       # isolate failure per job
        rec.update(status="failed", error=f"{type(e).__name__}: {e}", trace=traceback.format_exc()[-800:])
    rec["finished"] = time.time()
    atomic_write_json(job_dir / "job_record.json", rec)
    return rec
