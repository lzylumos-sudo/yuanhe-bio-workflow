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
    raise ValueError('Legacy PDBQT-only execution disabled: use vina_official.run_vina_job with pre-docking evidence')
