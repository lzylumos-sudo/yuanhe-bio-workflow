"""Project / run records: run_manifest, environment capture, state log (plan §25-26)."""
from __future__ import annotations
import json, platform, sys, time, uuid
from pathlib import Path

from . import SCHEMA_VERSION, __version__
from .hashing import atomic_write_json, sha256_file, write_checksums
from .paths import RuntimePaths, init_project, init_run
from .states import check_transition


def new_project(paths: RuntimePaths, name: str) -> Path:
    pid = f"{time.strftime('%Y%m%d')}_{uuid.uuid4().hex[:8]}"
    p = init_project(paths, pid)
    atomic_write_json(p / "project.json", {"schema_version": SCHEMA_VERSION, "project_id": pid, "name": name,
                                           "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "sbvs2_version": __version__})
    return p


def capture_environment() -> dict:
    return {"python": sys.version.split()[0], "platform": platform.platform(), "machine": platform.machine(),
            "sbvs2": __version__, "schema_version": SCHEMA_VERSION}


def new_run(paths: RuntimePaths, project_id: str, protocol_hash: str, inputs: dict, seed_plan: list) -> Path:
    rid = f"{time.strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:10]}"
    r = init_run(paths, project_id, rid)
    atomic_write_json(r / "00_run_record" / "run_manifest.json", {
        "schema_version": SCHEMA_VERSION, "run_id": rid, "protocol_hash": protocol_hash,
        "input_hashes": inputs, "seeds": seed_plan, "environment": capture_environment(), "state": "IMPORTED"})
    log_state(r, None, "IMPORTED")
    return r


def log_state(run: Path, cur: str | None, new: str, note: str = "") -> None:
    if cur is not None:
        check_transition(cur, new)
    with open(run / "00_run_record" / "state_log.jsonl", "a", encoding="utf-8") as f:
        f.write(json.dumps({"t": time.time(), "from": cur, "to": new, "note": note}, ensure_ascii=False) + "\n")


def seal(run: Path) -> int:
    """Final checksum manifest for the run directory."""
    return write_checksums(run)
