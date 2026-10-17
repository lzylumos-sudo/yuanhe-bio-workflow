#!/usr/bin/env python3
"""Wait for this module's completed run, then independently audit it. No MD steps."""
from __future__ import annotations
import argparse
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

PROJECT_ROOT = Path(__file__).resolve().parents[1]

def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def read(path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}

def write(path, record):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(record, indent=2), encoding="utf-8")
    temporary.replace(path)

def file_hash(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()

def isolate_current_candidate(folder, run_id, evidence):
    archived = []
    for index, candidate in enumerate([folder / "MD_TEST_PASS.md", PROJECT_ROOT / "MD_TEST_PASS.md"]):
        if candidate.is_file() and f"Run ID: {run_id}" in candidate.read_text(encoding="utf-8").splitlines():
            destination = evidence / f"candidate_PASS_before_failed_audit_{index}.md"
            if destination.exists():
                destination = evidence / f"candidate_PASS_before_failed_audit_{index}_{time.time_ns()}.md"
            digest = file_hash(candidate)
            shutil.copyfile(candidate, destination)
            if file_hash(destination) != digest:
                raise RuntimeError("Candidate evidence preservation failed; original left intact")
            candidate.unlink()
            archived.append({"file": destination.name, "sha256": digest})
    return archived

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--poll-seconds", type=float, default=60)
    parser.add_argument("--timeout-hours", type=float, default=72)
    parser.add_argument("--delivery-dir", type=Path)
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.run_id):
        parser.error("run-id must be a plain safe directory name")
    if args.poll_seconds < 5 or args.timeout_hours <= 0:
        parser.error("poll must be >=5 seconds and timeout positive")
    folder = PROJECT_ROOT / "results" / args.run_id
    if not (folder / "manifest.json").is_file():
        parser.error("The requested run does not exist")
    evidence = PROJECT_ROOT / "installation_evidence" / ("finalizer_" + args.run_id)
    evidence.mkdir(parents=True, exist_ok=True)
    lock_path = evidence / "observer.lock"
    # This observer's lock is separate from the simulation's lock.
    lock = lock_path.open("a")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        print("An observer for this run is already active", file=sys.stderr)
        return 2
    started = time.monotonic()
    record = {"run_id": args.run_id, "pid": os.getpid(), "started": now(),
              "status": "WAITING", "observer_sha256": file_hash(Path(__file__)),
              "scope": "Reads only the requested run; never advances, stops or restarts a simulation."}
    status_path = evidence / "status.json"
    def publish():
        record["updated"] = now()
        write(status_path, record)
        if args.delivery_dir:
            args.delivery_dir.mkdir(parents=True, exist_ok=True)
            write(args.delivery_dir / "10ns_Background acceptance status.json", record)
    publish()
    validator = PROJECT_ROOT / "scripts/validate_10ns.py"
    while True:
        manifest = read(folder / "manifest.json")
        record["manifest_status"] = manifest.get("status", "UNREADABLE")
        progress = read(folder / "md/replica_01/progress.json")
        record["checkpointed_production_steps"] = progress.get("production_steps", 0)
        record["phase"] = progress.get("phase")
        if manifest.get("status") in {"failed", "error"}:
            record.update(status="SIMULATION_FAILED", reason="The scientific workflow reported failure; no automatic retry.")
            break
        if manifest.get("status") == "complete" and not (folder / ".running").exists():
            record["status"] = "INDEPENDENT_AUDIT_RUNNING"
            publish()
            with (evidence / "validator.log").open("w") as output:
                result = subprocess.run([sys.executable, str(validator), "--run-id", args.run_id],
                                        cwd=PROJECT_ROOT, stdout=output, stderr=subprocess.STDOUT)
            record["validator_exit_code"] = result.returncode
            final = read(PROJECT_ROOT / "installation_evidence/10ns_final_audit.json")
            if result.returncode == 0 and final.get("passed") is True and final.get("run_id") == args.run_id:
                record.update(status="PASS", completed=now(), final_audit_sha256=file_hash(PROJECT_ROOT / "installation_evidence/10ns_final_audit.json"))
                if args.delivery_dir:
                    shutil.copyfile(PROJECT_ROOT / "MD_TEST_PASS.md", args.delivery_dir / "MD_TEST_PASS.md")
                    shutil.copyfile(PROJECT_ROOT / "installation_evidence/10ns_final_audit.json", args.delivery_dir / "10ns_final_audit.json")
                publish()
                print(json.dumps(record), flush=True)
                return 0
            record.update(status="INDEPENDENT_AUDIT_FAILED", reason="See this observer's validator.log; no PASS is certified.")
            record["preserved_candidate_reports"] = isolate_current_candidate(folder, args.run_id, evidence)
            break
        if time.monotonic() - started >= args.timeout_hours * 3600:
            record.update(status="TIMED_OUT", reason="Observer timeout; simulation was not interrupted and is not accepted.")
            break
        publish()
        time.sleep(args.poll_seconds)
    record["completed"] = now()
    if "preserved_candidate_reports" not in record:
        record["preserved_candidate_reports"] = isolate_current_candidate(folder, args.run_id, evidence)
    publish()
    previous = read(PROJECT_ROOT / "validation_status.json")
    if previous.get("status") != "PASS" or previous.get("run_id") == args.run_id:
        write(PROJECT_ROOT / "validation_status.json", {
            "status": record["status"], "run_id": args.run_id, "checked_at": now(),
            "scope": "10ns independent acceptance not passed; no scientific completion claim.",
            "evidence": str(status_path.relative_to(PROJECT_ROOT)),
        })
    print(json.dumps(record), flush=True)
    return 2

if __name__ == "__main__":
    raise SystemExit(main())
