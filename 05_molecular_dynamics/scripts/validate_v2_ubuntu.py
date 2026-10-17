#!/usr/bin/env python3
"""Run a bounded native Ubuntu/WSL engineering test, including real continuation.

Run AFTER ``bash setup.sh``:
    python3 scripts/validate_v2_ubuntu.py --max-seconds 1200

This deliberately confirms only the supplied 1fkj engineering fixture's construct.
It never approves a user's structure, starts a long production job, or changes a
saved run's configuration. Mac execution records BLOCKED without launching MD.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
import uuid

ROOT = Path(__file__).resolve().parents[1]
TOTAL_STEPS = 5000  # 10 ps / 2 fs; guarded against the generated effective config.


def stamp():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    checksum = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def read_json(path, default=None):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        return default


def snapshot(job, directories):
    records = {}
    for directory in directories:
        for path in sorted((job / directory).rglob("*")):
            if path.is_file():
                stat = path.stat()
                records[path.relative_to(job).as_posix()] = {
                    "sha256": digest(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    return records


def validate_committed(job, progress):
    replica = job / "md/replica_01"
    if not progress or not progress.get("segments"):
        raise ValueError("No committed production segment exists.")
    if sum(segment["steps"] for segment in progress["segments"]) != progress["production_steps"]:
        raise ValueError("Committed segment step sum differs from progress.")
    for segment in progress["segments"]:
        for kind in ("dcd", "csv"):
            filename = segment[kind]
            if Path(filename).name != filename or digest(replica / filename) != segment[kind + "_sha256"]:
                raise ValueError("Invalid committed segment path/checksum: " + filename)
    checkpoint = progress["checkpoint"]
    if Path(checkpoint).name != checkpoint or digest(replica / checkpoint) != progress["checkpoint_sha256"]:
        raise ValueError("Invalid committed binary checkpoint path/checksum.")


class BudgetExceeded(RuntimeError):
    pass


class Driver:
    def __init__(self, directory, max_seconds):
        self.directory = directory
        self.deadline = time.monotonic() + max_seconds
        self.summary = {
            "schema_version": 2, "started": stamp(), "status": "RUNNING",
            "scope": "Native engineering validation only; no scientific sampling/convergence claim.",
            "host": {"system": platform.system(), "release": platform.release(),
                     "machine": platform.machine(), "python": sys.version,
                     "driver_python": sys.executable},
            "limits": {"max_seconds": max_seconds, "production_ps": 10, "nvt_ps": 1,
                       "npt_ps": 1, "replicates": 1, "platform": "CPU", "mmgbsa_enabled": True},
            "commands": [], "checks": {name: {"status": "NOT_RUN"} for name in (
                "native_environment", "fixture_unchanged", "workflow_complete",
                "interrupted_checkpoint_reload", "interrupted_resume", "completed_resume_no_recompute")},
        }
        self.save()

    def save(self):
        write_json(self.directory / "summary.json", self.summary)

    def check(self, name, status, **details):
        self.summary["checks"][name] = {"status": status, **details}
        self.save()

    def remaining(self):
        remaining = self.deadline - time.monotonic()
        if remaining <= 0:
            raise BudgetExceeded("Native validation exceeded the explicit wall-clock budget.")
        return remaining

    @staticmethod
    def stop(process, first_signal=signal.SIGTERM, grace=10):
        """Signal only this driver's isolated child process group, never other jobs."""
        if process.poll() is not None:
            return {"already_exited": True, "exit_code": process.returncode}
        escalation = []
        for sig, allowance in ((first_signal, grace), (signal.SIGTERM, 5), (signal.SIGKILL, 5)):
            if process.poll() is not None:
                break
            try:
                os.killpg(process.pid, sig)
                escalation.append(signal.Signals(sig).name)
            except ProcessLookupError:
                break
            try:
                process.wait(timeout=max(0.1, allowance))
                break
            except subprocess.TimeoutExpired:
                pass
        return {"signals": escalation, "exit_code": process.poll()}

    def command(self, args, label, interrupt_job=None):
        log = self.directory / (label + ".log")
        record = {"argv": [str(arg) for arg in args], "cwd": str(ROOT), "started": stamp(), "log": str(log)}
        self.summary["commands"].append(record)
        self.save()
        interrupted = None
        with log.open("w", encoding="utf-8") as handle:
            process = subprocess.Popen(record["argv"], cwd=ROOT, stdin=subprocess.DEVNULL,
                                       stdout=handle, stderr=subprocess.STDOUT, start_new_session=True)
            record["pid"] = process.pid
            try:
                while process.poll() is None:
                    remaining = self.remaining()
                    if interrupt_job is not None:
                        progress = read_json(interrupt_job / "md/replica_01/progress.json")
                        if progress and progress.get("phase") == "production" and 0 < progress.get("production_steps", 0) < TOTAL_STEPS:
                            validate_committed(interrupt_job, progress)
                            write_json(self.directory / "observed_before_interrupt.json", progress)
                            interrupted = self.stop(process, signal.SIGINT, grace=min(30, remaining))
                            record["intentional_interrupt"] = interrupted
                            break
                    time.sleep(min(0.1, remaining))
            except BaseException:
                record["cleanup"] = self.stop(process, signal.SIGTERM, grace=5)
                record["exit_code"] = process.poll()
                record["ended"] = stamp()
                self.save()
                raise
        record.update(exit_code=process.poll(), ended=stamp())
        self.save()
        return record, interrupted


def probe_checkpoint(job, output):
    """Executed only by the restored Linux runtime. Loads state but takes zero steps."""
    if platform.system() != "Linux" or platform.machine().lower() not in ("x86_64", "amd64"):
        raise ValueError("The native checkpoint probe requires x86-64 Linux.")
    job = Path(job).resolve()
    if not job.is_relative_to((ROOT / "runs").resolve()):
        raise ValueError("Checkpoint probe job must belong to this project's runs directory.")
    sys.path.insert(0, str(ROOT))
    from modules.automation.backend import _reload_final_checkpoint
    config = read_json(job / "effective_config.json")
    progress = read_json(job / "md/replica_01/progress.json")
    validate_committed(job, progress)
    if progress["phase"] != "production" or not 0 < progress["production_steps"] < TOTAL_STEPS:
        raise ValueError("Checkpoint probe requires a genuinely interrupted, incomplete production run.")
    checkpoint_before = digest(job / "md/replica_01" / progress["checkpoint"])
    ps = _reload_final_checkpoint(config, job / "md/system", job / "md/replica_01", progress, progress["production_steps"])
    if checkpoint_before != digest(job / "md/replica_01" / progress["checkpoint"]):
        raise ValueError("Zero-step probe changed the checkpoint.")
    result = {"status": "PASS", "native_binary_checkpoint_reloaded": True, "integration_steps_taken": 0,
              "production_steps": progress["production_steps"], "checkpoint_time_ps": ps,
              "checkpoint_sha256": checkpoint_before, "runtime_python": sys.executable, "host": platform.platform()}
    write_json(output, result)
    print(json.dumps(result, indent=2))


def validate_acceptance(job):
    state = read_json(job / "state.json", {})
    acceptance = read_json(job / "acceptance.json", {})
    config = read_json(job / "effective_config.json", {})
    md = config.get("MD", {})
    if state.get("status") != "COMPLETE" or acceptance.get("passed") is not True:
        raise ValueError("Workflow did not complete its native audit.")
    if (md.get("production"), md.get("nvt"), md.get("npt"), md.get("timestep"), config.get("replicates")) != ("10ps", "1ps", "1ps", "2fs", 1):
        raise ValueError("Effective configuration differs from the bounded 10 ps engineering test.")
    if not config.get("mmgbsa", {}).get("enabled"):
        raise ValueError("MM/GBSA was disabled; the full requested chain was not tested.")
    if acceptance.get("verification_level") != "NATIVE_CHECKPOINT_AND_TRAJECTORY":
        raise ValueError("Final acceptance is not a real native checkpoint-and-trajectory audit.")
    replicas = acceptance.get("replicas", [])
    if len(replicas) != 1 or replicas[0].get("binary_checkpoint_reloaded") is not True or replicas[0].get("steps") != TOTAL_STEPS:
        raise ValueError("Final binary checkpoint reload or expected step count is missing.")
    if replicas[0].get("mmgbsa") != "COMPUTED":
        raise ValueError("Actual MM/GBSA computation was not verified.")
    artifacts = acceptance.get("artifact_sha256", {})
    if not artifacts:
        raise ValueError("Native audit did not record artifact checksums.")
    for relative, expected in artifacts.items():
        path = (job / relative).resolve()
        if not path.is_relative_to(job.resolve()) or digest(path) != expected:
            raise ValueError("Accepted artifact missing or changed: " + relative)
    return acceptance


def run_validation(driver, job_id):
    if platform.system() != "Linux" or platform.machine().lower() not in ("x86_64", "amd64"):
        driver.check("native_environment", "BLOCKED", reason="Requires x86-64 Ubuntu/WSL Linux; Mac is not a native MD validation host.")
        driver.summary["status"] = "BLOCKED"
        return 4
    driver.check("native_environment", "PASS", note="Host architecture only; actual dependency/forcefield checks remain in the workflow.")
    fixture = ROOT / "examples/1fkj/complex.pdb"
    template = ROOT / "examples/1fkj/ligand.sdf"
    fixture_hashes = {str(path.relative_to(ROOT)): digest(path) for path in (fixture, template)}
    job = ROOT / "runs" / job_id
    driver.summary.update(job=str(job), engineering_fixture={"input": str(fixture), "hashes": fixture_hashes,
        "construct_confirmation_scope": "Bundled 1fkj engineering fixture only, not a scientific approval of user structures."})
    configuration = driver.directory / "engineering_config.yaml"
    # JSON is a YAML subset and avoids requiring PyYAML in the driver interpreter.
    write_json(configuration, {"MD": {"platform": "CPU", "cpu_threads": 2, "timestep": "2fs",
        "production": "10ps", "nvt": "1ps", "npt": "1ps", "equilibration": "2ps",
        "report_interval": "1ps", "checkpoint_interval": "2ps"}, "replicates": 1,
        "analysis": {"stride": 1, "discard_ns": 0, "sasa_stride": 1},
        "mmgbsa": {"enabled": True, "max_frames": 4, "start_ns": 0},
        "decisions": {"construct_confirmed": True, "water_policy": "replace_bulk"}})
    driver.save()
    first, signal_record = driver.command(["bash", ROOT / "run_md.sh", "--input", fixture,
        "--config", configuration, "--run-id", job_id, "--smoke-test", "--non-interactive"],
        "01_start_and_interrupt", interrupt_job=job)
    stopped = read_json(job / "md/replica_01/progress.json", {})
    genuine_interruption = (signal_record is not None and signal_record.get("signals") == ["SIGINT"]
        and first["exit_code"] in (130, -signal.SIGINT) and stopped.get("phase") == "production"
        and 0 < stopped.get("production_steps", 0) < TOTAL_STEPS)
    if genuine_interruption:
        validate_committed(job, stopped)
        write_json(driver.directory / "interrupted_progress.json", stopped)
        probe_output = driver.directory / "interrupted_checkpoint_probe.json"
        # Static shell program; every project/job path travels as a quoted positional argument.
        shell = 'source "$1/runtime.sh"; exec "$NEXUS_ENV_PREFIX/bin/python" "$2" --probe-checkpoint "$3" "$4"'
        probe, _ = driver.command(["bash", "-c", shell, "nexus-checkpoint-probe", ROOT,
            Path(__file__).resolve(), job, probe_output], "02_zero_step_checkpoint_probe")
        evidence = read_json(probe_output, {})
        if probe["exit_code"] != 0 or evidence.get("native_binary_checkpoint_reloaded") is not True:
            driver.check("interrupted_checkpoint_reload", "FAIL", result=evidence)
            raise ValueError("Interrupted binary checkpoint could not be loaded by the native runtime.")
        driver.check("interrupted_checkpoint_reload", "PASS", **{k: v for k, v in evidence.items() if k != "status"})
        resumed, _ = driver.command(["bash", ROOT / "run_md.sh", "--resume", job, "--non-interactive"], "03_resume")
        if resumed["exit_code"] != 0:
            raise ValueError("Actual continuation failed; inspect 03_resume.log and the job state.")
        finished = read_json(job / "md/replica_01/progress.json", {})
        if finished.get("segments", [])[:len(stopped["segments"])] != stopped["segments"]:
            raise ValueError("Continuation changed its already committed segment history.")
        validate_committed(job, stopped)
        if finished.get("production_steps", 0) <= stopped["production_steps"]:
            raise ValueError("Continuation did not advance beyond the interrupted checkpoint.")
        driver.check("interrupted_resume", "PASS", resumed_from_steps=stopped["production_steps"],
                     final_steps=finished["production_steps"], committed_prefix_unchanged=True,
                     configuration_changed=False, checkpoint_reload_probe=str(probe_output))
    elif first["exit_code"] == 0:
        driver.check("interrupted_resume", "NOT_RUN", reason="Workflow completed before a genuine incomplete production checkpoint could be interrupted.")
        driver.check("interrupted_checkpoint_reload", "NOT_RUN", reason="No incomplete interrupted checkpoint was available.")
    else:
        state = read_json(job / "state.json", {})
        blocked = state.get("status") in ("DEPENDENCY_BLOCKED", "WAITING_FOR_DECISION") or first["exit_code"] == 2 and not state
        driver.check("workflow_complete", "BLOCKED" if blocked else "FAIL", state=state,
                     reason="Initial workflow did not complete or produce a verified graceful interruption.")
        driver.summary["status"] = "BLOCKED" if blocked else "FAIL"
        return 4 if blocked else 1
    acceptance = validate_acceptance(job)
    driver.check("workflow_complete", "PASS", acceptance=str(job / "acceptance.json"),
                 artifacts_verified=len(acceptance["artifact_sha256"]), actual_mmgbsa="COMPUTED",
                 final_binary_checkpoint_reloaded=True)
    before = snapshot(job, ["md", "analysis"])
    prior_state = read_json(job / "state.json")
    repeat, _ = driver.command(["bash", ROOT / "run_md.sh", "--resume", job, "--non-interactive"], "04_completed_resume")
    after = snapshot(job, ["md", "analysis"])
    after_state = read_json(job / "state.json")
    if repeat["exit_code"] != 0 or before != after or prior_state.get("stages") != after_state.get("stages"):
        raise ValueError("Completed-run resume recomputed/changed MD or analysis artifacts, or failed verification.")
    validate_acceptance(job)
    driver.check("completed_resume_no_recompute", "PASS", files_unchanged=len(before),
                 compared=["sha256", "size", "mtime_ns", "stage_records"])
    for relative, expected in fixture_hashes.items():
        if digest(ROOT / relative) != expected:
            raise ValueError("Validation changed its original bundled input fixture.")
    driver.check("fixture_unchanged", "PASS", sha256=fixture_hashes)
    statuses = [item["status"] for item in driver.summary["checks"].values()]
    driver.summary["status"] = "PASS" if all(value == "PASS" for value in statuses) else "NOT_RUN"
    return 0 if driver.summary["status"] == "PASS" else 3


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--max-seconds", type=float, default=1200, help="Whole validation wall-clock budget; default 1200 seconds.")
    parser.add_argument("--probe-checkpoint", nargs=2, metavar=("JOB", "OUTPUT"), help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.probe_checkpoint:
        probe_checkpoint(*args.probe_checkpoint)
        return 0
    if args.max_seconds <= 0 or args.max_seconds > 7200:
        parser.error("--max-seconds must be > 0 and <= 7200 (two hours); this never increases simulation length.")
    tag = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ_") + uuid.uuid4().hex[:8]
    directory = ROOT / "runs" / ("validation_" + tag)
    directory.mkdir(parents=True, exist_ok=False)
    driver = Driver(directory, args.max_seconds)
    try:
        code = run_validation(driver, "validation_md_" + tag)
    except BudgetExceeded as exc:
        driver.summary.update(status="BLOCKED", reason=str(exc))
        code = 4
    except KeyboardInterrupt:
        driver.summary.update(status="BLOCKED", reason="Validation interrupted by the user; no passing result is claimed.")
        code = 130
    except Exception as exc:
        driver.summary.update(status="FAIL", reason=f"{type(exc).__name__}: {exc}")
        code = 1
    driver.summary.update(ended=stamp(), summary_file=str(directory / "summary.json"))
    driver.save()
    print(json.dumps(driver.summary, ensure_ascii=False, indent=2))
    return code


if __name__ == "__main__":
    raise SystemExit(main())
