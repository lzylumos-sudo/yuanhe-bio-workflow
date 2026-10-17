"""Subprocess runner shared by all engines (plan §17/§21): argv arrays only, own process group,
timeout, cancel (flag file or callable), stdout/stderr to files, optional GPU sampling. Never invents results."""
from __future__ import annotations
import os, signal, subprocess, time
from pathlib import Path

from ..hashing import atomic_write_json


class Cancelled(Exception):
    pass


def gpu_sample(nvidia_smi: str | None, pid: int) -> dict:
    if not nvidia_smi:
        return {"t": time.time(), "error": "no nvidia-smi"}
    out = {"t": time.time()}
    try:
        q = subprocess.run([nvidia_smi, "--query-gpu=utilization.gpu,memory.used", "--format=csv,noheader,nounits"],
                           capture_output=True, text=True, timeout=5)
        out["gpu"] = q.stdout.strip()
        a = subprocess.run([nvidia_smi, "--query-compute-apps=pid,used_gpu_memory", "--format=csv,noheader,nounits"],
                           capture_output=True, text=True, timeout=5)
        out["our_process_rows"] = [l.strip() for l in a.stdout.splitlines() if l.split(",")[0].strip() == str(pid)]
    except Exception as e:   # noqa
        out["error"] = str(e)
    return out


def run_argv(argv, folder, *, timeout: float, env_extra: dict | None = None, cancel=None, cancel_file=None,
             sample_gpu: str | None = None, sample_every: float = 1.0) -> dict:
    """Run argv in `folder`. Returns audit dict (also written to execution.json). Raises Cancelled / TimeoutError
    only AFTER the audit is written, so an interrupted job still leaves evidence."""
    folder = Path(folder)
    env = {**os.environ, **(env_extra or {})}
    start, reason, samples, last = time.time(), None, [], 0.0
    argv = [str(a) for a in argv]
    with open(folder / "stdout.log", "wb") as so, open(folder / "stderr.log", "wb") as se:
        p = subprocess.Popen(argv, stdout=so, stderr=se, env=env, cwd=str(folder), start_new_session=True, shell=False)
        audit = {"argv": argv, "pid": p.pid, "started": start, "timeout_s": timeout}
        atomic_write_json(folder / "execution.json", audit)
        while p.poll() is None:
            if sample_gpu and time.time() - last >= sample_every:
                last = time.time(); samples.append(gpu_sample(sample_gpu, p.pid))
            if (cancel and cancel()) or (cancel_file and Path(cancel_file).exists()):
                reason = "cancelled"
            elif time.time() - start > timeout:
                reason = "timeout"
            if reason:
                os.killpg(p.pid, signal.SIGTERM)
                try:
                    p.wait(10)
                except subprocess.TimeoutExpired:
                    os.killpg(p.pid, signal.SIGKILL); p.wait()
                break
            time.sleep(0.1)
    audit.update(returncode=p.returncode, finished=time.time(), elapsed_s=time.time() - start, termination_reason=reason)
    atomic_write_json(folder / "execution.json", audit)
    if samples:
        atomic_write_json(folder / "gpu_samples.json", samples)
    if reason == "cancelled":
        raise Cancelled("cancelled")
    if reason == "timeout":
        raise TimeoutError(f"timeout after {timeout}s")
    return audit
