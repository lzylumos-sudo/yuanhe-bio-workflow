"""Run scheduler (plan §21, §25). Standard library only.

* submit() is idempotent: the run id is the hash of (stage manifest, protocol, device, seeds, operation).
  Re-submitting the same spec (e.g. after a browser reload or an SSH reconnect) returns the existing run, never a second one.
* The worker is a detached process (own session). Closing SSH / the browser does not signal it.
* Compute slots: an exclusive flock per device slot (GPU: 1 job at a time by default) so runs cannot oversubscribe.
* Each job gets its own directory per attempt (job.aN). Failed/corrupt attempts are kept (renamed *.quarantine_*), never deleted.
* Before every job: free-space check; if low the run goes WAITING_STORAGE (nothing is deleted, nothing is marked complete).
* A run is COMPLETED only when every planned job has a hash-verified status=ok record. Otherwise it is PARTIAL / FAILED.
"""
from __future__ import annotations
import fcntl, json, os, shutil, signal, sqlite3, subprocess, sys, time
from pathlib import Path

from .hashing import hash_obj, sha256_file, atomic_write_json, write_checksums
from .docking.box import Box
from .docking.protocol import Protocol
from .docking.vina_official import run_vina_job, validate_resources, VERSION
from .docking.execute import Cancelled

TERMINAL = {"COMPLETED", "FAILED", "CANCELLED", "PARTIAL"}


# ---------------------------------------------------------------- registry
class Registry:
    def __init__(self, db: Path):
        db = Path(db); db.parent.mkdir(parents=True, exist_ok=True)
        # SQLite's busy timeout does not reliably serialize simultaneous changes
        # of journal mode on a newly created database. Serialize initialization
        # across processes, and never reset journal mode for an existing WAL DB.
        with db.with_name(db.name + ".init.lock").open("a") as init_lock:
            fcntl.flock(init_lock, fcntl.LOCK_EX)
            self.db = sqlite3.connect(str(db), timeout=30, isolation_level=None)
            try:
                mode = self.db.execute("PRAGMA journal_mode").fetchone()[0]
                if mode.lower() != "wal":
                    self.db.execute("PRAGMA journal_mode=WAL")
                self.db.execute("CREATE TABLE IF NOT EXISTS runs(id TEXT PRIMARY KEY, path TEXT, state TEXT, pid INTEGER, created REAL, updated REAL)")
            except Exception:
                self.db.close()
                raise

    def get(self, rid):
        r = self.db.execute("SELECT id,path,state,pid,created,updated FROM runs WHERE id=?", (rid,)).fetchone()
        return dict(zip(("id", "path", "state", "pid", "created", "updated"), r)) if r else None

    def put(self, rid, path, state, pid=None):
        now = time.time()
        self.db.execute("INSERT INTO runs VALUES(?,?,?,?,?,?) ON CONFLICT(id) DO UPDATE SET state=excluded.state,pid=excluded.pid,updated=excluded.updated",
                        (rid, str(path), state, pid, now, now))

    def all(self):
        return [dict(zip(("id", "path", "state", "pid", "created", "updated"), r)) for r in
                self.db.execute("SELECT id,path,state,pid,created,updated FROM runs ORDER BY created")]


def _alive(pid, started_token_file: Path | None = None) -> bool:
    if not pid:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


# ---------------------------------------------------------------- spec / planning
def plan_jobs(stage: Path, seeds, limit_confs=0):
    ligs = sorted((Path(stage) / "ligands").glob("*.sdf"))
    if limit_confs:
        ligs = ligs[:limit_confs]
    return [{"name": f"{l.stem}_seed{s}", "ligand": l.name, "seed": int(s)} for l in ligs for s in seeds]


def make_spec(stage, *, engine_path, device, seeds, operation="redock", exhaustiveness=32, num_modes=9, cpu=4,
              timeout=1800, cnn_model="crossdock_default2018_1_3", lib_dirs=(), limit_confs=0, max_attempts=2,
              min_free_gb=5.0, slot="cpu", tag="", energy_range=3.0):
    if device != "cpu": raise ValueError("Official Vina requires device=cpu")
    validate_resources(cpu)
    if operation!="redock":raise ValueError("Only Vina redock is supported for new tasks")
    stage = Path(stage)
    from yuanhe_preflight.gate import require_ligand
    for ligand in (stage/'ligands').glob('*.sdf'): require_ligand(ligand)
    man = stage / "STAGE_MANIFEST.json"
    spec = {"engine":"vina", "engine_version":VERSION, "energy_range":energy_range, "scoring":"vina", "stage_dir": str(stage), "stage_manifest_sha256": sha256_file(man), "engine_path": str(engine_path),
            "device": device, "seeds": [int(s) for s in seeds], "operation": operation, "exhaustiveness": exhaustiveness,
            "num_modes": num_modes, "cpu": cpu, "timeout": timeout, "cnn_status":"not_applicable", "lib_dirs": [],
            "jobs": plan_jobs(stage, seeds, limit_confs), "max_attempts": max_attempts, "min_free_gb": min_free_gb,
            "slot": slot if device == "gpu" else "cpu", "tag": tag}
    spec["run_id"] = "r_" + hash_obj({k: v for k, v in spec.items() if k not in ("min_free_gb", "timeout", "max_attempts")})[:16]
    return spec


def submit(spec: dict, runs_root: Path, registry: Registry, *, python=sys.executable, detach=True) -> dict:
    """Idempotent. Returns {'run_id','created':bool,'state'}."""
    rid = spec["run_id"]; runs_root = Path(runs_root)
    rd = runs_root / rid
    ex = registry.get(rid)
    if ex or rd.exists():
        state = (ex or {}).get("state", "UNKNOWN")
        if state not in TERMINAL and _alive((ex or {}).get("pid")):
            return {"run_id": rid, "created": False, "state": state, "note": "already submitted; worker alive"}
        return {"run_id": rid, "created": False, "state": state, "note": "already submitted; use resume() to continue"}
    (rd / "jobs").mkdir(parents=True)
    atomic_write_json(rd / "run_spec.json", spec)
    registry.put(rid, rd, "QUEUED")
    _spawn(rd, python, detach); return {"run_id": rid, "created": True, "state": "QUEUED"}


def resume(rid: str, registry: Registry, *, python=sys.executable) -> dict:
    ex = registry.get(rid)
    if not ex:
        raise KeyError(rid)
    if ex["state"] not in TERMINAL and _alive(ex["pid"]):
        return {"run_id": rid, "note": "worker alive; nothing to do", "state": ex["state"]}
    if (Path(ex["path"]) / "CANCEL").exists():
        (Path(ex["path"]) / "CANCEL").rename(Path(ex["path"]) / f"CANCEL.consumed_{int(time.time())}")
    registry.put(rid, ex["path"], "QUEUED")
    _spawn(Path(ex["path"]), python, True); return {"run_id": rid, "state": "QUEUED", "note": "worker restarted; verified jobs will be reused"}


def cancel(rid: str, registry: Registry) -> dict:
    ex = registry.get(rid)
    if not ex:
        raise KeyError(rid)
    (Path(ex["path"]) / "CANCEL").write_text(str(time.time()))
    return {"run_id": rid, "note": "cancel flag written; the running engine process group is terminated by its worker"}


def _spawn(rd: Path, python, detach):
    """Start the worker in its own session so an SSH/browser disconnect cannot kill it.

    NOTE: `-I` must NOT be used here - it implies `-E`, which drops PYTHONPATH, and a dead worker looks exactly
    like a queued one. The package root is therefore injected on the command line instead of via the environment.
    """
    pkg_root = str(Path(__file__).resolve().parents[1])
    bootstrap = f"import sys; sys.path.insert(0, {pkg_root!r}); from sbvs2 import scheduler; sys.exit(scheduler.worker({str(rd)!r}))"
    cmd = [python, "-B", "-c", bootstrap]
    with open(rd / "worker.log", "ab") as log:
        return subprocess.Popen(cmd, stdout=log, stderr=log, stdin=subprocess.DEVNULL, start_new_session=detach, close_fds=True)


# ---------------------------------------------------------------- worker
def _free_gb(path: Path) -> float:
    return shutil.disk_usage(path).free / 1e9


def verify_job(jd: Path, proto_hash: str) -> bool:
    try:
        r = json.loads((jd / "job_record.json").read_text())
        return r["status"] == "ok" and r["protocol_hash"] == proto_hash and sha256_file(jd / "poses.sdf") == r["output_sha256"] and bool(r["poses"])
    except Exception:
        return False


def _quarantine(jd: Path) -> Path:
    q = jd.with_name(jd.name + f".quarantine_{int(time.time()*1000)}")
    jd.rename(q); return q


class _Slot:
    def __init__(self, path: Path):
        self.f = open(path, "a+")

    def acquire(self, cancelled, poll=0.5):
        while True:
            try:
                fcntl.flock(self.f, fcntl.LOCK_EX | fcntl.LOCK_NB); return True
            except BlockingIOError:
                if cancelled():
                    return False
                time.sleep(poll)

    def release(self):
        fcntl.flock(self.f, fcntl.LOCK_UN)


def worker(rd: Path, registry: Registry | None = None):
    rd = Path(rd); spec = json.loads((rd / "run_spec.json").read_text()); rid = spec["run_id"]
    reg = registry or Registry(Path(os.environ.get("SBVS_DB") or rd.parent / "registry.sqlite3"))
    stage = Path(spec["stage_dir"]); st_manifest = json.loads((stage / "STAGE_MANIFEST.json").read_text())
    log = lambda *a: print(time.strftime("%F %T"), *a, flush=True)
    state = lambda s: (reg.put(rid, rd, s, os.getpid()), atomic_write_json(rd / "state.json", {"state": s, "t": time.time(), "pid": os.getpid()}))
    cancelled = lambda: (rd / "CANCEL").exists()
    state("RUNNING")
    # inputs must still be exactly what was staged / approved
    if sha256_file(stage / "STAGE_MANIFEST.json") != spec["stage_manifest_sha256"]:
        state("FAILED"); log("STAGE_MANIFEST_CHANGED"); return 2
    for rel, h in st_manifest.items():
        if rel != "STAGE_MANIFEST.json" and sha256_file(stage / rel) != h:
            state("FAILED"); log("STAGE_FILE_CHANGED", rel); return 2
    meta = json.loads((stage / "stage_meta.json").read_text()); b = meta["box"]
    box = Box(tuple(b["center"]), tuple(b["size"]), b["source"], b["revision"], b["approved_reference_pose"])
    if box.hash() != b["hash"]:
        state("FAILED"); log("BOX_HASH_MISMATCH"); return 2
    proto = Protocol("vina_" + spec["operation"], "vina", spec["exhaustiveness"], spec["num_modes"], tuple(spec["seeds"]), "none", energy_range=spec["energy_range"])
    env = {"LD_LIBRARY_PATH": os.pathsep.join(list(spec["lib_dirs"]) + [os.environ.get("LD_LIBRARY_PATH", "")])} if spec["lib_dirs"] else {}
    slot = _Slot(rd.parent / f".slot_{spec['slot']}.lock")
    out = {"reused": 0, "ran": 0, "failed": [], "waiting_storage": False}
    for j in spec["jobs"]:
        from yuanhe_preflight.gate import require_ligand
        require_ligand(stage/"ligands"/j["ligand"])
        if cancelled():
            state("CANCELLED"); log("cancelled"); return 3
        done = None
        for a in range(1, spec["max_attempts"] + 1):
            jd = rd / "jobs" / f"{j['name']}.a{a}"
            if jd.exists() and verify_job(jd, proto.hash()):
                done = jd; out["reused"] += 1; break
            if jd.exists():
                log("quarantine unverifiable attempt", jd.name); _quarantine(jd)
        if done:
            continue
        for a in range(1, spec["max_attempts"] + 1):
            jd = rd / "jobs" / f"{j['name']}.a{a}"
            if jd.exists():
                continue
            while _free_gb(rd) < spec["min_free_gb"]:
                state("WAITING_STORAGE"); out["waiting_storage"] = True; log("low disk, waiting"); time.sleep(5)
                if cancelled():
                    state("CANCELLED"); return 3
            state("RUNNING")
            if not slot.acquire(cancelled):
                state("CANCELLED"); return 3
            try:
                r = run_vina_job(spec["engine_path"], stage / "receptor_prepared.pdb", stage / "ligands" / j["ligand"], box, proto, j["seed"], jd,
                                  device=spec["device"], cpu=spec["cpu"], timeout=spec["timeout"], operation=spec["operation"], env_extra=env,
                                  nvidia_smi=shutil.which("nvidia-smi") if spec["device"] == "gpu" else None, cancel_file=rd / "CANCEL",
                                  meta={"ligand": j["ligand"], "attempt": a})
            finally:
                slot.release()
            log(j["name"], f"a{a}", r["status"], r.get("error") or "")
            if r["status"] == "cancelled":
                state("CANCELLED"); return 3
            if r["status"] == "ok" and verify_job(jd, proto.hash()):
                out["ran"] += 1; done = jd; break
            if r["status"] == "ok":                  # engine said ok but artefact does not verify
                _quarantine(jd)
        if not done:
            out["failed"].append(j["name"])
    ok_jobs = sum(any(verify_job(rd / "jobs" / f"{j['name']}.a{a}", proto.hash()) for a in range(1, spec["max_attempts"] + 1)) for j in spec["jobs"])
    complete = ok_jobs == len(spec["jobs"])
    summary = {"run_id": rid, "planned": len(spec["jobs"]), "verified_ok": ok_jobs, "failed_jobs": out["failed"], "reused": out["reused"], "ran_now": out["ran"], "complete": complete}
    atomic_write_json(rd / "summary.json", summary)
    log("finished", summary)                          # all logging BEFORE the seal so the checksums stay valid
    if complete:                                      # seal only a truly complete run, and seal LAST
        state("COMPLETED"); write_checksums(rd, "RUN_CHECKSUMS.sha256")
    else:
        state("PARTIAL")
    return 0 if complete else 1


if __name__ == "__main__":
    sys.exit(worker(Path(sys.argv[1])))
