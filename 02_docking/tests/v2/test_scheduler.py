"""Scheduler behaviour tests. The engine here is a STUB script (tests only): it proves queue/recovery logic, NOT docking science."""
import json, os, stat, sys, time, pathlib, textwrap
import pytest
from sbvs2.hashing import sha256_file, atomic_write_json
from sbvs2 import scheduler as S

STUB = textwrap.dedent('''\
    #!{py}
    import sys, os, time, pathlib
    a = sys.argv[1:]
    if "--version" in a: print("stub-gnina 0.0 (TEST ONLY)"); sys.exit(0)
    g = lambda k: a[a.index(k) + 1]
    mode = os.environ.get("STUB_MODE", "ok") or "ok"
    marker = pathlib.Path(os.environ.get("STUB_DIR", ".")) / "mode"
    if marker.exists(): mode = marker.read_text().strip()
    lig = pathlib.Path(g("-l")).read_text()
    out = pathlib.Path(g("-o"))
    if mode == "sleep": time.sleep(60)
    if mode == "fail": sys.stderr.write("boom\\n"); sys.exit(7)
    if mode == "empty": out.write_text(""); sys.exit(0)
    body = lig.replace("$$$$", "") .rstrip("\\n")
    out.write_text(body + "\\n> <minimizedAffinity>\\n-7.5\\n\\n> <CNNscore>\\n0.5\\n\\n> <CNNaffinity>\\n5.5\\n\\n$$$$\\n")
''')


@pytest.fixture
def env(tmp_path):
    stage = tmp_path / "stage"; (stage / "ligands").mkdir(parents=True)
    from rdkit import Chem
    from rdkit.Chem import AllChem
    for i, smi in enumerate(["CCO", "CCN"]):
        m = Chem.AddHs(Chem.MolFromSmiles(smi)); AllChem.EmbedMolecule(m, randomSeed=1)
        w = Chem.SDWriter(str(stage / "ligands" / f"L{i}_c0.sdf")); w.write(m); w.close()
    (stage / "receptor_prepared.pdb").write_text("END\n")
    from sbvs2.docking.box import Box
    b = Box((0.0, 0.0, 0.0), (10.0, 10.0, 10.0), "manual")
    atomic_write_json(stage / "stage_meta.json", {"box": {"center": list(b.center), "size": list(b.size), "source": b.source,
                      "revision": b.revision, "approved_reference_pose": b.approved_reference_pose, "hash": b.hash()}})
    files = {str(p.relative_to(stage)): sha256_file(p) for p in stage.rglob("*") if p.is_file()}
    atomic_write_json(stage / "STAGE_MANIFEST.json", files)
    stub = tmp_path / "stub_engine"; stub.write_text(STUB.format(py=sys.executable)); stub.chmod(stub.stat().st_mode | stat.S_IEXEC)
    os.environ["STUB_DIR"] = str(tmp_path)
    reg = S.Registry(tmp_path / "reg.sqlite3")
    yield dict(stage=stage, stub=stub, tmp=tmp_path, reg=reg, runs=tmp_path / "runs")
    os.environ.pop("STUB_DIR", None)
    reg.db.close()


def spec_for(e, **kw):
    kw.setdefault("seeds", [1, 2]); kw.setdefault("cpu", 1); kw.setdefault("timeout", 20)
    return S.make_spec(e["stage"], engine_path=e["stub"], device="cpu", **kw)


def run_inline(e, spec):
    S.submit(spec, e["runs"], e["reg"], detach=False) if False else None


def start(e, spec):
    rd = e["runs"] / spec["run_id"]
    if not rd.exists():
        (rd / "jobs").mkdir(parents=True); atomic_write_json(rd / "run_spec.json", spec)
    return S.worker(rd, e["reg"]), rd


def test_idempotent_submit_returns_same_run(env, monkeypatch):
    monkeypatch.setattr(S, "_spawn", lambda *a, **k: None)
    sp = spec_for(env)
    r1 = S.submit(sp, env["runs"], env["reg"]); r2 = S.submit(sp, env["runs"], env["reg"])
    assert r1["created"] and not r2["created"] and r1["run_id"] == r2["run_id"]
    assert len(env["reg"].all()) == 1
    sp2 = spec_for(env, seeds=[1, 2, 3])
    assert sp2["run_id"] != sp["run_id"]                  # different seeds => different run


def test_complete_run_is_sealed_and_verifiable(env):
    sp = spec_for(env); rc, rd = start(env, sp)
    assert rc == 0 and env["reg"].get(sp["run_id"])["state"] == "COMPLETED"
    summ = json.loads((rd / "summary.json").read_text())
    assert summ["planned"] == 4 == summ["verified_ok"] and summ["complete"]
    from sbvs2.hashing import verify_checksums
    assert verify_checksums(rd, "RUN_CHECKSUMS.sha256") == []


def test_corrupt_output_is_quarantined_and_rerun(env):
    sp = spec_for(env); rc, rd = start(env, sp); assert rc == 0
    victim = next((rd / "jobs").glob("L0_c0_seed1.a1")); (victim / "poses.sdf").write_text("corrupted")
    (rd / "RUN_CHECKSUMS.sha256").unlink()
    rc, _ = start(env, sp)
    assert rc == 0
    assert list((rd / "jobs").glob("L0_c0_seed1.a1.quarantine_*")), "corrupt attempt must be kept, not deleted"
    assert json.loads((rd / "summary.json").read_text())["reused"] == 3     # the other three verified jobs reused untouched


def test_failed_engine_is_never_complete(env):
    (env["tmp"] / "mode").write_text("fail")
    sp = spec_for(env); rc, rd = start(env, sp)
    assert rc == 1 and env["reg"].get(sp["run_id"])["state"] == "PARTIAL"
    assert not (rd / "RUN_CHECKSUMS.sha256").exists()                       # no "all complete" seal
    assert len(json.loads((rd / "summary.json").read_text())["failed_jobs"]) == 4
    assert len(list((rd / "jobs").glob("*.a2"))) == 4                       # independent retry directory


def test_empty_output_rejected(env):
    (env["tmp"] / "mode").write_text("empty")
    rc, rd = start(env, spec_for(env)); assert rc == 1
    rec = json.loads(next((rd / "jobs").glob("*.a1")).joinpath("job_record.json").read_text())
    assert rec["status"] == "failed" and any(k in rec["error"] for k in ("EMPTY_OUTPUT", "MISSING_", "INVALID_SDF"))
    assert "poses" not in rec                                                # no score was invented


def test_timeout_recorded(env):
    (env["tmp"] / "mode").write_text("sleep")
    rc, rd = start(env, spec_for(env, timeout=1, seeds=[1])); assert rc == 1
    rec = json.loads(next((rd / "jobs").glob("*.a1")).joinpath("job_record.json").read_text())
    assert rec["status"] == "timeout"


def test_recovery_after_failure_reuses_and_does_not_overwrite(env):
    (env["tmp"] / "mode").write_text("fail")
    sp = spec_for(env, seeds=[1]); rc, rd = start(env, sp); assert rc == 1
    (env["tmp"] / "mode").write_text("ok")
    rc2, _ = start(env, sp)
    assert rc2 == 0
    # failed attempts a1/a2 stay on disk as evidence; the success is in a new attempt dir or quarantined predecessors
    names = sorted(p.name for p in (rd / "jobs").glob("L0_c0_seed1.*"))
    q = [n for n in names if ".quarantine_" in n]
    assert len(q) == 2 and "L0_c0_seed1.a1" in names, names                   # both failed attempts preserved as evidence
    for n in q:                                                               # and their failure records are intact
        assert json.loads((rd / "jobs" / n / "job_record.json").read_text())["status"] == "failed"


def test_stage_tamper_blocks_run(env):
    sp = spec_for(env)
    (env["stage"] / "receptor_prepared.pdb").write_text("TAMPERED\n")
    rc, rd = start(env, sp)
    assert rc == 2 and env["reg"].get(sp["run_id"])["state"] == "FAILED"
    assert not list((rd / "jobs").iterdir())                                # nothing ran


def test_cancel_flag_stops_before_next_job(env):
    sp = spec_for(env); rd = env["runs"] / sp["run_id"]; (rd / "jobs").mkdir(parents=True)
    atomic_write_json(rd / "run_spec.json", sp); (rd / "CANCEL").write_text("x")
    assert S.worker(rd, env["reg"]) == 3 and env["reg"].get(sp["run_id"])["state"] == "CANCELLED"


def test_low_disk_waits_and_never_completes(env, monkeypatch):
    sp = spec_for(env, min_free_gb=1e9); rd = env["runs"] / sp["run_id"]; (rd / "jobs").mkdir(parents=True)
    atomic_write_json(rd / "run_spec.json", sp)
    seen = []
    monkeypatch.setattr(S.time, "sleep", lambda s: (seen.append(env["reg"].get(sp["run_id"])["state"]), (rd / "CANCEL").write_text("x")))
    assert S.worker(rd, env["reg"]) == 3
    assert "WAITING_STORAGE" in seen
    assert not (rd / "RUN_CHECKSUMS.sha256").exists()


def test_spawned_worker_can_import_the_package(env, monkeypatch):
    """Exercise the actual spawn once, and reap that child instead of duplicating it."""
    import subprocess
    sp = spec_for(env); rd = env["runs"] / sp["run_id"]
    (rd / "jobs").mkdir(parents=True)
    atomic_write_json(rd / "run_spec.json", sp)
    (rd / "CANCEL").write_text("x")
    captured = {}
    real_popen = subprocess.Popen
    def capture(cmd, **kwargs):
        captured["cmd"] = cmd
        return real_popen(cmd, **kwargs)
    monkeypatch.setattr(S.subprocess, "Popen", capture)
    proc = S._spawn(rd, sys.executable, detach=False)
    try:
        code = proc.wait(timeout=60)
    finally:
        if proc.poll() is None:
            proc.kill(); proc.wait()
    output = (rd / "worker.log").read_text()
    assert "-I" not in captured["cmd"]
    assert "ModuleNotFoundError" not in output
    assert code == 3, (code, output)


def test_concurrent_registry_initialization(tmp_path):
    """Independent workers can initialize one new database without losing records."""
    import subprocess
    db = tmp_path / "shared.sqlite3"
    gate = tmp_path / "start"
    source = str(pathlib.Path(S.__file__).resolve().parents[1])
    script = """
import sys, time
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from sbvs2.scheduler import Registry
deadline = time.monotonic() + 20
while not Path(sys.argv[3]).exists():
    if time.monotonic() > deadline: raise RuntimeError('Start barrier timed out')
    time.sleep(.01)
r = Registry(Path(sys.argv[2]))
r.put(sys.argv[4], sys.argv[4], 'QUEUED')
r.db.close()
"""
    children = [subprocess.Popen([sys.executable, '-c', script, source, str(db), str(gate), str(i)],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                for i in range(6)]
    try:
        gate.touch()
        for proc in children:
            out, err = proc.communicate(timeout=30)
            assert proc.returncode == 0, (out, err)
        registry = S.Registry(db)
        try:
            assert {row['id'] for row in registry.all()} == {str(i) for i in range(6)}
        finally:
            registry.db.close()
    finally:
        for proc in children:
            if proc.poll() is None:
                proc.kill(); proc.wait()
