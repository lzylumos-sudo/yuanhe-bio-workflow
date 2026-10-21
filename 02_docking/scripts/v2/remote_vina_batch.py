"""Run a staged, approved input set through REAL official AutoDock Vina CPU on the compute host. Resume-safe and idempotent.

usage: remote_vina_batch.py --stage DIR --out RUN_DIR --vina EXE --device cpu --seeds 11,23,47
         [--operation redock] [--exhaustiveness N] [--timeout S] [--cpu N] [--limit-confs K]
A job dir that already holds a verifiable job_record (status ok + output hash match) is reused; anything else is
renamed to <job>.quarantine_<ts> and re-run. Existing run dirs are never wiped; --out may be reused only for resume.
"""
import argparse, json, os, sys, time
from pathlib import Path
_here = Path(__file__).resolve()
for _c in (_here.parents[2] / "src", _here.parents[1]):      # repo layout, or deployed layout code/{scripts,sbvs2}
    if (_c / "sbvs2").is_dir():
        sys.path.insert(0, str(_c)); break
from sbvs2.hashing import sha256_file, atomic_write_json, write_checksums
from sbvs2.docking.box import Box
from sbvs2.docking.protocol import Protocol
from sbvs2.docking.vina_official import run_vina_job

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True); ap.add_argument("--out", required=True); ap.add_argument("--vina", required=True)
ap.add_argument("--device", choices=["cpu"], default="cpu"); ap.add_argument("--seeds", default="11,23,47")
ap.add_argument("--operation", choices=["redock"], default="redock"); ap.add_argument("--exhaustiveness", type=int, default=32)
ap.add_argument("--num-modes", type=int, default=9); ap.add_argument("--timeout", type=float, default=1800)
ap.add_argument("--cpu", type=int, default=4); ap.add_argument("--limit-confs", type=int, default=0)
ap.add_argument("--lib-dirs", default="")
ap.add_argument("--cancel-file", default="")
a = ap.parse_args()
stage, out = Path(a.stage), Path(a.out)
meta = json.loads((stage / "stage_meta.json").read_text())
man = json.loads((stage / "STAGE_MANIFEST.json").read_text())
for rel, h in man.items():                     # inputs must still be what was staged
    if rel != "STAGE_MANIFEST.json" and sha256_file(stage / rel) != h:
        sys.exit(f"STAGE_HASH_MISMATCH {rel}")
b = meta["box"]; box = Box(tuple(b["center"]), tuple(b["size"]), b["source"], b["revision"], b["approved_reference_pose"])
assert box.hash() == meta["box"]["hash"], "BOX_HASH_MISMATCH"
seeds = tuple(int(x) for x in a.seeds.split(","))
proto = Protocol("vina_" + a.operation, "vina", a.exhaustiveness, a.num_modes, seeds, "none")
env = {"LD_LIBRARY_PATH": os.pathsep.join([x for x in a.lib_dirs.split(os.pathsep) if x] + [os.environ.get("LD_LIBRARY_PATH", "")])} if a.lib_dirs else {}
(out / "jobs").mkdir(parents=True, exist_ok=True)
smi = "nvidia-smi"
rm = out / "run_manifest.json"
if not rm.exists():
    atomic_write_json(rm, {"created": time.time(), "stage_manifest_sha256": sha256_file(stage / "STAGE_MANIFEST.json"),
        "protocol_hash": proto.hash(), "device": a.device, "operation": a.operation, "seeds": seeds, "vina": a.vina,
        "exhaustiveness": a.exhaustiveness, "num_modes": a.num_modes, "cpu": a.cpu, "timeout": a.timeout, "host": os.uname().nodename})
else:
    old = json.loads(rm.read_text())
    assert old["protocol_hash"] == proto.hash() and old["device"] == a.device, "RESUME_PROTOCOL_MISMATCH"
ligs = sorted((stage / "ligands").glob("*.sdf"))
if a.limit_confs:
    ligs = ligs[:a.limit_confs]
results, n_reuse = [], 0
for lig in ligs:
    for seed in seeds:
        name = f"{lig.stem}_seed{seed}"; jd = out / "jobs" / name
        if (out / "CANCEL").exists() or (a.cancel_file and Path(a.cancel_file).exists()):
            print("CANCEL flag seen; stopping before", name); sys.exit(3)
        if jd.exists():
            try:
                r = json.loads((jd / "job_record.json").read_text())
                if r["status"] == "ok" and sha256_file(jd / "poses.sdf") == r["output_sha256"] and r["protocol_hash"] == proto.hash():
                    results.append(r); n_reuse += 1; continue
            except Exception:
                pass
            jd.rename(jd.with_name(name + f".quarantine_{int(time.time())}"))
        r = run_vina_job(a.vina, stage / "receptor_prepared.pdb", lig, box, proto, seed, jd, device=a.device, cpu=a.cpu,
                          timeout=a.timeout, operation=a.operation, env_extra=env, nvidia_smi=smi,
                          cancel_file=a.cancel_file or out / "CANCEL", meta={"ligand": lig.name})
        print(name, r["status"], r.get("error") or "", flush=True); results.append(r)
summ = [{"job": r.get("job", ""), "status": r["status"], "device": "cpu", "seed": r["seed"], "elapsed_s": (r.get("audit") or {}).get("elapsed_s"),
         "n_poses": len(r.get("poses", [])), "top": (r.get("poses") or [None])[0], "error": r.get("error")} for r in results]
atomic_write_json(out / "summary.json", {"jobs": summ, "reused": n_reuse, "ok": sum(r["status"] == "ok" for r in results), "total": len(results)})
write_checksums(out, "RUN_CHECKSUMS.sha256")
print(f"done total={len(results)} ok={sum(r['status']=='ok' for r in results)} reused={n_reuse}")
