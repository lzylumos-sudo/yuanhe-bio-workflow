"""Thin CLI over sbvs2.scheduler (runs on the compute host). usage: sched_cli.py <submit|resume|cancel|status> ..."""
import argparse, json, os, sys
from pathlib import Path
_h = Path(__file__).resolve()
for _c in (_h.parents[2] / "src", _h.parents[1]):
    if (_c / "sbvs2").is_dir():
        sys.path.insert(0, str(_c)); break
from sbvs2 import scheduler as S
ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest="cmd", required=True)
s = sp.add_parser("submit"); s.add_argument("--stage", required=True); s.add_argument("--engine", required=True)
s.add_argument("--device", choices=["cpu", "gpu"], required=True); s.add_argument("--seeds", default="11,23,47")
s.add_argument("--operation", default="redock"); s.add_argument("--exhaustiveness", type=int, default=16)
s.add_argument("--cpu", type=int, default=4); s.add_argument("--timeout", type=float, default=1800)
s.add_argument("--limit-confs", type=int, default=0); s.add_argument("--lib-dirs", default=""); s.add_argument("--tag", default="")
s.add_argument("--min-free-gb", type=float, default=5.0)
for n in ("resume", "cancel", "status"):
    p = sp.add_parser(n); p.add_argument("run_id", nargs="?")
ap.add_argument("--runs", default=os.path.expanduser("~/sbvs_platform/runs_sched"))
a = ap.parse_args()
runs = Path(a.runs); runs.mkdir(parents=True, exist_ok=True)
reg = S.Registry(runs / "registry.sqlite3"); os.environ["SBVS_DB"] = str(runs / "registry.sqlite3")
if a.cmd == "submit":
    spec = S.make_spec(a.stage, engine_path=a.engine, device=a.device, seeds=[int(x) for x in a.seeds.split(",")], operation=a.operation,
                       exhaustiveness=a.exhaustiveness, cpu=a.cpu, timeout=a.timeout, limit_confs=a.limit_confs, tag=a.tag,
                       lib_dirs=[x for x in a.lib_dirs.split(os.pathsep) if x], min_free_gb=a.min_free_gb)
    print(json.dumps(S.submit(spec, runs, reg, python=sys.executable)))
elif a.cmd == "resume": print(json.dumps(S.resume(a.run_id, reg, python=sys.executable)))
elif a.cmd == "cancel": print(json.dumps(S.cancel(a.run_id, reg)))
else:
    rows = [r for r in reg.all() if not a.run_id or r["id"] == a.run_id]
    for r in rows:
        rd = Path(r["path"]); jobs = list((rd / "jobs").glob("*.a*")) if rd.exists() else []
        done = sum(1 for j in jobs if (j / "job_record.json").exists() and json.loads((j / "job_record.json").read_text()).get("status") == "ok")
        print(json.dumps({**r, "jobs_dirs": len(jobs), "jobs_ok": done}))
