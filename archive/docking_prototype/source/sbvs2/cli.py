"""sbvs2 command line. Every command prints plain JSON so results are scriptable and auditable."""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path

from .hashing import atomic_write_json, sha256_file, verify_checksums
from .paths import RuntimePaths
from .project import new_project


def _p(o): print(json.dumps(o, ensure_ascii=False, indent=1, default=str))


def cmd_doctor(a):
    from .docking.engines import engine_status, gpu_status
    import platform
    rp = RuntimePaths.resolve(a.home)
    try:
        import rdkit; rd = rdkit.__version__
    except Exception:
        rd = None
    _p({"paths": {k: str(v) for k, v in rp.__dict__.items()}, "arch": platform.machine(), "rdkit": rd,
        "engines": {n: engine_status(n) for n in ("gnina", "vina")}, "gpu": gpu_status()})


def cmd_new(a):
    _p({"project": str(new_project(RuntimePaths.resolve(a.home), a.name))})


def cmd_receptor(a):
    from .receptor.ingest import ingest_pdb
    from dataclasses import asdict
    r = ingest_pdb(a.pdb)
    d = asdict(r)
    if a.out:
        atomic_write_json(a.out, d)
    _p({"sha256": r.sha256, "chains": {c: {k: v for k, v in x.items() if k != "sequence_from_atoms"} for c, x in r.chains.items()},
        "issues": r.issues, "plan": r.plan, "hetero_counts": {k: len(v) for k, v in r.hetero.items()}})


def cmd_ligands(a):
    from .ligands.csv_io import ingest_ligands
    res = ingest_ligands(Path(a.csv).read_bytes(), a.id_col, a.smiles_col)
    out = {"parents": len(res.parents), "source_rows": len(res.source_map), "failures": [f.__dict__ for f in res.failures], "flags": res.flags}
    if a.out:
        atomic_write_json(a.out, {"parents": res.parents, "source_map": res.source_map, **out})
    _p(out)


def cmd_estimate(a):
    from .ligands.prepare import estimate_jobs
    _p(estimate_jobs(a.states, a.conformers, a.receptors, a.seeds))


def cmd_verify(a):
    probs = verify_checksums(a.dir)
    _p({"dir": a.dir, "ok": not probs, "problems": probs}); sys.exit(1 if probs else 0)


def main(argv=None):
    ap = argparse.ArgumentParser("sbvs2"); ap.add_argument("--home")
    sp = ap.add_subparsers(dest="cmd", required=True)
    sp.add_parser("doctor").set_defaults(f=cmd_doctor)
    s = sp.add_parser("new-project"); s.add_argument("name"); s.set_defaults(f=cmd_new)
    s = sp.add_parser("ingest-receptor"); s.add_argument("pdb"); s.add_argument("--out"); s.set_defaults(f=cmd_receptor)
    s = sp.add_parser("ingest-ligands"); s.add_argument("csv"); s.add_argument("--id-col"); s.add_argument("--smiles-col"); s.add_argument("--out"); s.set_defaults(f=cmd_ligands)
    s = sp.add_parser("estimate"); [s.add_argument(f"--{k}", type=int, required=True) for k in ("states", "conformers", "receptors", "seeds")]; s.set_defaults(f=cmd_estimate)
    s = sp.add_parser("verify"); s.add_argument("dir"); s.set_defaults(f=cmd_verify)
    a = ap.parse_args(argv); a.f(a)


if __name__ == "__main__":
    main()
