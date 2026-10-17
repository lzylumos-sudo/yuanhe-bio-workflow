"""ZLIS command line:  python -m zlis <command> ..."""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")


def run_all(pdb: str, out_root: str, ligand=None, template=None, charge=None, model=None,
            offline=False, no_png=False, layout=None, png_size=(2400, 1800), plip=False) -> list[dict]:
    from .diagram2d import write_diagrams
    from .pipeline import analyze
    from .pymol_scene import build_scene, write_scene_files
    from .tables import write_tables

    summaries = []
    for res in analyze(pdb, out_root, ligand=ligand, template=template, net_charge=charge,
                       model=model, allow_network=not offline):
        out = Path(res["out_dir"])
        files = {}
        files.update(write_tables(res, out))
        files.update(write_diagrams(res, out, layout_path=layout))
        files.update(write_scene_files(res, out, width=png_size[0], height=png_size[1]))
        if plip:
            from .pipeline import PROJECT_ROOT
            from .plip_check import run_plip
            pc = run_plip(res, out, PROJECT_ROOT)
            (out / "plip_crosscheck.json").write_text(json.dumps(pc, indent=1, ensure_ascii=False), encoding="utf-8")
            files["plip"] = pc["status"]
        scene = build_scene(out, no_png=no_png)
        files["pse_ok"] = scene["ok"]
        if not scene["ok"]:
            files["pymol_error"] = scene["stderr"]
        # canonical interactions.json (without private objects) was written by write_scene_files
        summaries.append({"ligand": res["ligand"]["label"], "out_dir": str(out),
                          "n_events": len(res["events"]), "chemistry": res["chemistry"]["status"],
                          "files": files})
    return summaries


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="zlis", description="Zhiyu Ligand Interaction Studio")
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("ligands", help="list ligand candidates in a PDB")
    s.add_argument("pdb")

    a = sub.add_parser("analyze", help="run tables + 2D + PyMOL package")
    a.add_argument("pdb")
    a.add_argument("--out", default=None, help="results folder (default: <project>/results)")
    a.add_argument("--ligand", help="e.g. A:BEN1 or BEN")
    a.add_argument("--template", help="SDF/MOL/SMILES file giving bond orders and protonation")
    a.add_argument("--charge", type=int, help="net formal charge if no template is known")
    a.add_argument("--model", type=int)
    a.add_argument("--offline", action="store_true", help="never contact RCSB for ligand templates")
    a.add_argument("--no-png", action="store_true", help="skip the ray-traced PyMOL PNG")
    a.add_argument("--plip", action="store_true", help="optional PLIP cross-check (needs a separate PLIP interpreter)")
    a.add_argument("--layout", help="layout json saved from the 2D editor")

    w = sub.add_parser("web", help="start the local web UI (127.0.0.1 only)")
    w.add_argument("--port", type=int, default=8765)
    w.add_argument("--no-browser", action="store_true")

    args = ap.parse_args(argv)
    from .pipeline import PROJECT_ROOT, list_candidates

    if args.cmd == "ligands":
        print(json.dumps(list_candidates(args.pdb), indent=1, ensure_ascii=False))
        return 0
    if args.cmd == "analyze":
        out = args.out or str(PROJECT_ROOT / "results")
        res = run_all(args.pdb, out, args.ligand, args.template, args.charge, args.model,
                      args.offline, args.no_png, args.layout, plip=args.plip)
        for r in res:
            print(f'{r["ligand"]}: {r["n_events"]} events, chemistry={r["chemistry"]} -> {r["out_dir"]}')
        return 0
    if args.cmd == "web":
        from .web import serve
        serve(args.port, open_browser=not args.no_browser)
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
