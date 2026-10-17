"""Optional PLIP cross-check (separate environment, run as a subprocess).

PLIP is NOT part of the core analysis. It is an independent second opinion that
perceives interactions with its own rules (OpenBabel). Differences are reported,
never silently reconciled: ZLIS/ProLIF remains the source of truth.
Configure the interpreter with ZLIS_PLIP_PYTHON or env/plip_python.txt.
"""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

_RUNNER = r'''
import sys, json
from plip.structure.preparation import PDBComplex
pdb, het, chain, num = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
m = PDBComplex(); m.load_pdb(pdb); m.analyze()
out = {"plip_version": None, "hbonds": [], "salt_bridges": [], "hydrophobic": [], "pi_stacking": [], "pi_cation": [], "halogen": [], "metal": [], "water_bridges": []}
try:
    import plip; out["plip_version"] = getattr(plip, "__version__", None)
except Exception: pass
for k, v in m.interaction_sets.items():
    h, c, n = k.split(":")
    if h != het or c != chain or int(n) != num: continue
    for b in v.hbonds_ldon + v.hbonds_pdon:
        out["hbonds"].append({"res": f"{b.reschain}:{b.restype}{b.resnr}", "distance_ad": round(b.distance_ad, 2), "angle": round(b.angle, 1), "protein_is_donor": bool(b.protisdon)})
    for s in v.saltbridge_lneg + v.saltbridge_pneg:
        out["salt_bridges"].append({"res": f"{s.reschain}:{s.restype}{s.resnr}", "distance": round(s.distance, 2)})
    for x in v.hydrophobic_contacts:
        out["hydrophobic"].append({"res": f"{x.reschain}:{x.restype}{x.resnr}", "distance": round(x.distance, 2)})
    for x in v.pistacking:
        out["pi_stacking"].append({"res": f"{x.reschain}:{x.restype}{x.resnr}", "distance": round(x.distance, 2)})
    for x in v.pication_laro + v.pication_paro:
        out["pi_cation"].append({"res": f"{x.reschain}:{x.restype}{x.resnr}", "distance": round(x.distance, 2)})
    for x in v.halogen_bonds:
        out["halogen"].append({"res": f"{x.reschain}:{x.restype}{x.resnr}", "distance": round(x.distance, 2)})
    for x in v.metal_complexes:
        out["metal"].append({"res": f"{x.reschain}:{x.restype}{x.resnr}", "distance": round(x.distance, 2)})
    for x in v.water_bridges:
        out["water_bridges"].append({"res": f"{x.reschain}:{x.restype}{x.resnr}", "distance": round(x.distance_aw, 2)})
print(json.dumps(out), flush=True)
'''


def find_plip_python(project_root: Path) -> str | None:
    cand = os.environ.get("ZLIS_PLIP_PYTHON")
    f = project_root / "env" / "plip_python.txt"
    if not cand and f.exists():
        cand = f.read_text(encoding="utf-8").strip()
    return cand if cand and Path(cand).exists() else None


def run_plip(result: dict, out_dir: Path, project_root: Path, timeout: int = 300) -> dict:
    py = find_plip_python(project_root)
    if not py:
        return {"status": "skipped", "reason": "PLIP interpreter not configured (set ZLIS_PLIP_PYTHON or env/plip_python.txt)"}
    lig = result["ligand"]
    try:
        p = subprocess.run([py, "-c", _RUNNER, str(Path(out_dir) / "complex_analysis.pdb"), lig["resname"], lig["chain"], str(lig["resseq"])],
                           capture_output=True, text=True, timeout=timeout)
        line = [l for l in p.stdout.splitlines() if l.startswith("{")]
        if p.returncode != 0 or not line:
            return {"status": "failed", "reason": (p.stderr or p.stdout)[-400:]}
        plip = json.loads(line[-1])
    except Exception as ex:  # timeout, missing binary, ...
        return {"status": "failed", "reason": str(ex)}
    return {"status": "ok", "plip": plip, "comparison": compare(result, plip)}


def compare(result: dict, plip: dict) -> dict:
    zl = [e for e in result["events"] if e["type"] == "hydrogen_bond"]
    zl_res = sorted(f'{e["residue"]["label"].replace(":", ":", 1)}' for e in zl)
    pl_res = sorted(h["res"].replace(":", ":") for h in plip["hbonds"])
    norm = lambda s: s  # ZLIS label A:ASP189 and PLIP A:ASP189 share the same format
    z_set, p_set = {}, {}
    for e in zl:
        z_set.setdefault(e["residue"]["label"], []).append(e["distance"])
    for h in plip["hbonds"]:
        p_set.setdefault(h["res"], []).append(h["distance_ad"])
    only_z = sorted(set(z_set) - set(p_set))
    only_p = sorted(set(p_set) - set(z_set))
    both = sorted(set(z_set) & set(p_set))
    z_ion = sorted({e["residue"]["label"] for e in result["events"] if e["type"] == "ionic"})
    p_ion = sorted({s["res"] for s in plip["salt_bridges"]})
    return {
        "hbond": {"zlis_count": len(zl), "plip_count": len(plip["hbonds"]), "residues_in_both": both,
                  "only_zlis": only_z, "only_plip": only_p,
                  "note": "ZLIS (ProLIF) counts every donor-acceptor pair; PLIP keeps its own set. Differences reflect different hydrogen/angle rules, not an error flag by themselves."},
        "salt_bridge": {"zlis_residues": z_ion, "plip_residues": p_ion, "only_zlis": sorted(set(z_ion) - set(p_ion)), "only_plip": sorted(set(p_ion) - set(z_ion))},
        "hydrophobic_residues_plip": sorted({x["res"] for x in plip["hydrophobic"]}),
    }
