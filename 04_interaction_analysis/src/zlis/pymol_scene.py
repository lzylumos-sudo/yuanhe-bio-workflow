"""Write the PyMOL deliverables for one analysed ligand and (optionally) build the scene."""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

from rdkit import Chem

from .diagram2d import COLORS
from .structure import Atom

TEMPLATE = Path(__file__).with_name("scene_template.py")


def write_complex_pdb(result: dict, out_dir: Path) -> Path:
    """Real atoms only: prepared protein + ligand (input pose) + input waters/ions. No pseudo-atoms."""
    lig = result["_lig"]
    prot_lines = [l for l in Path(result["_pmol_path"]).read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
    g = result["_group"]
    c, n, ic, rn = g.key
    lines = list(prot_lines)
    lines.append("TER")
    m = lig.mol
    conf = m.GetConformer()
    for i, a in enumerate(m.GetAtoms()):
        p = conf.GetAtomPosition(i)
        info = lig.atoms[i]
        name = info["name"][:4]
        name = f" {name:<3s}" if len(name) < 4 else name
        el = a.GetSymbol()
        lines.append(f"HETATM{90000 + i:5d} {name} {rn:>3s} {c[:1]:1s}{n:4d}{ic[:1] or ' '}   "
                     f"{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}{1.0:6.2f}{0.0:6.2f}          {el:>2s}")
    # connect records so PyMOL draws the correct ligand bonds
    for a in m.GetAtoms():
        nb = [90000 + n.GetIdx() for n in a.GetNeighbors()]
        for k in range(0, len(nb), 4):
            lines.append("CONECT" + f"{90000 + a.GetIdx():5d}" + "".join(f"{x:5d}" for x in nb[k:k + 4]))
    wat = [a for a in result.get("_waters", [])]
    for k, w in enumerate(wat):
        lines.append(f"HETATM{95000 + k:5d}  O   HOH {w.chain[:1]:1s}{w.resseq:4d}{w.icode[:1] or ' '}   "
                     f"{w.x:8.3f}{w.y:8.3f}{w.z:8.3f}{1.0:6.2f}{w.bfactor:6.2f}           O")
    lines.append("END")
    p = out_dir / "complex_analysis.pdb"
    p.write_text("\n".join(lines) + "\n")
    sdf = out_dir / "ligand.sdf"
    mm = Chem.Mol(m)
    mm.SetProp("_Name", result["ligand"]["label"])
    w = Chem.SDWriter(str(sdf))
    w.write(mm)
    w.close()
    return p


def write_scene_files(result: dict, out_dir: Path, hidden_types=("hydrophobic", "vdw_contact"),
                      width: int = 2400, height: int = 1800, background: str = "white") -> dict:
    out_dir = Path(out_dir)
    write_complex_pdb(result, out_dir)
    # interactions.json in the output folder is the single source of truth
    clean = {k: v for k, v in result.items() if not k.startswith("_")}
    (out_dir / "interactions.json").write_text(json.dumps(clean, indent=1, ensure_ascii=False), encoding="utf-8")
    cfg = {"colors": COLORS, "hidden_types": list(hidden_types), "background": background,
           "png_width": width, "png_height": height, "ligand_carbon": "cyan", "show_ligand_h": False,
           "hbond_arrow": "donor -> acceptor (role marker, not electron flow)"}
    (out_dir / "scene_config.json").write_text(json.dumps(cfg, indent=1), encoding="utf-8")
    shutil.copyfile(TEMPLATE, out_dir / "rebuild_scene.py")
    (out_dir / "scene.pml").write_text(
        "# Rebuilds the ZLIS scene. Run from this folder:  pymol scene.pml\n"
        "run rebuild_scene.py\n", encoding="utf-8")
    return {"pdb": str(out_dir / "complex_analysis.pdb"), "sdf": str(out_dir / "ligand.sdf")}


def build_scene(out_dir: Path, no_png: bool = False, timeout: int = 600) -> dict:
    """Run rebuild_scene.py in a fresh headless PyMOL (same interpreter that has pymol installed)."""
    import os
    out_dir = Path(out_dir)
    env = dict(os.environ)
    if no_png:
        env["ZLIS_NO_PNG"] = "1"
    code = (f"import pymol,sys,os; pymol.finish_launching(['pymol','-qc']); os.chdir({str(out_dir)!r}); "
            f"exec(compile(open('rebuild_scene.py').read(),'rebuild_scene.py','exec'))")
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=timeout,
                          env=env, cwd=str(out_dir))
    ok = (out_dir / "scene.pse").exists() and proc.returncode == 0
    return {"ok": ok, "stdout": proc.stdout[-400:], "stderr": proc.stderr[-1200:],
            "pse": str(out_dir / "scene.pse"), "png": str(out_dir / "binding_site_main.png")}
