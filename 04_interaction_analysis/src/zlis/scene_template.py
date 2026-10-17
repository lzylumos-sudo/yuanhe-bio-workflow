"""Standalone PyMOL scene rebuilder (copied next to each result as rebuild_scene.py).

Reads complex_analysis.pdb, interactions.json and scene_config.json from its own
folder (relative paths only, so the folder can be moved or copied between
machines) and rebuilds the scene. PyMOL never re-decides what is an
interaction: every dash, arrow and label comes from interactions.json.

Run:   pymol -cq rebuild_scene.py          (writes scene.pse and PNGs)
   or: pymol rebuild_scene.py              (interactive)
H-bond arrows point donor -> acceptor. They show donor/acceptor roles,
not electron flow or force.
"""
import json
import math
import os
import sys

from pymol import cmd, cgo

HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
if not os.path.exists(os.path.join(HERE, "scene_config.json")):
    HERE = os.getcwd()
CFG = json.load(open(os.path.join(HERE, "scene_config.json"), encoding="utf-8"))
DATA = json.load(open(os.path.join(HERE, "interactions.json"), encoding="utf-8"))
COL = CFG["colors"]


def rgb(name):
    h = COL[name].lstrip("#")
    return [int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4)]


def sub(a, b): return [a[i] - b[i] for i in range(3)]
def add(a, b): return [a[i] + b[i] for i in range(3)]
def mul(a, s): return [a[i] * s for i in range(3)]
def norm(a): return math.sqrt(sum(x * x for x in a))
def unit(a):
    n = norm(a)
    return [x / n for x in a] if n else [0, 0, 1]


def dashes(p, q, color, radius, dash=0.28, gap=0.18):
    v = sub(q, p)
    L = norm(v)
    u = unit(v)
    out = []
    t = 0.0
    while t < L:
        e = min(t + dash, L)
        a, b = add(p, mul(u, t)), add(p, mul(u, e))
        out += [cgo.CYLINDER, *a, *b, radius, *color, *color]
        t = e + gap
    return out


def dashdot(p, q, color, radius):
    v = sub(q, p)
    L = norm(v)
    u = unit(v)
    out, t = [], 0.0
    while t < L:
        for seg, gap in ((0.30, 0.12), (0.05, 0.12)):
            if t >= L:
                break
            e = min(t + seg, L)
            out += [cgo.CYLINDER, *add(p, mul(u, t)), *add(p, mul(u, e)), radius, *color, *color]
            t = e + gap
    return out


def solid(p, q, color, radius):
    return [cgo.CYLINDER, *p, *q, radius, *color, *color]


def arrow_head(tip, direction, color, length=0.55, radius=0.2):
    u = unit(direction)
    base = sub(tip, mul(u, length))
    return [cgo.CONE, *base, *tip, radius, 0.0, *color, *color, 1.0, 0.0]


def sphere(p, color, r):
    return [cgo.COLOR, *color, cgo.SPHERE, *p, r]


def add_label(name, pos, text, color):
    cmd.pseudoatom(name, pos=pos, label=text)
    cmd.set("label_color", color, name)


def load_scene():
    cmd.reinitialize()
    cmd.load(os.path.join(HERE, "complex_analysis.pdb"), "complex")
    lig = DATA["ligand"]
    sel_lig = f"complex and chain {lig['chain']} and resn {lig['resname']} and resi {lig['resseq']}{lig['icode']}"
    cmd.select("ligand", sel_lig)
    cmd.select("protein_atoms", "complex and polymer.protein")
    cmd.bg_color(CFG["background"])
    cmd.set("antialias", 2)
    cmd.set("ray_opaque_background", 1 if CFG["background"] != "transparent" else 0)
    cmd.hide("everything")
    cmd.show("cartoon", "protein_atoms")
    cmd.color("grey85", "protein_atoms")
    cmd.set("cartoon_transparency", 0.35)
    cmd.show("sticks", "ligand")
    cmd.util.cbaw("ligand")
    cmd.color(CFG["ligand_carbon"], "ligand and elem C")
    cmd.set("stick_radius", 0.17, "ligand")

    # residues involved in any drawn interaction
    shown = set()
    for ev in DATA["events"]:
        if ev["type"] in CFG["hidden_types"]:
            continue
        r = ev["residue"]
        shown.add((r["chain"], r["resseq"], r["icode"], r["resname"]))
    if shown:
        expr = " or ".join(f"(chain {c} and resi {n}{i} and resn {rn})" for c, n, i, rn in shown)
        cmd.select("site_residues", f"complex and ({expr})")
        cmd.show("sticks", "site_residues")
        cmd.set("stick_radius", 0.13, "site_residues")
        cmd.util.cbag("site_residues")
        cmd.hide("everything", "site_residues and hydro and not (neighbor (ligand))")
        cmd.label("site_residues and name CA", '"%s%s" % (resn, resi)')
        cmd.set("label_size", 16)
        cmd.set("label_color", "black", "site_residues")
    cmd.show("sticks", "ligand")
    if CFG.get("show_ligand_h", True) is False:
        cmd.hide("everything", "ligand and hydro")
    cmd.show("nb_spheres", "complex and (resn HOH) within 6 of ligand")
    cmd.set("nb_spheres_size", 0.25)
    cmd.color("red", "complex and resn HOH")


def draw_events():
    objs = {}
    labs = {}

    def bucket(t):
        objs.setdefault(t, [])
        return objs[t]

    for ev in DATA["events"]:
        t = ev["type"]
        col = rgb(t)
        b = bucket(t)
        lig_xyz = [a["xyz"] for a in ev["ligand_atoms"]]
        prot_xyz = [a["xyz"] for a in ev["protein_atoms"]]
        eid = ev["event_id"]
        if t == "hydrogen_bond":
            d, a = ev["donor"]["xyz"], ev["acceptor"]["xyz"]
            h = ev["hydrogen"]["xyz"]
            # dashes along H...A, solid stub D-H is the covalent bond already drawn by sticks
            b += dashes(h, a, col, 0.045)
            tip = sub(a, mul(unit(sub(a, h)), 0.35))   # stop at acceptor surface
            b += arrow_head(tip, sub(a, h), col)
            labs[eid] = (mul(add(d, a), 0.5), f'{ev["distance"]:.2f}', col)
        elif t == "ionic":
            p, q = lig_xyz[0], prot_xyz[0]
            b += dashdot(p, q, col, 0.04)
            b += sphere(p, col, 0.16) + sphere(q, col, 0.16)
            labs[eid] = (mul(add(p, q), 0.5), f'{ev["distance"]:.2f}', col)
        elif t in ("pi_pi", "cation_pi"):
            p, q = ev.get("ligand_anchor", lig_xyz[0]), ev.get("protein_anchor", prot_xyz[0])
            b += dashes(p, q, col, 0.05, 0.35, 0.2) + sphere(p, col, 0.25) + sphere(q, col, 0.25)
            labs[eid] = (mul(add(p, q), 0.5), f'{ev["distance"]:.2f}', col)
        elif t == "halogen_bond":
            p, q = lig_xyz[-1], prot_xyz[0]
            b += dashdot(p, q, col, 0.04)
            labs[eid] = (mul(add(p, q), 0.5), f'{ev["distance"]:.2f}', col)
        elif t == "metal_coordination":
            p, q = lig_xyz[0], prot_xyz[0]
            b += solid(p, q, col, 0.05) + sphere(q, col, 0.55)
            labs[eid] = (mul(add(p, q), 0.5), f'{ev["distance"]:.2f}', col)
        elif t == "water_bridge":
            w, p, q = ev["water"]["xyz"], lig_xyz[0], prot_xyz[0]
            b += dashes(p, w, col, 0.04) + dashes(w, q, col, 0.04) + sphere(w, col, 0.3)
            labs[eid] = (w, ev["water"]["label"].split(":")[-1], col)
        elif t == "steric_clash":
            p, q = lig_xyz[0], prot_xyz[0]
            b += solid(p, q, col, 0.07) + sphere(mul(add(p, q), 0.5), col, 0.22)
            labs[eid] = (mul(add(p, q), 0.5), f'clash {ev["distance"]:.2f}', col)
        else:  # hydrophobic, vdw_contact
            p, q = lig_xyz[0], prot_xyz[0]
            b += dashes(p, q, col, 0.025, 0.12, 0.2)
    for t, obj in objs.items():
        if obj:
            cmd.load_cgo(obj, f"ix_{t}")
    # distance labels (pseudoatoms live in their own objects, never in the exported PDB)
    cmd.set("label_size", 14)
    for eid, (pos, text, col) in labs.items():
        t = next(e["type"] for e in DATA["events"] if e["event_id"] == eid)
        if t in CFG["hidden_types"]:
            continue
        name = f"lbl_{t}_{eid}"
        cmd.pseudoatom(name, pos=[float(x) for x in pos], label=text)
        cname = f"c_{name}"
        cmd.set_color(cname, col)
        cmd.set("label_color", cname, name)
        cmd.set("label_position", (0, 0.4, 2.0), name)
    for t in CFG["hidden_types"]:
        for o in cmd.get_object_list("all"):
            if o == f"ix_{t}" or o.startswith(f"lbl_{t}_"):
                cmd.disable(o)
    cmd.group("interactions", "ix_* lbl_*")


def finish(png=True):
    cmd.orient("ligand or site_residues")
    cmd.zoom("ligand or site_residues", 6)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_shadow", 0)
    cmd.save(os.path.join(HERE, "scene.pse"))
    if png and not CFG.get("skip_png"):
        cmd.png(os.path.join(HERE, "binding_site_main.png"), width=CFG["png_width"], height=CFG["png_height"], dpi=300, ray=1)


if __name__ in ("__main__", "pymol", "builtins") or True:
    load_scene()
    draw_events()
    finish(png=os.environ.get("ZLIS_NO_PNG") != "1")
    print("ZLIS scene rebuilt:", os.path.join(HERE, "scene.pse"))
