"""2D interaction diagrams: RDKit ligand depiction + ZLIS atom-anchored overlay.

Views
  publication : residue-level, hydrophobic / ionic / water merged per residue, vdW hidden
  detailed    : one line per atom-level event (vdW hidden), atom names, distances, angles
  audit       : every event, plus a full table (HTML)

Hidden items are counted in the footer; nothing is removed from the data.
Label positions can be edited (HTML editor) and stored in layout.json.
"""
from __future__ import annotations

import html
import json
import math
import textwrap
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from .tables import TYPE_LABEL

COLORS = {
    "hydrogen_bond": "#1f5fbf", "ionic": "#d55e00", "hydrophobic": "#8a8a8a", "pi_pi": "#7b2d9b",
    "cation_pi": "#1b9e3e", "halogen_bond": "#e69f00", "metal_coordination": "#8b5a2b",
    "water_bridge": "#00a5b5", "steric_clash": "#d62728", "vdw_contact": "#c0c0c0",
}
DASH = {
    "hydrogen_bond": "7,4", "ionic": "9,3,2,3", "hydrophobic": "2,4", "pi_pi": "6,3", "cation_pi": "6,3",
    "halogen_bond": "9,3,2,3", "metal_coordination": "", "water_bridge": "5,3", "steric_clash": "3,3",
    "vdw_contact": "2,3",
}
W, H = 1180, 900
FONT = "Helvetica, Arial, sans-serif"


def ligand_2d(result: dict):
    """Return (mol2d, {rdkit_idx: (x, y)}, drawer_svg_body) for the heavy-atom ligand."""
    lig = result["_lig"]
    m = Chem.RemoveHs(lig.mol)
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(m)
    d = rdMolDraw2D.MolDraw2DSVG(W, H)
    o = d.drawOptions()
    o.padding = 0.30
    o.clearBackground = False
    o.bondLineWidth = 2
    o.baseFontSize = 0.55
    o.addStereoAnnotation = False
    d.DrawMolecule(m)
    d.FinishDrawing()
    svg = d.GetDrawingText()
    body = svg[svg.index("<!-- END OF HEADER -->") + len("<!-- END OF HEADER -->"):svg.rindex("</svg>")]
    pos = {i: (d.GetDrawCoords(i).x, d.GetDrawCoords(i).y) for i in range(m.GetNumAtoms())}
    return m, pos, body


def _view_events(events: list[dict], view: str):
    hidden = defaultdict(int)
    keep = []
    for e in events:
        t = e["type"]
        if t == "vdw_contact" and view != "audit":
            hidden["van der Waals contacts"] += 1
            continue
        keep.append(e)
    if view != "publication":
        return keep, hidden, []
    # publication: merge by (residue, type) for ionic / hydrophobic / water_bridge
    merged, out = defaultdict(list), []
    for e in keep:
        if e["type"] in ("ionic", "hydrophobic", "water_bridge"):
            merged[(e["residue"]["label"], e["type"])].append(e)
        else:
            out.append(e)
    groups = []
    for (res, t), evs in merged.items():
        best = min(evs, key=lambda x: x["distance"])
        m = dict(best)
        seen, atoms = set(), []
        for x in evs:
            for a in x["ligand_atoms"]:
                if a["uid"] not in seen:
                    seen.add(a["uid"])
                    atoms.append(a)
        m["ligand_atoms"] = atoms
        m["protein_atoms"] = sorted({a["uid"]: a for x in evs for a in x["protein_atoms"]}.values(), key=lambda a: a["name"])
        m["merged_n"] = len(evs)
        m["event_ids"] = [x["event_id"] for x in evs]
        out.append(m)
        if len(evs) > 1:
            hidden[f"{TYPE_LABEL.get(t, t)} events merged into residue-level lines"] += len(evs) - 1
    return out, hidden, groups


def _anchor(e: dict, pos: dict) -> tuple[float, float]:
    pts = [pos[a["rdkit_idx"]] for a in e["ligand_atoms"] if a.get("rdkit_idx") in pos]
    if not pts:
        raise ValueError(f"event {e['event_id']}: ligand atom has no 2D position")
    return (sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts))


def _clip(cx, cy, w, h, tx, ty):
    """Point on the rectangle (centre cx,cy size w,h) in the direction of (tx,ty)."""
    dx, dy = tx - cx, ty - cy
    if dx == 0 and dy == 0:
        return cx, cy
    sx = (w / 2) / abs(dx) if dx else float("inf")
    sy = (h / 2) / abs(dy) if dy else float("inf")
    s = min(sx, sy)
    return cx + dx * s, cy + dy * s


def _label_text(res: dict, evs: list[dict], view: str) -> list[str]:
    head = f'{res["resname"]}{res["resseq"]}{res["icode"]} ({res["chain"]})'
    if view == "publication":
        parts = sorted({e.get("protein_part", "") for e in evs} - {""})
        return [head, "/".join(parts)] if parts else [head]
    atoms = sorted({a["name"] for e in evs for a in e["protein_atoms"]})
    parts = sorted({e.get("protein_part", "") for e in evs} - {""})
    return [head, f'{"/".join(parts)}: {",".join(atoms)}'[:34]]


def _layout(groups: dict, pos: dict, saved: dict | None, centre, radius):
    # Include space for RDKit atom labels, not just the atom anchor points.
    left = min(p[0] for p in pos.values()) - 38
    right = max(p[0] for p in pos.values()) + 38
    top = min(p[1] for p in pos.values()) - 38
    bottom = max(p[1] for p in pos.values()) + 38
    boxes = {}
    for res, evs in groups.items():
        lines = evs["lines"]
        w = max(len(t) for t in lines) * 7.4 + 18
        h = 18 * len(lines) + 10
        if saved and res in saved.get("labels", {}):
            cx, cy = saved["labels"][res]
        else:
            pts = [_anchor(e, pos) for e in evs["events"]]
            ax, ay = np.mean([p[0] for p in pts]), np.mean([p[1] for p in pts])
            ang = math.atan2(ay - centre[1], ax - centre[0])
            cx, cy = centre[0] + math.cos(ang) * (radius + 70), centre[1] + math.sin(ang) * (radius + 70)
        boxes[res] = [cx, cy, w, h, saved is not None and res in saved.get("labels", {})]
    keys = list(boxes)
    for _ in range(400):
        moved = False
        for i, a in enumerate(keys):
            for b in keys[i + 1:]:
                A, B = boxes[a], boxes[b]
                ox = (A[2] + B[2]) / 2 + 10 - abs(A[0] - B[0])
                oy = (A[3] + B[3]) / 2 + 10 - abs(A[1] - B[1])
                if ox > 0 and oy > 0:
                    moved = True
                    if ox < oy:
                        s = 1 if A[0] >= B[0] else -1
                        if not A[4]: A[0] += s * ox / 2
                        if not B[4]: B[0] -= s * ox / 2
                    else:
                        s = 1 if A[1] >= B[1] else -1
                        if not A[4]: A[1] += s * oy / 2
                        if not B[4]: B[1] -= s * oy / 2
        for k in keys:
            A = boxes[k]
            if not A[4]:
                dx, dy = A[0] - centre[0], A[1] - centre[1]
                dist = math.hypot(dx, dy) or 1
                if dist < radius + 40:
                    A[0] = centre[0] + dx / dist * (radius + 40)
                    A[1] = centre[1] + dy / dist * (radius + 40)
                if (A[0]+A[2]/2 > left and A[0]-A[2]/2 < right and
                        A[1]+A[3]/2 > top and A[1]-A[3]/2 < bottom):
                    choices = [(left-A[2]/2-10, A[1]), (right+A[2]/2+10, A[1]),
                               (A[0], top-A[3]/2-10), (A[0], bottom+A[3]/2+10)]
                    A[0], A[1] = min(choices, key=lambda p: math.hypot(p[0]-A[0], p[1]-A[1]))
                    moved = True
                A[0] = min(max(A[0], A[2] / 2 + 8), W - A[2] / 2 - 8)
                A[1] = min(max(A[1], A[3] / 2 + 40), H - A[3] / 2 - 60)
        if not moved:
            break
    return boxes


def render_svg(result: dict, view: str = "detailed", layout: dict | None = None, interactive: bool = False):
    mol2d, pos, body = ligand_2d(result)
    events, hidden, _ = _view_events(result["events"], view)
    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    centre = ((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2)
    radius = max(max(xs) - min(xs), max(ys) - min(ys)) / 2 + 20

    groups: dict[str, dict] = {}
    for e in events:
        r = e["residue"]
        g = groups.setdefault(r["label"], {"res": r, "events": []})
        g["events"].append(e)
    for g in groups.values():
        g["lines"] = _label_text(g["res"], g["events"], view)
    boxes = _layout(groups, pos, layout, centre, radius)

    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="{FONT}">',
           '<defs>']
    for t, c in COLORS.items():
        out.append(f'<marker id="arr_{t}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="9" markerHeight="9" orient="auto">'
                   f'<path d="M0,0 L10,5 L0,10 z" fill="{c}"/></marker>')
        out.append(f'<marker id="arrs_{t}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="9" markerHeight="9" orient="auto-start-reverse">'
                   f'<path d="M0,0 L10,5 L0,10 z" fill="{c}"/></marker>')
    out.append('</defs>')
    background_index = len(out)
    out.append(f'<rect width="{W}" height="{H}" fill="white"/>')
    out.append(f'<text x="18" y="28" font-size="17" font-weight="bold">{html.escape(result["ligand"]["label"])} '
               f'&#8211; {view} view</text>')
    out.append(f'<g id="ligand">{body}</g>')

    used_atoms: dict[int, tuple] = {}
    line_els, label_els, dist_els = [], [], []
    for res, g in groups.items():
        cx, cy, bw, bh, _ = boxes[res]
        for e in g["events"]:
            t = e["type"]
            ax, ay = _anchor(e, pos)
            for a in e["ligand_atoms"]:
                if a.get("rdkit_idx") in pos:
                    used_atoms[a["rdkit_idx"]] = (a["name"], COLORS.get(t, "#444"), pos[a["rdkit_idx"]])
            x1, y1 = _clip(cx, cy, bw, bh, ax, ay)
            direction = ""
            if t == "hydrogen_bond":
                direction = "toligand" if e["donor_side"] == "protein" else "toprotein"
            col = COLORS.get(t, "#444")
            dash = DASH.get(t, "")
            sw = 1.2 if t in ("hydrophobic", "vdw_contact") else 2.0
            attrs = (f'class="ev" data-event="{e["event_id"]}" data-res="{html.escape(res)}" data-tx="{ax:.1f}" data-ty="{ay:.1f}" '
                     f'data-dir="{direction}" data-type="{t}"')
            marker = ""
            if direction == "toligand":
                marker = f'marker-end="url(#arr_{t})"'
            elif direction == "toprotein":
                marker = f'marker-start="url(#arrs_{t})"'
            line_els.append(f'<line {attrs} x1="{x1:.1f}" y1="{y1:.1f}" x2="{ax:.1f}" y2="{ay:.1f}" stroke="{col}" '
                            f'stroke-width="{sw}" stroke-dasharray="{dash}" {marker}/>')
            bits = []
            dist = e.get("distance")
            if dist is not None and dist == dist:
                bits.append(f"{dist:.2f} Å")
            if view != "publication":
                if e["type"] == "hydrogen_bond":
                    bits.append(f"∠DHA {e['angle_DHA']:.0f}°")
                if e["type"] == "pi_pi" and "plane_angle" in e:
                    bits.append(f"{e.get('geometry', '')} {e['plane_angle']:.0f}°".strip())
                if e["type"] == "water_bridge":
                    bits.append(e["water"]["label"])
            if e.get("merged_n", 1) > 1:
                bits.append(f"×{e['merged_n']}")
            if t == "water_bridge" and view == "publication":
                bits.append(e["water"]["label"])
            label = " ".join(bits)
            mx, my = (x1 + ax) / 2, (y1 + ay) / 2
            dx, dy = (layout or {}).get("dist", {}).get(e["event_id"], [0, 0])
            dist_els.append(f'<text class="dist" data-event="{e["event_id"]}" x="{mx + dx:.1f}" y="{my + dy - 4:.1f}" font-size="11" '
                            f'fill="{col}" text-anchor="middle" paint-order="stroke" stroke="white" stroke-width="3">{html.escape(label)}</text>')
    for res, g in groups.items():
        cx, cy, bw, bh, _ = boxes[res]
        hot = [e for e in g["events"] if e["type"] in ("hydrogen_bond", "ionic", "metal_coordination")]
        edge = COLORS.get(hot[0]["type"], "#555") if hot else "#555"
        tx = "".join(f'<tspan x="{cx:.1f}" dy="{0 if i == 0 else 15}" font-weight="{"bold" if i == 0 else "normal"}" '
                     f'font-size="{13 if i == 0 else 10.5}">{html.escape(t)}</tspan>' for i, t in enumerate(g["lines"]))
        label_els.append(
            f'<g class="lab" data-res="{html.escape(res)}" data-w="{bw:.1f}" data-h="{bh:.1f}" data-cx="{cx:.1f}" data-cy="{cy:.1f}" '
            f'style="cursor:move">'
            f'<rect x="{cx - bw / 2:.1f}" y="{cy - bh / 2:.1f}" width="{bw:.1f}" height="{bh:.1f}" rx="7" fill="#fbfbf4" stroke="{edge}" stroke-width="1.8"/>'
            f'<text x="{cx:.1f}" y="{cy - bh / 2 + 19:.1f}" text-anchor="middle">{tx}</text></g>')
    out += line_els
    # involved ligand atoms: coloured dot + atom name pushed away from ligand centre
    for idx, (name, col, (px, py)) in used_atoms.items():
        out.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="5" fill="none" stroke="{col}" stroke-width="2"/>')
        if view != "publication":
            vx, vy = px - centre[0], py - centre[1]
            n = math.hypot(vx, vy) or 1
            out.append(f'<text x="{px + vx / n * 20:.1f}" y="{py + vy / n * 20 + 4:.1f}" font-size="10" fill="#222" '
                       f'text-anchor="middle" font-style="italic">{html.escape(name)}</text>')
    out += dist_els + label_els

    # legend + footer
    present = sorted({e["type"] for e in events}, key=lambda t: list(COLORS).index(t) if t in COLORS else 99)
    ly = H - 36
    lx = 18
    for t in present:
        entry_width = 36 + len(TYPE_LABEL.get(t, t)) * 6.4 + 16
        if lx > 18 and lx + entry_width > W - 18:
            lx = 18
            ly += 18
        out.append(f'<line x1="{lx}" y1="{ly - 4}" x2="{lx + 30}" y2="{ly - 4}" stroke="{COLORS[t]}" stroke-width="2.4" stroke-dasharray="{DASH[t]}"/>')
        out.append(f'<text x="{lx + 36}" y="{ly}" font-size="11">{html.escape(TYPE_LABEL.get(t, t))}</text>')
        lx += entry_width
    foot = "H-bond arrow: donor → acceptor (not electron flow). "
    chem = result["chemistry"]
    foot += f'Ligand chemistry: {chem["source"]} ({chem["status"]}). '
    if hidden:
        foot += "Hidden in this view: " + "; ".join(f"{v} {k}" for k, v in hidden.items()) + ". "
    foot += "Static-structure geometry, not binding energy."
    footer_lines = textwrap.wrap(foot, width=max(20, int((W-36)/6.2)))
    footer_y = ly + 26
    render_height = max(H, footer_y + 14*len(footer_lines) + 4)
    out[0] = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {render_height}" width="{W}" height="{render_height}" font-family="{FONT}">'
    out[background_index] = f'<rect width="{W}" height="{render_height}" fill="white"/>'
    for index, line in enumerate(footer_lines):
        out.append(f'<text x="18" y="{footer_y+index*14}" font-size="10" fill="#555">{html.escape(line)}</text>')
    out.append('</svg>')
    svg = "\n".join(out)
    layout_out = {"view": view, "labels": {k: [round(v[0], 1), round(v[1], 1)] for k, v in boxes.items()},
                  "dist": (layout or {}).get("dist", {})}
    return svg, layout_out


_JS = r"""
<script>
(function(){
const svg=document.querySelector('svg');let drag=null,off=[0,0];
function clip(cx,cy,w,h,tx,ty){const dx=tx-cx,dy=ty-cy;if(!dx&&!dy)return[cx,cy];
 const sx=dx?(w/2)/Math.abs(dx):1e9,sy=dy?(h/2)/Math.abs(dy):1e9,s=Math.min(sx,sy);return[cx+dx*s,cy+dy*s];}
function pt(e){const p=svg.createSVGPoint();p.x=e.clientX;p.y=e.clientY;return p.matrixTransform(svg.getScreenCTM().inverse());}
function update(g){const cx=+g.dataset.cx,cy=+g.dataset.cy,w=+g.dataset.w,h=+g.dataset.h;
 g.setAttribute('transform',`translate(${cx-g.dataset.ox},${cy-g.dataset.oy})`);
 document.querySelectorAll(`.ev[data-res="${CSS.escape(g.dataset.res)}"]`).forEach(l=>{
  const [x,y]=clip(cx,cy,w,h,+l.dataset.tx,+l.dataset.ty);l.setAttribute('x1',x);l.setAttribute('y1',y);
  const t=document.querySelector(`.dist[data-event="${l.dataset.event}"]`);
  if(t){const dx=+(t.dataset.dx||0),dy=+(t.dataset.dy||0);t.setAttribute('x',(x+ +l.dataset.tx)/2+dx);t.setAttribute('y',(y+ +l.dataset.ty)/2+dy-4);}});}
document.querySelectorAll('.lab').forEach(g=>{g.dataset.ox=g.dataset.cx;g.dataset.oy=g.dataset.cy;
 g.addEventListener('pointerdown',e=>{drag=g;const p=pt(e);off=[+g.dataset.cx-p.x,+g.dataset.cy-p.y];g.setPointerCapture(e.pointerId);});
 g.addEventListener('pointermove',e=>{if(drag!==g)return;const p=pt(e);g.dataset.cx=p.x+off[0];g.dataset.cy=p.y+off[1];update(g);});
 g.addEventListener('pointerup',()=>{drag=null;});});
document.querySelectorAll('.dist').forEach(t=>{t.style.cursor='move';t.dataset.bx=t.getAttribute('x');t.dataset.by=t.getAttribute('y');
 let d=false,o=[0,0];t.addEventListener('pointerdown',e=>{d=true;const p=pt(e);o=[+(t.dataset.dx||0)-p.x,+(t.dataset.dy||0)-p.y];t.setPointerCapture(e.pointerId);e.stopPropagation();});
 t.addEventListener('pointermove',e=>{if(!d)return;const p=pt(e);t.dataset.dx=p.x+o[0];t.dataset.dy=p.y+o[1];
  const l=document.querySelector(`.ev[data-event="${t.dataset.event}"]`);const x1=+l.getAttribute('x1'),y1=+l.getAttribute('y1');
  t.setAttribute('x',(x1+ +l.dataset.tx)/2+ +t.dataset.dx);t.setAttribute('y',(y1+ +l.dataset.ty)/2+ +t.dataset.dy-4);});
 t.addEventListener('pointerup',()=>{d=false;});});
window.zlisLayout=function(){const L={labels:{},dist:{}};
 document.querySelectorAll('.lab').forEach(g=>L.labels[g.dataset.res]=[+(+g.dataset.cx).toFixed(1),+(+g.dataset.cy).toFixed(1)]);
 document.querySelectorAll('.dist').forEach(t=>{if(t.dataset.dx||t.dataset.dy)L.dist[t.dataset.event]=[+(+(t.dataset.dx||0)).toFixed(1),+(+(t.dataset.dy||0)).toFixed(1)];});
 return L;};
document.getElementById('save').onclick=()=>{const b=new Blob([JSON.stringify(window.zlisLayout(),null,1)],{type:'application/json'});
 const a=document.createElement('a');a.href=URL.createObjectURL(b);a.download='layout.json';a.click();};
document.getElementById('svgdl').onclick=()=>{const s=new XMLSerializer().serializeToString(svg);
 const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([s],{type:'image/svg+xml'}));a.download='edited.svg';a.click();};
})();
</script>"""


def editor_html(svg: str, result: dict, table_html: str = "") -> str:
    return (f"<!doctype html><html><head><meta charset='utf-8'><title>ZLIS 2D editor</title>"
            "<style>body{font-family:Helvetica,Arial,sans-serif;margin:12px}button{margin-right:8px;padding:6px 12px}"
            "table{border-collapse:collapse;font-size:12px}td,th{border:1px solid #ccc;padding:3px 6px}</style></head><body>"
            "<p>Drag residue boxes and distance labels. <button id='save'>Save layout.json</button>"
            "<button id='svgdl'>Download edited SVG</button></p>" + svg + table_html + _JS + "</body></html>")


def audit_table_html(result: dict) -> str:
    from .tables import events_to_frame
    df = events_to_frame(result)
    cols = ["event_id", "type", "chain", "residue", "resseq", "protein_part", "protein_atoms", "ligand_atoms",
            "distance_A", "angle_DHA_deg", "donor_side", "water", "confidence", "flags"]
    return "<h3>All events</h3>" + df[cols].to_html(index=False, escape=True, na_rep="")


def write_diagrams(result: dict, out_dir: Path, layout_path: Path | None = None) -> dict:
    out_dir = Path(out_dir)
    saved = None
    if layout_path and Path(layout_path).exists():
        saved = json.loads(Path(layout_path).read_text(encoding="utf-8"))
    files = {}
    for view in ("publication", "detailed", "audit"):
        svg, lay = render_svg(result, view, saved if saved and saved.get("view", view) == view else None)
        p = out_dir / f"interaction_2d_{view}.svg"
        p.write_text(svg, encoding="utf-8")
        files[f"{view}_svg"] = str(p)
        (out_dir / f"layout_{view}.json").write_text(json.dumps(lay, indent=1), encoding="utf-8")
        if view == "publication" or view == "detailed":
            try:
                import cairosvg
                cairosvg.svg2png(bytestring=svg.encode(), write_to=str(out_dir / f"interaction_2d_{view}.png"), scale=2.5)
                files[f"{view}_png"] = str(out_dir / f"interaction_2d_{view}.png")
                cairosvg.svg2pdf(bytestring=svg.encode(), write_to=str(out_dir / f"interaction_2d_{view}.pdf"))
                files[f"{view}_pdf"] = str(out_dir / f"interaction_2d_{view}.pdf")
            except Exception as ex:  # cairo missing: SVG is still produced
                files[f"{view}_raster_error"] = str(ex)
    svg, lay = render_svg(result, "detailed", saved)
    audit_svg, _ = render_svg(result, "audit", None)
    (out_dir / "interaction_2d_editor.html").write_text(editor_html(svg, result), encoding="utf-8")
    (out_dir / "interaction_2d_audit.html").write_text(editor_html(audit_svg, result, audit_table_html(result)), encoding="utf-8")
    files["editor_html"] = str(out_dir / "interaction_2d_editor.html")
    files["audit_html"] = str(out_dir / "interaction_2d_audit.html")
    return files
