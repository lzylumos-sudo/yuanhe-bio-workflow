"""Programmatic checks of the 2D diagram (labels in bounds, no overlaps, arrow orientation, layout round-trip)."""
import itertools
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from zlis.diagram2d import render_svg, write_diagrams
from zlis.pipeline import analyze

DATA = Path(__file__).parent / "data"
NS = "{http://www.w3.org/2000/svg}"


@pytest.fixture(scope="module")
def res(tmp_path_factory):
    return analyze(DATA / "3PTB.pdb", tmp_path_factory.mktemp("d"), allow_network=False)[0]


def boxes(root):
    out = {}
    for g in root.iter(NS + "g"):
        if g.get("class") == "lab":
            cx, cy, w, h = (float(g.get(k)) for k in ("data-cx", "data-cy", "data-w", "data-h"))
            out[g.get("data-res")] = (cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)
    return out


@pytest.mark.parametrize("view", ["publication", "detailed", "audit"])
def test_labels_inside_canvas_and_not_overlapping(res, view):
    svg = render_svg(res, view)[0] if isinstance(render_svg(res, view), tuple) else render_svg(res, view)
    root = ET.fromstring(svg)
    W, H = (float(x) for x in root.get("viewBox").split()[2:])
    b = boxes(root)
    assert b, "no residue labels drawn"
    for k, (x0, y0, x1, y1) in b.items():
        assert x0 >= 0 and y0 >= 0 and x1 <= W and y1 <= H, k
    for (k1, a), (k2, c) in itertools.combinations(b.items(), 2):
        overlap = a[0] < c[2] and c[0] < a[2] and a[1] < c[3] and c[1] < a[3]
        assert not overlap, f"{k1} overlaps {k2}"


def test_hbond_arrow_points_donor_to_acceptor(res):
    root = ET.fromstring(render_svg(res, "detailed")[0] if isinstance(render_svg(res, "detailed"), tuple) else render_svg(res, "detailed"))
    hb = {e["event_id"]: e for e in res["events"] if e["type"] == "hydrogen_bond"}
    lines = {l.get("data-event"): l for l in root.iter(NS + "line") if l.get("data-type") == "hydrogen_bond"}
    assert set(lines) == set(hb)
    for eid, l in lines.items():
        e = hb[eid]
        if e["donor_side"] == "ligand":     # arrow head at the protein end
            assert l.get("marker-start") and not l.get("marker-end"), eid
        else:                               # arrow head at the ligand end
            assert l.get("marker-end") and not l.get("marker-start"), eid


def test_no_arrows_on_non_hbond_lines(res):
    root = ET.fromstring(render_svg(res, "detailed")[0] if isinstance(render_svg(res, "detailed"), tuple) else render_svg(res, "detailed"))
    for l in root.iter(NS + "line"):
        if l.get("data-type") not in ("hydrogen_bond",):
            assert not l.get("marker-start") and not l.get("marker-end"), l.get("data-type")


def test_layout_json_roundtrip(res, tmp_path):
    write_diagrams(res, tmp_path)
    lay = json.loads((tmp_path / "layout_detailed.json").read_text())
    key = next(iter(lay["labels"]))
    lay["labels"][key] = [333.0, 777.0]
    lp = tmp_path / "edited.json"
    lp.write_text(json.dumps(lay))
    out2 = tmp_path / "again"
    out2.mkdir()
    write_diagrams(res, out2, layout_path=lp)
    root = ET.fromstring((out2 / "interaction_2d_detailed.svg").read_text(encoding="utf-8"))
    g = [g for g in root.iter(NS + "g") if g.get("data-res") == key][0]
    assert float(g.get("data-cx")) == pytest.approx(333.0, abs=0.5) and float(g.get("data-cy")) == pytest.approx(777.0, abs=0.5)


def test_all_formats_written(res, tmp_path):
    write_diagrams(res, tmp_path)
    for v in ("publication", "detailed", "audit"):
        assert (tmp_path / f"interaction_2d_{v}.svg").stat().st_size > 1000
    for v in ("publication", "detailed"):
        assert (tmp_path / f"interaction_2d_{v}.png").read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
        assert (tmp_path / f"interaction_2d_{v}.pdf").read_bytes()[:4] == b"%PDF"
    assert "layout" in (tmp_path / "interaction_2d_editor.html").read_text().lower()
    audit = (tmp_path / "interaction_2d_audit.html").read_text()
    assert all(e["event_id"] in audit for e in res["events"])
