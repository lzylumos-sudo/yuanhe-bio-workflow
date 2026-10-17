from sbvs2.receptor.ingest import ingest_pdb


def atom(i, name, res, chain, n, x, el, rec="ATOM", alt=" "):
    return f"{rec:<6}{i:5d} {name:<4}{alt}{res:>3} {chain}{n:4d}    {x:8.3f}{0:8.3f}{0:8.3f}  1.00 10.00          {el:>2}"


def build(tmp_path, with_seqres=True, gap=3, near_chain=False):
    L = []
    if with_seqres:
        L.append("SEQRES   1 A   10  ALA GLY SER ALA GLY SER ALA GLY SER ALA")
    i = 1
    for n in range(1, 11):
        if 3 < n <= 3 + gap:
            continue
        for nm, el, dx in (("N", "N", 0), ("CA", "C", 1.46), ("C", "C", 2.5), ("O", "O", 3)):
            L.append(atom(i, nm, "ALA", "A", n, n * 3.8 + dx * 0.0 + (0 if nm == "CA" else 0.1), el)); i += 1
        L.append(atom(i, "CB", "ALA", "A", n, n * 3.8, "C")); i += 1
    L.append(atom(i, "O", "HOH", "A", 101, 1.0, "O", "HETATM")); i += 1
    L.append(atom(i, "ZN", "ZN", "A", 102, 2.0, "ZN", "HETATM")); i += 1
    p = tmp_path / "t.pdb"; p.write_text("\n".join(L) + "\nEND\n"); return p


def test_gap_with_seqres_is_short_loop(tmp_path):
    r = ingest_pdb(build(tmp_path, gap=3))
    assert r.gaps and r.gaps[0]["n_missing_by_numbering"] == 3
    assert any(p["case"] == "SHORT_LOOP" and p["auto"] for p in r.plan)


def test_long_gap_requires_review(tmp_path):
    r = ingest_pdb(build(tmp_path, gap=6))
    assert any(i["code"] == "LOOP_REVIEW" for i in r.issues)


def test_no_seqres_never_guessed(tmp_path):
    r = ingest_pdb(build(tmp_path, with_seqres=False))
    assert any(i["code"] == "MISSING_SEQUENCE" for i in r.issues)
    assert not any(p["auto"] for p in r.plan if p["case"] == "MISSING_SEQUENCE")


def test_hetero_not_blanket_deleted(tmp_path):
    r = ingest_pdb(build(tmp_path))
    assert any(i["code"] == "HETERO_METAL" for i in r.issues)
    assert r.hetero["water"][0]["resname"] == "HOH"
