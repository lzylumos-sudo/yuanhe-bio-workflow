"""Acceptance tests on real structures (3PTB trypsin-benzamidine, 1HSG protease-MK1) and edited copies."""
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from zlis.chemistry import build_ligand
from zlis.pipeline import analyze
from zlis.structure import classify_components, find_ligand_candidates, read_pdb, resolve_altlocs

DATA = Path(__file__).parent / "data"
ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def r3ptb(tmp_path_factory):
    out = tmp_path_factory.mktemp("r3ptb")
    return analyze(DATA / "3PTB.pdb", out, allow_network=False)[0]


def hb(res):
    return [e for e in res["events"] if e["type"] == "hydrogen_bond"]


def edit_pdb(src: Path, dst: Path, fn):
    dst.write_text("".join(fn(l) for l in src.read_text().splitlines(keepends=True)))
    return dst


# --- same residue, two H-bonds ------------------------------------------------------
def test_two_hbonds_from_same_residue_both_listed(r3ptb):
    asp = [e for e in hb(r3ptb) if e["residue"]["label"] == "A:ASP189"]
    assert len(asp) == 2
    assert {e["acceptor"]["name"] for e in asp} | {e["donor"]["name"] for e in asp} >= {"OD1", "OD2"}
    assert len({e["event_id"] for e in asp}) == 2


def test_donor_side_and_backbone_vs_sidechain_on_real_structures(r3ptb, tmp_path):
    # 3PTB: benzamidinium N-H donates to ASP189 OD1/OD2 (side chain) and GLY219 O (backbone carbonyl)
    assert {e["donor_side"] for e in hb(r3ptb)} == {"ligand"}
    gly = [e for e in hb(r3ptb) if e["residue"]["label"] == "A:GLY219"][0]
    assert gly["protein_part"] == "backbone" and gly["acceptor"]["name"] == "O"
    assert gly["donor"]["uid"].startswith("M1:A:BEN1")
    asp = [e for e in hb(r3ptb) if e["residue"]["label"] == "A:ASP189"]
    assert all(e["protein_part"] == "sidechain" for e in asp)
    # 1HSG: both directions occur (protein-donor ARG8 NH2 / ASP29 backbone N, ligand-donor ASP25 / GLY27)
    res = analyze(DATA / "1HSG.pdb", tmp_path, allow_network=False)[0]
    assert {e["donor_side"] for e in hb(res)} == {"ligand", "protein"}
    pd = [e for e in hb(res) if e["donor_side"] == "protein"]
    assert any(e["protein_part"] == "backbone" and e["donor"]["name"] == "N" for e in pd)
    assert all(e["acceptor"]["uid"].startswith("M1:B:MK1902") for e in pd)


def test_heavy_atoms_never_move(r3ptb):
    assert r3ptb["protein_prep"]["max_heavy_atom_shift_A"] == 0.0
    lig = r3ptb["_lig"]
    st = read_pdb(DATA / "3PTB.pdb")
    atoms, _ = resolve_altlocs(st.models[1])
    ben = [a for a in atoms if a.resname == "BEN"]
    conf = lig.mol.GetConformer()
    for a in ben:
        idx = next(i for i, info in enumerate(lig.atoms) if info["name"] == a.name.strip() and info["source"] == "input")
        p = conf.GetAtomPosition(idx)
        assert (p.x, p.y, p.z) == pytest.approx((a.x, a.y, a.z), abs=1e-6)


# --- atom_uid stability across stages -------------------------------------------------
def test_atom_uid_exists_in_input_pdb_and_exports(r3ptb, tmp_path):
    from zlis.pymol_scene import write_scene_files
    from zlis.tables import write_tables
    write_scene_files(r3ptb, tmp_path)
    write_tables(r3ptb, tmp_path)
    exported = json.loads((tmp_path / "interactions.json").read_text())
    st = read_pdb(DATA / "3PTB.pdb")
    uids = set()
    for a in st.models[1]:
        uids.add(f"M1:{a.chain}:{a.resname}{a.resseq}{a.icode}:{a.name.strip()}")
    for e in exported["events"]:
        for a in e["protein_atoms"] + e["ligand_atoms"]:
            base = a["uid"].split(":H")[0] if a["element"] == "H" else a["uid"]
            assert base in uids or a["element"] == "H", a["uid"]
    csv = (tmp_path / "interactions.csv").read_text(encoding="utf-8-sig")
    for e in exported["events"]:
        assert e["event_id"] in csv


# --- exported PDB: real atoms only ---------------------------------------------------
def test_exported_pdb_contains_only_real_atoms(r3ptb, tmp_path):
    from zlis.pymol_scene import write_complex_pdb
    p = write_complex_pdb(r3ptb, tmp_path)
    lines = [l for l in p.read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
    assert not any(l[17:20].strip() in ("PSE", "PSD", "XXX", "DUM") for l in lines)
    waters_out = sum(1 for l in lines if l[17:20] == "HOH")
    assert waters_out == r3ptb["n_waters"] == 62
    assert 'ligand' in r3ptb and (tmp_path / "ligand.sdf").exists()
    assert not any("pseudo" in l.lower() for l in p.read_text().splitlines())


def test_ligand_sdf_has_input_pose(r3ptb, tmp_path):
    from zlis.pymol_scene import write_complex_pdb
    write_complex_pdb(r3ptb, tmp_path)
    m = Chem.SDMolSupplier(str(tmp_path / "ligand.sdf"), removeHs=False)[0]
    assert m.GetNumAtoms() == r3ptb["_lig"].mol.GetNumAtoms()
    a = np.array(m.GetConformer().GetPositions())
    b = np.array(r3ptb["_lig"].mol.GetConformer().GetPositions())
    assert np.abs(a - b).max() < 1e-3


# --- no water ------------------------------------------------------------------------
def test_no_waters_means_no_water_bridge_and_no_water_in_pdb(tmp_path):
    nowat = edit_pdb(DATA / "3PTB.pdb", tmp_path / "nowat.pdb", lambda l: "" if l[17:20] == "HOH" else l)
    res = analyze(nowat, tmp_path / "o", allow_network=False)[0]
    assert res["n_waters"] == 0
    assert not [e for e in res["events"] if e["type"] == "water_bridge"]
    from zlis.pymol_scene import write_complex_pdb
    p = write_complex_pdb(res, tmp_path)
    assert "HOH" not in p.read_text()
    assert len(hb(res)) == 3        # direct H-bonds are unaffected


# --- chemistry: SDF atom-order mismatch / template coordinates differ ----------------
def test_sdf_template_with_different_atom_order_and_coordinates(tmp_path):
    st = read_pdb(DATA / "3PTB.pdb")
    atoms, _ = resolve_altlocs(st.models[1])
    poly, w, ions, hets = classify_components(atoms)
    cand = [c for c in find_ligand_candidates(atoms) if c.group.key[3] == "BEN"][0]
    ref = build_ligand(cand.group, library_path=ROOT / "configs" / "ligand_templates.json", allow_network=False)
    # template: same chemistry, scrambled atom order, unrelated coordinates
    m = Chem.AddHs(Chem.MolFromSmiles("NC(=[NH2+])c1ccccc1"))
    AllChem.EmbedMolecule(m, randomSeed=3)
    order = list(range(m.GetNumAtoms()))[::-1]
    m = Chem.RenumberAtoms(m, order)
    conf = m.GetConformer()
    for i in range(m.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x + 50, p.y - 20, p.z + 7))
    sdf = tmp_path / "tmpl.sdf"
    w_ = Chem.SDWriter(str(sdf)); w_.write(m); w_.close()
    lig = build_ligand(cand.group, template=str(sdf), allow_network=False)
    assert lig.chem_source == "template" and lig.formal_charge == 1
    # heavy-atom coordinates equal the PDB pose, not the template
    c0, c1 = ref.mol.GetConformer(), lig.mol.GetConformer()
    names0 = {i["name"]: k for k, i in enumerate(ref.atoms)}
    for k, info in enumerate(lig.atoms):
        if info["source"] == "input" and info["element"] != "H":
            a, b = c0.GetAtomPosition(names0[info["name"]]), c1.GetAtomPosition(k)
            assert a.Distance(b) < 1e-6
    # atom names/uids still match PDB atoms
    src = {a.name.strip() for a in cand.group.atoms}
    assert {i["name"] for i in lig.atoms if i["source"] == "input"} == src


def test_wrong_template_is_reported_not_silently_used(tmp_path):
    st = read_pdb(DATA / "3PTB.pdb")
    atoms, _ = resolve_altlocs(st.models[1])
    cand = [c for c in find_ligand_candidates(atoms) if c.group.key[3] == "BEN"][0]
    lig = build_ligand(cand.group, template="CCO", allow_network=False)
    assert lig.status == "needs_confirmation"
    assert any("did not match" in n for n in lig.notes)


# --- missing hydrogens / uncertain bond orders ----------------------------------------
def test_unknown_ligand_without_template_is_flagged(tmp_path):
    f = edit_pdb(DATA / "3PTB.pdb", tmp_path / "unk.pdb",
                 lambda l: l[:17] + "XYZ" + l[20:] if l.startswith("HETATM") and l[17:20] == "BEN" else l)
    res = analyze(f, tmp_path / "o", allow_network=False)[0]
    assert res["chemistry"]["status"] == "needs_confirmation"
    assert res["chemistry"]["source"] in ("inferred", "unverified")
    assert any("NOT confirmed" in n or "Provide an SDF" in n for n in res["chemistry"]["notes"])
    assert all("ligand_chemistry_unconfirmed" in e["flags"] for e in res["events"] if e["type"] == "hydrogen_bond")


def test_missing_hydrogens_are_added_and_reported(r3ptb):
    assert r3ptb["protein_prep"]["input_h_atoms"] == 0 and r3ptb["protein_prep"]["hydrogens_added"]
    assert any("pH-model" in n for n in r3ptb["protein_prep"]["notes"])
    assert any("hydrogens were absent" in n for n in r3ptb["chemistry"]["notes"])


# --- insertion codes, MODEL, multi-chain --------------------------------------------------
def test_insertion_code_is_carried_through(tmp_path):
    def f(l):
        if l.startswith(("ATOM", "HETATM")) and l[21] == "A" and l[22:26].strip() == "189":
            return l[:26] + "A" + l[27:]
        return l
    p = edit_pdb(DATA / "3PTB.pdb", tmp_path / "ic.pdb", f)
    res = analyze(p, tmp_path / "o", allow_network=False)[0]
    asp = [e for e in hb(res) if e["residue"]["resseq"] == 189]
    assert len(asp) == 2 and all(e["residue"]["icode"] == "A" for e in asp)
    uids = [a["uid"] for e in asp for a in e["protein_atoms"] + e["ligand_atoms"]]
    assert any("ASP189A:" in u for u in uids)


def test_multi_model_models_are_not_mixed(tmp_path):
    body = [l for l in (DATA / "3PTB.pdb").read_text().splitlines(keepends=True) if l.startswith(("ATOM", "HETATM", "TER"))]
    def shifted(l):
        if l.startswith(("ATOM", "HETATM")):
            return l[:30] + f"{float(l[30:38]) + 100:8.3f}" + l[38:]
        return l
    txt = "MODEL        1\n" + "".join(body) + "ENDMDL\nMODEL        2\n" + "".join(shifted(l) for l in body) + "ENDMDL\nEND\n"
    p = tmp_path / "mm.pdb"; p.write_text(txt)
    res = analyze(p, tmp_path / "o", allow_network=False)
    assert sorted(r["model"] for r in res) == [1, 2]
    for r in res:
        assert len(hb(r)) == 3
        assert all(a["uid"].startswith(f"M{r['model']}:") for e in r["events"] for a in e["protein_atoms"])
    assert len({r["out_dir"] for r in res}) == 2


def test_multichain_same_resseq_distinguished(tmp_path):
    res = analyze(DATA / "1HSG.pdb", tmp_path, allow_network=False)[0]
    assert res["ligand"]["label"] == "B:MK1902"
    chains = {e["residue"]["chain"] for e in res["events"]}
    assert chains == {"A", "B"}
    labels = {e["residue"]["label"] for e in res["events"]}
    assert "A:ASP25" in labels      # catalytic Asp of chain A
    assert res["chemistry"]["status"] == "needs_confirmation"    # CCD default protonation is not verified
    assert res["protein_prep"]["temporary_OXT_placeholders_removed"] == 2


# --- CLI / launcher: moved project directory, spaces and Chinese in the path ------------------
def test_moved_project_dir_with_spaces_and_chinese(tmp_path):
    dst = tmp_path / "Relocated \u7814\u7a76 ZLIS Contents"
    shutil.copytree(ROOT, dst, ignore=shutil.ignore_patterns("__pycache__", "results", ".pytest_cache", "._*"))
    out = tmp_path / "Output \u7ed3\u679c Directory"
    env_py = sys.executable
    proc = subprocess.run([env_py, "-m", "zlis", "analyze", str(DATA / "3PTB.pdb"), "--out", str(out), "--offline", "--no-png"],
                          cwd=dst, capture_output=True, text=True, env={**__import__("os").environ, "PYTHONPATH": str(dst / "src")})
    assert proc.returncode == 0, proc.stderr[-800:]
    d = json.loads((out / "3PTB__A_BEN1" / "interactions.json").read_text())
    assert len([e for e in d["events"] if e["type"] == "hydrogen_bond"]) == 3
    # template library is found relative to the moved project, so BEN keeps its protonated template
    assert d["chemistry"]["source"] == "library"
    assert (out / "3PTB__A_BEN1" / "scene.pse").exists()


# --- PyMOL: scene reopens with objects --------------------------------------------------------
def test_pse_reopens_with_objects_and_arrows(r3ptb, tmp_path):
    from zlis.pymol_scene import build_scene, write_scene_files, pymol_python
    write_scene_files(r3ptb, tmp_path)
    b = build_scene(tmp_path, no_png=True)
    assert b["ok"], b["stderr"]
    code = ("import pymol,json;pymol.finish_launching(['pymol','-qc']);from pymol import cmd;"
            f"cmd.load({str(tmp_path / 'scene.pse')!r});"
            "print(json.dumps({'objs':cmd.get_names('objects'),'complex':cmd.count_atoms('complex'),'ligand':cmd.count_atoms('ligand')}), flush=True)")
    p = subprocess.run([pymol_python(), "-c", code], capture_output=True, text=True)
    info = json.loads([l for l in p.stdout.splitlines() if l.startswith("{")][-1])
    assert "ix_hydrogen_bond" in info["objs"] and "ix_ionic" in info["objs"]
    assert info["ligand"] == r3ptb["_lig"].mol.GetNumAtoms()
    # exported PDB must not contain any pseudo-atom used for labels
    pdb_atoms = [l for l in (tmp_path / "complex_analysis.pdb").read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
    assert info["complex"] == len(pdb_atoms)


def test_scene_reads_same_json_and_has_no_interaction_logic(tmp_path):
    src = (ROOT / "src" / "zlis" / "scene_template.py").read_text()
    assert "interactions.json" in src
    for banned in ("find_pairs", "cmd.distance", "cmd.dist(", "h_bond", "polar_contacts"):
        assert banned not in src, f"PyMOL must not recompute interactions ({banned})"
