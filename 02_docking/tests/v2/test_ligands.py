import math
from sbvs2.ligands.csv_io import ingest_ligands
from sbvs2.ligands.prepare import PrepSettings, prepare_parent, charge_check, estimate_jobs
from sbvs2.ligands.review import state_svg, bulk_approvable

CSV = ("compound_id,smiles,name\n"
       "A1,CCO,ethanol\nA2,CC(=O)O,acetic\nA3,CCO,dup-structure\nA3,CCN,dupid\nB1,not_a_smiles,bad\nB2,,empty\n"
       "C1,C[C@H](N)C(=O)O,L-ala\nC2,C[C@@H](N)C(=O)O,D-ala\nC3,CC(N)C(=O)O,unk-ala\nD1,CC(=O)O.[Na+],salt\nE1,C[Si](C)(C)C,si\n").encode()


def test_csv_traceability_and_failures():
    r = ingest_ligands(CSV)
    codes = {f.compound_id + ":" + f.reason_code for f in r.failures}
    assert {"B1:INVALID_SMILES", "B2:EMPTY_SMILES", "E1:UNSUPPORTED_ELEMENT"} <= codes
    assert any(f.reason_code == "DUPLICATE_ID" for f in r.failures)
    assert len(r.source_map) == 11                                         # every row accounted for
    a1 = next(s for s in r.source_map if s["compound_id"] == "A1"); a3 = [s for s in r.source_map if s["compound_id"] == "A3" and s["status"] == "ok"]
    assert a1["parent_id"] == a3[0]["parent_id"]                           # same structure merged, rows retained
    sm = {s["compound_id"]: s["parent_id"] for s in r.source_map}
    assert sm["C1"] != sm["C2"]                                            # L/D never merged
    assert any(f["flag"] == "UNKNOWN_STEREO" and f["compound_id"] == "C3" for f in r.flags)
    assert any(f["flag"] == "SALT_STRIPPED" for f in r.flags)


def test_prepare_states_charges_and_hash():
    r = ingest_ligands(b"id,smiles\nX,CC(=O)O\nY,CCN\n")
    s = PrepSettings(ph_confirmed=True, state_method="rules")      # deterministic fallback path
    states = [st for p in r.parents for st in prepare_parent(p, s)]
    acid = [x for x in states if x["parent_id"] == "P000001"]
    assert {x["formal_charge"] for x in acid} == {0, -1}
    amine = [x for x in states if x["parent_id"] == "P000002"]
    assert {x["formal_charge"] for x in amine} == {0, 1}
    for x in states:
        assert x["conformer_status"] == "OK" and x["charges"]["check"]["pass"], x["state_id"]
        assert x["charges"]["check"]["abs_diff"] <= 1e-3 and x["method"] == "rdkit_rules_v0"
    again = prepare_parent(r.parents[0], s)
    assert [x["state_hash"] for x in again] == [x["state_hash"] for x in acid]       # deterministic
    s2 = PrepSettings(ph_confirmed=True, ph_lo=2.0, ph_hi=3.0, state_method="rules")
    assert prepare_parent(r.parents[0], s2)[0]["state_hash"] != acid[0]["state_hash"]  # pH change -> new hash


def test_unknown_stereo_enumerated_and_truncation_blocks_bulk():
    r = ingest_ligands(b"id,smiles\nZ,CC(N)C(O)C(F)C(Cl)C(Br)CC(I)C(S)C\n")
    s = PrepSettings(ph_confirmed=True, max_stereo=4, max_states=2, n_conf_generate=4)
    st = prepare_parent(r.parents[0], s)
    assert any(x["truncated"] for x in st)
    ok, why = bulk_approvable(next(x for x in st if x["truncated"]), {"max_abs_formal_charge": 2})
    assert not ok and why in ("truncated_enumeration", "unknown_stereo")


def test_charge_check_rejects_nan_and_svg():
    from rdkit import Chem
    mh = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert not charge_check(mh, [float("nan")] * mh.GetNumAtoms())["pass"]
    assert not charge_check(mh, [0.5] * mh.GetNumAtoms())["pass"]
    assert "<svg" in state_svg("CC(=O)[O-]")
    assert estimate_jobs(10, 3, 2, 3)["n_jobs"] == 180


def test_molscrub_states_are_pka_based_and_tagged():
    import pytest
    from sbvs2.ligands import prepare as P
    if not P.HAVE_MOLSCRUB:
        pytest.skip("molscrub not installed")
    r = ingest_ligands(b"id,smiles\nX,CC(=O)O\nY,CCN\n")
    s = PrepSettings(ph_confirmed=True, n_conf_generate=6)
    acid = prepare_parent(r.parents[0], s); amine = prepare_parent(r.parents[1], s)
    assert {x["formal_charge"] for x in acid} == {-1}          # acetic acid, pH 6.4-8.4: carboxylate
    assert {x["formal_charge"] for x in amine} == {1}          # ethylamine: ammonium
    assert all(x["method"].startswith("molscrub_") and x["charges"]["check"]["pass"] for x in acid + amine)
    s_acid_ph = PrepSettings(ph_confirmed=True, ph_lo=1.0, ph_hi=2.0, n_conf_generate=6)
    assert {x["formal_charge"] for x in prepare_parent(r.parents[0], s_acid_ph)} == {0}   # low pH: neutral
    assert acid[0]["state_hash"] != prepare_parent(r.parents[0], s_acid_ph)[0]["state_hash"]
