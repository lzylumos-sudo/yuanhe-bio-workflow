"""Positive / negative / boundary tests for every interaction type (synthetic geometry)."""
import pytest

import helpers as H


def of(events, t):
    return [e for e in events if e["type"] == t]


# ---- hydrogen bonds ----------------------------------------------------
def test_hbond_ligand_donor_positive():
    ev = of(H.hbond_ligand_donor(3.0), "hydrogen_bond")
    assert len(ev) == 1
    e = ev[0]
    assert e["donor_side"] == "ligand"
    assert e["donor"]["uid"].startswith("M1:B:LIG1") and e["acceptor"]["uid"].startswith("M1:A:ASN1")
    assert e["distance"] == pytest.approx(3.0, abs=0.01)
    assert e["angle_DHA"] == pytest.approx(180.0, abs=1.0)


def test_hbond_protein_donor_positive():
    ev = of(H.hbond_protein_donor(3.0), "hydrogen_bond")
    assert len(ev) == 1
    assert ev[0]["donor_side"] == "protein"
    assert ev[0]["donor"]["uid"].startswith("M1:A:LYS1") and ev[0]["acceptor"]["uid"].startswith("M1:B:LIG1")


@pytest.mark.parametrize("fn", [H.hbond_ligand_donor, H.hbond_protein_donor])
def test_hbond_distance_boundary(fn):
    assert of(fn(3.49), "hydrogen_bond"), "3.49 A is inside the 3.5 A cutoff"
    assert not of(fn(3.55), "hydrogen_bond"), "3.55 A is outside"
    assert not of(fn(4.2), "hydrogen_bond")


def test_hbond_angle_boundary():
    assert of(H.hbond_ligand_donor(3.0, 135.0), "hydrogen_bond")
    assert not of(H.hbond_ligand_donor(3.0, 125.0), "hydrogen_bond"), "D-H...A below 130 degrees is rejected"


# ---- ionic -------------------------------------------------------------
def test_ionic_positive_both_orientations():
    assert of(H.ionic(4.0), "ionic")
    assert of(H.ionic(4.0, ligand_cation=False), "ionic")


def test_ionic_boundary_and_negative():
    assert of(H.ionic(4.49), "ionic")
    assert not of(H.ionic(4.6), "ionic")
    assert not of(H.ionic(4.0, same_sign=True), "ionic"), "like charges must not produce a salt bridge"


# ---- hydrophobic --------------------------------------------------------
def test_hydrophobic_boundary():
    assert of(H.hydrophobic(3.8), "hydrophobic")
    assert of(H.hydrophobic(4.4), "hydrophobic")
    assert not of(H.hydrophobic(4.6), "hydrophobic")


# ---- pi-pi / cation-pi ---------------------------------------------------
def test_pi_pi_boundary():
    assert of(H.pi_stack(3.8), "pi_pi")
    assert of(H.pi_stack(5.4), "pi_pi")
    assert not of(H.pi_stack(5.8), "pi_pi")


def test_cation_pi_boundary():
    assert of(H.cation_pi(4.0), "cation_pi")
    assert not of(H.cation_pi(4.8), "cation_pi")


def test_no_cation_pi_for_neutral_amine():
    # a neutral ligand amine over a ring must not be called cation-pi
    import numpy as np
    from rdkit import Chem
    pro = H.embed("c1ccccc1")
    p = H.xyz(pro) - H.xyz(pro)[:6].mean(axis=0)
    _, _, vt = np.linalg.svd(p[:6])
    H.set_xyz(pro, p @ H.rot_between(vt[2], np.array([0, 0, 1.0])).T)
    lig = H.embed("CN")
    n = [a.GetIdx() for a in lig.GetAtoms() if a.GetSymbol() == "N"][0]
    H.place(lig, n, (0, 0, 4.0))
    ev = H.analyse(H.make_ligand(lig, {n: "N1"}), H.make_protein(pro, "PHE", {}))
    assert not of(ev, "cation_pi")


# ---- halogen ------------------------------------------------------------
def test_halogen_boundary():
    assert of(H.halogen(3.2), "halogen_bond")
    assert of(H.halogen(3.45), "halogen_bond")
    assert not of(H.halogen(3.7), "halogen_bond")


# ---- metal --------------------------------------------------------------
def test_metal_boundary():
    assert of(H.metal(2.2), "metal_coordination")
    assert of(H.metal(2.75), "metal_coordination")
    assert not of(H.metal(3.0), "metal_coordination")


# ---- water bridge -------------------------------------------------------
def test_water_bridge_present_only_with_a_real_water():
    ev = of(H.water_bridge(2.8, 2.8), "water_bridge")
    assert len(ev) == 1
    assert ev[0]["confidence"] == "candidate"
    assert ev[0]["water"]["uid"].startswith("M1:A:HOH500")
    assert ev[0]["distance"] == pytest.approx(2.8, abs=0.01), "reported distance is one segment, not a sum"
    assert "donor" not in ev[0] or ev[0].get("donor") is None, "roles are not invented without water hydrogens"


def test_water_bridge_boundary_and_no_water():
    assert not of(H.water_bridge(3.6, 2.8), "water_bridge")
    assert not of(H.water_bridge(2.8, 3.6), "water_bridge")
    assert not of(H.hbond_ligand_donor(3.0), "water_bridge")  # no water in the input at all


# ---- event ids ------------------------------------------------------------
def test_event_ids_unique_and_ordered():
    ev = H.hbond_ligand_donor(3.0)
    ids = [e["event_id"] for e in ev]
    assert len(ids) == len(set(ids))
    assert ids == sorted(ids)
