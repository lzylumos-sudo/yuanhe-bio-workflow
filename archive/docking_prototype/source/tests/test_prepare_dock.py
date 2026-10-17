import json, pathlib, pytest
from rdkit import Chem
from rdkit.Chem import AllChem

FIX = pathlib.Path(__file__).parent / "fixtures"
pytest.importorskip("openmm"); pytest.importorskip("pdbfixer"); pytest.importorskip("vina"); pytest.importorskip("meeko")

from sbvs2.receptor.prepare import ReceptorPolicy, prepare_receptor, PolicyError
from sbvs2.analysis.rmsd import symmetric_rmsd
from sbvs2.docking.run_vina import parse_vina_pdbqt, run_vina_job
from sbvs2.docking.pdbqt import receptor_to_pdbqt, ligand_to_pdbqt
from sbvs2.docking.box import box_from_reference
from sbvs2.docking.protocol import Protocol
from sbvs2.hashing import write_checksums, verify_checksums

POL = dict(ph_confirmed=True, remove_hetero=("HOH", "FOL", "SO4"))


@pytest.fixture(scope="module")
def receptor(tmp_path_factory):
    out = tmp_path_factory.mktemp("rec")
    return out, prepare_receptor(FIX / "1DRF.pdb", ReceptorPolicy(**POL), out)


def test_policy_blocks(tmp_path):
    with pytest.raises(PolicyError, match="PH_NOT_CONFIRMED"): prepare_receptor(FIX / "1DRF.pdb", ReceptorPolicy(), tmp_path)
    with pytest.raises(PolicyError, match="HETERO_UNRESOLVED"): prepare_receptor(FIX / "1DRF.pdb", ReceptorPolicy(ph_confirmed=True), tmp_path)
    with pytest.raises(PolicyError, match="HETERO_KEEP_UNSUPPORTED"):
        prepare_receptor(FIX / "1DRF.pdb", ReceptorPolicy(ph_confirmed=True, remove_hetero=("HOH", "SO4"), keep_hetero=("FOL",)), tmp_path)


def test_receptor_prepare_real_1drf(receptor):
    _, m = receptor
    assert m["removed"]["water"] == 249 and m["removed"]["hetero"] == {"SO4": 10, "FOL": 32}
    assert list(m["missing_atoms_completed"]) == ["A:GLN12"]
    assert m["n_hydrogens"] > 1000 and m["minimization"]["energy_after"] < m["minimization"]["energy_before"]
    assert not m["drift"]["warn_backbone"] and m["drift"]["backbone_rmsd_A"] < 0.5


def test_revision_hash_changes_with_policy(receptor, tmp_path):
    _, m = receptor
    m2 = prepare_receptor(FIX / "1DRF.pdb", ReceptorPolicy(ph=6.0, **POL), tmp_path)
    assert m2["revision_hash"] != m["revision_hash"] and m2["policy_hash"] != m["policy_hash"]


def test_rmsd_symmetry_and_protonation_independent():
    m = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1C(=O)O")); AllChem.EmbedMolecule(m, randomSeed=1)
    assert symmetric_rmsd(m, m) == pytest.approx(0.0, abs=1e-6)
    # swap the two equivalent ring-carbon (ortho) positions: symmetric RMSD must still be ~0
    anion = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1C(=O)[O-]")); AllChem.EmbedMolecule(anion, randomSeed=1)
    assert symmetric_rmsd(anion, anion) == pytest.approx(0.0, abs=1e-6)
    shifted = Chem.Mol(m); c = shifted.GetConformer()
    for i in range(m.GetNumAtoms()):
        p = c.GetAtomPosition(i); c.SetAtomPosition(i, (p.x + 3.0, p.y, p.z))
    assert symmetric_rmsd(shifted, m) == pytest.approx(3.0, abs=1e-6)          # bad pose is detected, not re-aligned
    with pytest.raises(ValueError): symmetric_rmsd(Chem.AddHs(Chem.MolFromSmiles("CCO")), m)


def test_parse_vina_output():
    txt = "MODEL 1\nREMARK VINA RESULT:     -9.8      0.000      0.000\nENDMDL\nMODEL 2\nREMARK VINA RESULT:     -9.1      1.500      2.200\nENDMDL\n"
    p = parse_vina_pdbqt(txt)
    assert [x["vina_affinity"] for x in p] == [-9.8, -9.1] and p[1]["rmsd_ub"] == 2.2
    assert parse_vina_pdbqt("garbage") == []


def test_real_vina_job_and_no_overwrite(receptor, tmp_path):
    out, m = receptor
    pdbqt = receptor_to_pdbqt(out / "receptor_prepared.pdb", out / "rec")
    ref = Chem.MolFromMolFile(str(FIX / "FOL_crystal.sdf"), removeHs=False, sanitize=False)
    box = box_from_reference([tuple(p) for p in ref.GetConformer().GetPositions()])
    lig = Chem.AddHs(Chem.MolFromSmiles("CC(=O)Nc1ccccc1")); AllChem.EmbedMolecule(lig, randomSeed=3)
    pq, _ = ligand_to_pdbqt(lig, 0)
    proto = Protocol.preset("quick", exhaustiveness=2, num_modes=3)
    rec = run_vina_job(pdbqt, pq, box, proto, 11, tmp_path / "job")
    assert rec["status"] == "ok" and rec["poses"] and rec["poses"][0]["vina_affinity"] < 0
    assert json.loads((tmp_path / "job" / "job_record.json").read_text())["output_sha256"] == rec["output_sha256"]
    with pytest.raises(FileExistsError): run_vina_job(pdbqt, pq, box, proto, 11, tmp_path / "job")   # job dir never reused
    with pytest.raises(ValueError, match="CHARGE_MAP_MISMATCH"): ligand_to_pdbqt(lig, 0, [0.0])


def test_checksums_ignore_os_junk(tmp_path):
    (tmp_path / "a.txt").write_text("x"); (tmp_path / "._a.txt").write_bytes(b"\0"); write_checksums(tmp_path)
    (tmp_path / "._b").write_bytes(b"junk")
    assert verify_checksums(tmp_path) == [] and "._a.txt" not in (tmp_path / "checksums.sha256").read_text()
