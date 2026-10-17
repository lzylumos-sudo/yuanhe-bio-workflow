"""Local chemical identity tests; these do not constitute Ubuntu MD validation."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from modules.automation.ligand_chemistry import recover_component

HAS_RDKIT = importlib.util.find_spec("rdkit") is not None
if HAS_RDKIT:
    from rdkit import Chem
    from rdkit.Chem import AllChem


def component_from_smiles(smiles, explicit_h=True):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=4321) != 0:
        raise AssertionError("Fixture embedding failed")
    if not explicit_h:
        mol = Chem.RemoveHs(mol)
    atoms = []
    for atom in mol.GetAtoms():
        index = atom.GetIdx()
        point = mol.GetConformer().GetAtomPosition(index)
        atoms.append({"id": f"a{index + 1:06d}", "element": atom.GetSymbol(),
                      "serial": index + 1, "name": atom.GetSymbol(),
                      "xyz": list(point), "formal_charge": atom.GetFormalCharge(),
                      "charge_explicit": atom.GetFormalCharge() != 0})
    bonds = [{"atom_ids": [atoms[b.GetBeginAtomIdx()]["id"], atoms[b.GetEndAtomIdx()]["id"]],
              "multiplicity": 1, "source": "CONECT"} for b in mol.GetBonds()]
    return {"id": "c001", "atoms": atoms, "bonds": bonds, "conect_complete": True}, mol


class EvidenceValidationTests(unittest.TestCase):
    def test_missing_atoms_never_ready(self):
        with tempfile.TemporaryDirectory() as directory:
            result = recover_component({"atoms": [], "bonds": []}, directory)
            self.assertEqual(result["status"], "WAITING")
            self.assertTrue((Path(directory) / "chemistry.json").exists())

    def test_non_finite_coordinates_rejected_before_chemistry(self):
        with tempfile.TemporaryDirectory() as directory:
            result = recover_component({"atoms": [{"id": "a1", "element": "C", "xyz": [0, float("nan"), 0]}]}, directory)
            self.assertEqual(result["reason"], "invalid_chemical_evidence")

    def test_duplicate_stable_ids_rejected(self):
        atom = {"id": "a1", "element": "C", "xyz": [0, 0, 0]}
        with tempfile.TemporaryDirectory() as directory:
            result = recover_component({"atoms": [atom, atom]}, directory)
            self.assertEqual(result["status"], "WAITING")


@unittest.skipUnless(HAS_RDKIT, "RDKit is not installed in this interpreter; chemistry dynamic tests BLOCKED")
class ChemicalRecoveryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.output = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def recover(self, component, name="job", **kwargs):
        return recover_component(component, self.output / name, **kwargs)

    def assertCoordinatesPreserved(self, component, result):
        mol = Chem.SDMolSupplier(result["sdf"], removeHs=False)[0]
        self.assertIsNotNone(mol)
        for index, atom in enumerate(component["atoms"]):
            point = mol.GetConformer().GetAtomPosition(index)
            self.assertLess(max(abs(point[axis] - atom["xyz"][axis]) for axis in range(3)), 0.00011)
        return mol

    def test_pdb_candidate_needs_exact_confirmation_and_keeps_every_atom(self):
        component, reference = component_from_smiles("CCO")
        first = self.recover(component)
        self.assertEqual(first["status"], "WAITING", first)
        self.assertEqual(first["canonical_isomeric_smiles"], "CCO")
        self.assertCoordinatesPreserved(component, first)
        wrong = self.recover(component, decision={"accept_candidate_sha256": "not-the-hash"})
        self.assertEqual(wrong["status"], "WAITING")
        confirmed = self.recover(component, decision={"accept_candidate_sha256": first["candidate_sha256"]})
        self.assertEqual(confirmed["status"], "READY", confirmed)
        self.assertTrue(confirmed["cache_hit"])
        self.assertEqual(confirmed["atom_count"], reference.GetNumAtoms())
        self.assertCoordinatesPreserved(component, confirmed)

    def test_trusted_template_checks_charge_sites(self):
        component, _ = component_from_smiles("C[NH3+]")
        result = self.recover(component, template={"smiles": "C[NH3+]"})
        self.assertEqual(result["status"], "READY", result)
        self.assertEqual(result["formal_charge"], 1)
        bad = self.recover(component, name="neutral", template={"smiles": "CN"})
        self.assertEqual(bad["status"], "WAITING", bad)

    def test_wrong_element_or_connectivity_template_rejected(self):
        component, _ = component_from_smiles("CCO")
        bad = self.recover(component, template={"smiles": "COC"})
        self.assertEqual(bad["status"], "WAITING", bad)

    def test_missing_h_needs_definition_then_template_adds_h_preserving_pose(self):
        component, _ = component_from_smiles("CCO", explicit_h=False)
        no_template = self.recover(component)
        self.assertEqual(no_template["status"], "WAITING", no_template)
        supplied = self.recover(component, template={"smiles": "CCO"})
        self.assertEqual(supplied["status"], "READY", supplied)
        self.assertEqual(supplied["added_hydrogens"], 6)
        mol = self.assertCoordinatesPreserved(component, supplied)
        self.assertEqual(Chem.MolToSmiles(Chem.RemoveHs(mol)), "CCO")

    def test_stereochemistry_positive_and_conflicting_template(self):
        component, _ = component_from_smiles("N[C@@H](C)C(=O)O")
        accepted = self.recover(component, template={"smiles": "N[C@@H](C)C(=O)O"})
        self.assertEqual(accepted["status"], "READY", accepted)
        rejected = self.recover(component, name="wrong_stereo", template={"smiles": "N[C@H](C)C(=O)O"})
        self.assertEqual(rejected["status"], "WAITING", rejected)
        self.assertIn("stereochemistry", rejected["message"])

    def test_double_bond_stereochemistry_cannot_be_relabelled(self):
        component, _ = component_from_smiles("C/C=C/C")
        accepted = self.recover(component, template={"smiles": "C/C=C/C"})
        self.assertEqual(accepted["status"], "READY", accepted)
        rejected = self.recover(component, name="cis", template={"smiles": "C/C=C\\C"})
        self.assertEqual(rejected["status"], "WAITING", rejected)

    def test_fixed_charge_fallback_does_not_move_protons_or_charge(self):
        from modules.automation.ligand_chemistry import _fixed_site_valence_candidates, _input_molecule, _normalise_component
        component, reference = component_from_smiles("C[NH+]1CCN(C)CC1")
        data = _normalise_component(component)
        candidate, checks = _fixed_site_valence_candidates(_input_molecule(data, Chem), data, Chem)
        self.assertEqual(Chem.GetFormalCharge(candidate), 1)
        self.assertEqual(Chem.MolToSmiles(Chem.RemoveHs(candidate)), Chem.MolToSmiles(Chem.RemoveHs(reference)))
        self.assertEqual(candidate.GetNumAtoms(), reference.GetNumAtoms())
        self.assertEqual([a.GetFormalCharge() for a in candidate.GetAtoms()], [a.GetFormalCharge() for a in reference.GetAtoms()])
        self.assertGreater(checks["fixed_site_valence_assignments"], 0)

    def test_cache_key_tracks_evidence_not_filename(self):
        component, _ = component_from_smiles("CCO")
        first = self.recover(component, name="Name (1)", template={"smiles": "CCO"})
        renamed = copy.deepcopy(component)
        renamed["id"] = "another_filename"
        for atom in renamed["atoms"]:
            atom["name"] = "renamed"
        second = self.recover(renamed, name="different", template={"smiles": "CCO"})
        self.assertEqual(first["cache_key"], second["cache_key"])
        self.assertTrue(second["cache_hit"])
        changed = copy.deepcopy(component)
        changed["atoms"][0]["xyz"][0] += 0.01
        third = self.recover(changed, name="moved", template={"smiles": "CCO"})
        self.assertNotEqual(first["cache_key"], third["cache_key"])

    def test_incomplete_or_disconnected_graph_cannot_infer(self):
        component, _ = component_from_smiles("CCO")
        component["bonds"] = []
        result = self.recover(component)
        self.assertEqual(result["status"], "WAITING", result)
        self.assertNotIn("sdf", result)

    def test_subprocess_timeout_is_bounded_waiting(self):
        component, _ = component_from_smiles("CCO")
        result = self.recover(component, timeout=0.0001)
        self.assertIn(result["status"], ["WAITING", "DEPENDENCY_BLOCKED"])
        self.assertNotEqual(result["status"], "READY")

    def test_modified_cache_payload_is_rejected(self):
        component, _ = component_from_smiles("CCO")
        first = self.recover(component, template={"smiles": "CCO"})
        cache = self.output / ".chemistry_cache" / (first["cache_key"] + ".json")
        value = json.loads(cache.read_text())
        value["molblock"] = value["molblock"].replace("0.0000", "9.9999", 1) + "\n"
        cache.write_text(json.dumps(value))
        second = self.recover(component, template={"smiles": "CCO"})
        self.assertEqual(second["status"], "WAITING", second)
        self.assertIn("checksum", second["message"])

    def test_original_r2_regression_when_explicit_fixture_available(self):
        fixture = os.environ.get("NEXUS_R2_REGRESSION_PDB", str(Path(__file__).parent / "fixtures" / "complex_R2.pdb"))
        if not Path(fixture).is_file():
            self.skipTest("Local user-provided R2 fixture missing; set NEXUS_R2_REGRESSION_PDB.")
        from modules.automation.pdb_inventory import inventory_pdb
        original = Path(fixture).read_bytes()
        inv = inventory_pdb(fixture)
        component = next(c for c in inv["components"] if c["kind"] == "ligand")
        result = self.recover(component)
        self.assertEqual(result["status"], "WAITING", result)
        self.assertEqual(result["input_atom_count"], 92)
        self.assertEqual(result["atom_count"], 92)
        self.assertEqual(result["formal_charge"], 2)
        molecule = self.assertCoordinatesPreserved(component, result)
        charged = [component["atoms"][a.GetIdx()]["serial"] for a in molecule.GetAtoms() if a.GetFormalCharge() == 1]
        self.assertEqual(sorted(charged), [3022, 3028])
        self.assertEqual(Path(fixture).read_bytes(), original)


if __name__ == "__main__":
    unittest.main(verbosity=2)
