"""Standard-library regression tests; no MD engine or third-party parser needed."""
from pathlib import Path
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from modules.automation.pdb_inventory import InventoryDecisionRequired, inventory_pdb, normalize_pdb


def atom(serial, name="C", resname="UNK", chain=" ", resid=1, element="C", charge="", xyz=(1., 2., 3.),
         record="HETATM", alt="", icode=""):
    aligned = (" " + name.ljust(3)) if len(element) == 1 and len(name) < 4 else name.ljust(4)
    return (f"{record:<6}{serial:5d} {aligned}{alt or ' '}{resname:>3} {chain}{resid:4d}{icode or ' '}   "
            + "".join(f"{v:8.3f}" for v in xyz) + f"{1:6.2f}{0:6.2f}          {element:>2}{charge:>2}")


def residue(start=1, resname="ALA", chain="A", resid=1, xyz=(0., 0., 0.), icode=""):
    return [atom(start + i, name=name, resname=resname, chain=chain, resid=resid, icode=icode,
                 element=element, record="ATOM", xyz=(xyz[0] + i * .4, xyz[1], xyz[2]))
            for i, (name, element) in enumerate((("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")))]


class InventoryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def inventory(self, lines, name="input.pdb"):
        path = self.root / name
        path.write_text("\n".join(lines + ["END"]) + "\n")
        return inventory_pdb(path)

    def test_blank_names_repeated_atoms_and_charges_are_lossless(self):
        p = self.inventory([atom(1, resname="", element="N", name="N", charge="1+"),
                            atom(2, resname="", element="N", name="N", charge="1+", xyz=(2., 2., 3.)),
                            "CONECT    1    2", "CONECT    2    1"])
        self.assertEqual(p["atom_count"], 2)
        self.assertEqual(p["components"][0]["total_formal_charge"], 2)
        r = normalize_pdb(p, {}, self.root / "out")
        self.assertEqual([a["atom_name"] for a in r["atoms"]], ["N1", "N2"])
        self.assertEqual([a["formal_charge"] for a in r["atoms"]], [1, 1])
        self.assertEqual([a["xyz"] for a in r["atoms"]], [a["xyz"] for a in p["atoms"]])
        self.assertEqual(r["retained_atom_count"], 2)
        saved = json.loads(Path(r["atom_map_path"]).read_text())
        self.assertEqual(saved["inverse_serial_to_id"], {"1": "a000001", "2": "a000002"})

    def test_filename_unicode_spaces_parentheses_uppercase_extension(self):
        p = self.inventory([atom(1), atom(2)], name="my_complex R2 (new).v1.PDB")
        self.assertEqual(p["source"]["filename"], "my_complex R2 (new).v1.PDB")
        r = normalize_pdb(p, {}, self.root / "Output (1)")
        self.assertTrue(Path(r["normalized_path"]).is_file())

    def test_bond_multiplicity_does_not_add_forward_reverse(self):
        p = self.inventory([atom(1), atom(2, element="O", name="O"), "CONECT    1    2    2", "CONECT    2    1    1"])
        self.assertEqual(len(p["bonds"]), 1)
        self.assertEqual(p["bonds"][0]["multiplicity"], 2)
        r = normalize_pdb(p, {}, self.root / "out")
        self.assertEqual(inventory_pdb(r["normalized_path"])["bonds"][0]["multiplicity"], 2)

    def test_duplicate_serial_is_not_silently_overwritten(self):
        p = self.inventory([atom(1), atom(1, element="O", name="O"), "CONECT    1    1"])
        self.assertEqual(p["atom_count"], 2)
        self.assertIn("duplicate_serial", [d["type"] for d in p["required_decisions"]])
        with self.assertRaises(InventoryDecisionRequired):
            normalize_pdb(p, {}, self.root / "out")

    def test_multiple_models_must_be_selected_never_combined(self):
        p = self.inventory(["MODEL        1", atom(1), "ENDMDL", "MODEL        2", atom(1, xyz=(9., 8., 7.)), "ENDMDL"])
        self.assertEqual(p["atom_count"], 2)
        with self.assertRaises(InventoryDecisionRequired):
            normalize_pdb(p, {}, self.root / "out")
        r = normalize_pdb(p, {"model": 2}, self.root / "selected")
        self.assertEqual(r["retained_atom_count"], 1)
        self.assertEqual(r["atoms"][0]["xyz"], [9., 8., 7.])
        mapping = json.loads(Path(r["atom_map_path"]).read_text())["atoms"]
        self.assertEqual(mapping[0]["excluded_reason"], "unselected_model")

    def test_altloc_is_residue_wide_and_all_exclusions_are_mapped(self):
        p = self.inventory([atom(1, name="C", alt="A"), atom(2, name="C", alt="B"),
                            atom(3, name="O", element="O", alt="A"), atom(4, name="O", element="O", alt="B")])
        with self.assertRaises(InventoryDecisionRequired):
            normalize_pdb(p, {}, self.root / "out")
        r = normalize_pdb(p, {"altloc": {"r0001": "A"}}, self.root / "selected")
        self.assertEqual(r["selected_atom_ids"], ["a000001", "a000003"])
        self.assertEqual(r["excluded_atom_count"], 2)
        self.assertTrue(all(a["altloc"] == "" for a in r["atoms"]))
        self.assertEqual(len(json.loads(Path(r["atom_map_path"]).read_text())["atoms"]), 4)

    def test_aliases_identified_from_residues_not_atom_record(self):
        lines = residue(resname="HIE")
        lines = ["HETATM" + line[6:] for line in lines]
        p = self.inventory(lines)
        self.assertEqual(p["components"][0]["kind"], "protein")
        r = normalize_pdb(p, {}, self.root / "out")
        self.assertEqual(r["variants"][0]["variant"], "HIE")
        self.assertTrue(all(a["resname"] == "HIS" for a in r["atoms"]))

    def test_hxt_is_preserved_and_requires_chemical_decision(self):
        p = self.inventory(residue() + [atom(5, name="HXT", element="H", resname="ALA", chain="A", record="ATOM")])
        self.assertIn("terminus_hxt_without_oxt", [d["type"] for d in p["required_decisions"]])
        r = normalize_pdb(p, {}, self.root / "out")
        self.assertIn("HXT", [a["atom_name"] for a in r["atoms"]])
        self.assertEqual(r["retained_atom_count"], 5)
        self.assertEqual(r["pending_decisions"][0]["type"], "terminus_hxt_without_oxt")

    def test_histidine_ring_hydrogens_preserve_tautomer(self):
        for label, hydrogen_names in (("HID", ["HD1"]), ("HIE", ["HE2"]), ("HIP", ["HD1", "HE2"])):
            p = self.inventory(residue(resname="HIS") + [atom(5 + i, name=name, element="H", resname="HIS", chain="A", record="ATOM") for i, name in enumerate(hydrogen_names)])
            r = normalize_pdb(p, {}, self.root / label)
            self.assertEqual(r["variants"][0]["variant"], label)
            self.assertEqual(r["variants"][0]["source"], "explicit_ring_hydrogens")

    def test_histidine_name_hydrogen_conflict_is_explicit(self):
        p = self.inventory(residue(resname="HIE") + [atom(5, name="HD1", element="H", resname="HIE", chain="A", record="ATOM")])
        r = normalize_pdb(p, {}, self.root / "conflict")
        self.assertIn("histidine_variant_conflict", [d["type"] for d in r["pending_decisions"]])
        r2 = normalize_pdb(p, {"histidine_variants": {"r0001": "HIE"}}, self.root / "resolved")
        self.assertNotIn("histidine_variant_conflict", [d["type"] for d in r2["pending_decisions"]])
        self.assertEqual(r2["variants"][0]["source"], "explicit_decision")
        self.assertEqual(r2["retained_atom_count"], 5)

    def test_unsupported_protonation_variant_is_not_reset(self):
        p = self.inventory(residue(resname="ASH"))
        r = normalize_pdb(p, {}, self.root / "out")
        self.assertTrue(all(a["resname"] == "ASH" for a in r["atoms"]))
        self.assertIn("unsupported_protonation_variant", [d["type"] for d in r["pending_decisions"]])

    def test_water_ions_multiligand_and_ter_are_classified(self):
        p = self.inventory(residue() + ["TER"] + residue(start=5, resid=2) + [
            atom(9, name="O", resname="HOH", resid=3, element="O"),
            atom(10, name="ZN", resname="ZN", resid=4, element="ZN", charge="2+"),
            atom(11, name="C", resname="UNL", resid=5), atom(12, name="C", resname="UNL", resid=6)])
        self.assertEqual([c["kind"] for c in p["components"]], ["protein", "protein", "water", "ion", "ligand", "ligand"])
        r = normalize_pdb(p, {}, self.root / "out")
        self.assertEqual(len(r["ligand_paths"]), 2)
        self.assertEqual(r["retained_atom_count"], 12)

    def test_component_exclusion_requires_explicit_selection_and_is_mapped(self):
        p = self.inventory(residue() + [atom(5, resid=2)])
        r = normalize_pdb(p, {"keep_components": ["c001"]}, self.root / "out")
        self.assertEqual(r["retained_atom_count"], 4)
        rows = json.loads(Path(r["atom_map_path"]).read_text())["atoms"]
        self.assertEqual(rows[-1]["excluded_reason"], "excluded_component")

    def test_explicit_covalent_link_and_ssbond_are_gated(self):
        p = self.inventory(residue() + [atom(5, resid=2), "CONECT    3    5", "SSBOND   1 CYS A    1    CYS A    2"])
        types = [d["type"] for d in p["required_decisions"]]
        self.assertIn("covalent_or_metal_link", types)
        self.assertIn("explicit_link_record", types)

    def test_numbering_gap_does_not_build_missing_residues(self):
        p = self.inventory(residue(resid=93) + residue(start=5, resid=97, xyz=(1.73, 0., 0.)))
        self.assertEqual(p["residue_count"], 2)
        self.assertFalse(p["has_seqres"])
        self.assertFalse(p["sequence_verified"])
        self.assertIn("no_full_sequence", [w["type"] for w in p["warnings"]])

    def test_nonfinite_and_malformed_coordinates_remain_inventory_evidence(self):
        bad = atom(1)
        bad = bad[:30] + "     nan" + bad[38:]
        p = self.inventory([bad])
        self.assertEqual(p["atom_count"], 1)
        self.assertIsNone(p["atoms"][0]["xyz"])
        with self.assertRaises(InventoryDecisionRequired):
            normalize_pdb(p, {}, self.root / "out")

    def test_no_output_overwrite(self):
        p = self.inventory([atom(1)])
        normalize_pdb(p, {}, self.root / "out")
        with self.assertRaises(FileExistsError):
            normalize_pdb(p, {}, self.root / "out")

    def test_coordinate_text_precision_is_not_rounded(self):
        line = atom(1)
        line = line[:30] + "  1.2345" + line[38:]
        p = self.inventory([line])
        r = normalize_pdb(p, {}, self.root / "out")
        self.assertEqual(inventory_pdb(r["normalized_path"])["atoms"][0]["xyz"][0], 1.2345)

    def test_many_blank_chain_waters_do_not_exhaust_chain_identifiers(self):
        p = self.inventory([atom(i, name="O", element="O", resname="HOH", resid=i) for i in range(1, 101)])
        r = normalize_pdb(p, {}, self.root / "out")
        self.assertEqual(r["retained_atom_count"], 100)
        self.assertEqual(len({a["chain"] for a in r["atoms"]}), 1)

    def test_ter_fragments_with_reused_source_chain_get_distinct_internal_chains(self):
        p = self.inventory(residue() + ["TER"] + residue(start=5, resid=2))
        r = normalize_pdb(p, {}, self.root / "out")
        self.assertEqual([res["chain"] for res in r["residues"]], ["A", "B"])
        self.assertEqual([res["original"]["chain"] for res in r["residues"]], ["A", "A"])


class RealInputRegression(unittest.TestCase):
    def test_original_complex_R2_without_preediting(self):
        # Local user-provided fixture; never uploaded to any external service.
        candidate = os.environ.get("NEXUS_R2_REGRESSION_PDB")
        path = Path(candidate) if candidate else Path(__file__).parent / "fixtures" / "complex_R2.pdb"
        if not path.exists():
            self.skipTest("Set NEXUS_R2_REGRESSION_PDB to the original complex_R2.pdb")
        original = path.read_bytes()
        p = inventory_pdb(path)
        self.assertEqual(p["atom_count"], 3094)
        protein, ligand = p["components"]
        self.assertEqual((protein["kind"], protein["atom_count"]), ("protein", 3002))
        self.assertEqual((ligand["atom_count"], ligand["heavy_atom_count"], ligand["hydrogen_count"]), (92, 45, 47))
        self.assertEqual(ligand["total_formal_charge"], 2)
        self.assertEqual(len(ligand["bonds"]), 98)
        self.assertEqual(len(p["conect_records"]), 92)
        self.assertTrue(ligand["connectivity_complete"])
        charged = [(a["serial"], a["formal_charge"]) for a in ligand["atoms"] if a["charge_explicit"]]
        self.assertEqual(charged, [(3022, 1), (3028, 1)])
        self.assertEqual([r["resseq"] for r in p["residues"] if r["resname"] == "HIE"], [151, 268])
        self.assertEqual([d["type"] for d in p["required_decisions"]], ["terminus_hxt_without_oxt"])
        with tempfile.TemporaryDirectory() as tmp:
            r = normalize_pdb(p, {}, tmp)
            self.assertEqual(r["retained_atom_count"], 3094)
            self.assertTrue(r["coordinates_unchanged"])
            self.assertEqual([a["xyz"] for a in r["atoms"]], [a["xyz"] for a in p["atoms"]])
            self.assertEqual([(v["resid"], v["variant"]) for v in r["variants"]],
                             [("95", "HID"), ("120", "HID"), ("151", "HIE"), ("181", "HID"), ("227", "HID"), ("268", "HIE")])
            ligand_out = inventory_pdb(r["ligand_paths"][ligand["id"]])
            self.assertEqual(ligand_out["atom_count"], 92)
            self.assertEqual(len({a["atom_name"] for a in ligand_out["atoms"]}), 92)
            self.assertEqual(len(ligand_out["bonds"]), 98)
            self.assertEqual(sum(a["formal_charge"] for a in ligand_out["atoms"]), 2)
        self.assertEqual(path.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
