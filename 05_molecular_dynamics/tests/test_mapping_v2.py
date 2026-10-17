"""Synthetic topology fixtures verify mappings without native Linux MD imports."""
from pathlib import Path
import csv
import json
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from modules.automation.pdb_inventory import inventory_pdb, normalize_pdb
from modules.automation.analysis_mapping import MappingError, add_original_identifiers


def atom(serial, name, resname="ALA", chain="A", resid=93, element="C", xyz=(1., 2., 3.), record="ATOM"):
    aligned = " " + name.ljust(3) if len(element) == 1 and len(name) < 4 else name.ljust(4)
    return (f"{record:<6}{serial:5d} {aligned} {resname:>3} {chain}{resid:4d}    "
            + "".join(f"{v:8.3f}" for v in xyz) + f"{1:6.2f}{0:6.2f}          {element:>2}  ")


class MappingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.job = Path(self.temp.name)
        self.paths = {name: self.job / name for name in ("md", "analysis", "reports")}
        for path in self.paths.values():
            path.mkdir()
        (self.paths["md"] / "system").mkdir()
        (self.paths["md"] / "prepared").mkdir()
        lines = [atom(1, "N", element="N", xyz=(0., 0., 0.)), atom(2, "CA", xyz=(1., 0., 0.)),
                 atom(3, "C", xyz=(2., 0., 0.)), atom(4, "O", element="O", xyz=(3., 0., 0.)),
                 atom(5, "H", element="H", xyz=(0., 1., 0.)),
                 atom(6, "C", resname="", chain=" ", resid=900, xyz=(2., 4., 1.), record="HETATM"),
                 atom(7, "O", resname="", chain=" ", resid=900, element="O", xyz=(3., 4., 1.), record="HETATM"),
                 "CONECT    6    7", "CONECT    7    6", "END"]
        raw = self.job / "source.pdb"
        raw.write_text("\n".join(lines) + "\n")
        self.inv = inventory_pdb(raw)
        self.plan = normalize_pdb(self.inv, {}, self.job / "normalized")
        self.plan["decisions"] = {"construct_confirmed": True}
        chemistry = self.job / "chemical_atom_map.json"
        chemistry.write_text(json.dumps([{"original_id": "a000006", "sdf_atom_index": 0, "element": "C", "added": False},
                                        {"original_id": "a000007", "sdf_atom_index": 1, "element": "O", "added": False}]))
        self.plan["chemistry"] = {"atom_map": str(chemistry)}
        self.plan["ligand_sdf"] = "confirmed.sdf"
        prep = {"explicit_chemical_removals": [], "added_heavy_atom_identities": [["A", "93", "", "OXT"]]}
        (self.paths["md"] / "prepared" / "protein_preparation.json").write_text(json.dumps(prep))
        # Uniformly shift the entire complex +10 A, retaining original protein
        # heavy coordinates and ligand pose. Protein H has been regenerated.
        self.final = [atom(1, "N", element="N", xyz=(10., 0., 0.)), atom(2, "CA", xyz=(11., 0., 0.)),
                      atom(3, "C", xyz=(12., 0., 0.)), atom(4, "O", element="O", xyz=(13., 0., 0.)),
                      atom(5, "H", element="H", xyz=(10., 1.1, .2)), atom(6, "OXT", element="O", xyz=(12., -1., 0.)),
                      atom(7, "C1", resname="LIG", chain="Z", resid=1, xyz=(12., 4., 1.), record="HETATM"),
                      atom(8, "O2", resname="LIG", chain="Z", resid=1, element="O", xyz=(13., 4., 1.), record="HETATM"),
                      atom(9, "O", resname="HOH", chain="W", resid=1, element="O", xyz=(20., 20., 20.), record="HETATM")]
        self.meta = {"atoms": 9, "elements": ["N", "C", "C", "O", "H", "O", "C", "O", "O"],
                     "residues": ["A:ALA93"] * 6 + ["Z:LIG1"] * 2 + ["W:HOH1"],
                     "protein_indices": list(range(6)), "ligand_indices": [6, 7]}
        self.write_topology()

    def tearDown(self):
        self.temp.cleanup()

    def write_topology(self):
        base = self.paths["md"] / "system"
        (base / "topology.pdb").write_text("\n".join(self.final + ["END"]) + "\n")
        (base / "topology.json").write_text(json.dumps(self.meta))

    def test_all_originals_accounted_for_without_false_hydrogen_identity(self):
        result = add_original_identifiers(self.job, self.plan, self.paths)
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["mapped_original_atoms"], 6)
        report = json.loads(Path(result["mapping_path"]).read_text())
        self.assertEqual(len(report["source_atoms"]), 7)
        self.assertEqual(report["uniform_translation_angstrom"], [10., 0., 0.])
        self.assertEqual(report["source_atoms"][4]["status"], "removed_before_protein_hydrogen_regeneration")
        hydrogen = report["solute_atoms"][4]
        self.assertIsNone(hydrogen["original_atom_id"])
        self.assertEqual(hydrogen["name_correspondence_original_id"], "a000005")
        self.assertEqual(report["solute_atoms"][5]["identity_status"], "audited_added_protein_heavy_atom")
        self.assertEqual(report["generated_bulk_index_ranges_inclusive"], [[8, 8]])
        self.assertEqual(report["residue_mapping"]["Z:LIG1"]["original_labels"], ["<blank>:<blank>900"])

    def test_relative_ligand_pose_change_is_rejected(self):
        self.final[6] = atom(7, "C1", resname="LIG", chain="Z", resid=1, xyz=(13., 4., 1.), record="HETATM")
        self.write_topology()
        with self.assertRaisesRegex(MappingError, "relative coordinates"):
            add_original_identifiers(self.job, self.plan, self.paths)
        self.assertFalse((self.job / "mappings" / "topology_atom_map.json").exists())

    def test_reordered_ligand_atoms_are_rejected(self):
        self.meta["ligand_indices"] = [7, 6]
        self.write_topology()
        with self.assertRaisesRegex(MappingError, "element order"):
            add_original_identifiers(self.job, self.plan, self.paths)

    def test_unknown_heavy_atom_cannot_be_called_new_without_audit(self):
        path = self.paths["md"] / "prepared" / "protein_preparation.json"
        path.write_text(json.dumps({"explicit_chemical_removals": [], "added_heavy_atom_identities": []}))
        with self.assertRaisesRegex(MappingError, "without preparation audit"):
            add_original_identifiers(self.job, self.plan, self.paths)

    def test_analysis_tables_add_original_labels_and_atom_identifiers(self):
        out = self.paths["analysis"] / "replica_01"
        out.mkdir()
        (out / "rmsf.csv").write_text("residue,CA_RMSF_A\nA:ALA93,0.5\n")
        (out / "contact_map.csv").write_text("interface,residue_A,residue_B,frames,occupancy\nprotein_ligand,A:ALA93,Z:LIG1,2,0.5\n")
        (out / "hbond_occupancy.csv").write_text("donor_index,acceptor_index,donor_residue,acceptor_residue,frames,occupancy\n0,7,A:ALA93,Z:LIG1,2,0.5\n")
        original = (out / "rmsf.csv").read_bytes()
        result = add_original_identifiers(self.job, self.plan, self.paths)
        self.assertEqual(result["analysis_tables_written"], 3)
        with (out / "rmsf_original.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(json.loads(rows[0]["residue_original_residue_ids"]), ["r0001"])
        self.assertEqual(rows[0]["CA_RMSF_A"], "0.5")
        with (out / "hbond_occupancy_original.csv").open() as handle:
            row = next(csv.DictReader(handle))
        self.assertEqual(row["donor_original_atom_id"], "a000001")
        self.assertEqual(row["acceptor_original_atom_id"], "a000007")
        self.assertEqual(row["acceptor_original_atom_name"], "O")
        self.assertEqual((out / "rmsf.csv").read_bytes(), original)

    def test_unknown_analysis_residue_is_rejected(self):
        out = self.paths["analysis"] / "replica_01"
        out.mkdir()
        (out / "rmsf.csv").write_text("residue,CA_RMSF_A\nB:ALA99,0.5\n")
        with self.assertRaisesRegex(MappingError, "no verified original mapping"):
            add_original_identifiers(self.job, self.plan, self.paths)

    def test_identical_second_run_is_idempotent(self):
        result = add_original_identifiers(self.job, self.plan, self.paths)
        first = Path(result["mapping_path"]).read_bytes()
        add_original_identifiers(self.job, self.plan, self.paths)
        self.assertEqual(Path(result["mapping_path"]).read_bytes(), first)


if __name__ == "__main__":
    unittest.main()
