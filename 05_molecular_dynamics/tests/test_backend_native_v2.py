"""Opt-in Ubuntu preparation API regression; never launches production MD.

Run on the target environment with:
  NEXUS_NATIVE_TESTS=1 python -m unittest discover -s tests -p test_backend_native_v2.py
A six-residue piece of the local fixture is used only to exercise missing-side-
chain repair, explicit histidine tautomer preservation and a reviewed terminus.
"""
from __future__ import annotations
import os
from pathlib import Path
import platform
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from modules.automation.backend import prepare


@unittest.skipUnless(os.environ.get("NEXUS_NATIVE_TESTS") == "1" and platform.system() == "Linux", "Opt-in target Ubuntu native API test")
class NativePreparationTests(unittest.TestCase):
    def test_sidechain_histidine_and_explicit_hxt_terminus(self):
        import yaml
        from openmm import app
        fixture = ROOT / "tests/fixtures/complex_R2.pdb"
        if not fixture.is_file():
            self.skipTest("Local complex_R2 fixture not installed")
        config = yaml.safe_load((ROOT / "config/default.yaml").read_text())
        config.update(ligand=None, system_type="protein_protein", replicates=1)
        config["MD"].update(production="4ps", equilibration="2ps", nvt="1ps", npt="1ps", report_interval="1ps", checkpoint_interval="2ps")
        lines = []
        terminal = None
        for line in fixture.read_text().splitlines():
            if not line.startswith("ATOM  ") or line[21:22] != "A" or not 148 <= int(line[22:26]) <= 153:
                continue
            resid = int(line[22:26])
            if resid == 148 and line[12:16].strip() == "CB":
                continue
            if line[17:20] == "HIE":
                line = line[:17] + "HIS" + line[20:]
            lines.append(line)
            if resid == 153 and line[12:16].strip() == "O":
                terminal = line
        self.assertIsNotNone(terminal)
        # Deliberately model the anomalous O-element HXT record. Its position is
        # never reused as OXT: the explicit decision deletes it before parsing.
        hxt = terminal[:6] + f"{9999:5d}" + terminal[11:12] + " HXT" + terminal[16:]
        lines.extend([hxt, "TER", "END"])
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / "protein.pdb"
            source.write_text("\n".join(lines) + "\n")
            paths = {"md": Path(folder) / "md"}
            plan = {"protein_path": str(source), "variants": [{"chain": "A", "resid": "151", "insertion_code": "", "variant": "HIE"}],
                    "decisions": {"construct_confirmed": True, "termini": "standard_charged"}}
            record = prepare(config, paths, plan)
            topology = app.PDBFile(str(paths["md"] / "prepared/protein.pdb")).topology
            by_residue = {str(r.id): {a.name for a in r.atoms()} for r in topology.residues()}
            self.assertIn("CB", by_residue["148"])
            self.assertIn("HE2", by_residue["151"])
            self.assertNotIn("HD1", by_residue["151"])
            self.assertIn("OXT", by_residue["153"])
            self.assertNotIn("HXT", by_residue["153"])
            self.assertLessEqual(record["maximum_existing_heavy_displacement_angstrom"], 1e-4)
            self.assertEqual(len(record["explicit_chemical_removals"]), 1)
            self.assertTrue(any(a["atom"] == "CB" for a in record["added_missing_heavy_atoms"]))


if __name__ == "__main__":
    unittest.main()
