"""Regression tests for early failure and optional-output acceptance.

These tests require only Python/PyYAML. Full native regression is the Ubuntu
smoke procedure; no mock test is presented as a molecular-dynamics validation.
"""
from __future__ import annotations
import copy
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from modules.automation.forcefields import validate_config, select_forcefield
from modules.automation.backend import required_analysis_files, prepare, normalized_executable_sha256, build
from modules.utils.core import WorkflowError


def configuration():
    import yaml
    config = yaml.safe_load((ROOT / "config/default.yaml").read_text())
    config.update(workflow="md_from_complex_v2", project_name="v2_test", complex={"file": "example.pdb"},
                  ligand={"template_sdf": "ligand.sdf", "resname": "LIG"}, system_type="protein_ligand", replicates=1)
    config["MD"].update(production="4ps", equilibration="2ps", nvt="1ps", npt="1ps",
                        report_interval="1ps", checkpoint_interval="2ps", cpu_threads=2)
    return config


class EarlyValidationTests(unittest.TestCase):
    def test_baseline_remains_accepted(self):
        c = configuration()
        self.assertIs(validate_config(c), c)
        self.assertEqual(select_forcefield(c)[0], "amber14-gaff211-tip3p")

    def reject(self, section, key, value, fragment):
        c = configuration()
        c[section][key] = value
        with self.assertRaisesRegex(WorkflowError, fragment):
            validate_config(c)

    def test_equilibration_cannot_be_silently_ignored(self):
        self.reject("MD", "equilibration", "2ns", "must equal")

    def test_retained_water_rejected_before_build(self):
        self.reject("protein", "remove_water", False, "crystallographic water")

    def test_entropy_not_silently_ignored(self):
        self.reject("mmgbsa", "entropy", True, "not implemented")

    def test_forcefield_combinations_fail_closed(self):
        self.reject("MD", "water_forcefield", "charmm36/water.xml", "Unsupported force-field")
        self.reject("MD", "ligand_forcefield", "openff-2.2.1", "Unsupported force-field")

    def test_unknown_solvent_and_misspelled_option_rejected(self):
        self.reject("MD", "solvent", "ethanol", "Unimplemented")
        self.reject("MD", "temprature", 300, "misspelled")

    def test_phase_and_output_intervals_checked(self):
        self.reject("MD", "checkpoint_interval", "1.5ps", "whole reporting")
        self.reject("MD", "production", "4.5ps", "whole reporting")
        self.reject("MD", "timestep", "4fs", "2 fs")

    def test_analysis_does_not_fail_after_md_due_to_zero_frames(self):
        self.reject("analysis", "discard_ns", .0039, "fewer than two")
        self.reject("analysis", "stride", 10, "fewer than two")
        self.reject("mmgbsa", "start_ns", .0039, "fewer than two")

    def test_analysis_boundaries_exactly_two_frames(self):
        c = configuration()
        c["analysis"]["discard_ns"] = .003
        c["mmgbsa"]["start_ns"] = .003
        validate_config(c)

    def test_analysis_group_overlap_checked(self):
        self.reject("analysis", "groups", {"p1": {"chains": ["A"]}, "p2": {"chains": ["A"]}}, "overlap")

    def test_sasa_and_mmgbsa_options_checked_early(self):
        self.reject("analysis", "sasa_points", 10, ">= 20")
        self.reject("analysis", "sasa_probe_angstrom", 0, "supported range")
        self.reject("mmgbsa", "igb", 8, "igb=5")
        self.reject("mmgbsa", "max_frames", 1, ">= 2")

    def test_wrong_types_rejected(self):
        self.reject("MD", "temperature", True, "finite number")
        self.reject("analysis", "sasa", "false", "true or false")
        self.reject("protein", "ph", "7.0", "finite number")
        self.reject("MD", "cpu_threads", 2.5, "integer")

    def test_missing_construct_confirmation_stops_before_native_imports(self):
        with self.assertRaisesRegex(WorkflowError, "Confirm the intended"):
            prepare(configuration(), {"md": ROOT / "runs/unused/md"}, {"protein_path": "does-not-exist.pdb"})

    def test_hxt_requires_explicit_chemical_decision_before_native_imports(self):
        c = configuration()
        c["ligand"] = None
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder) / "protein.pdb"
            p.write_text("ATOM      1  HXT ALA A   3       1.000   2.000   3.000  1.00  0.00           O  \nEND\n")
            with self.assertRaisesRegex(WorkflowError, "never renamed to OXT"):
                prepare(c, {"md": Path(folder) / "md"}, {"protein_path": str(p), "decisions": {"construct_confirmed": True}})
            self.assertIn(" HXT ", p.read_text())
            self.assertNotIn(" OXT ", p.read_text())

    def test_multiligand_model_rejected_before_native_imports(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder) / "protein.pdb"
            p.write_text("END\n")
            plan = {"protein_path": str(p), "ligand_sdf": "ligand.sdf", "ligand_component_count": 2,
                    "decisions": {"construct_confirmed": True}}
            with self.assertRaisesRegex(WorkflowError, "exactly one noncovalent"):
                prepare(configuration(), {"md": Path(folder) / "md"}, plan)

    def test_temporary_legacy_factory_is_restored_if_build_fails(self):
        from modules.md import openmm_engine
        original_builder = openmm_engine.build_system
        original_factory = openmm_engine.ligand_template_generator
        def failed_builder(config, paths):
            self.assertIsNot(openmm_engine.ligand_template_generator, original_factory)
            raise RuntimeError("deliberate build failure")
        openmm_engine.build_system = failed_builder
        try:
            with self.assertRaisesRegex(RuntimeError, "deliberate build failure"):
                build(configuration(), {})
            self.assertIs(openmm_engine.ligand_template_generator, original_factory)
        finally:
            openmm_engine.build_system = original_builder
            openmm_engine.ligand_template_generator = original_factory

    def test_huge_duration_validation_is_constant_time(self):
        c = configuration()
        c["MD"]["production"] = "1000000000ns"
        validate_config(c)


class OptionalOutputTests(unittest.TestCase):
    def test_disabled_rmsd_rmsf_do_not_require_png(self):
        c = configuration()
        c["analysis"]["rmsd"] = c["analysis"]["rmsf"] = False
        validate_config(c)
        files = required_analysis_files(c)
        self.assertNotIn("RMSD.png", files)
        self.assertNotIn("Ligand_RMSD.png", files)
        self.assertNotIn("RMSF.png", files)
        self.assertIn("quality.json", files)

    def test_disabled_optional_analyses_not_required(self):
        c = configuration()
        for key in ("hbond", "contacts", "distance", "sasa"):
            c["analysis"][key] = False
        files = required_analysis_files(c)
        for name in ("Hbond.png", "Contact_map.png", "distance.png", "SASA.csv"):
            self.assertNotIn(name, files)
        self.assertIn("MMGBSA.csv", files)

    def test_no_ligand_plot_required_for_protein_protein(self):
        c = configuration()
        c["ligand"] = None
        c["system_type"] = "protein_protein"
        self.assertNotIn("Ligand_RMSD.png", required_analysis_files(c))
        validate_config(c)


class EnvironmentFingerprintTests(unittest.TestCase):
    def test_script_prefix_relocation_has_same_compatibility_hash(self):
        with tempfile.TemporaryDirectory() as folder:
            first, second = Path(folder) / "a", Path(folder) / "b"
            first.write_bytes(b"#!/old/runtime/bin/python\nprint(123)\n")
            second.write_bytes(b"#!/new/location/runtime/bin/python\nprint(123)\n")
            self.assertEqual(normalized_executable_sha256(first, "/old/runtime"), normalized_executable_sha256(second, "/new/location/runtime"))
            second.write_bytes(b"#!/new/location/runtime/bin/python\nprint(456)\n")
            self.assertNotEqual(normalized_executable_sha256(first, "/old/runtime"), normalized_executable_sha256(second, "/new/location/runtime"))

    def test_fixed_binary_prefix_slot_relocation_preserves_other_bytes(self):
        def binary(prefix):
            slot = prefix.encode() + b"/lib:val"
            return b"\x7fELF\x00" + slot + b"\x00" * (128-len(slot)) + b"executable instructions"
        with tempfile.TemporaryDirectory() as folder:
            first, second = Path(folder) / "a", Path(folder) / "b"
            first.write_bytes(binary("/old/runtime"))
            second.write_bytes(binary("/new/location/runtime"))
            self.assertEqual(normalized_executable_sha256(first, "/old/runtime"), normalized_executable_sha256(second, "/new/location/runtime"))


if __name__ == "__main__":
    unittest.main()
