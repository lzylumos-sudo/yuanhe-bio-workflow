"""Identity survives receptor/path and SDF-title handoffs."""
from types import SimpleNamespace
from pathlib import Path
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from zvss.ids import sha256_file
from zvss.ingest.import_run import resolve_receptor
from zvss.ingest.poses import read_sdf_poses
from zvss.select.qualify import _target_region
from zvss.config import RegionConfig


def test_exact_receptor_relative_path_and_hash(tmp_path):
    path = tmp_path / "docked.pdb"
    path.write_text("REMARK exact docking frame\n")
    ctx = SimpleNamespace(archive_input=lambda p: str(p))
    bundle = SimpleNamespace(root=str(tmp_path), receptor_file=None)
    hint = {"receptor_file": "docked.pdb", "receptor_sha256": sha256_file(path),
            "coordinate_frame": "dock_frame"}
    result = resolve_receptor(ctx, bundle, hint)
    assert result["path"] == str(path)
    assert result["sha256"] == hint["receptor_sha256"]
    assert result["coordinate_frame"] == "dock_frame"
    path.write_text("REMARK changed coordinates\n")
    with pytest.raises(ValueError, match="sha256"):
        resolve_receptor(ctx, bundle, hint)


def test_manifest_title_preserves_parent_label_without_changing_graph(tmp_path):
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=17)
    mol.SetProp("_Name", "engine output generic title")
    path = tmp_path / "poses.sdf"
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)
    discovered = SimpleNamespace(path=str(path), relpath="poses.sdf", sha256=sha256_file(path))
    records, _ = read_sdf_poses(discovered, title="input-parent-123")
    assert records[0].parent_label == "input-parent-123"
    assert records[0].mol.GetNumAtoms() == mol.GetNumAtoms()


@pytest.mark.parametrize("contacts,status", [(3, "pass"), (0, "fail"), (None, "not_evaluated")])
def test_configured_region_evaluates_actual_contacts(contacts, status):
    region = RegionConfig.from_dict({"sides": [{"name": "site", "native_residues": [50]}]})
    rows = [] if contacts is None else [{"metric": "region_contacts_total", "value": contacts}]
    ctx = SimpleNamespace(config=SimpleNamespace(target=SimpleNamespace(region=region)),
                          db=SimpleNamespace(run_id="current", query=lambda *a: rows))
    assert _target_region(ctx, {"pose_id": "pose"})["status"] == status


@pytest.mark.parametrize("payload", [{"sides": []}, {"sides": [{"name": "site"}]},
    {"contact_cutoff_A": float("nan"), "sides": [{"name": "site", "native_residues": [50]}]}])
def test_empty_or_invalid_region_is_not_usable(payload):
    assert not RegionConfig.from_dict(payload).usable
