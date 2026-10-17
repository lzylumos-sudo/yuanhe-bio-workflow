"""Export an *approved* redocking input set for a remote engine host (plan Mode B).
Everything is derived from the approved receptor revision + approval ledger of an existing project. The remote host never
sees /Volumes paths. Crystal ligand coordinates are shipped ONLY as reference.sdf (box + RMSD), never as a docking start.

usage: export_stage.py <project_dir> <out_dir>
"""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from rdkit import Chem
from sbvs2.hashing import sha256_file, atomic_write_json
from sbvs2.approvals import ApprovalLedger
from sbvs2.ligands.csv_io import ingest_ligands
from sbvs2.ligands.prepare import PrepSettings, prepare_parent
from sbvs2.docking.box import box_from_reference

proj, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=False)
led = ApprovalLedger(proj / "01_preparation" / "approvals.jsonl")
rman = json.loads((proj / "01_preparation/receptor/receptor_manifest.json").read_text())
rpdb = proj / "01_preparation/receptor/receptor_prepared.pdb"
assert sha256_file(rpdb) == rman["prepared_sha256"], "RECEPTOR_HASH_MISMATCH"
assert led.is_valid("receptor", "R1", {"prepared": rman["prepared_sha256"]}, rman["policy_hash"]), "RECEPTOR_NOT_APPROVED"
(out / "receptor_prepared.pdb").write_bytes(rpdb.read_bytes())
ref = Chem.MolFromMolFile(str(proj / "00_original/FOL_crystal.sdf"), removeHs=False, sanitize=False)
box = box_from_reference([tuple(p) for p in ref.GetConformer().GetPositions()], 5.0)
assert led.is_valid("box", "B1", {"box": box.hash()}, box.hash()), "BOX_NOT_APPROVED"
(out / "reference_FOL.sdf").write_bytes((proj / "00_original/FOL_crystal.sdf").read_bytes())
lr = ingest_ligands((proj / "00_original/ligands.csv").read_bytes())
ps = PrepSettings(ph_confirmed=True, n_conf_generate=30, n_conf_keep=3, unspecified_stereo="keep")
(out / "ligands").mkdir()
states = []
for s in prepare_parent(lr.parents[0], ps):
    assert led.is_valid("ligand_state", s["state_id"], {"state_hash": s["state_hash"]}, ps.policy_hash()), f"STATE_NOT_APPROVED {s['state_id']}"
    for ci, cid in enumerate(s["_keep"]):
        w = Chem.SDWriter(str(out / "ligands" / f"{s['state_id']}_c{ci}.sdf")); w.write(s["_mol"], confId=cid); w.close()
    states.append({"state_id": s["state_id"], "smiles": s["smiles"], "state_hash": s["state_hash"], "n_conf": len(s["_keep"])})
atomic_write_json(out / "stage_meta.json", {"box": {"center": list(box.center), "size": list(box.size), "source": box.source, "revision": box.revision, "approved_reference_pose": box.approved_reference_pose, "hash": box.hash()},
    "receptor_revision_hash": rman["revision_hash"], "receptor_prepared_sha256": rman["prepared_sha256"], "policy_hash": ps.policy_hash(),
    "states": states, "approvals": [r for r in led.records() if r["kind"] in ("receptor", "box", "ligand_state")],
    "note": "approvals here are redock_test_harness (automated), NOT a human 2D review"})
print("staged", len(states), "states ->", out)
