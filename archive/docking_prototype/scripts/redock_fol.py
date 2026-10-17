"""T06/T07-style real redocking of 1DRF/FOL through the v2 modules.

Start 3D comes from SMILES (never the crystal pose). Crystal coordinates are used ONLY for the box and the RMSD.
Reviewer is recorded as 'redock_test_harness' - this is an automated test approval, not a human 2D review.
"""
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "source"))
from rdkit import Chem
from sbvs2.paths import RuntimePaths
from sbvs2.project import new_project, new_run, log_state, seal
from sbvs2.hashing import atomic_write_json, sha256_file, hash_obj
from sbvs2.approvals import ApprovalLedger, require_queue_preconditions
from sbvs2.receptor.prepare import ReceptorPolicy, prepare_receptor
from sbvs2.ligands.csv_io import ingest_ligands
from sbvs2.ligands.prepare import PrepSettings, prepare_parent
from sbvs2.docking.box import box_from_reference
from sbvs2.docking.protocol import Protocol
from sbvs2.docking.pdbqt import receptor_to_pdbqt, ligand_to_pdbqt
from sbvs2.docking.run_vina import run_vina_job
from sbvs2.analysis.rmsd import symmetric_rmsd
from sbvs2.analysis.qc import seed_stability

L = Path("/Volumes/Zhiyu Li/02_docking/data")
rp = RuntimePaths.resolve()
import shutil
RESUME = sys.argv[1:3] if len(sys.argv) >= 3 else None      # usage: redock_fol.py [<project_id> <run_id>]
if RESUME:
    pid, rid = RESUME; proj = rp.project(pid); run = rp.run_dir(pid, rid)
else:
    proj = new_project(rp, "FOL_redock_1DRF"); pid = proj.name
orig = proj / "00_original"; prep = proj / "01_preparation"
if not RESUME:
    shutil.copy(L / "raw/1DRF.pdb", orig / "1DRF.pdb"); shutil.copy(L / "reference/FOL_crystal.sdf", orig / "FOL_crystal.sdf")
    (orig / "ligands.csv").write_text("compound_id,smiles\n" + "FOL_redock,Nc1nc(=O)c2nc(CNc3ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc3)cnc2[nH]1\n")
ledger = ApprovalLedger(proj / "01_preparation" / "approvals.jsonl")

# receptor
pol = ReceptorPolicy(ph=7.4, ph_confirmed=True, remove_hetero=("HOH", "FOL", "SO4"))
if RESUME:   # never re-prepare: verify the approved revision by hash
    rman = json.loads((prep / "receptor/receptor_manifest.json").read_text())
    assert sha256_file(prep / "receptor/receptor_prepared.pdb") == rman["prepared_sha256"], "RECEPTOR_HASH_MISMATCH"
    pdbqt = prep / "receptor/receptor.pdbqt"; rec_hash = {"prepared": rman["prepared_sha256"]}
else:
    rman = prepare_receptor(orig / "1DRF.pdb", pol, prep / "receptor")
    rec_hash = {"prepared": rman["prepared_sha256"]}
    ledger.record("receptor", "R1", "approved", "redock_test_harness", rec_hash, rman["policy_hash"], "automated test approval")
    pdbqt = receptor_to_pdbqt(prep / "receptor/receptor_prepared.pdb", prep / "receptor/receptor")

# box from crystal ligand
ref = Chem.MolFromMolFile(str(orig / "FOL_crystal.sdf"), removeHs=False, sanitize=False)
coords = [tuple(p) for p in ref.GetConformer().GetPositions()]
box = box_from_reference(coords, 5.0); assert not box.validate()
if not ledger.is_valid("box", "B1", {"box": box.hash()}, box.hash()):
    ledger.record("box", "B1", "approved", "redock_test_harness", {"box": box.hash()}, box.hash(), "reference ligand +5A")

# ligand
lr = ingest_ligands((orig / "ligands.csv").read_bytes())
ps = PrepSettings(ph_confirmed=True, n_conf_generate=30, n_conf_keep=3, unspecified_stereo="keep")
states = prepare_parent(lr.parents[0], ps)
proto = Protocol.preset("quick", seeds=(11, 23, 47))
ok_states = []
for s in states:
    h = {"state_hash": s["state_hash"]}
    if not ledger.is_valid("ligand_state", s["state_id"], h, ps.policy_hash()):
        ledger.record("ligand_state", s["state_id"], "approved", "redock_test_harness", h, ps.policy_hash(), "automated test approval")
    ok_states.append((s["state_id"], h, ps.policy_hash()))
ok = require_queue_preconditions(ledger, receptor=("R1", rec_hash, rman["policy_hash"]), box=("B1", {"box": box.hash()}, box.hash()),
                                 states=ok_states, device_ok=True, budget_ok=True)
inputs = {"receptor": rman["prepared_sha256"], "box": box.hash(), "states": [s["state_hash"] for s in states]}
if RESUME:
    man = json.loads((run / "00_run_record/run_manifest.json").read_text())
    assert man["protocol_hash"] == proto.hash() and man["input_hashes"] == inputs, "RESUME_PROTOCOL_OR_INPUT_MISMATCH"
    log_state(run, "RUNNING", "INTERRUPTED", "session/drive interruption"); log_state(run, "INTERRUPTED", "RUNNING", "resumed after protocol+input hash check")
else:
    run = new_run(rp, pid, proto.hash(), inputs, list(proto.seeds))
    for x, y in (("IMPORTED","VALIDATING"),("VALIDATING","PREPARING"),("PREPARING","REVIEW_REQUIRED"),("REVIEW_REQUIRED","APPROVED"),("APPROVED","BUDGET_CONFIRMED"),("BUDGET_CONFIRMED","QUEUED"),("QUEUED","RUNNING")):
        log_state(run, x, y)
results = []
for s in states:
    if s["state_id"] not in {o[0] for o in ok}: continue
    for ci, cid in enumerate(s["_keep"]):
        pq, _ = ligand_to_pdbqt(s["_mol"], cid, s["charges"]["values"])
        for seed in proto.seeds:
            t = time.time()
            jd = run / "02_work_products" / f"{s['state_id']}_c{ci}_seed{seed}"
            reused = False
            if jd.exists():   # resume: reuse only a record whose hashes all still verify; otherwise quarantine, never overwrite
                try:
                    old = json.loads((jd / "job_record.json").read_text())
                    good = (old["status"] == "ok" and old["protocol_hash"] == proto.hash() and old["box_hash"] == box.hash()
                            and old["receptor_sha256"] == sha256_file(pdbqt) and old["ligand_sha256"] == sha256_file(jd / "ligand.pdbqt")
                            and old["output_sha256"] == sha256_file(jd / "poses.pdbqt") and old["seed"] == seed
                            and (jd / "ligand.pdbqt").read_text() == pq)
                except Exception:
                    good = False
                if good:
                    rec, reused = old, True
                else:
                    q = jd.with_name(jd.name + f".quarantine_{int(time.time())}"); jd.rename(q); print("quarantined incomplete/invalid:", q.name)
            if not reused:
                rec = run_vina_job(pdbqt, pq, box, proto, seed, jd)
            row = {"state_id": s["state_id"], "formal_charge": s["formal_charge"], "conf": ci, "seed": seed, "status": rec["status"], "sec": round(time.time() - t, 1), "reused_verified": reused}
            if rec["status"] == "ok":
                from meeko import PDBQTMolecule, RDKitMolCreate
                pm = PDBQTMolecule.from_file(str(jd / "poses.pdbqt"), skip_typing=True)
                mols = RDKitMolCreate.from_pdbqt_mol(pm)[0]
                rm = [symmetric_rmsd(mols, ref, conf_pose=k) for k in range(mols.GetNumConformers())]
                row.update(poses=rec["poses"], rmsd_per_pose=rm, top1_rmsd=rm[0], best_rmsd=min(rm), best_aff=rec["poses"][0]["vina_affinity"])
            else: row["error"] = rec["error"]
            results.append(row); print(row["state_id"], ci, seed, row["status"], row.get("best_aff"), row.get("top1_rmsd"), flush=True)
log_state(run, "RUNNING", "QC")
atomic_write_json(run / "03_results" / "redock_results.json", results)
okr = [r for r in results if r["status"] == "ok"]
summary = {"n_jobs": len(results), "n_ok": len(okr), "n_failed": len(results) - len(okr),
           "top1_le_2A": sum(r["top1_rmsd"] <= 2.0 for r in okr), "best_any_pose_le_2A": sum(r["best_rmsd"] <= 2.0 for r in okr),
           "min_top1_rmsd": min(r["top1_rmsd"] for r in okr) if okr else None,
           "note": "Vina empirical score only; single target; NOT an enrichment validation."}
atomic_write_json(run / "04_analysis" / "summary.json", summary)
log_state(run, "QC", "ANALYZING"); log_state(run, "ANALYZING", "ARCHIVING"); log_state(run, "ARCHIVING", "COMPLETED", "sealed with checksums.sha256"); n = seal(run)   # seal LAST so the final log line is covered
print(json.dumps(summary, indent=1)); print("RUN", run)
