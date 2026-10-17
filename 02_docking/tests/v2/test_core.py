import json, math, pytest
from sbvs2.paths import RuntimePaths, init_project, init_run, safe_join
from sbvs2.hashing import hash_obj, write_checksums, verify_checksums, atomic_write_json
from sbvs2.states import check_transition, IllegalTransition
from sbvs2.approvals import ApprovalLedger, require_queue_preconditions, GateError
from sbvs2.docking.box import box_from_reference, Box
from sbvs2.docking.protocol import Protocol
from sbvs2.analysis.qc import enrichment_factor, seed_stability, hard_checks
from sbvs2.analysis.ranking import rank_states, parent_summary


def test_paths_chinese_space(tmp_path):
    rp = RuntimePaths.resolve(tmp_path / "containing Space Contents")
    p = init_project(rp, "proj1"); r = init_run(rp, "proj1", "run1")
    assert r.is_dir() and (r / "00_run_record").is_dir()
    with pytest.raises(FileExistsError): init_run(rp, "proj1", "run1")     # never overwrite
    with pytest.raises(ValueError): rp.project("../x")
    with pytest.raises(ValueError): safe_join(p, "../../etc")


def test_checksums_detect_tamper(tmp_path):
    (tmp_path / "a.txt").write_text("x"); write_checksums(tmp_path)
    assert verify_checksums(tmp_path) == []
    (tmp_path / "a.txt").write_text("y")
    assert verify_checksums(tmp_path)


def test_state_machine():
    check_transition("IMPORTED", "VALIDATING")
    with pytest.raises(IllegalTransition): check_transition("IMPORTED", "RUNNING")      # no skipping approval
    with pytest.raises(IllegalTransition): check_transition("PREPARING", "QUEUED")
    check_transition("INTERRUPTED", "QUEUED")


def test_approval_invalidated_by_hash_change(tmp_path):
    L = ApprovalLedger(tmp_path / "ap.jsonl")
    h = {"smiles": "aaa"}; pol = hash_obj({"pH": 7.4})
    L.record("ligand_state", "S1", "approved", "alice", h, pol)
    assert L.is_valid("ligand_state", "S1", h, pol)
    assert not L.is_valid("ligand_state", "S1", h, hash_obj({"pH": 6.0}))     # pH changed
    assert not L.is_valid("ligand_state", "S1", {"smiles": "bbb"}, pol)       # structure changed
    L.record("ligand_state", "S1", "rejected", "alice", h, pol)
    assert not L.is_valid("ligand_state", "S1", h, pol)                        # latest wins
    with pytest.raises(ValueError): L.record("box", "b", "approved", "", {}, "x")


def test_queue_gate(tmp_path):
    L = ApprovalLedger(tmp_path / "ap.jsonl"); pol = "p"
    rec, box = ("R1", {"f": "1"}, pol), ("B1", {"b": "1"}, pol)
    st = [("S1", {"s": "1"}, pol), ("S2", {"s": "2"}, pol)]
    with pytest.raises(GateError) as e: require_queue_preconditions(L, receptor=rec, box=box, states=st, device_ok=True, budget_ok=True)
    assert e.value.code == "RECEPTOR_NOT_APPROVED"
    L.record("receptor", "R1", "approved", "u", rec[1], pol); L.record("box", "B1", "approved", "u", box[1], pol)
    with pytest.raises(GateError) as e: require_queue_preconditions(L, receptor=rec, box=box, states=st, device_ok=True, budget_ok=True)
    assert e.value.code == "STATE_REVIEW"
    L.record("ligand_state", "S1", "approved", "u", st[0][1], pol)
    ok = require_queue_preconditions(L, receptor=rec, box=box, states=st, device_ok=True, budget_ok=True)
    assert [s[0] for s in ok] == ["S1"]                                         # S2 never passes
    with pytest.raises(GateError): require_queue_preconditions(L, receptor=rec, box=box, states=st, device_ok=True, budget_ok=False)


def test_box():
    pts = [(0, 0, 0), (4, 2, 6)]
    b = box_from_reference(pts)
    assert b.size == (14.0, 12.0, 16.0) and b.center == (2.0, 1.0, 3.0) and b.contains(pts)
    assert "NON_POSITIVE_SIZE" in Box((0, 0, 0), (0, 1, 1), "manual").validate()
    assert "NON_FINITE" in Box((float("nan"), 0, 0), (1, 1, 1), "manual").validate()


def test_protocol_hash_and_argv():
    p = Protocol.preset("standard"); q = Protocol.preset("standard", exhaustiveness=17)
    assert p.hash() != q.hash() and p.hash() == Protocol.preset("standard").hash()
    argv = p.gnina_args("r.pdb", "l.sdf", Box((1, 2, 3), (20, 20, 20), "manual"), 11, "o.sdf", "o.log")
    assert all(isinstance(x, str) for x in argv) and "--cnn_scoring" in argv


def test_ef_and_na():
    assert enrichment_factor([1, 0, 0, 0, 0, 0, 0, 0, 0, 0], 0.1) == pytest.approx(10.0)
    assert enrichment_factor([0, 0, 0]) is None and enrichment_factor([1, 1]) is None     # NA, not 0
    assert seed_stability([[0, 1, 5], [1, 0, 6], [5, 6, 0]]) == "stable"
    assert seed_stability([[0, 1], [1, 0]]) == "insufficient_seeds"
    assert "H_APPROVAL_OR_PROTOCOL_MISMATCH" in hard_checks({"approval_valid": False})
    assert "H_PARTIAL_CHARGE" in hard_checks({"approval_valid": True, "charges": [0.1, float("nan")]})


def test_ranking_vina_and_gnina():
    rows = []
    for st, aff in (("S1", -7.0), ("S2", -8.0)):
        for seed in (11, 23, 47):
            rows += [{"parent_id": "P1", "state_id": st, "seed": seed, "pose_id": 1, "vina_affinity": aff + 0.1, "cnn_score": .5, "cnn_affinity": 5.0},
                     {"parent_id": "P1", "state_id": st, "seed": seed, "pose_id": 2, "vina_affinity": aff, "cnn_score": .9, "cnn_affinity": 6.0 if st == "S1" else 7.0}]
    v = rank_states(rows, "vina"); assert v[0]["state_id"] == "S2" and v[0]["robustness"] == "complete"
    g = rank_states(rows, "gnina"); assert g[0]["state_id"] == "S2" and g[0]["score"] == 7.0
    ps = parent_summary(v); assert ps[0]["tested_state_count"] == 2
    two = rank_states([r for r in rows if r["seed"] != 47], "vina"); assert two[0]["robustness"] == "incomplete_2_of_3"


def test_approval_origin_separates_human_from_automated(tmp_path):
    L = ApprovalLedger(tmp_path / "a.jsonl"); h = {"x": "1"}
    L.record("ligand_state", "S1", "approved", "redock_test_harness", h, "p", "auto", origin="automated_test_harness")
    assert L.is_valid("ligand_state", "S1", h, "p")                         # fine for harness runs
    assert not L.is_valid("ligand_state", "S1", h, "p", require_human=True)  # but never counts as a human 2D review
    with pytest.raises(ValueError): L.record("box", "b", "approved", "redock_test_harness", h, "p", origin="human_2d_review")
    with pytest.raises(ValueError): L.record("box", "b", "approved", "alice", h, "p", origin="made_up")
    L.record("ligand_state", "S2", "approved", "alice", h, "p")             # default origin = human_2d_review
    assert L.is_valid("ligand_state", "S2", h, "p", require_human=True)
