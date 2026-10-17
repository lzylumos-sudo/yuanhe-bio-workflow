"""Append-only approval ledger + hard gates (plan §15, §25, App. B).

A 2D approval binds (artifact hashes + policy hash). Changing pH / structure / charge method changes
those hashes, so the old approval silently stops matching: it is *invalidated by construction*.
"""
from __future__ import annotations
import json, time, uuid
from dataclasses import dataclass, field
from pathlib import Path

from .hashing import hash_obj


class GateError(Exception):
    def __init__(self, code: str, msg: str):
        super().__init__(f"{code}: {msg}")
        self.code = code


@dataclass
class Approval:
    review_id: str
    kind: str                    # receptor | box | ligand_state
    subject_id: str              # e.g. state_id
    decision: str                # approved | rejected
    reviewer: str
    artifact_hashes: dict
    policy_hash: str
    reason: str
    timestamp: float
    batch_rule: str | None = None
    origin: str = "legacy_unlabelled"   # human_2d_review | automated_test_harness | human_bulk_rule | legacy_unlabelled


ORIGINS = ("human_2d_review", "automated_test_harness", "human_bulk_rule")
HUMAN_ORIGINS = ("human_2d_review", "human_bulk_rule")


class ApprovalLedger:
    def __init__(self, path: str | Path):
        self.path = Path(path)

    def _all(self) -> list[dict]:
        if not self.path.exists():
            return []
        return [json.loads(l) for l in self.path.read_text(encoding="utf-8").splitlines() if l.strip()]

    def records(self) -> list[dict]:
        return self._all()

    def record(self, kind, subject_id, decision, reviewer, artifact_hashes: dict, policy_hash: str,
               reason: str = "", batch_rule: str | None = None, origin: str = "human_2d_review") -> Approval:
        if decision not in ("approved", "rejected"):
            raise ValueError("decision must be approved|rejected")
        if origin not in ORIGINS:
            raise ValueError(f"origin must be one of {ORIGINS}")
        if reviewer == "redock_test_harness" and origin != "automated_test_harness":
            raise ValueError("the test harness reviewer may only record origin=automated_test_harness")
        if not reviewer:
            raise ValueError("reviewer is required")
        a = Approval(uuid.uuid4().hex, kind, subject_id, decision, reviewer, dict(artifact_hashes),
                     policy_hash, reason, time.time(), batch_rule, origin)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.path, "a", encoding="utf-8") as f:       # append-only
            f.write(json.dumps(a.__dict__, ensure_ascii=False, sort_keys=True) + "\n")
        return a

    def latest(self, kind: str, subject_id: str) -> dict | None:
        hits = [r for r in self._all() if r["kind"] == kind and r["subject_id"] == subject_id]
        return hits[-1] if hits else None

    def is_valid(self, kind, subject_id, artifact_hashes: dict, policy_hash: str, *, require_human: bool = False) -> bool:
        r = self.latest(kind, subject_id)
        ok = bool(r and r["decision"] == "approved" and r["artifact_hashes"] == dict(artifact_hashes)
                  and r["policy_hash"] == policy_hash)
        if ok and require_human:
            ok = r.get("origin") in HUMAN_ORIGINS
        return ok


def require_queue_preconditions(ledger: ApprovalLedger, *, receptor: tuple, box: tuple,
                                states: list[tuple], device_ok: bool, budget_ok: bool,
                                require_human: bool = False) -> list[tuple]:
    """Hard gate to enter QUEUED. Each tuple = (subject_id, artifact_hashes, policy_hash).

    Returns only the states whose approval is currently valid. Raises GateError otherwise.
    Unapproved states are never returned, hence can never reach an engine.
    """
    if not ledger.is_valid("receptor", receptor[0], receptor[1], receptor[2], require_human=require_human):
        raise GateError("RECEPTOR_NOT_APPROVED", "receptor revision has no valid approval")
    if not ledger.is_valid("box", box[0], box[1], box[2], require_human=require_human):
        raise GateError("BOX_NOT_APPROVED", "box revision has no valid approval")
    ok = [s for s in states if ledger.is_valid("ligand_state", s[0], s[1], s[2], require_human=require_human)]
    if not ok:
        raise GateError("STATE_REVIEW", "no ligand state has a valid approval")
    if not device_ok:
        raise GateError("DEVICE_UNVERIFIED", "engine/device check failed")
    if not budget_ok:
        raise GateError("BUDGET_NOT_CONFIRMED", "budget not confirmed by user")
    return ok
