"""Run state machine (plan §25). Approval, compute and archive are independent states."""
from __future__ import annotations

STATES = ["IMPORTED", "VALIDATING", "PREPARING", "REVIEW_REQUIRED", "APPROVED", "BUDGET_CONFIRMED",
          "QUEUED", "RUNNING", "QC", "ANALYZING", "ARCHIVING", "COMPLETED",
          "BLOCKED", "FAILED", "CANCELLED", "INTERRUPTED", "PENDING_ARCHIVE", "WAITING_STORAGE"]

_MAIN = ["IMPORTED", "VALIDATING", "PREPARING", "REVIEW_REQUIRED", "APPROVED", "BUDGET_CONFIRMED",
         "QUEUED", "RUNNING", "QC", "ANALYZING", "ARCHIVING", "COMPLETED"]
_ABORT = {"BLOCKED", "FAILED", "CANCELLED", "INTERRUPTED"}

TRANSITIONS: dict[str, set[str]] = {s: set() for s in STATES}
for a, b in zip(_MAIN, _MAIN[1:]):
    TRANSITIONS[a].add(b)
for s in _MAIN[:-1]:
    TRANSITIONS[s] |= _ABORT
# preparation may pause for review, or return to PREPARING after an edit invalidates approval
TRANSITIONS["PREPARING"] |= {"REVIEW_REQUIRED"}
TRANSITIONS["APPROVED"] |= {"REVIEW_REQUIRED", "PREPARING"}
TRANSITIONS["BUDGET_CONFIRMED"] |= {"REVIEW_REQUIRED"}
TRANSITIONS["BLOCKED"] |= {"REVIEW_REQUIRED", "PREPARING"}
TRANSITIONS["INTERRUPTED"] |= {"QUEUED", "RUNNING"}          # resume after hash re-check
TRANSITIONS["FAILED"] |= {"PREPARING", "QUEUED"}             # new attempt, never overwrite
TRANSITIONS["ARCHIVING"] |= {"PENDING_ARCHIVE", "WAITING_STORAGE"}
TRANSITIONS["RUNNING"] |= {"WAITING_STORAGE"}
TRANSITIONS["PENDING_ARCHIVE"] |= {"ARCHIVING", "FAILED"}
TRANSITIONS["WAITING_STORAGE"] |= {"RUNNING", "ARCHIVING", "QUEUED"}


class IllegalTransition(Exception):
    pass


def check_transition(cur: str, new: str) -> None:
    if cur not in TRANSITIONS or new not in TRANSITIONS:
        raise IllegalTransition(f"unknown state: {cur}->{new}")
    if new not in TRANSITIONS[cur]:
        raise IllegalTransition(f"illegal transition {cur} -> {new}")
