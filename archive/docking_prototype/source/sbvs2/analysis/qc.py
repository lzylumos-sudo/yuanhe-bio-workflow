"""Hard QC (plan §19). H = hard, L = literature, P = tunable; the three are reported separately."""
from __future__ import annotations
import math

DEFAULT_P = {"backbone_drift_A": 0.5, "pocket_drift_A": 1.0, "stable_rmsd_A": 2.0}


def hard_checks(row: dict) -> list[str]:
    fails = []
    if not row.get("approval_valid"):
        fails.append("H_APPROVAL_OR_PROTOCOL_MISMATCH")
    if not row.get("identity_ok", True):
        fails.append("H_CHEMICAL_IDENTITY")
    if not row.get("parsable", True):
        fails.append("H_UNPARSABLE_STRUCTURE")
    q = row.get("charges")
    if q is not None and not all(math.isfinite(x) for x in q):
        fails.append("H_PARTIAL_CHARGE")
    return fails


def seed_stability(rmsd_matrix, thr=2.0):
    """P-class: >=2 of 3 valid seeds within thr of each other -> 'stable'. <3 valid seeds -> no strong claim."""
    n = len(rmsd_matrix)
    if n < 3:
        return "insufficient_seeds" if n else "no_data"
    best = 0
    for i in range(n):
        best = max(best, 1 + sum(1 for j in range(n) if j != i and rmsd_matrix[i][j] <= thr))
    return "stable" if best >= 2 else "unstable"


def enrichment_factor(labels_ranked, x=0.01):
    """EF_x = (A_top/N_top)/(A_total/N_total); N_top = ceil(x*N). Undefined -> None (NA, never 0)."""
    n = len(labels_ranked); a = sum(labels_ranked)
    if n == 0 or a == 0 or a == n:
        return None
    ntop = math.ceil(x * n)
    return (sum(labels_ranked[:ntop]) / ntop) / (a / n)
