"""Ranking (plan §18). Different-unit scores are never summed into a total score."""
from __future__ import annotations
import statistics


def _med(v): return statistics.median(v) if v else None


def rank_states(rows, engine: str):
    """rows: dicts with parent_id,state_id,seed,pose_id,vina_affinity|cnn_score,cnn_affinity, valid(bool).

    Representative pose per (state,seed): Vina -> lowest affinity; GNINA -> highest CNNscore.
    State order: Vina by median affinity asc; GNINA by median representative-pose CNNaffinity desc.
    Ties broken by (parent,state,pose) ids; raw values never truncated.
    """
    by = {}
    for r in rows:
        if r.get("valid", True):
            by.setdefault((r["parent_id"], r["state_id"], r["seed"]), []).append(r)
    reps = {}
    for k, v in by.items():
        reps[k] = min(v, key=lambda r: (r["vina_affinity"], r["pose_id"])) if engine == "vina" \
            else sorted(v, key=lambda r: (-r["cnn_score"], str(r["pose_id"])))[0]
    states = {}
    for (p, s, seed), r in reps.items():
        states.setdefault((p, s), []).append(r)
    out = []
    for (p, s), rs in states.items():
        n = len(rs)
        score = _med([x["vina_affinity"] for x in rs]) if engine == "vina" else _med([x["cnn_affinity"] for x in rs])
        out.append({"parent_id": p, "state_id": s, "n_valid_seeds": n, "score": score,
                    "robustness": "complete" if n >= 3 else ("incomplete_2_of_3" if n == 2 else "no_stability_claim")})
    out.sort(key=lambda r: ((r["score"] if engine == "vina" else -r["score"]), r["parent_id"], r["state_id"]))
    for i, r in enumerate(out, 1):
        r["state_rank"] = i
    return out


def parent_summary(state_rows):
    """Parent = best valid state; also reports tested_state_count and rank range (no fake population weighting)."""
    by = {}
    for r in state_rows:
        by.setdefault(r["parent_id"], []).append(r)
    res = []
    for p, rs in by.items():
        rs.sort(key=lambda r: r["state_rank"])
        res.append({"parent_id": p, "best_state": rs[0]["state_id"], "best_score": rs[0]["score"],
                    "tested_state_count": len(rs), "state_rank_range": [rs[0]["state_rank"], rs[-1]["state_rank"]],
                    "robustness": rs[0]["robustness"]})
    res.sort(key=lambda r: (r["state_rank_range"][0], r["parent_id"]))
    return res
