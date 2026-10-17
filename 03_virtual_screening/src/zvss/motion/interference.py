"""Conformational interference: do repeated poses agree with each other?

For one chemical state the available poses are compared pairwise with a
symmetry-corrected heavy-atom RMSD and clustered.  A compound whose poses split
into unrelated modes is *not* a confident hit even if one mode scores well, and
the number of poses is never treated as the number of hits.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ..chem import rdio
from ..constants import AssessmentStatus, Reason
from ..util import pairwise_distances

CLUSTER_RMSD = 2.0


def _rmsd(a, b) -> float | None:
    mapping = rdio.heavy_atom_mapping(a, b)
    if not mapping:
        return None
    return rdio.rmsd_fixed_frame(a, b, mapping)


def run(ctx) -> dict[str, Any]:
    from .. import data

    poses = data.pose_rows(ctx)
    by_state: dict[str, list[dict[str, Any]]] = {}
    for pose in poses:
        by_state.setdefault(pose["state_id"], []).append(pose)

    summary: dict[str, Any] = {"states": len(by_state), "evaluated": 0, "single_pose_states": 0, "clusters": {}}
    for state_id in sorted(by_state):
        group = sorted(by_state[state_id], key=lambda p: p["pose_id"])
        mols = []
        for pose in group:
            mol = data.pose_mol(ctx, pose)
            mols.append(mol)
        usable = [(pose, mol) for pose, mol in zip(group, mols) if mol is not None and mol.GetNumAtoms() > 2]
        if len(usable) < 2:
            summary["single_pose_states"] += 1
            for pose in group:
                ctx.db.record_assessment("state", state_id, "motion", "pose_cluster_size", 1, AssessmentStatus.FLAG,
                                         Reason.SINGLE_POSE)
            continue

        n = len(usable)
        matrix = np.full((n, n), np.nan)
        np.fill_diagonal(matrix, 0.0)
        for i in range(n):
            for j in range(i + 1, n):
                value = _rmsd(usable[i][1], usable[j][1])
                if value is not None:
                    matrix[i, j] = matrix[j, i] = value
        clusters = _cluster(matrix, CLUSTER_RMSD)
        sizes = [len(c) for c in clusters]
        largest = max(sizes)
        summary["clusters"][state_id] = {
            "poses": n,
            "clusters": [sorted(usable[i][0]["pose_id"] for i in c) for c in clusters],
            "largest": largest,
            "spread": float(np.nanmax(matrix)) if n > 1 else 0.0,
        }
        summary["evaluated"] += 1
        for index, cluster in enumerate(clusters):
            for member in cluster:
                pose = usable[member][0]
                others = [matrix[member, k] for k in range(n) if k != member and not np.isnan(matrix[member, k])]
                ctx.db.record_assessment("pose", pose["pose_id"], "motion", "pose_cluster_size", len(cluster), AssessmentStatus.PASS)
                ctx.db.record_assessment("pose", pose["pose_id"], "motion", "pose_cluster_id", index, AssessmentStatus.PASS)
                ctx.db.record_assessment("pose", pose["pose_id"], "motion", "rmsd_to_nearest_pose",
                                         round(float(min(others)), 3) if others else None, AssessmentStatus.PASS)
                ctx.db.record_assessment("pose", pose["pose_id"], "motion", "pose_consistency",
                                         round(len(cluster) / n, 4), AssessmentStatus.PASS if len(cluster) == largest else AssessmentStatus.FLAG)
        for pose, mol in zip(group, mols):
            if mol is None:
                ctx.db.record_assessment("pose", pose["pose_id"], "motion", "pose_consistency", None,
                                         AssessmentStatus.NOT_EVALUATED, Reason.MISSING_COORDINATES)
    return summary


def _cluster(matrix: np.ndarray, threshold: float) -> list[list[int]]:
    """Deterministic single-linkage clustering over the RMSD matrix."""
    n = matrix.shape[0]
    parent = list(range(n))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)

    for i in range(n):
        for j in range(i + 1, n):
            value = matrix[i, j]
            if not np.isnan(value) and value <= threshold:
                union(i, j)
    groups: dict[int, list[int]] = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return [sorted(v) for _, v in sorted(groups.items())]
