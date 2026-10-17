"""Docking box (plan §16). Reference-ligand boxes start at +5 A per side [P]."""
from __future__ import annotations
import math
from dataclasses import dataclass, asdict
from ..hashing import hash_obj


@dataclass
class Box:
    center: tuple
    size: tuple
    source: str                  # reference_ligand | residues | manual | pocket_detector
    revision: int = 1
    approved_reference_pose: bool = False

    def validate(self, max_edge=40.0):
        out = []
        for v in (*self.center, *self.size):
            if not math.isfinite(v):
                out.append("NON_FINITE")
        if any(s <= 0 for s in self.size):
            out.append("NON_POSITIVE_SIZE")
        if math.prod(self.size) > max_edge ** 3:
            out.append("LARGE_BOX_COST_WARNING")
        return out

    def contains(self, pts) -> bool:
        return all(abs(p[i] - self.center[i]) <= self.size[i] / 2 for p in pts for i in range(3))

    def hash(self):
        return hash_obj(asdict(self))


def box_from_reference(coords, pad=5.0) -> Box:
    if not coords:
        raise ValueError("empty reference ligand")
    lo = [min(c[i] for c in coords) for i in range(3)]
    hi = [max(c[i] for c in coords) for i in range(3)]
    b = Box(tuple(round((l + h) / 2, 3) for l, h in zip(lo, hi)),
            tuple(round(h - l + 2 * pad, 3) for l, h in zip(lo, hi)), "reference_ligand")
    assert b.contains(coords)
    return b
