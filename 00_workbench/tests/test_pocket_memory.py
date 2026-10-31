"""Bounded geometry must preserve the original cutoffs, ordering and tie rules."""
import itertools

import numpy as np
import pytest

from pepwedge.core import Complex, Residue
from pepwedge import pocket


def synthetic_complex():
    return Complex('synthetic', {
        1: Residue(1, 'LYS', {'NZ': np.array([0., 0., 0.])}),
        2: Residue(2, 'PHE', {'CG': np.array([3., 0., 0.])}),
        3: Residue(3, 'GLY', {'CA': np.array([0., 3., 0.])}),
    }, {1: Residue(1, 'GLY', {'CA': np.array([1., 1., 0.])})})


@pytest.mark.parametrize('chunk', [1, 7, 200000])
def test_blocked_voxels_match_scalar_reference_and_boundary_rules(monkeypatch, chunk):
    cx = synthetic_complex()
    # Independent per-point oracle for this integral grid, including boundaries.
    expected = []
    for point in itertools.product(range(-1, 4), range(-1, 4), range(-2, 3)):
        p = np.asarray(point, dtype=float)
        inside = min(np.linalg.norm(p - atom) for atom in cx.l_xyz) <= 2.0
        clash = any(np.sum((p - atom) ** 2) < 1.0 for atom in cx.r_xyz)
        if inside and not clash:
            expected.append(point)
    monkeypatch.setattr(pocket, '_DISTANCE_WORK_BYTES', 160)
    points, step = pocket.pocket_voxels(cx, step=1.0, envelope=2.0, clearance=1.0, chunk=chunk)
    np.testing.assert_array_equal(points, np.asarray(expected, dtype=float))
    assert step == 1.0
    assert (0., 1., 0.) in map(tuple, points)  # exactly at receptor clearance
    assert (1., 1., 2.) in map(tuple, points)  # exactly at envelope distance
    assert (0., 0., 0.) not in map(tuple, points)


def test_blocked_chemistry_keeps_first_atom_ties_and_input_order(monkeypatch):
    cx = synthetic_complex()
    points = np.tile([[1.5, 0., 0.], [3., 0., 0.], [0., 3., 0.], [1.5, 1.5, 0.]], (9, 1))
    monkeypatch.setattr(pocket, '_DISTANCE_WORK_BYTES', 160)
    assert pocket.label_chemistry(cx, points) == ['pos', 'aromatic', 'hydrophobic', 'pos'] * 9
    assert pocket.label_chemistry(cx, np.empty((0, 3))) == []


@pytest.mark.parametrize('chunk', [0, -1])
def test_invalid_chunk_cannot_silently_skip_receptor_clashes(chunk):
    with pytest.raises(ValueError, match='positive'):
        pocket.pocket_voxels(synthetic_complex(), step=1., envelope=2., clearance=1., chunk=chunk)
