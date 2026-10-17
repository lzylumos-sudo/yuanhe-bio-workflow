"""Deterministic identity: no dependence on process hash seed or input order."""

from __future__ import annotations

from zvss import ids


def test_sha256_text_is_stable() -> None:
    assert ids.sha256_text("abc") == ids.sha256_text("abc")
    assert ids.sha256_text("abc") != ids.sha256_text("abd")


def test_compound_id_depends_only_on_canonical_smiles() -> None:
    first = ids.compound_id("CCO")
    again = ids.compound_id("CCO")
    other = ids.compound_id("CCN")
    assert first == again
    assert first.startswith("CMP")
    assert len(first) == 19
    assert first != other


def test_chemical_state_id_separates_charge_states() -> None:
    neutral = ids.chemical_state_id("CCO", 0)
    charged = ids.chemical_state_id("CCO", -1)
    assert neutral.startswith("STA")
    assert neutral != charged
    assert ids.chemical_state_id("CCO", 0) == neutral


def test_pose_id_is_reproducible_and_distinct_per_record() -> None:
    coords = [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0)]
    digest = ids.coords_digest(coords)
    one = ids.pose_id("JOBabc", 0, digest, 1)
    two = ids.pose_id("JOBabc", 1, digest, 2)
    assert one == ids.pose_id("JOBabc", 0, digest, 1)
    assert one != two
    assert one.startswith("POS")


def test_coords_digest_ignores_sub_milliangstrom_noise() -> None:
    base = [(1.0, 2.0, 3.0)]
    nudged = [(1.0000001, 1.9999999, 3.0000002)]
    assert ids.coords_digest(base) == ids.coords_digest(nudged)


def test_selection_version_is_order_independent() -> None:
    rules = ids.sha256_text("rules")
    a = ids.selection_version(["CMP1", "CMP2", "CMP3"], rules)
    b = ids.selection_version(["CMP3", "CMP1", "CMP2"], rules)
    assert a == b
    assert a.startswith("SEL")
    assert a != ids.selection_version(["CMP1", "CMP2"], rules)


def test_cache_key_changes_with_every_component() -> None:
    base = ids.cache_key("scoring", "in", "cfg", "0.9.0", {"rdkit": "1"})
    assert base == ids.cache_key("scoring", "in", "cfg", "0.9.0", {"rdkit": "1"})
    assert base != ids.cache_key("scoring", "in2", "cfg", "0.9.0", {"rdkit": "1"})
    assert base != ids.cache_key("scoring", "in", "cfg2", "0.9.0", {"rdkit": "1"})
    assert base != ids.cache_key("scoring", "in", "cfg", "0.9.1", {"rdkit": "1"})
