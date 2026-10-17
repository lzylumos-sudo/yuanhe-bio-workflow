"""Configuration loading: explicit failures, stable digests, honest defaults."""

from __future__ import annotations

import os

import pytest

from zvss.config import ConfigError, load_budget, load_target_profile, validate_target_profile

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIGS = os.path.join(os.path.dirname(HERE), "configs")


def test_generic_profile_declares_no_local_region() -> None:
    profile = load_target_profile(os.path.join(CONFIGS, "targets", "generic.yaml"))
    assert profile.target_id == "generic"
    assert profile.primary.pdb_id == "1DRF"
    assert profile.primary.chain == "A"
    assert profile.is_generic is True
    assert profile.region is None


def test_murj_profile_declares_construct_and_regions() -> None:
    profile = load_target_profile(os.path.join(CONFIGS, "targets", "ecmurj.yaml"))
    assert profile.uniprot == "P0AF16"
    assert profile.is_generic is False
    assert profile.construct.native_range == (3, 500)
    assert {m.auth_resseq for m in profile.construct.mutations} == {5, 12, 13}
    assert profile.region is not None
    assert [side.name for side in profile.region.sides] == ["TM2", "TM7"]


def test_profile_digest_is_reproducible() -> None:
    path = os.path.join(CONFIGS, "targets", "generic.yaml")
    assert load_target_profile(path).digest() == load_target_profile(path).digest()


def test_budget_defaults_are_bounded() -> None:
    budget = load_budget(os.path.join(CONFIGS, "budgets", "quick.yaml"))
    assert budget.md_primary_slots >= 1
    assert budget.md_primary_slots + budget.md_exploratory_slots <= 12
    assert budget.poses_per_parent_max >= 1
    assert budget.seeds


def test_missing_required_key_fails_loudly(tmp_path) -> None:
    bad = tmp_path / "bad.yaml"
    bad.write_text("target_id: broken\ntitle: x\n", encoding="utf-8")
    with pytest.raises(ConfigError):
        load_target_profile(str(bad))


def test_validation_rejects_sequence_shorter_than_construct() -> None:
    profile = load_target_profile(os.path.join(CONFIGS, "targets", "ecmurj.yaml"))
    profile.native_sequence = "MK" * 5
    assert len(profile.native_sequence) < profile.construct.native_range[1]
    with pytest.raises(ConfigError):
        validate_target_profile(profile)

def test_input_cache_detects_same_path_content_changes(tmp_path):
    from zvss.stages import _input_digest
    p = tmp_path / 'poses.sdf'
    p.write_text('first')
    inputs = {'input_dirs': [str(tmp_path)]}
    before = _input_digest(inputs)
    p.write_text('second')
    assert _input_digest(inputs) != before
