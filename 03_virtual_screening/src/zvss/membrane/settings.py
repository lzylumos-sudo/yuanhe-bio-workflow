"""Resolve Y-MAVS settings from the target profile and stage defaults.

Scientific thresholds live in YAML.  The target profile wins over the shared
``stages.yaml`` defaults so a target can tighten a rule without changing other
targets.  This module only merges dictionaries; it never invents a measured
value.
"""

from __future__ import annotations

import copy
import hashlib
import json
from typing import Any


def deep_merge(base: dict[str, Any] | None, override: dict[str, Any] | None) -> dict[str, Any]:
    result: dict[str, Any] = copy.deepcopy(base or {})
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def settings_for(ctx) -> dict[str, Any]:
    """Return the effective membrane-screening settings for this run.

    ``screening_mode`` is an explicit target-level switch.  A filename, a PDB
    code, or the word "MurJ" is never used to infer the mode.
    """
    target = ctx.config.target
    mode = str(getattr(target, "screening_mode", "generic") or "generic")
    base = ctx.config.stages.get("membrane_screening", default={}) or {}
    target_block = getattr(target, "membrane_screening", None) or {}
    merged = deep_merge(base, target_block)
    enabled = bool(merged.get("enabled", False)) and mode == "membrane_aware"
    merged["enabled"] = enabled
    merged["screening_mode"] = mode
    return merged


def is_enabled(ctx) -> bool:
    return bool(settings_for(ctx).get("enabled", False))


def orientation_rules(settings: dict[str, Any]) -> dict[str, Any]:
    block = settings.get("orientation") or {}
    return {
        "require_verified_source": bool(block.get("require_verified_source", True)),
        "allow_proxy_for_exploration": bool(block.get("allow_proxy_for_exploration", True)),
    }


def geometry_rules(settings: dict[str, Any]) -> dict[str, Any]:
    block = settings.get("geometry") or {}
    try:
        buffer_A = float(block.get("interface_buffer_A", 5.0))
        proxy_thickness_A = float(block.get("proxy_thickness_A", 30.0))
    except (TypeError, ValueError):
        buffer_A = 5.0
        proxy_thickness_A = 30.0
    return {"interface_buffer_A": buffer_A, "proxy_thickness_A": proxy_thickness_A}


def pocket_rules(settings: dict[str, Any]) -> dict[str, Any]:
    block = settings.get("pocket") or {}
    return {
        "contact_cutoff_A": float(block.get("contact_cutoff_A", 5.0)),
        "embedded_core_fraction": float(block.get("embedded_core_fraction", 0.60)),
        "interface_fraction": float(block.get("interface_fraction", 0.50)),
        "aqueous_outside_fraction": float(block.get("aqueous_outside_fraction", 0.60)),
        "mixed_min_fraction": float(block.get("mixed_min_fraction", 0.20)),
    }


def accessibility_rules(settings: dict[str, Any]) -> dict[str, Any]:
    block = settings.get("accessibility") or {}
    return {
        "enabled": bool(block.get("enabled", True)),
        "grid_spacing_A": float(block.get("grid_spacing_A", 1.0)),
        "probe_radius_A": float(block.get("probe_radius_A", 1.4)),
        "margin_A": float(block.get("margin_A", 8.0)),
        "seed_radius_A": float(block.get("seed_radius_A", 3.0)),
        "max_grid_points": int(block.get("max_grid_points", 1_500_000)),
        "method": str(block.get("method", "grid_connectivity")),
    }


def compatibility_rules(settings: dict[str, Any]) -> dict[str, Any]:
    block = settings.get("compatibility") or {}
    return {
        "enabled": bool(block.get("enabled", True)),
        "explicit_lipid_check": str(block.get("explicit_lipid_check", "optional")),
        "lipid_coordinates": block.get("lipid_coordinates"),
        "lipid_coordinate_frame": block.get("lipid_coordinate_frame"),
        "lipid_clash_cutoff_A": float(block.get("lipid_clash_cutoff_A", 2.0)),
        "contact_cutoff_A": float(block.get("contact_cutoff_A", 4.5)),
    }


def scoring_rules(settings: dict[str, Any]) -> dict[str, Any]:
    block = settings.get("scoring") or {}
    return {
        "enabled": bool(block.get("enabled", True)),
        "strategy": str(block.get("strategy", "evidence_first")),
        "preserve_baseline": bool(block.get("preserve_baseline", True)),
        "allow_proxy_for_exploration": bool(block.get("allow_proxy_for_exploration", True)),
        "require_verified_orientation_for_formal": bool(
            block.get("require_verified_orientation_for_formal", True)
        ),
        "allow_pose_reselection": bool(block.get("allow_pose_reselection", False)),
    }


def ensemble_rules(settings: dict[str, Any]) -> dict[str, Any]:
    block = settings.get("ensemble") or {}
    return {
        "enabled": bool(block.get("enabled", False)),
        "receptors": list(block.get("receptors") or []),
        "require_state_evidence_for_formal": bool(block.get("require_state_evidence_for_formal", False)),
        "alignment_min_ca": int(block.get("alignment_min_ca", 12)),
        "alignment_max_rmsd_A": float(block.get("alignment_max_rmsd_A", 2.5)),
    }


def digest(settings: dict[str, Any]) -> str:
    """Stable short digest of the effective scientific configuration."""
    payload = json.dumps(settings, sort_keys=True, default=str, ensure_ascii=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]
