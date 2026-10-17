"""Configuration loading: host paths, target profiles, budgets, thresholds.

Design rule: every number that a scientific decision depends on lives in a YAML
file, carries a ``source`` string, and is echoed into the run's
``resolved_config.yaml``.  Nothing is hard-coded per target, and no MurJ residue
number appears in this module.
"""

from __future__ import annotations

import copy
import hashlib
import json
import os
from dataclasses import dataclass, field
from typing import Any

import yaml

from .constants import ConfigError
from .util import read_text

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _require(mapping: dict, key: str, where: str) -> Any:
    if key not in mapping:
        raise ConfigError(f"missing required key '{key}' in {where}")
    return mapping[key]


def _as_list(value: Any) -> list:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


def _digest(payload: Any) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode("utf-8")).hexdigest()[:16]


# --------------------------------------------------------------------------- #
# host configuration
# --------------------------------------------------------------------------- #


@dataclass
class ExecutorConfig:
    name: str
    mode: str = "disabled"
    settings: dict[str, Any] = field(default_factory=dict)

    def get(self, key: str, default: Any = None) -> Any:
        return self.settings.get(key, default)


@dataclass
class HostConfig:
    raw: dict[str, Any]
    source_path: str
    paths: dict[str, str] = field(default_factory=dict)
    executors: dict[str, ExecutorConfig] = field(default_factory=dict)
    defaults: dict[str, Any] = field(default_factory=dict)

    def resolve(self, logical: str) -> str:
        """Resolve a logical path name (``md_project``, ``gnina_project``, ...)."""
        if logical in self.paths:
            value = self.paths[logical]
        elif os.path.isabs(logical):
            return logical
        else:
            raise ConfigError(f"unknown logical path '{logical}'; known: {sorted(self.paths)}")
        value = os.path.expandvars(os.path.expanduser(value))
        if not os.path.isabs(value):
            value = os.path.join(os.path.dirname(self.source_path), value)
        return os.path.normpath(value)

    def maybe(self, logical: str) -> str | None:
        try:
            path = self.resolve(logical)
        except ConfigError:
            return None
        return path if os.path.exists(path) else None

    def executor(self, name: str) -> ExecutorConfig:
        return self.executors.get(name, ExecutorConfig(name=name, mode="disabled"))

    def digest(self) -> str:
        return _digest({"paths": self.paths, "executors": {k: v.settings for k, v in self.executors.items()}})


def load_host_config(path: str | None = None) -> HostConfig:
    path = path or os.environ.get("ZVSS_HOST_CONFIG") or os.path.join(PROJECT_ROOT, "configs", "host.yaml")
    if not os.path.exists(path):
        raise ConfigError(f"host config not found: {path}")
    raw = yaml.safe_load(read_text(path)) or {}
    paths = {k: str(v) for k, v in (raw.get("paths") or {}).items() if v is not None}
    executors: dict[str, ExecutorConfig] = {}
    for name, block in (raw.get("executors") or {}).items():
        block = block or {}
        executors[name] = ExecutorConfig(
            name=name, mode=str(block.get("mode", "disabled")), settings=block
        )
    return HostConfig(
        raw=raw,
        source_path=os.path.abspath(path),
        paths=paths,
        executors=executors,
        defaults=raw.get("defaults") or {},
    )


# --------------------------------------------------------------------------- #
# target profile
# --------------------------------------------------------------------------- #


@dataclass
class StructureRef:
    pdb_id: str
    chain: str
    role: str
    species: str = ""
    source_url: str = ""
    sha256: str = ""
    conformational_state: str = ""
    notes: str = ""

    @classmethod
    def from_dict(cls, payload: dict[str, Any], where: str) -> "StructureRef":
        pdb_id = str(_require(payload, "pdb_id", where)).upper()
        chain = str(payload.get("chain", "A")).strip()
        return cls(
            pdb_id=pdb_id,
            chain=chain,
            role=str(payload.get("role", "")),
            species=str(payload.get("species", "")),
            source_url=str(payload.get("source_url", "")),
            sha256=str(payload.get("sha256", "")),
            conformational_state=str(payload.get("conformational_state", "")),
            notes=str(payload.get("notes", "")),
        )


@dataclass
class Mutation:
    auth_resseq: int
    construct_resname: str
    native_resname: str

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "Mutation":
        return cls(
            auth_resseq=int(_require(payload, "auth_resseq", "construct.mutations")),
            construct_resname=str(_require(payload, "construct_resname", "construct.mutations")).upper(),
            native_resname=str(_require(payload, "native_resname", "construct.mutations")).upper(),
        )


@dataclass
class ConstructInfo:
    declared_length: int | None
    native_range: tuple[int, int] | None
    label_native_offset: int | None
    mutations: list[Mutation]
    fusion_nterm: tuple[int, int] | None
    ctag: dict[str, Any]
    raw: dict[str, Any]

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "ConstructInfo":
        native_range = payload.get("native_range")
        fusion = payload.get("fusion_nterm")
        return cls(
            declared_length=payload.get("declared_length"),
            native_range=(int(native_range[0]), int(native_range[1])) if native_range else None,
            label_native_offset=payload.get("label_native_offset"),
            mutations=[Mutation.from_dict(m) for m in _as_list(payload.get("mutations"))],
            fusion_nterm=(int(fusion[0]), int(fusion[1])) if fusion else None,
            ctag=payload.get("c_terminal_tag") or {},
            raw=payload,
        )


@dataclass
class RegionSide:
    name: str
    native_residues: list[int]
    note: str = ""

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "RegionSide":
        return cls(
            name=str(_require(payload, "name", "region.sides")),
            native_residues=[int(x) for x in _as_list(payload.get("native_residues"))],
            note=str(payload.get("note", "")),
        )


@dataclass
class ReferencePose:
    pdb_id: str
    chain: str
    resname: str
    auth_resseq: int
    native_resseq: int | None
    note: str = ""

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "ReferencePose":
        return cls(
            pdb_id=str(_require(payload, "pdb_id", "region.reference_poses")).upper(),
            chain=str(payload.get("chain", "B")),
            resname=str(payload.get("resname", "")).upper(),
            auth_resseq=int(_require(payload, "auth_resseq", "region.reference_poses")),
            native_resseq=int(payload["native_resseq"]) if payload.get("native_resseq") is not None else None,
            note=str(payload.get("note", "")),
        )


@dataclass
class RegionConfig:
    contact_cutoff_A: float
    sides: list[RegionSide]
    reference_poses: list[ReferencePose]
    box_padding_A: float
    box_min_size_A: float
    box_max_size_A: float
    raw: dict[str, Any]

    @classmethod
    def from_dict(cls, payload: dict[str, Any] | None) -> "RegionConfig | None":
        if not payload:
            return None
        return cls(
            contact_cutoff_A=float(payload.get("contact_cutoff_A", 4.5)),
            sides=[RegionSide.from_dict(side) for side in _as_list(payload.get("sides"))],
            reference_poses=[ReferencePose.from_dict(r) for r in _as_list(payload.get("reference_poses"))],
            box_padding_A=float(payload.get("box_padding_A", 5.0)),
            box_min_size_A=float(payload.get("box_min_size_A", 14.0)),
            box_max_size_A=float(payload.get("box_max_size_A", 30.0)),
            raw=payload,
        )


@dataclass
class MembraneConfig:
    status: str
    normal: list[float] | None
    center: list[float] | None
    thickness_A: float | None
    source: str
    method: str
    coordinate_frame: str
    uncertainty: str = ""
    raw: dict[str, Any] = field(default_factory=dict)

    @property
    def usable(self) -> bool:
        return bool(
            self.status in {"sourced", "explicit"}
            and self.normal
            and self.center
            and self.thickness_A
        )

    @classmethod
    def from_dict(cls, payload: dict[str, Any] | None) -> "MembraneConfig":
        payload = payload or {}
        return cls(
            status=str(payload.get("status", "unavailable")),
            normal=payload.get("normal"),
            center=payload.get("center"),
            thickness_A=(float(payload["thickness_A"]) if payload.get("thickness_A") is not None else None),
            source=str(payload.get("source", "")),
            method=str(payload.get("method", "")),
            coordinate_frame=str(payload.get("coordinate_frame", "")),
            uncertainty=str(payload.get("uncertainty", "")),
            raw=payload,
        )


@dataclass
class TargetProfile:
    target_id: str
    display_name: str
    organism: str
    uniprot: str
    native_length: int
    native_sequence: str
    primary: StructureRef
    comparisons: list[StructureRef]
    construct: ConstructInfo
    region: RegionConfig | None
    membrane: MembraneConfig
    md_protocol: dict[str, Any]
    analysis_selections: dict[str, Any]
    structures_dir: str
    raw: dict[str, Any]
    source_path: str = ""

    @property
    def is_generic(self) -> bool:
        return self.target_id == "generic"

    def side_names(self) -> list[str]:
        return [side.name for side in (self.region.sides if self.region else [])]

    def native_to_auth(self, native_resseq: int) -> int | None:
        offset = self.construct.label_native_offset
        if offset is None:
            return None
        return native_resseq + offset

    def digest(self) -> str:
        return _digest(
            {
                "target_id": self.target_id,
                "uniprot": self.uniprot,
                "native_length": self.native_length,
                "native_sequence": self.native_sequence,
                "primary": self.primary.__dict__,
                "comparisons": [c.__dict__ for c in self.comparisons],
                "construct": self.construct.raw,
                "region": self.region.raw if self.region else None,
                "membrane": self.membrane.raw,
                "md_protocol": self.md_protocol,
                "analysis_selections": self.analysis_selections,
            }
        )


def load_target_profile(path: str, *, structures_dir: str | None = None) -> TargetProfile:
    if not os.path.exists(path):
        raise ConfigError(f"target profile not found: {path}")
    raw = yaml.safe_load(read_text(path)) or {}
    where = os.path.basename(path)
    structures = raw.get("structures") or {}
    primary_payload = _require(structures, "primary", f"{where}:structures")
    comparisons = [
        StructureRef.from_dict(item, f"{where}:structures.comparisons")
        for item in _as_list(structures.get("comparisons"))
    ]
    region = RegionConfig.from_dict(raw.get("region"))
    profile = TargetProfile(
        target_id=str(_require(raw, "target_id", where)),
        display_name=str(raw.get("display_name", raw["target_id"])),
        organism=str(raw.get("organism", "")),
        uniprot=str(raw.get("uniprot", "")),
        native_length=int(raw.get("native_length", 0)),
        native_sequence=str(raw.get("native_sequence", "")).strip(),
        primary=StructureRef.from_dict(primary_payload, f"{where}:structures.primary"),
        comparisons=comparisons,
        construct=ConstructInfo.from_dict(raw.get("construct") or {}),
        region=region,
        membrane=MembraneConfig.from_dict(raw.get("membrane")),
        md_protocol=raw.get("md_protocol") or {},
        analysis_selections=raw.get("analysis_selections") or {},
        structures_dir=structures_dir or str(raw.get("structures_dir", "")),
        raw=raw,
        source_path=os.path.abspath(path),
    )
    validate_target_profile(profile)
    return profile


def validate_target_profile(profile: TargetProfile) -> None:
    errors: list[str] = []
    if profile.native_length and profile.native_sequence and len(profile.native_sequence) != profile.native_length:
        errors.append(
            f"native_sequence length {len(profile.native_sequence)} != native_length {profile.native_length}"
        )
    if profile.construct.native_range:
        start, end = profile.construct.native_range
        if start < 1 or end < start:
            errors.append(f"construct.native_range invalid: {profile.construct.native_range}")
        if profile.native_length and end > profile.native_length:
            errors.append("construct.native_range exceeds native_length")
    if profile.construct.label_native_offset is None and not profile.is_generic:
        errors.append("construct.label_native_offset missing; native numbering cannot be resolved")
    if profile.region:
        if not profile.region.sides:
            errors.append("region.sides empty")
        for side in profile.region.sides:
            if not side.native_residues:
                errors.append(f"region side '{side.name}' has no native_residues")
    if profile.membrane.status in {"sourced", "explicit"}:
        if not profile.membrane.normal or not profile.membrane.center or not profile.membrane.thickness_A:
            errors.append("membrane declared sourced/explicit but normal/center/thickness incomplete")
        if profile.membrane.normal:
            norm = sum(float(c) ** 2 for c in profile.membrane.normal) ** 0.5
            if abs(norm - 1.0) > 1e-3:
                errors.append(f"membrane.normal must be a unit vector (|n|={norm:.4f})")
    if errors:
        raise ConfigError("; ".join(errors))


# --------------------------------------------------------------------------- #
# budget / rules
# --------------------------------------------------------------------------- #


@dataclass
class BudgetConfig:
    budget_id: str
    description: str
    redock_max_parents: int
    deep_max_parents: int
    md_max_parents: int
    md_primary_slots: int
    md_exploratory_slots: int
    poses_per_parent_max: int
    seeds: list[int]
    cluster_rmsd_A: float
    redock_explore_fraction: float
    deep_explore_fraction: float
    fingerprint_radius: int
    fingerprint_bits: int
    scaffold_max_per_group: int
    similarity_ceiling: float
    pareto_objectives: list[dict[str, Any]]
    tier_requirements: dict[str, Any]
    membrane_alert_burden_max: float
    strain_max_kcal: float | None
    not_run: list[str]
    raw: dict[str, Any]

    def digest(self) -> str:
        return _digest(self.raw)

    @classmethod
    def from_dict(cls, payload: dict[str, Any], source_path: str = "") -> "BudgetConfig":
        where = f"budget:{source_path or payload.get('budget_id', '?')}"
        return cls(
            budget_id=str(_require(payload, "budget_id", where)),
            description=str(payload.get("description", "")),
            redock_max_parents=int(_require(payload, "redock_max_parents", where)),
            deep_max_parents=int(_require(payload, "deep_max_parents", where)),
            md_max_parents=int(_require(payload, "md_max_parents", where)),
            md_primary_slots=int(payload.get("md_primary_slots", 8)),
            md_exploratory_slots=int(payload.get("md_exploratory_slots", 4)),
            poses_per_parent_max=int(payload.get("poses_per_parent_max", 2)),
            seeds=[int(s) for s in _as_list(payload.get("seeds"))] or [17, 41, 97],
            cluster_rmsd_A=float(payload.get("cluster_rmsd_A", 2.0)),
            redock_explore_fraction=float(payload.get("redock_explore_fraction", 0.2)),
            deep_explore_fraction=float(payload.get("deep_explore_fraction", 0.2)),
            fingerprint_radius=int(payload.get("fingerprint_radius", 2)),
            fingerprint_bits=int(payload.get("fingerprint_bits", 2048)),
            scaffold_max_per_group=int(payload.get("scaffold_max_per_group", 2)),
            similarity_ceiling=float(payload.get("similarity_ceiling", 0.7)),
            pareto_objectives=_as_list(payload.get("pareto_objectives")),
            tier_requirements=payload.get("tier_requirements") or {},
            membrane_alert_burden_max=float(payload.get("membrane_alert_burden_max", 3)),
            strain_max_kcal=(
                float(payload["strain_max_kcal"]) if payload.get("strain_max_kcal") is not None else None
            ),
            not_run=[str(x) for x in _as_list(payload.get("not_run"))],
            raw=payload,
        )


def load_budget(path: str) -> BudgetConfig:
    if not os.path.exists(path):
        raise ConfigError(f"budget config not found: {path}")
    payload = yaml.safe_load(read_text(path)) or {}
    return BudgetConfig.from_dict(payload, os.path.basename(path))


@dataclass
class StageConfig:
    """Thresholds for the seven assessment families."""

    raw: dict[str, Any]
    source_path: str = ""

    def get(self, *keys: str, default: Any = None) -> Any:
        node: Any = self.raw
        for key in keys:
            if not isinstance(node, dict) or key not in node:
                return default
            node = node[key]
        return node

    def digest(self) -> str:
        return _digest(self.raw)


def load_stages(path: str | None = None) -> StageConfig:
    path = path or os.path.join(PROJECT_ROOT, "configs", "stages.yaml")
    if not os.path.exists(path):
        raise ConfigError(f"stage config not found: {path}")
    return StageConfig(yaml.safe_load(read_text(path)) or {}, os.path.abspath(path))


# --------------------------------------------------------------------------- #
# resolved config record
# --------------------------------------------------------------------------- #


@dataclass
class ResolvedConfig:
    host: HostConfig
    target: TargetProfile
    budget: BudgetConfig
    stages: StageConfig
    code_version: str
    overrides: dict[str, Any] = field(default_factory=dict)

    def digest(self) -> str:
        return _digest(
            {
                "host": self.host.digest(),
                "target": self.target.digest(),
                "budget": self.budget.digest(),
                "stages": self.stages.digest(),
                "code_version": self.code_version,
                "overrides": self.overrides,
            }
        )

    def to_yaml(self) -> str:
        payload = {
            "code_version": self.code_version,
            "digest": self.digest(),
            "host_config": self.host.source_path,
            "host_paths": {k: v for k, v in self.host.paths.items()},
            "executors": {name: exec_.settings for name, exec_ in self.host.executors.items()},
            "target": {
                "target_id": self.target.target_id,
                "source": self.target.source_path,
                "digest": self.target.digest(),
                "uniprot": self.target.uniprot,
                "primary": self.target.primary.__dict__,
                "comparisons": [c.__dict__ for c in self.target.comparisons],
                "construct": self.target.construct.raw,
                "region": self.target.region.raw if self.target.region else None,
                "membrane": self.target.membrane.raw,
                "md_protocol": self.target.md_protocol,
            },
            "budget": self.budget.raw,
            "stages": self.stages.raw,
            "overrides": self.overrides,
        }
        return yaml.safe_dump(copy.deepcopy(payload), allow_unicode=True, sort_keys=False)

    def redacted(self) -> dict[str, Any]:
        """Config safe for logs: no credentials, host paths reduced to basenames."""
        payload = json.loads(json.dumps(self.to_yaml(), default=str))
        return payload

    def get(self, key: str, default: Any = None) -> Any:
        """Dotted lookup across stage settings, overrides, budget and host settings."""
        for source in (
            getattr(self.stages, "settings", {}),
            self.overrides,
            getattr(self.budget, "raw", {}),
            getattr(self.host, "settings", {}),
        ):
            if not isinstance(source, dict):
                continue
            current: Any = source
            found = True
            for part in str(key).split("."):
                if isinstance(current, dict) and part in current:
                    current = current[part]
                else:
                    found = False
                    break
            if found and current is not None:
                return current
        return default


def build_resolved(
    *,
    host: HostConfig,
    target: TargetProfile,
    budget: BudgetConfig,
    stages: StageConfig,
    code_version: str,
    overrides: dict[str, Any] | None = None,
) -> ResolvedConfig:
    return ResolvedConfig(
        host=host, target=target, budget=budget, stages=stages,
        code_version=code_version, overrides=overrides or {},
    )
