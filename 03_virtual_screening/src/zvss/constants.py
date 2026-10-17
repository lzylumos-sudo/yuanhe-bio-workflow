"""Typed constants and exceptions shared by the whole pipeline.

The three status families are deliberately separate.  A tool exiting 0 is an
*execution* fact; whether a scientific check passed is an *assessment* fact;
whether a molecule ended up in a list is a *selection* fact.  Mixing them is how
screening pipelines end up reporting un-run checks as PASS.
"""

from __future__ import annotations

from enum import StrEnum


class ExecutionStatus(StrEnum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CACHED = "cached"
    SKIPPED = "skipped"


class AssessmentStatus(StrEnum):
    PASS = "pass"
    FLAG = "flag"
    FAIL = "fail"
    NEEDS_INPUT = "needs_input"
    NOT_EVALUATED = "not_evaluated"
    UNSUPPORTED = "unsupported"


class SelectionStatus(StrEnum):
    ELIGIBLE = "eligible"
    EXCLUDED = "excluded"
    EXPLORATORY = "exploratory"
    SELECTED = "selected"
    NOT_SELECTED = "not_selected"


EXECUTION_ORDER = [
    ExecutionStatus.PENDING,
    ExecutionStatus.RUNNING,
    ExecutionStatus.COMPLETED,
    ExecutionStatus.FAILED,
    ExecutionStatus.CACHED,
    ExecutionStatus.SKIPPED,
]

# Assessment statuses that mean "a real measurement happened and it was fine".
PASSING = {AssessmentStatus.PASS}
# Statuses that represent an absent measurement rather than a good or bad one.
MISSING = {
    AssessmentStatus.NOT_EVALUATED,
    AssessmentStatus.NEEDS_INPUT,
    AssessmentStatus.UNSUPPORTED,
}
# Statuses that block a molecule from the primary candidate set.
BLOCKING = {AssessmentStatus.FAIL}


class ZvssError(Exception):
    """Base class for recoverable pipeline errors."""

    code = "ZVSS_ERROR"

    def __init__(self, message: str, *, code: str | None = None, detail: dict | None = None):
        super().__init__(message)
        self.message = message
        if code:
            self.code = code
        self.detail = detail or {}


class ConfigError(ZvssError):
    code = "CONFIG_ERROR"


class InputError(ZvssError):
    code = "INPUT_ERROR"


class ChemistryError(ZvssError):
    code = "CHEMISTRY_ERROR"


class ReceptorError(ZvssError):
    code = "RECEPTOR_ERROR"


class ExternalToolError(ZvssError):
    code = "EXTERNAL_TOOL_ERROR"


class CapabilityError(ZvssError):
    code = "CAPABILITY_MISSING"


class IdentityConflict(ZvssError):
    code = "IDENTITY_CONFLICT"


# Reason codes used across the pipeline.  Keeping them in one place makes the
# exclusion report grep-able and stable across versions.
class Reason(StrEnum):
    OK = "ok"
    CHEM_TEMPLATE_MISSING = "chem_template_missing"
    CHEM_ATOM_COUNT_MISMATCH = "chem_atom_count_mismatch"
    CHEM_ELEMENT_MISMATCH = "chem_element_mismatch"
    CHEM_BOND_ORDER_UNRECOVERED = "chem_bond_order_unrecovered"
    CHEM_IDENTITY_CONFLICT = "chem_identity_conflict"
    CHEM_STEREO_CONFLICT = "chem_stereo_conflict"
    POSE_CLASH_PROTEIN = "pose_clash_protein"
    POSE_CLASH_INTERNAL = "pose_clash_internal"
    POSE_BOND_GEOMETRY = "pose_bond_geometry"
    POSE_RING_PLANARITY = "pose_ring_planarity"
    POSE_OUT_OF_REGION = "pose_out_of_region"
    POSE_NOT_IN_BOX = "pose_out_of_box"
    SCORE_MISSING = "score_missing"
    SCORE_GROUP_MISMATCH = "score_group_mismatch"
    SCORE_MAPPING_UNVERIFIED = "score_mapping_unverified"
    LIGAND_AMBIGUOUS = "ligand_ambiguity_unresolved"
    TEMPLATE_NOT_FOUND = "chemical_template_not_found"
    RECEPTOR_UNBOUND = "receptor_not_bound_to_pose"
    RECEPTOR_VERSION_UNKNOWN = "receptor_version_unknown"
    RECEPTOR_FLEXIBLE_GEOMETRY = "receptor_flexible_geometry_unknown"
    REGION_NOT_CONFIGURED = "region_not_configured"
    INTERACTION_ENGINE_FAILED = "interaction_engine_failed"
    MOTION_REFERENCE_MISSING = "motion_reference_missing"
    MOTION_ALIGNMENT_FAILED = "motion_alignment_failed"
    REDOCK_TOOL_UNAVAILABLE = "redock_tool_unavailable"
    REDOCK_RUN_FAILED = "redock_run_failed"
    REDOCK_SEED_INCOMPLETE = "redock_seed_incomplete"
    MEMBRANE_ORIENTATION_MISSING = "membrane_orientation_missing"
    MEMBRANE_ORIENTATION_UNSOURCED = "membrane_orientation_unsourced"
    FORCEFIELD_UNSUPPORTED = "forcefield_unsupported"
    STRAIN_NOT_CONVERGED = "strain_not_converged"
    EVIDENCE_INCOMPLETE = "evidence_incomplete"
    BUDGET_EXHAUSTED = "budget_exhausted"
    POSE_GEOMETRY_INVALID = "pose_geometry_invalid"
    SEVERE_CLASH = "severe_clash"
    MISSING_COORDINATES = "missing_coordinates"
    NO_CONTACTS = "no_contacts"
    RESIDUE_MAPPING_MISSING = "residue_mapping_missing"
    REFERENCE_POSE_MISSING = "reference_pose_missing"
    MD_MEMBRANE_CAPABILITY_MISSING = "blocked_missing_membrane_capability"
    MD_PLATFORM_LINUX_ONLY = "blocked_platform_linux_required"
    EXCESSIVE_STRAIN = "excessive_strain"
    PROPERTY_VIOLATION = "property_violation"
    SINGLE_POSE = "single_pose_only"
    SKELETON_DUPLICATE = "skeleton_duplicate"
    TIER_BELOW_CUTOFF = "tier_below_cutoff"
    PARETO_DOMINATED = "pareto_dominated"
    FAILED_POSE_EXCLUDED = "failed_pose_excluded"
    NEEDS_INPUT = "needs_input"
