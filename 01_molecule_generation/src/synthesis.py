"""Minimal, versioned and auditable synthesis support built on RDKit."""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from datetime import date, datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from .cores import CoreError, attach_substituents, canonical_smiles, validate_product


class ContractError(ValueError):
    pass


class RouteStatus(str, Enum):
    SUPPORTED_FORWARD = "supported_forward"
    SUPPORTED_RETRO = "supported_retro"
    UNRESOLVED = "unresolved"
    TIMEOUT = "timeout"
    ERROR = "error"


@dataclass(frozen=True)
class BuildingBlock:
    block_id: str
    smiles: str
    version: str
    in_stock: bool
    checked_at: date
    vendor: str = ""

    @property
    def stale(self) -> bool:
        return False  # staleness is relative to an audit date; use stock_is_stale


@dataclass(frozen=True)
class ReactionRule:
    reaction_id: str
    smarts: str
    version: str
    direction: str = "forward"


@dataclass(frozen=True)
class IdentityCheck:
    matches: bool
    reason: str
    expected: str
    observed: str


@dataclass(frozen=True)
class RouteStep:
    reaction_id: str
    status: RouteStatus
    product_smiles: str | None = None
    block_ids: tuple[str, ...] = ()
    identity: IdentityCheck | None = None
    note: str = ""


@dataclass(frozen=True)
class RouteAudit:
    steps: tuple[RouteStep, ...]
    catalog_version: str
    reaction_version: str
    needs_review: bool = False
    review_reasons: tuple[str, ...] = ()

    @property
    def complete(self) -> bool:
        supported = {RouteStatus.SUPPORTED_FORWARD, RouteStatus.SUPPORTED_RETRO}
        return bool(self.steps) and all(step.status in supported for step in self.steps)

    @property
    def formal(self) -> bool:
        # A route with stale/unknown stock is evidence, but not a formal route.
        return self.complete and not self.needs_review

    @property
    def status(self) -> RouteStatus:
        if not self.steps:
            return RouteStatus.UNRESOLVED
        for status in (RouteStatus.ERROR, RouteStatus.TIMEOUT, RouteStatus.UNRESOLVED):
            if any(step.status == status for step in self.steps):
                return status
        if all(step.status == RouteStatus.SUPPORTED_FORWARD for step in self.steps):
            return RouteStatus.SUPPORTED_FORWARD
        return RouteStatus.SUPPORTED_RETRO


_TRUE = {"1", "true", "yes", "y", "in_stock", "available"}
_FALSE = {"0", "false", "no", "n", "out_of_stock", "unavailable"}


def _read_csv(source: str | Path) -> tuple[list[dict[str, str]], set[str]]:
    path = Path(source)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ContractError(f"{path.name}: missing CSV header")
        rows = [{str(k).strip(): (v or "").strip() for k, v in row.items()} for row in reader]
        return rows, {x.strip() for x in reader.fieldnames}


def _column(row: Mapping[str, str], *names: str) -> str:
    for name in names:
        if name in row and row[name] != "":
            return row[name]
    return ""


def _version(rows: Sequence[Mapping[str, str]], kind: str, expected: str | None) -> str:
    versions = {_column(row, "contract_version", "version", f"{kind}_version") for row in rows}
    versions.discard("")
    if len(versions) != 1:
        raise ContractError(f"{kind} CSV must contain exactly one non-empty contract version")
    value = next(iter(versions))
    if expected is not None and value != expected:
        raise ContractError(f"{kind} contract version {value!r} does not match {expected!r}")
    return value


def _parse_date(text: str) -> date:
    if not text:
        raise ContractError("stock checked_at/snapshot_date is required")
    try:
        return date.fromisoformat(text[:10])
    except ValueError as exc:
        raise ContractError(f"invalid ISO stock date {text!r}") from exc


def load_building_blocks_csv(source: str | Path, *, expected_version: str | None = None) -> tuple[BuildingBlock, ...]:
    rows, _ = _read_csv(source)
    version = _version(rows, "building_blocks", expected_version)
    blocks: list[BuildingBlock] = []
    seen: set[str] = set()
    for line, row in enumerate(rows, 2):
        block_id = _column(row, "block_id", "id", "building_block_id")
        smiles = _column(row, "smiles", "mapped_smiles")
        stock = _column(row, "in_stock", "stock_status", "available").lower()
        if not block_id or block_id in seen:
            raise ContractError(f"building_blocks line {line}: missing or duplicate block_id")
        if Chem.MolFromSmiles(smiles) is None:
            raise ContractError(f"building_blocks line {line}: invalid SMILES")
        if stock not in _TRUE | _FALSE:
            raise ContractError(f"building_blocks line {line}: invalid stock status {stock!r}")
        seen.add(block_id)
        blocks.append(BuildingBlock(block_id, smiles, version, stock in _TRUE, _parse_date(_column(row, "checked_at", "snapshot_date", "stock_checked_at")), _column(row, "vendor")))
    if not blocks:
        raise ContractError("building_blocks CSV is empty")
    return tuple(blocks)


def load_reactions_csv(source: str | Path, *, expected_version: str | None = None) -> tuple[ReactionRule, ...]:
    rows, _ = _read_csv(source)
    version = _version(rows, "reactions", expected_version)
    rules: list[ReactionRule] = []
    seen: set[str] = set()
    for line, row in enumerate(rows, 2):
        rid = _column(row, "reaction_id", "id", "name")
        smarts = _column(row, "reaction_smarts", "smarts", "mapped_reaction_smarts")
        direction = (_column(row, "direction") or "forward").lower()
        if not rid or rid in seen:
            raise ContractError(f"reactions line {line}: missing or duplicate reaction_id")
        try:
            reaction = rdChemReactions.ReactionFromSmarts(smarts)
        except Exception:
            reaction = None
        if reaction is None or ">>" not in smarts:
            raise ContractError(f"reactions line {line}: invalid reaction SMARTS")
        # Auditing requires maps on both sides, not merely syntactically valid SMARTS.
        left, right = smarts.split(">>", 1)
        if ":" not in left or ":" not in right:
            raise ContractError(f"reactions line {line}: reaction SMARTS must be atom-mapped")
        if direction not in {"forward", "retro", "both"}:
            raise ContractError(f"reactions line {line}: invalid direction")
        seen.add(rid)
        rules.append(ReactionRule(rid, smarts, version, direction))
    if not rules:
        raise ContractError("reactions CSV is empty")
    return tuple(rules)


def stock_is_stale(block: BuildingBlock, *, as_of: date | None = None, max_age_days: int = 30) -> bool:
    as_of = as_of or datetime.now(timezone.utc).date()
    return (as_of - block.checked_at).days > max_age_days or (as_of - block.checked_at).days < 0


def _identity_smiles(value: str | Chem.Mol, *, isomeric: bool) -> str:
    """Canonical chemical identity, ignoring reaction bookkeeping maps."""
    mol = Chem.Mol(value) if isinstance(value, Chem.Mol) else Chem.MolFromSmiles(value)
    if mol is None:
        raise CoreError("invalid molecule")
    Chem.SanitizeMol(mol)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=isomeric)


def check_identity(expected: str | Chem.Mol, observed: str | Chem.Mol) -> IdentityCheck:
    """Compare constitution and stereochemistry, reporting stereo separately."""
    try:
        exp_iso = _identity_smiles(expected, isomeric=True)
        obs_iso = _identity_smiles(observed, isomeric=True)
        exp_flat = _identity_smiles(expected, isomeric=False)
        obs_flat = _identity_smiles(observed, isomeric=False)
    except CoreError as exc:
        return IdentityCheck(False, f"invalid molecule: {exc}", "", "")
    if exp_flat != obs_flat:
        return IdentityCheck(False, "constitutional identity mismatch", exp_iso, obs_iso)
    if exp_iso != obs_iso:
        return IdentityCheck(False, "stereochemistry mismatch", exp_iso, obs_iso)
    return IdentityCheck(True, "identity and stereochemistry match", exp_iso, obs_iso)


def run_mapped_reaction(rule: ReactionRule | str, reactants: Iterable[str | Chem.Mol]) -> tuple[Chem.Mol, ...]:
    smarts = rule.smarts if isinstance(rule, ReactionRule) else rule
    reaction = rdChemReactions.ReactionFromSmarts(smarts)
    if reaction is None:
        raise ContractError("invalid mapped reaction SMARTS")
    mols = tuple(Chem.Mol(x) if isinstance(x, Chem.Mol) else Chem.MolFromSmiles(x) for x in reactants)
    if any(mol is None for mol in mols):
        raise ContractError("invalid reaction reactant")
    if len(mols) != reaction.GetNumReactantTemplates():
        raise ContractError(f"reaction expects {reaction.GetNumReactantTemplates()} reactants, got {len(mols)}")
    unique: dict[str, Chem.Mol] = {}
    for outcome in reaction.RunReactants(mols):
        for product in outcome:
            try:
                Chem.SanitizeMol(product)
                smiles = Chem.MolToSmiles(product, canonical=True, isomericSmiles=True)
                unique[smiles] = product
            except Exception:
                continue
    return tuple(unique[key] for key in sorted(unique))


def forward_reconstruct(
    proposal: str | Chem.Mol,
    *,
    core: str | Chem.Mol | None = None,
    substituents: Mapping[int, str | Chem.Mol] | Sequence[str | Chem.Mol] | None = None,
    reaction: ReactionRule | str | None = None,
    reactants: Iterable[str | Chem.Mol] | None = None,
    key_atom_maps: Iterable[int] = (),
) -> IdentityCheck:
    """Forward-reconstruct a proposal, then require identity and stereo equality."""
    if reaction is not None:
        if reactants is None:
            raise ContractError("reactants are required with reaction SMARTS")
        products = run_mapped_reaction(reaction, reactants)
    else:
        if core is None or substituents is None:
            raise ContractError("core and substituents are required for star attachment")
        product = attach_substituents(core, substituents)
        validate_product(core, product, key_atom_maps)
        products = (product,)
    checks = [check_identity(proposal, product) for product in products]
    for check in checks:
        if check.matches:
            return check
    if checks:
        stereo = next((x for x in checks if x.reason == "stereochemistry mismatch"), None)
        return stereo or checks[0]
    return IdentityCheck(False, "reaction produced no sanitized product", canonical_smiles(proposal), "")


def audit_route(
    steps: Iterable[RouteStep],
    blocks: Iterable[BuildingBlock],
    reactions: Iterable[ReactionRule],
    *,
    as_of: date | None = None,
    max_stock_age_days: int = 30,
) -> RouteAudit:
    steps = tuple(steps)
    block_by_id = {block.block_id: block for block in blocks}
    reaction_by_id = {rule.reaction_id: rule for rule in reactions}
    b_versions = {block.version for block in block_by_id.values()}
    r_versions = {rule.version for rule in reaction_by_id.values()}
    if len(b_versions) != 1 or len(r_versions) != 1:
        raise ContractError("route catalogs must each use exactly one version")
    reasons: list[str] = []
    for i, step in enumerate(steps, 1):
        if step.reaction_id not in reaction_by_id:
            reasons.append(f"step {i}: unknown reaction {step.reaction_id}")
        if step.status == RouteStatus.SUPPORTED_FORWARD and (step.identity is None or not step.identity.matches):
            reasons.append(f"step {i}: forward support lacks an identity/stereo match")
        for bid in step.block_ids:
            block = block_by_id.get(bid)
            if block is None:
                reasons.append(f"step {i}: unknown building block {bid}")
            elif not block.in_stock:
                reasons.append(f"step {i}: building block {bid} is out of stock")
            elif stock_is_stale(block, as_of=as_of, max_age_days=max_stock_age_days):
                reasons.append(f"step {i}: building block {bid} stock snapshot is stale")
    return RouteAudit(steps, next(iter(b_versions), ""), next(iter(r_versions), ""), bool(reasons), tuple(reasons))


# Short aliases used by pipeline callers.
load_building_blocks = load_building_blocks_csv
load_reactions = load_reactions_csv
validate_forward_identity = check_identity
