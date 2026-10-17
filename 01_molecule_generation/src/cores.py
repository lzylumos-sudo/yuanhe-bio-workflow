"""Auditable selection of actual design cores from mapped Peptide2Mol seeds.

A *statistical Murcko scaffold* is useful for clustering but is not a synthesis
object.  A :class:`CoreSelection` therefore keeps it separate from the actual
(mapped, decorated) design core and from the uncut Peptide2Mol precursor.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


class CoreError(ValueError):
    """Raised when a proposed core is not an auditable design core."""


@dataclass(frozen=True)
class CoreSelection:
    seed_id: str
    statistical_murcko: str
    design_core: str
    precursor: str
    key_atom_maps: tuple[int, ...]
    attachment_ids: tuple[int, ...]

    @property
    def murcko_smiles(self) -> str:  # explicit compatibility spelling
        return self.statistical_murcko

    @property
    def core_smiles(self) -> str:
        return self.design_core

    @property
    def precursor_smiles(self) -> str:
        return self.precursor


def _mol(value: str | Chem.Mol, what: str = "molecule") -> Chem.Mol:
    mol = Chem.Mol(value) if isinstance(value, Chem.Mol) else Chem.MolFromSmiles(value)
    if mol is None:
        raise CoreError(f"invalid {what}")
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:
        raise CoreError(f"invalid {what}: {exc}") from exc
    return mol


def canonical_smiles(value: str | Chem.Mol, *, isomeric: bool = True) -> str:
    return Chem.MolToSmiles(_mol(value), canonical=True, isomericSmiles=isomeric)


def statistical_murcko(seed: str | Chem.Mol) -> str:
    """Return the unmapped Bemis--Murcko scaffold used only for statistics."""
    mol = _mol(seed, "seed")
    clean = Chem.Mol(mol)
    for atom in clean.GetAtoms():
        atom.SetAtomMapNum(0)
    scaffold = MurckoScaffold.GetScaffoldForMol(clean)
    if scaffold.GetNumAtoms() == 0:
        raise CoreError("seed has no Murcko scaffold")
    return Chem.MolToSmiles(scaffold, canonical=True, isomericSmiles=True)


def _mapped_atoms(mol: Chem.Mol) -> dict[int, int]:
    result: dict[int, int] = {}
    for atom in mol.GetAtoms():
        number = atom.GetAtomMapNum()
        if number:
            if number in result:
                raise CoreError(f"duplicate atom map {number}")
            result[number] = atom.GetIdx()
    return result


def _core_atom_indices(mol: Chem.Mol, key_maps: Sequence[int]) -> set[int]:
    maps = _mapped_atoms(mol)
    missing = sorted(set(key_maps) - maps.keys())
    if missing:
        raise CoreError(f"seed is missing mapped key atoms {missing}")
    murcko = MurckoScaffold.GetScaffoldForMol(mol)
    match = mol.GetSubstructMatch(murcko)
    keep = set(match)
    if not keep:
        raise CoreError("seed has no ring/linker Murcko scaffold")
    # Preserve every key atom and the shortest chemically real path linking it
    # to the scaffold.  This is intentionally not Murcko extraction.
    for map_no in key_maps:
        start = maps[map_no]
        if start in keep:
            continue
        paths = (Chem.GetShortestPath(mol, start, target) for target in keep)
        path = min(paths, key=len)
        keep.update(path)
    return keep


def _boundary_bonds(mol: Chem.Mol, keep: set[int]) -> list[Chem.Bond]:
    bonds = [b for b in mol.GetBonds() if (b.GetBeginAtomIdx() in keep) != (b.GetEndAtomIdx() in keep)]
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=True, includeChirality=True))
    return sorted(
        bonds,
        key=lambda b: (
            ranks[b.GetBeginAtomIdx() if b.GetBeginAtomIdx() in keep else b.GetEndAtomIdx()],
            ranks[b.GetBeginAtomIdx() if b.GetBeginAtomIdx() not in keep else b.GetEndAtomIdx()],
        ),
    )


def _fragment_core(mol: Chem.Mol, keep: set[int], boundaries: Sequence[Chem.Bond]) -> Chem.Mol:
    rw = Chem.RWMol(mol)
    # Remove atoms outside the selected core, then add stable mapped dummies.
    for idx in sorted(set(range(mol.GetNumAtoms())) - keep, reverse=True):
        rw.RemoveAtom(idx)
    core = rw.GetMol()
    old_to_new = {old: new for new, old in enumerate(sorted(keep))}
    rw = Chem.RWMol(core)
    for label, bond in enumerate(boundaries, 1):
        inside = bond.GetBeginAtomIdx() if bond.GetBeginAtomIdx() in keep else bond.GetEndAtomIdx()
        dummy = Chem.Atom(0)
        dummy.SetAtomMapNum(label)
        dummy.SetIsotope(label)  # survives software which strips atom maps
        d_idx = rw.AddAtom(dummy)
        rw.AddBond(old_to_new[inside], d_idx, bond.GetBondType())
    result = rw.GetMol()
    Chem.SanitizeMol(result)
    return result


def validate_attachment_labels(core: str | Chem.Mol, *, expected_count: int | None = None) -> tuple[int, ...]:
    """Validate 1--2 unique, continuous ``[*:n]`` labels and return them."""
    mol = _mol(core, "design core")
    dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    labels = [a.GetAtomMapNum() for a in dummies]
    if any(label <= 0 for label in labels):
        raise CoreError("every dummy atom must have a positive atom-map label")
    if len(labels) != len(set(labels)):
        raise CoreError("dummy atom-map labels must be unique")
    if sorted(labels) != list(range(1, len(labels) + 1)):
        raise CoreError("dummy atom-map labels must be continuous from 1")
    if not 1 <= len(labels) <= 2:
        raise CoreError("a design core must have one or two attachment points")
    if expected_count is not None and len(labels) != expected_count:
        raise CoreError(f"expected {expected_count} attachment points, found {len(labels)}")
    for atom in dummies:
        if atom.GetDegree() != 1:
            raise CoreError(f"attachment [*:{atom.GetAtomMapNum()}] must have degree one")
    return tuple(sorted(labels))


def derive_design_core(
    seed: str | Chem.Mol,
    key_atom_maps: Iterable[int],
    *,
    seed_id: str = "seed",
    max_attachment_points: int = 2,
) -> CoreSelection:
    """Derive a design core while retaining all mapped pharmacophore atoms.

    The seed must be a valid mapped Peptide2Mol molecule.  Natural bonds from
    the retained ring/linker/key-atom region to discarded substituents become
    attachment points.  Seeds with zero or more than two such cuts are rejected
    rather than silently inventing synthetic handles.
    """
    mol = _mol(seed, "Peptide2Mol seed")
    keys = tuple(sorted(set(int(x) for x in key_atom_maps)))
    if not keys or any(x <= 0 for x in keys):
        raise CoreError("at least one positive mapped key atom is required")
    keep = _core_atom_indices(mol, keys)
    boundaries = _boundary_bonds(mol, keep)
    if not 1 <= len(boundaries) <= min(2, max_attachment_points):
        raise CoreError(f"seed yields {len(boundaries)} attachment points; expected 1-{min(2, max_attachment_points)}")
    core = _fragment_core(mol, keep, boundaries)
    labels = validate_attachment_labels(core)
    core_maps = {a.GetAtomMapNum() for a in core.GetAtoms() if a.GetAtomicNum() != 0 and a.GetAtomMapNum()}
    missing = sorted(set(keys) - core_maps)
    if missing:
        raise CoreError(f"design core lost key atom maps {missing}")
    return CoreSelection(
        seed_id=seed_id,
        statistical_murcko=statistical_murcko(mol),
        design_core=Chem.MolToSmiles(core, canonical=True, isomericSmiles=True),
        precursor=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        key_atom_maps=keys,
        attachment_ids=labels,
    )


def select_core(
    seeds: Iterable[str | Chem.Mol],
    key_atom_maps: Iterable[int],
    *,
    seed_ids: Iterable[str] | None = None,
) -> CoreSelection:
    """Select the deterministic least-cut valid core from candidate seeds."""
    seed_list = list(seeds)
    ids = list(seed_ids) if seed_ids is not None else [f"seed-{i + 1}" for i in range(len(seed_list))]
    if len(ids) != len(seed_list):
        raise CoreError("seed_ids and seeds must have equal length")
    valid: list[CoreSelection] = []
    errors: list[str] = []
    for sid, seed in zip(ids, seed_list):
        try:
            valid.append(derive_design_core(seed, key_atom_maps, seed_id=sid))
        except CoreError as exc:
            errors.append(f"{sid}: {exc}")
    if not valid:
        raise CoreError("no valid Peptide2Mol seed: " + "; ".join(errors))
    return min(valid, key=lambda x: (len(x.attachment_ids), x.design_core, x.seed_id))


def _substituent_anchor(mol: Chem.Mol, label: int) -> tuple[int, int]:
    matches = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetAtomMapNum() == label]
    if len(matches) != 1 or matches[0].GetDegree() != 1:
        raise CoreError(f"substituent must contain exactly one terminal [*:{label}]")
    return matches[0].GetIdx(), matches[0].GetNeighbors()[0].GetIdx()


def attach_substituents(
    core: str | Chem.Mol,
    substituents: Mapping[int, str | Chem.Mol] | Sequence[str | Chem.Mol],
) -> Chem.Mol:
    """Consume exactly the expected mapped stars and return a sanitized product."""
    base = _mol(core, "design core")
    labels = validate_attachment_labels(base)
    supplied = dict(enumerate(substituents, 1)) if not isinstance(substituents, Mapping) else {int(k): v for k, v in substituents.items()}
    if set(supplied) != set(labels):
        raise CoreError(f"substituent labels {sorted(supplied)} do not match core labels {list(labels)}")
    product = Chem.Mol(base)
    # One label at a time makes atom-index renumbering harmless.
    for label in labels:
        sub = _mol(supplied[label], f"substituent {label}")
        c_matches = [a for a in product.GetAtoms() if a.GetAtomicNum() == 0 and a.GetAtomMapNum() == label]
        if len(c_matches) != 1 or c_matches[0].GetDegree() != 1:
            raise CoreError(f"core attachment [*:{label}] is missing or ambiguous")
        c_dummy = c_matches[0]
        c_neighbor = c_dummy.GetNeighbors()[0]
        c_bond = product.GetBondBetweenAtoms(c_dummy.GetIdx(), c_neighbor.GetIdx()).GetBondType()
        s_dummy_idx, s_neighbor_idx = _substituent_anchor(sub, label)
        s_bond = sub.GetBondBetweenAtoms(s_dummy_idx, s_neighbor_idx).GetBondType()
        if c_bond != s_bond and Chem.BondType.SINGLE not in (c_bond, s_bond):
            raise CoreError(f"incompatible bond types at attachment {label}")
        combo = Chem.CombineMols(product, sub)
        offset = product.GetNumAtoms()
        rw = Chem.RWMol(combo)
        rw.AddBond(c_neighbor.GetIdx(), offset + s_neighbor_idx, c_bond)
        for idx in sorted((c_dummy.GetIdx(), offset + s_dummy_idx), reverse=True):
            rw.RemoveAtom(idx)
        product = rw.GetMol()
        try:
            Chem.SanitizeMol(product)
        except Exception as exc:
            raise CoreError(f"attachment {label} gives an invalid product: {exc}") from exc
    if any(a.GetAtomicNum() == 0 for a in product.GetAtoms()):
        raise CoreError("known substitution did not consume exactly the expected stars")
    return product


def validate_product(core: str | Chem.Mol, product: str | Chem.Mol, key_atom_maps: Iterable[int] = ()) -> bool:
    """Check sanitization plus mapped core/key-atom retention in a product."""
    core_mol = _mol(core, "design core")
    prod = _mol(product, "product")
    query = Chem.DeleteSubstructs(core_mol, Chem.MolFromSmarts("[#0]"))
    Chem.SanitizeMol(query)
    if not prod.HasSubstructMatch(query, useChirality=True):
        raise CoreError("product does not retain the design-core substructure/stereo")
    maps = {a.GetAtomMapNum() for a in prod.GetAtoms() if a.GetAtomicNum() != 0}
    missing = sorted(set(int(x) for x in key_atom_maps) - maps)
    if missing:
        raise CoreError(f"product lost mapped key atoms {missing}")
    return True


# Clear aliases for callers preferring noun phrases.
get_statistical_murcko = statistical_murcko
select_design_core = select_core
validate_dummy_labels = validate_attachment_labels
