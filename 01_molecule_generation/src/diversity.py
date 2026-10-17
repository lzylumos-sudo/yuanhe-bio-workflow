"""Graph-based scaffold and Morgan diversity with explicit acyclic topology.

Acyclic molecules never collapse into the empty Murcko scaffold: their labeled
heavy-atom graph and generic carbon topology are retained as separate keys.
Clustering uses deterministic Butina-style greedy sphere exclusion, with an
explicit Tanimoto similarity threshold (not a distance threshold).
"""
from __future__ import annotations

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

from .geometry import deterministic_id


def _mol(value):
    mol = Chem.MolFromSmiles(value) if isinstance(value, str) else Chem.Mol(value)
    if mol is None:
        raise ValueError("invalid molecular graph for diversity")
    return Chem.RemoveHs(mol)


def scaffold_keys(smiles_or_mol):
    mol = _mol(smiles_or_mol)
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold.GetNumAtoms():
        murcko = Chem.MolToSmiles(scaffold, isomericSmiles=False)
        generic = Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(scaffold), isomericSmiles=False)
        return dict(murcko="ring:" + murcko, generic_murcko="ring:" + generic,
                    acyclic_topology=None, is_acyclic=False)
    # Full tree/forest connectivity is the only honest acyclic analogue of a
    # framework: neither a shared empty string nor the token 'acyclic' suffices.
    topology = Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(mol), isomericSmiles=False)
    return dict(murcko="acyclic:" + Chem.MolToSmiles(mol, isomericSmiles=False),
                generic_murcko="acyclic:" + topology, acyclic_topology=topology, is_acyclic=True)


def morgan_clusters(smiles_or_mols, similarity_threshold=0.65, radius=2, n_bits=2048):
    """Return one deterministic cluster ID per input, independent of input order.

    Duplicate chemical graphs are clustered once, so multiple poses do not
    influence neighbor counts. A sphere member is similar to its center, not
    necessarily to every other member (this is not complete-link clustering).
    """
    if not 0 <= similarity_threshold <= 1:
        raise ValueError("similarity_threshold must be in [0, 1]")
    mols = [_mol(value) for value in smiles_or_mols]
    keys = [Chem.MolToSmiles(m, isomericSmiles=True) for m in mols]
    unique = sorted(set(keys))
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=int(radius), fpSize=int(n_bits))
    fps = {key: generator.GetFingerprint(Chem.MolFromSmiles(key)) for key in unique}
    neighbors = {key: {other for other in unique if DataStructs.TanimotoSimilarity(fps[key], fps[other]) >= similarity_threshold}
                 for key in unique}
    remaining, assignment = set(unique), {}
    while remaining:
        center = min(remaining, key=lambda key: (-len(neighbors[key] & remaining), key))
        members = neighbors[center] & remaining
        cid = deterministic_id("morgan", radius, n_bits, similarity_threshold, center)
        for member in members:
            assignment[member] = cid
        remaining -= members
    return [assignment[key] for key in keys]


def assign_diversity_groups(records, similarity_threshold=0.65, radius=2, n_bits=2048):
    """Copy records, preserving all compound/state/conformer/pose identities."""
    records = list(records)
    values = [r["mol"] if r.get("mol") is not None else r.get("smiles", r.get("canonical_smiles")) for r in records]
    clusters = morgan_clusters(values, similarity_threshold, radius, n_bits)
    result = []
    for record, value, cluster in zip(records, values, clusters):
        item = dict(record)
        keys = scaffold_keys(value)
        item.update(keys)
        item.update(murcko_group=deterministic_id("murcko", keys["murcko"]),
                    generic_murcko_group=deterministic_id("generic_murcko", keys["generic_murcko"]),
                    morgan_cluster=cluster)
        result.append(item)
    return result
