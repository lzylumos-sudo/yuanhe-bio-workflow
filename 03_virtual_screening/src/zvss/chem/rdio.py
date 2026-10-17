"""RDKit byte I/O that survives non-ASCII paths.

RDKit's C++ readers take narrow ``char*`` filenames.  On a host whose working
directory contains CJK characters that yields ``OSError: Bad input file`` even
though the file exists.  All molecule I/O in ZVSS therefore goes through this
module: Python does the bytes, RDKit parses/serialises blocks.
"""

from __future__ import annotations

import gzip
import io
import os
from typing import Any, Iterator

from rdkit import Chem
from rdkit.Chem import AllChem

from .. import compat

compat.bootstrap()

COMPRESSED_SUFFIXES = (".gz",)


def _read_bytes(path: str) -> bytes:
    if path.lower().endswith(COMPRESSED_SUFFIXES):
        with gzip.open(path, "rb") as handle:
            return handle.read()
    with open(path, "rb") as handle:
        return handle.read()


def read_text_maybe_gzip(path: str, encoding: str = "utf-8") -> str:
    return _read_bytes(path).decode(encoding, errors="replace")


def iter_sdf(path: str, *, sanitize: bool = True, strict: bool = False) -> Iterator[tuple[Chem.Mol, dict[str, str]]]:
    """Yield ``(mol, properties)`` for each record, preserving record order.

    ``record_index`` is the 0-based position inside this file.  It is part of
    the pose identity, so it must never be recomputed from a sorted name list.
    """
    data = _read_bytes(path)
    stream = io.BytesIO(data)
    supplier = Chem.ForwardSDMolSupplier(stream, sanitize=sanitize, removeHs=False, strictParsing=strict)
    index = 0
    for mol in supplier:
        if mol is None:
            yield None, {"_record_index": str(index), "_parse_error": "true"}
            index += 1
            continue
        props = {key: mol.GetProp(key) for key in mol.GetPropNames(includePrivate=False, includeComputed=False)}
        props["_record_index"] = str(index)
        index += 1
        yield mol, props


def mol_from_mol_block(block: str, *, sanitize: bool = True, remove_hs: bool = False) -> Chem.Mol | None:
    return Chem.MolFromMolBlock(block, sanitize=sanitize, removeHs=remove_hs, strictParsing=False)


def mol_from_pdb_block(block: str, *, sanitize: bool = True, remove_hs: bool = False) -> Chem.Mol | None:
    """Parse a PDB block into a molecule.

    Only use this for *ligand-sized* blocks: the PDB reader guesses bonds from
    distances, which is fine for a small molecule with complete heavy atoms and
    hopeless for a protein.
    """
    return Chem.MolFromPDBBlock(block, sanitize=sanitize, removeHs=remove_hs)


def mol_from_smiles(smiles: str) -> Chem.Mol | None:
    return Chem.MolFromSmiles(smiles)


def _add_props(mol: Chem.Mol, props: dict[str, Any]) -> None:
    for key, value in props.items():
        if value is None:
            continue
        if isinstance(value, (dict, list, tuple)):
            import json

            mol.SetProp(str(key), json.dumps(value, ensure_ascii=False, default=str))
        elif isinstance(value, float):
            mol.SetProp(str(key), f"{value:.6f}")
        else:
            mol.SetProp(str(key), str(value))


def mol_block_with_props(mol: Chem.Mol, props: dict[str, Any] | None = None, *, name: str | None = None) -> str:
    mol = Chem.Mol(mol)
    if name is not None:
        mol.SetProp("_Name", name)
    if props:
        _add_props(mol, props)
    return Chem.MolToMolBlock(mol)


def write_sdf(path: str, records: list[tuple[Chem.Mol, dict[str, Any]]], *, name_key: str = "_Name") -> str:
    """Serialise an SD file with Python doing the byte write."""
    chunks: list[str] = []
    for mol, props in records:
        combined = dict(props or {})
        name = combined.pop(name_key, None) or (mol.GetProp("_Name") if mol.HasProp("_Name") else "")
        chunks.append(mol_block_with_props(mol, combined, name=str(name)))
        for key, value in (props or {}).items():
            if key == name_key or value is None:
                continue
            if isinstance(value, (dict, list, tuple)):
                import json

                text = json.dumps(value, ensure_ascii=False, default=str)
            else:
                text = str(value)
            chunks.append(f">  <{key}>\n{text}\n")
        chunks.append("$$$$\n")
    from ..util import atomic_write_text

    return atomic_write_text(path, "".join(chunks))


def write_mol(path: str, mol: Chem.Mol, props: dict[str, Any] | None = None) -> str:
    from ..util import atomic_write_text

    return atomic_write_text(path, mol_block_with_props(mol, props))


def canonical_isomeric_smiles(mol: Chem.Mol) -> str | None:
    if mol is None:
        return None
    try:
        copy = Chem.Mol(mol)
        Chem.RemoveStereochemistry(copy)
        Chem.AssignStereochemistryFrom3D(copy)
        return Chem.MolToSmiles(copy, isomericSmiles=True, canonical=True)
    except Exception:
        try:
            return Chem.MolToSmiles(Chem.Mol(mol), isomericSmiles=True, canonical=True)
        except Exception:
            return None


def canonical_smiles_flat(mol: Chem.Mol) -> str | None:
    """Stereo-stripped canonical SMILES, used only for flat comparisons."""
    if mol is None:
        return None
    try:
        return Chem.MolToSmiles(Chem.Mol(mol), isomericSmiles=False, canonical=True)
    except Exception:
        return None


def mol_formula(mol: Chem.Mol) -> str:
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula

    return CalcMolFormula(mol)


def net_charge(mol: Chem.Mol) -> int:
    return int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms()))


def heavy_atom_count(mol: Chem.Mol) -> int:
    return int(sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1))


def coordinates(mol: Chem.Mol, conf_id: int = 0) -> list[tuple[float, float, float]]:
    conf = mol.GetConformer(conf_id)
    return [
        (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z)
        for i in range(mol.GetNumAtoms())
    ]


def set_coordinates(mol: Chem.Mol, coords: list[tuple[float, float, float]], conf_id: int = -1) -> Chem.Mol:
    mol = Chem.Mol(mol)
    if conf_id < 0 or conf_id >= mol.GetNumConformers():
        mol.RemoveAllConformers()
        conf = Chem.Conformer(mol.GetNumAtoms())
        conf.Set3D(True)
        mol.AddConformer(conf, assignId=True)
        conf_id = 0
    conf = mol.GetConformer(conf_id)
    for index, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(index, Chem.rdGeometry.Point3D(float(x), float(y), float(z)))
    return mol


def add_hydrogens(mol: Chem.Mol, *, add_coords: bool = True, ph: float | None = None) -> Chem.Mol:
    """Add explicit hydrogens without moving heavy atoms.

    ``add_coords=True`` places the new hydrogens from local geometry while every
    heavy-atom coordinate stays byte-identical to the input pose.
    """
    mol_h = Chem.Mol(mol)
    mol_h.UpdatePropertyCache(strict=False)
    if add_coords and mol_h.GetNumConformers():
        mol_h = Chem.AddHs(mol_h, addCoords=True)
    else:
        mol_h = Chem.AddHs(mol_h)
    return mol_h


def heavy_atom_mapping(mol_a: Chem.Mol, mol_b: Chem.Mol) -> list[tuple[int, int]]:
    """Best-effort element-ordered heavy-atom pairing used for sanity checks."""
    a = [i for i, atom in enumerate(mol_a.GetAtoms()) if atom.GetAtomicNum() > 1]
    b = [i for i, atom in enumerate(mol_b.GetAtoms()) if atom.GetAtomicNum() > 1]
    if len(a) != len(b):
        return []
    pairs = []
    for index, atom_index in enumerate(a):
        pairs.append((atom_index, b[index]))
    return pairs


def rmsd_fixed_frame(mol_a: Chem.Mol, mol_b: Chem.Mol, mapping: list[tuple[int, int]] | None = None) -> float | None:
    """Heavy-atom RMSD in the *receptor* frame (no ligand re-fitting).

    This is the only RMSD ZVSS is allowed to use for "same binding site or not".
    Fitting the ligand onto itself (GetBestRMS) would erase exactly the signal we
    want to measure.
    """
    import numpy as np

    if mol_a.GetNumConformers() == 0 or mol_b.GetNumConformers() == 0:
        return None
    pairs = mapping or heavy_atom_mapping(mol_a, mol_b)
    if not pairs:
        return None
    conf_a = mol_a.GetConformer()
    conf_b = mol_b.GetConformer()
    deltas = []
    for i, j in pairs:
        pa = conf_a.GetAtomPosition(i)
        pb = conf_b.GetAtomPosition(j)
        deltas.append((pa.x - pb.x, pa.y - pb.y, pa.z - pb.z))
    arr = np.asarray(deltas, dtype=float)
    return float(np.sqrt((arr * arr).sum(axis=1).mean()))


def symmetry_corrected_rmsd_numpy(coords_a, coords_b, elements) -> float | None:
    """RMSD with permutation over symmetry-equivalent atoms of the same element.

    Used for re-docked replicate comparison: it must not depend on the order in
    which GNINA happened to emit atoms.  Only same-element permutations inside
    the same molecule are considered, and the pairing is chosen to minimise the
    fixed-frame deviation, not to minimise a fitted RMSD.
    """
    import itertools

    import numpy as np

    a = np.asarray(coords_a, dtype=float)
    b = np.asarray(coords_b, dtype=float)
    if a.shape != b.shape or len(a) < 1:
        return None
    groups: dict[int, list[int]] = {}
    for index, element in enumerate(elements):
        groups.setdefault(int(element), []).append(index)
    best = None
    per_group_options = []
    for element, indices in groups.items():
        per_group_options.append((element, list(itertools.permutations(indices))))
    # Bound the search: at most 200k combinations.
    total = 1
    for _element, perms in per_group_options:
        total *= len(perms)
    if total > 200_000:
        per_group_options = [
            (element, [tuple(indices)]) for element, indices in groups.items()
        ]
    for combo in itertools.product(*[perms for _e, perms in per_group_options]):
        mapping = []
        for perm in combo:
            mapping.extend(perm)
        deltas = a - b[mapping]
        value = float(np.sqrt((deltas * deltas).sum(axis=1).mean()))
        if best is None or value < best:
            best = value
    return best


def murcko_scaffold(mol: Chem.Mol) -> str | None:
    from rdkit.Chem.Scaffolds import MurckoScaffold

    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=True)
    except Exception:
        try:
            return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        except Exception:
            return None


def morgan_fingerprint(mol: Chem.Mol, radius: int = 2, n_bits: int = 2048):
    from rdkit.Chem import rdFingerprintGenerator

    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    return generator.GetFingerprint(mol)


def morgan_similarity(fp_a, fp_b) -> float:
    from rdkit import DataStructs

    return float(DataStructs.TanimotoSimilarity(fp_a, fp_b))


def conformer_from_rdkit(mol: Chem.Mol, random_seed: int = 42, force_field: str = "MMFF94") -> tuple[Chem.Mol, dict]:
    """Generate one reference conformation for strain analysis.

    Records exactly what was used so the caller can mark UFF fallbacks and
    failures separately instead of silently mixing them.
    """
    work = Chem.AddHs(Chem.Mol(mol))
    params = AllChem.ETKDGv3()
    params.randomSeed = int(random_seed)
    status = AllChem.EmbedMolecule(work, params)
    report: dict[str, Any] = {"engine": "ETKDGv3", "seed": int(random_seed), "embed_status": int(status)}
    if status != 0:
        report["force_field"] = None
        report["converged"] = False
        return work, report
    used = None
    converged = False
    if force_field.upper() == "MMFF94" and AllChem.MMFFHasAllMoleculeParams(work):
        result = AllChem.MMFFOptimizeMolecule(work, maxIters=2000)
        used = "MMFF94"
        converged = result == 0
    elif AllChem.UFFHasAllMoleculeParams(work):
        result = AllChem.UFFOptimizeMolecule(work, maxIters=2000)
        used = "UFF"
        converged = result == 0
    report["force_field"] = used
    report["converged"] = converged
    return work, report


def energy_of_conformer(mol: Chem.Mol, force_field: str = "MMFF94") -> float | None:
    """Single-point energy in kcal/mol, or None when the force field cannot
    parameterise the molecule (never silently reported as zero strain)."""
    work = Chem.Mol(mol)
    try:
        if force_field.upper() == "MMFF94":
            if not AllChem.MMFFHasAllMoleculeParams(work):
                return None
            props = AllChem.MMFFGetMoleculeProperties(work)
            if props is None:
                return None
            field = AllChem.MMFFGetMoleculeForceField(work, props)
        else:
            if not AllChem.UFFHasAllMoleculeParams(work):
                return None
            field = AllChem.UFFGetMoleculeForceField(work)
        if field is None:
            return None
        return float(field.CalcEnergy())
    except Exception:
        return None


def force_field_supported(mol: Chem.Mol, force_field: str = "MMFF94") -> bool:
    try:
        if force_field.upper() == "MMFF94":
            return bool(AllChem.MMFFHasAllMoleculeParams(Chem.Mol(mol)))
        return bool(AllChem.UFFHasAllMoleculeParams(Chem.Mol(mol)))
    except Exception:
        return False
