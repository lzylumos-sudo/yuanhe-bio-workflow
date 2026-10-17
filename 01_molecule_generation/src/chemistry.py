"""Fail-closed RDKit chemistry. No bond inference or silent fragment removal.

Standardization preserves protonation, tautomer and stereochemistry; it only
sanitizes, removes ordinary explicit H and canonicalizes atom order/SMILES.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

from rdkit import Chem, rdBase
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

STANDARDIZATION_VERSION = "rdkit-v2-isomeric-preserve-charge-1"
RDKIT_VERSION = rdBase.rdkitVersion


class ChemistryError(ValueError):
    """Invalid molecular graph or unsupported coordinate serialization."""


def stable_id(kind: str, content: Any) -> str:
    payload = json.dumps(content, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=False, allow_nan=False)
    return kind + "_" + hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]


def validate_molecule(mol: Chem.Mol) -> Chem.Mol:
    if mol is None or mol.GetNumAtoms() == 0:
        raise ChemistryError("empty or unparseable molecule")
    mol = Chem.Mol(mol)
    if any(a.GetAtomicNum() == 0 or a.HasQuery() for a in mol.GetAtoms()):
        raise ChemistryError("dummy/query atoms are forbidden")
    if any(b.HasQuery() for b in mol.GetBonds()):
        raise ChemistryError("query bonds are forbidden")
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ChemistryError("disconnected molecules are forbidden")
    try:
        Chem.SanitizeMol(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except Exception as exc:
        raise ChemistryError(f"sanitization failed: {exc}") from exc
    return mol


def parse_smiles(smiles: str) -> Chem.Mol:
    if not isinstance(smiles, str) or not smiles.strip():
        raise ChemistryError("SMILES must be nonempty text")
    # Do not accept RDKit's optional trailing molecule name as chemistry.
    if any(ch.isspace() for ch in smiles.strip()):
        raise ChemistryError("whitespace/name suffix in SMILES")
    try:
        mol = Chem.MolFromSmiles(smiles.strip(), sanitize=False)
    except Exception as exc:
        raise ChemistryError(f"SMILES parse failed: {exc}") from exc
    return validate_molecule(mol)


def standardize_smiles(smiles: str | Chem.Mol) -> str:
    mol = parse_smiles(smiles) if isinstance(smiles, str) else validate_molecule(smiles)
    mol = Chem.RemoveHs(mol)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


canonicalize_smiles = standardize_smiles


def same_molecule(mol: Chem.Mol, canonical: str) -> bool:
    """Graph identity with stereo preserved where the canonical SMILES specifies it.

    3D coordinates make RDKit assign stereo to centres the SMILES left unspecified;
    that is not a mismatch, but any specified centre or double bond must agree.
    """
    observed = standardize_smiles(mol)
    if observed == canonical:
        return True
    expected = parse_smiles(canonical)
    flat = lambda m: Chem.MolToSmiles(m, canonical=True, isomericSmiles=False)
    if flat(Chem.RemoveHs(validate_molecule(mol))) != flat(expected):
        return False
    return Chem.RemoveHs(validate_molecule(mol)).HasSubstructMatch(expected, useChirality=True)


def compound_id(smiles: str | Chem.Mol) -> str:
    return stable_id("cmp", [STANDARDIZATION_VERSION, standardize_smiles(smiles)])


def state_id(smiles: str | Chem.Mol, compound: str | None = None) -> str:
    canonical = standardize_smiles(smiles)
    return stable_id("sta", [compound or compound_id(canonical), canonical, STANDARDIZATION_VERSION])


def provenance_id(record: Mapping[str, Any]) -> str:
    return stable_id("prv", {k: v for k, v in record.items() if k != "provenance_id"})


def pose_id(mol: Chem.Mol, state: str | None = None, coordinate_frame: str = "") -> str:
    mol = validate_molecule(mol)
    _validate_coordinates(mol)
    # Rank atoms by graph symmetry class (ties unbroken, stereo flags ignored so an SDF
    # round trip cannot reorder them), then break ties by the 4-decimal coordinates
    # that the SDF format stores. Atom order in the input is therefore immaterial.
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=False))
    positions = mol.GetConformer().GetPositions()
    coords = [[round(float(v), 4) + 0.0 for v in positions[i]] for i in range(mol.GetNumAtoms())]
    coords = [c for _, c in sorted(zip(ranks, coords), key=lambda rc: (rc[0], rc[1]))]
    return stable_id("pos", [state or state_id(mol), coordinate_frame, coords])


def calculate_properties(smiles: str | Chem.Mol) -> dict:
    mol = parse_smiles(smiles) if isinstance(smiles, str) else validate_molecule(smiles)
    return {"MW": Descriptors.MolWt(mol), "cLogP": Crippen.MolLogP(mol),
            "TPSA": rdMolDescriptors.CalcTPSA(mol), "HBD": Lipinski.NumHDonors(mol),
            "HBA": Lipinski.NumHAcceptors(mol), "RotB": Lipinski.NumRotatableBonds(mol),
            "formal_charge": Chem.GetFormalCharge(mol)}


properties = calculate_properties


def _validate_coordinates(mol: Chem.Mol) -> None:
    if not mol.GetNumConformers():
        raise ChemistryError("coordinates required; no conformer present")
    conf = mol.GetConformer()
    if any(not math.isfinite(float(v)) for p in conf.GetPositions() for v in p):
        raise ChemistryError("nonfinite coordinates")


def write_sdf(mol: Chem.Mol, path: str | Path, properties: Mapping[str, Any] | None = None) -> Path:
    """Write the supplied graph and coordinates, never infer edges from distances."""
    mol = validate_molecule(mol)
    _validate_coordinates(mol)
    canonical = standardize_smiles(mol)
    mol.SetProp("canonical_smiles", canonical)
    mol.SetProp("standardization_version", STANDARDIZATION_VERSION)
    for key, value in (properties or {}).items():
        if key in {"canonical_smiles", "standardization_version"}:
            if str(value) != mol.GetProp(key):
                raise ChemistryError(f"reserved SDF property mismatch: {key}")
        mol.SetProp(str(key), str(value))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(path))
    try:
        writer.write(mol)
    finally:
        writer.close()
    restored = read_sdf(path)
    if len(restored) != 1 or not same_molecule(restored[0], canonical):
        raise ChemistryError("SDF graph/stereochemistry roundtrip mismatch")
    return path


def read_sdf(path: str | Path) -> list[Chem.Mol]:
    """Strict graph-based SDF import; malformed legacy coordinate-only files fail."""
    path = Path(path)
    raw = path.read_bytes()
    if not raw.strip() or not raw.rstrip().endswith(b"$$$$"):
        raise ChemistryError("missing SDF record terminator")
    expected = raw.count(b"$$$$")
    supplier = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False, strictParsing=True)
    molecules = []
    try:
        for mol in supplier:
            mol = validate_molecule(mol)
            _validate_coordinates(mol)
            if mol.HasProp("canonical_smiles") and not same_molecule(mol, mol.GetProp("canonical_smiles")):
                raise ChemistryError("SDF canonical_smiles disagrees with molecular graph")
            if mol.HasProp("standardization_version") and mol.GetProp("standardization_version") != STANDARDIZATION_VERSION:
                raise ChemistryError("unsupported SDF standardization_version")
            molecules.append(mol)
    except ChemistryError:
        raise
    except Exception as exc:
        raise ChemistryError(f"invalid SDF: {exc}") from exc
    if len(molecules) != expected or not molecules:
        raise ChemistryError("malformed SDF records")
    return molecules
