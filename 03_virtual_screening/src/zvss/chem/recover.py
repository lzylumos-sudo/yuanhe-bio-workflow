"""Bound-ligand chemical recovery.

Situation this solves: a GNINA pose PDB (or a complex PDB) gives us coordinates
and element labels but not bond orders, formal charges, aromaticity or defined
stereochemistry.  Re-docking scores are meaningless without those, and an MD
handoff is unsafe without them.

Contract:
* the *template* (original ligand SDF / library SDF / explicit SMILES) supplies
  the chemical graph;
* the *PDB* supplies the bound coordinates;
* heavy-atom coordinates are never regenerated, moved or re-embedded;
* everything that could not be verified is reported as such rather than assumed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem

from ..constants import ChemistryError, Reason
from . import rdio


@dataclass
class RecoveredLigand:
    mol: Chem.Mol
    template: Chem.Mol
    atom_map: dict[str, int]              # PDB atom name -> 1-based atom index in `mol`
    heavy_atom_count: int
    atom_names: list[str]
    status: str = "pass"                  # pass | flag | fail | needs_input
    reason_code: str = str(Reason.OK)
    notes: list[str] = field(default_factory=list)
    evidence: dict[str, Any] = field(default_factory=dict)

    @property
    def smiles(self) -> str | None:
        return rdio.canonical_isomeric_smiles(self.mol)


def _heavy_indices(mol: Chem.Mol) -> list[int]:
    return [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]


def element_sequence(mol: Chem.Mol, indices: list[int] | None = None) -> list[int]:
    idx = indices or range(mol.GetNumAtoms())
    return [int(mol.GetAtomWithIdx(int(i)).GetAtomicNum()) for i in idx]


def recover_from_template(
    template: Chem.Mol,
    pdb_block: str,
    *,
    atom_names: list[str] | None = None,
    resname: str = "LIG",
    keep_hydrogens: bool = False,
    max_symmetry_matches: int = 64,
) -> RecoveredLigand:
    """Merge a template chemical graph with bound coordinates from a PDB block."""
    if template is None:
        raise ChemistryError("recovery needs a template molecule; refusing to guess bonds from distances")

    observed = rdio.mol_from_pdb_block(pdb_block, sanitize=False, remove_hs=False)
    if observed is None:
        raise ChemistryError("could not parse the ligand PDB block")

    # Remove hydrogens from the observed block: template supplies protonation,
    # and the pose hydrogens (if any) are usually auto-added generic names.
    editable = Chem.RWMol(observed)
    to_remove = [atom.GetIdx() for atom in editable.GetAtoms() if atom.GetAtomicNum() == 1]
    for index in sorted(to_remove, reverse=True):
        editable.RemoveAtom(index)
    observed_heavy = editable.GetMol()
    if observed_heavy.GetNumAtoms() == 0:
        raise ChemistryError("ligand block contains no heavy atoms")

    names = list(atom_names or [])
    if not names:
        names = [
            (atom.GetPDBResidueInfo().GetName().strip() if atom.GetPDBResidueInfo() else f"{atom.GetSymbol()}{atom.GetIdx() + 1}")
            for atom in observed_heavy.GetAtoms()
        ]
    if len(names) != observed_heavy.GetNumAtoms():
        raise ChemistryError(
            f"atom name list length {len(names)} != ligand heavy atoms {observed_heavy.GetNumAtoms()}"
        )

    template_heavy = Chem.RemoveHs(Chem.Mol(template), sanitize=False)
    try:
        Chem.SanitizeMol(template_heavy)
    except Exception as exc:  # pragma: no cover - template comes from validated input
        raise ChemistryError(f"template molecule is not chemically valid: {exc}") from exc

    if template_heavy.GetNumAtoms() != observed_heavy.GetNumAtoms():
        return RecoveredLigand(
            mol=observed_heavy,
            template=template_heavy,
            atom_map={name: index + 1 for index, name in enumerate(names)},
            heavy_atom_count=observed_heavy.GetNumAtoms(),
            atom_names=names,
            status="needs_input",
            reason_code=str(Reason.CHEM_ATOM_COUNT_MISMATCH),
            notes=[
                f"template heavy atoms {template_heavy.GetNumAtoms()} != observed {observed_heavy.GetNumAtoms()}; "
                "coordinates kept as-is, chemical graph NOT recovered"
            ],
            evidence={"template_atoms": template_heavy.GetNumAtoms(), "observed_atoms": observed_heavy.GetNumAtoms()},
        )

    template_elements = sorted(element_sequence(template_heavy))
    observed_elements = sorted(element_sequence(observed_heavy))
    if template_elements != observed_elements:
        return RecoveredLigand(
            mol=observed_heavy,
            template=template_heavy,
            atom_map={name: index + 1 for index, name in enumerate(names)},
            heavy_atom_count=observed_heavy.GetNumAtoms(),
            atom_names=names,
            status="needs_input",
            reason_code=str(Reason.CHEM_ELEMENT_MISMATCH),
            notes=["element composition differs between template and bound coordinates"],
            evidence={"template_elements": template_elements, "observed_elements": observed_elements},
        )

    try:
        recovered = AllChem.AssignBondOrdersFromTemplate(template_heavy, observed_heavy)
    except Exception as exc:
        return RecoveredLigand(
            mol=observed_heavy,
            template=template_heavy,
            atom_map={name: index + 1 for index, name in enumerate(names)},
            heavy_atom_count=observed_heavy.GetNumAtoms(),
            atom_names=names,
            status="needs_input",
            reason_code=str(Reason.CHEM_BOND_ORDER_UNRECOVERED),
            notes=[f"AssignBondOrdersFromTemplate failed: {exc}"],
            evidence={"exception": str(exc)},
        )

    notes: list[str] = []
    status = "pass"
    reason = str(Reason.OK)

    try:
        Chem.SanitizeMol(recovered)
    except Exception as exc:
        return RecoveredLigand(
            mol=observed_heavy,
            template=template_heavy,
            atom_map={name: index + 1 for index, name in enumerate(names)},
            heavy_atom_count=observed_heavy.GetNumAtoms(),
            atom_names=names,
            status="fail",
            reason_code=str(Reason.CHEM_BOND_ORDER_UNRECOVERED),
            notes=[f"recovered molecule fails sanitisation: {exc}"],
            evidence={"exception": str(exc)},
        )

    recovered = Chem.Mol(recovered)
    recovered.SetProp("_Name", resname)
    if recovered.GetNumConformers():
        Chem.AssignStereochemistryFrom3D(recovered)

    # Coordinates must be identical to the observed block (no re-embedding).
    original_coords = rdio.coordinates(observed_heavy) if observed_heavy.GetNumConformers() else None
    recovered_coords = rdio.coordinates(recovered) if recovered.GetNumConformers() else None
    coordinate_drift = None
    if original_coords and recovered_coords and len(original_coords) == len(recovered_coords):
        drift = max(
            max(abs(a[0] - b[0]), abs(a[1] - b[1]), abs(a[2] - b[2]))
            for a, b in zip(original_coords, recovered_coords)
        )
        coordinate_drift = float(drift)
        if drift > 1e-6:
            notes.append(f"coordinate drift {drift:.4f} A during recovery")
            status = "flag"

    # Chemical identity check against the template.
    template_smiles = rdio.canonical_isomeric_smiles(template_heavy)
    recovered_smiles = rdio.canonical_isomeric_smiles(recovered)
    flat_template = rdio.canonical_smiles_flat(template_heavy)
    flat_recovered = rdio.canonical_smiles_flat(recovered)
    if flat_template and flat_recovered and flat_template != flat_recovered:
        status = "fail"
        reason = str(Reason.CHEM_IDENTITY_CONFLICT)
        notes.append("recovered connectivity graph differs from the template")
    elif template_smiles and recovered_smiles and template_smiles != recovered_smiles:
        # Same graph, different stereo assignment from 3D: report, do not overwrite.
        if status == "pass":
            status = "flag"
        reason = str(Reason.CHEM_STEREO_CONFLICT)
        notes.append(
            "3D coordinates imply a different stereoisomer than the template; "
            "template stereo retained for identity, pose kept as observed"
        )

    # Atom-name -> index map: AssignBondOrdersFromTemplate preserves atom order.
    atom_map = {name: index + 1 for index, name in enumerate(names)}

    # Symmetry check: is the name<->atom pairing unique at the graph level?
    symmetry = _symmetry_report(template_heavy, recovered, max_symmetry_matches)
    if symmetry["ambiguous"] and status == "pass":
        status = "flag"
        notes.append(
            "template topology admits multiple equivalent atom assignments; "
            "canonical mapping kept and equivalence classes recorded"
        )

    # Stereo centres present but undefined in the template, defined by 3D.
    undefined_stereo = Chem.FindMolChiralCenters(template_heavy, includeUnassigned=True, useLegacyImplementation=False)
    unassigned = [centre for centre in undefined_stereo if centre[1] == "?"]

    evidence = {
        "template_smiles": template_smiles,
        "recovered_smiles": recovered_smiles,
        "template_formula": rdio.mol_formula(template_heavy),
        "recovered_formula": rdio.mol_formula(recovered),
        "heavy_atoms": recovered.GetNumAtoms(),
        "net_charge": rdio.net_charge(recovered),
        "coordinate_drift_A": coordinate_drift,
        "symmetry": symmetry,
        "unassigned_stereocentres": [centre[0] for centre in unassigned],
        "hydrogens_added_explicitly": keep_hydrogens,
    }
    if unassigned:
        notes.append(
            f"template leaves {len(unassigned)} stereocentre(s) unassigned: "
            f"{[c[0] for c in unassigned]}"
        )
        if status == "pass":
            status = "flag"

    return RecoveredLigand(
        mol=recovered,
        template=template_heavy,
        atom_map=atom_map,
        heavy_atom_count=recovered.GetNumAtoms(),
        atom_names=names,
        status=status,
        reason_code=reason,
        notes=notes,
        evidence=evidence,
    )


def _symmetry_report(template: Chem.Mol, recovered: Chem.Mol, limit: int) -> dict[str, Any]:
    """Detect multiple graph-consistent assignments of template atoms to pose atoms."""
    try:
        query = Chem.Mol(template)
        target = Chem.Mol(recovered)
        Chem.RemoveStereochemistry(query)
        matches = target.GetSubstructMatches(
            query, uniquify=False, useChirality=False, maxMatches=limit + 1
        )
    except Exception as exc:
        return {"ambiguous": False, "error": str(exc), "matches": 0}
    classes: dict[int, set[int]] = {}
    for match in matches[:limit]:
        for template_index, pose_index in enumerate(match):
            classes.setdefault(template_index, set()).add(pose_index)
    equivalent = {key: sorted(value) for key, value in classes.items() if len(value) > 1}
    return {
        "ambiguous": len(matches) > 1,
        "matches": len(matches),
        "capped": len(matches) > limit,
        "equivalent_atom_classes": equivalent,
    }


def recover_from_smiles(
    smiles: str,
    pdb_block: str,
    *,
    atom_names: list[str] | None = None,
    resname: str = "LIG",
) -> RecoveredLigand:
    template = rdio.mol_from_smiles(smiles)
    if template is None:
        raise ChemistryError(f"invalid template SMILES: {smiles}")
    return recover_from_template(template, pdb_block, atom_names=atom_names, resname=resname)


def validate_bound_against_template(
    bound: Chem.Mol, template: Chem.Mol
) -> dict[str, Any]:
    """Independent post-hoc validation used by QA and by the MD preflight."""
    flat_a = rdio.canonical_smiles_flat(bound)
    flat_b = rdio.canonical_smiles_flat(template)
    iso_a = rdio.canonical_isomeric_smiles(bound)
    iso_b = rdio.canonical_isomeric_smiles(template)
    return {
        "graph_match": bool(flat_a and flat_b and flat_a == flat_b),
        "stereo_match": bool(iso_a and iso_b and iso_a == iso_b),
        "formula_match": rdio.mol_formula(bound) == rdio.mol_formula(template),
        "charge_match": rdio.net_charge(bound) == rdio.net_charge(template),
        "bound_smiles": iso_a,
        "template_smiles": iso_b,
    }
