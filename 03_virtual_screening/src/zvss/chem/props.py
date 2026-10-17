"""Molecular property panel and transparent soft alerts.

Everything computed here is a *calculated descriptor* with a named method.  In
particular:

* ``MolLogP`` is Crippen cLogP, **not** logD at pH 7.4 -- there is no pKa model
  wired into ZVSS, and the report says so;
* Lipinski/eNTRy-style rules are emitted as soft alerts with their reason, and
  never used as a silent filter that empties a membrane-target candidate list;
* no descriptor is named "permeability", "membrane free energy" or "antibacterial
  probability".
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

from . import rdio


@dataclass
class PropertyPanel:
    values: dict[str, float | int | None]
    alerts: list[dict[str, str]]
    methods: dict[str, str]

    def to_dict(self) -> dict[str, Any]:
        return {"values": self.values, "alerts": self.alerts, "methods": self.methods}


def _pain_filters():
    from rdkit.Chem import FilterCatalog

    params = FilterCatalog.FilterCatalogParams()
    for name in ("PAINS", "PAINS_A", "PAINS_B", "PAINS_C"):
        try:
            params.AddCatalog(getattr(FilterCatalog.FilterCatalogParams.FilterCatalogs, name))
        except Exception:
            continue
    return FilterCatalog.FilterCatalog(params)


def compute_properties(mol: Chem.Mol, *, ph: float = 7.4) -> PropertyPanel:
    if mol is None:
        raise ValueError("compute_properties needs a molecule")
    mol_h = Chem.AddHs(Chem.Mol(mol))
    heavy = mol.GetNumHeavyAtoms()

    values: dict[str, float | int | None] = {
        "molecular_weight": round(Descriptors.MolWt(mol), 3),
        "exact_mass": round(Descriptors.ExactMolWt(mol), 4),
        "heavy_atoms": int(heavy),
        "clogp_crippen": round(float(Crippen.MolLogP(mol)), 3),
        "tpsa_A2": round(float(rdMolDescriptors.CalcTPSA(mol)), 3),
        "formal_charge": int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms())),
        "h_bond_donors": int(Lipinski.NumHDonors(mol)),
        "h_bond_acceptors": int(Lipinski.NumHAcceptors(mol)),
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
        "ring_count": int(rdMolDescriptors.CalcNumRings(mol)),
        "aromatic_rings": int(rdMolDescriptors.CalcNumAromaticRings(mol)),
        "aromatic_atom_fraction": round(
            sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic()) / max(heavy, 1), 4
        ),
        "fraction_csp3": round(float(rdMolDescriptors.CalcFractionCSP3(mol)), 4),
        "n_stereocentres": len(Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)),
        "labute_asa_A2": round(float(Descriptors.LabuteASA(mol)), 3),
        "n_acidic_groups": len(
            mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]"))
        ) if Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]") else None,
        "n_basic_groups": len(
            mol.GetSubstructMatches(Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(N-C=O);!$(N-S=O)]"))
        ) if Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(N-C=O);!$(N-S=O)]") else None,
    }

    alerts: list[dict[str, str]] = []

    def alert(code: str, level: str, message: str, source: str) -> None:
        alerts.append({"code": code, "level": level, "message": message, "source": source})

    mw = values["molecular_weight"] or 0
    clogp = values["clogp_crippen"] or 0
    tpsa = values["tpsa_A2"] or 0
    hbd = values["h_bond_donors"] or 0
    hba = values["h_bond_acceptors"] or 0

    if mw > 500:
        alert("MW_GT_500", "soft", f"MW {mw:.1f} > 500 (Lipinski-style soft alert)", "Lipinski rule of five (heuristic)")
    if clogp > 5:
        alert("CLOGP_GT_5", "soft", f"Crippen cLogP {clogp:.2f} > 5", "Lipinski rule of five (heuristic)")
    if hbd > 5:
        alert("HBD_GT_5", "soft", f"HBD {hbd} > 5", "Lipinski rule of five (heuristic)")
    if hba > 10:
        alert("HBA_GT_10", "soft", f"HBA {hba} > 10", "Lipinski rule of five (heuristic)")
    if tpsa > 140:
        alert("TPSA_GT_140", "soft", f"TPSA {tpsa:.1f} > 140", "generic physicochemical heuristic")

    # Gram-negative accumulation heuristics (eNTRy-style).  Explicitly soft:
    # the pocket in this project is at the inner/periplasmic face of a membrane
    # protein, so "must cross the inner membrane" is an unproven assumption.
    if mw and (mw > 600 or mw < 200):
        alert("ENTRY_MW_RANGE", "soft", f"MW {mw:.1f} outside the 200-600 heuristic window", "eNTRy-style heuristic (unvalidated)")
    if clogp < 1:
        alert("ENTRY_LOW_LOGP", "soft", f"cLogP {clogp:.2f} < 1", "eNTRy-style heuristic (unvalidated)")
    if not (values["n_basic_groups"] or 0):
        alert("ENTRY_NO_BASIC_N", "soft", "no basic nitrogen detected by SMARTS", "eNTRy-style heuristic (unvalidated)")
    if hbd > 3:
        alert("ENTRY_HBD_HIGH", "soft", f"HBD {hbd} > 3", "eNTRy-style heuristic (unvalidated)")

    if (values["rotatable_bonds"] or 0) > 10:
        alert("ROTABLE_GT_10", "soft", f"{values['rotatable_bonds']} rotatable bonds", "conformational-entropy heuristic")
    if (values["aromatic_rings"] or 0) <= 10 and (values["ring_count"] or 0) > 6:
        alert("MANY_RINGS", "soft", f"{values['ring_count']} rings", "complexity heuristic")

    try:
        catalog = _pain_filters()
        entry = catalog.GetFirstMatch(mol)
        if entry is not None:
            alert("PAINS", "soft", f"PAINS substructure: {entry.GetDescription()}", "RDKit FilterCatalog PAINS")
    except Exception as exc:  # pragma: no cover
        alert("PAINS_UNAVAILABLE", "info", f"PAINS catalog unavailable: {exc}", "RDKit FilterCatalog")

    # Aggregation / promiscuity proxy: very flat, very lipophilic, many aromatics.
    if (values["aromatic_atom_fraction"] or 0) > 0.6 and clogp > 4.5:
        alert("AGGREGATION_RISK_PROXY", "soft", "high aromatic fraction with high cLogP", "heuristic proxy, not a measurement")

    methods = {
        "molecular_weight": "RDKit Descriptors.MolWt",
        "clogp_crippen": "RDKit Crippen.MolLogP (NOT logD at pH 7.4; no pKa model in ZVSS)",
        "tpsa_A2": "RDKit rdMolDescriptors.CalcTPSA",
        "h_bond_donors": "RDKit Lipinski.NumHDonors",
        "h_bond_acceptors": "RDKit Lipinski.NumHAcceptors",
        "rotatable_bonds": "RDKit Lipinski.NumRotatableBonds",
        "labute_asa_A2": "RDKit Descriptors.LabuteASA",
        "n_acidic_groups": "SMARTS [CX3](=O)[OX2H1,OX1-]",
        "n_basic_groups": "SMARTS [NX3;H2,H1,H0;!$(N-C=O);!$(N-S=O)]",
        "ph_used_for_property_warnings": str(ph),
        "logd_available": "false",
        "pk_a_model": "none wired; charge state is whatever the input chemical state declares",
    }
    return PropertyPanel(values=values, alerts=alerts, methods=methods)


def alert_burden(panel: PropertyPanel, *, severe_weight: float = 2.0) -> float:
    """Scalar burden used by the selection objectives (lower is better)."""
    burden = 0.0
    for item in panel.alerts:
        if item["code"] in {"PAINS", "AGGREGATION_RISK_PROXY"}:
            burden += severe_weight
        else:
            burden += 1.0
    return burden


def smiles_or_none(mol: Chem.Mol) -> str | None:
    return rdio.canonical_isomeric_smiles(mol)
