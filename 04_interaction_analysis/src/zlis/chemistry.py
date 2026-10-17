"""Ligand chemistry recovery.

Priority: user template (SDF/SMILES) -> local template library (configs) ->
RCSB chemical-component SMILES (cached, only if online) -> constrained inference.

Templates only supply bond orders / formal charges. Coordinates always come
from the input PDB pose, never from the template.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDetermineBonds

from .structure import Atom, ResidueGroup


def atom_uid(a: Atom) -> str:
    """Stable identity that survives hydrogen addition and re-ordering."""
    return f"M{a.model}:{a.chain}:{a.resname}{a.resseq}{a.icode}:{a.name.strip()}"


@dataclass
class LigandChem:
    mol: Chem.Mol                       # explicit-H molecule, conformer = input pose
    atoms: list[dict]                   # per RDKit index: uid, serial, name, element, source
    chem_source: str                    # template | library | rcsb | inferred | unverified
    status: str                         # ok | needs_confirmation
    formal_charge: int
    smiles: str
    notes: list[str] = field(default_factory=list)


def _load_json(p: Path) -> dict:
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _template_from_sdf(path: Path):
    sup = Chem.SDMolSupplier(str(path), removeHs=True)
    for m in sup:
        if m is not None:
            return m
    raise ValueError(f"Cannot read a molecule from {path}")


def _fetch_rcsb_smiles(resname: str, timeout: float = 12.0) -> str | None:
    """Chemical Component Dictionary SMILES (neutral/default protonation)."""
    import urllib.request

    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{resname}"
    try:
        with urllib.request.urlopen(url, timeout=timeout) as r:
            data = json.loads(r.read().decode())
    except Exception:
        return None
    desc = data.get("rcsb_chem_comp_descriptor") or {}
    for key in ("SMILES_stereo", "SMILES"):
        if desc.get(key):
            return desc[key]
    for d in data.get("pdbx_chem_comp_descriptor") or []:
        if "SMILES" in str(d.get("type", "")):
            return d.get("descriptor")
    return None


def resolve_template(resname: str, template: str | None, library_path: Path | None,
                     cache_path: Path | None, allow_network: bool = True):
    """Return (template_mol_or_None, source, note)."""
    if template:
        p = Path(template)
        if p.exists() and p.suffix.lower() in (".sdf", ".mol"):
            return _template_from_sdf(p), "template", f"user SDF template {p.name}"
        m = Chem.MolFromSmiles(template)
        if m is None:
            raise ValueError(f"Template '{template}' is neither a readable SDF/MOL file nor valid SMILES")
        return m, "template", "user SMILES template"
    lib = _load_json(library_path) if library_path else {}
    if resname in lib:
        m = Chem.MolFromSmiles(lib[resname]["smiles"])
        if m is not None:
            return m, "library", lib[resname].get("note", "local template library")
    cache = _load_json(cache_path) if cache_path else {}
    if resname in cache:
        m = Chem.MolFromSmiles(cache[resname]["smiles"])
        if m is not None:
            return m, "rcsb", "RCSB chemical component (cached); protonation state is the CCD default"
    if allow_network:
        smi = _fetch_rcsb_smiles(resname)
        if smi:
            m = Chem.MolFromSmiles(smi)
            if m is not None:
                if cache_path:
                    cache[resname] = {"smiles": smi}
                    try:
                        cache_path.write_text(json.dumps(cache, indent=1), encoding="utf-8")
                    except OSError:
                        pass
                return m, "rcsb", "RCSB chemical component; protonation state is the CCD default"
    return None, "none", ""


def _pdb_block(atoms: list[Atom]) -> str:
    lines = []
    for i, a in enumerate(atoms, 1):
        el = a.element.upper() if len(a.element) == 1 else a.element.capitalize()
        lines.append(
            f"HETATM{i:5d} {a.name:<4s} {a.resname:>3s} {a.chain[:1]:1s}{a.resseq % 10000:4d}    "
            f"{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}{1.0:6.2f}{0.0:6.2f}          {el:>2s}"
        )
    lines.append("END")
    return "\n".join(lines) + "\n"


def build_ligand(group: ResidueGroup, template: str | None = None,
                 library_path: Path | None = None, cache_path: Path | None = None,
                 allow_network: bool = True, net_charge: int | None = None) -> LigandChem:
    atoms = list(group.atoms)
    heavy = [a for a in atoms if not a.is_h]
    hyd = [a for a in atoms if a.is_h]
    order = heavy + hyd  # heavy first so template mapping is simple
    notes: list[str] = []

    raw = Chem.MolFromPDBBlock(_pdb_block(order), removeHs=False, sanitize=False, proximityBonding=True)
    if raw is None:
        raise ValueError("RDKit could not read the ligand atoms")
    n_heavy = len(heavy)

    tmpl, source, tnote = resolve_template(group.key[3], template, library_path, cache_path, allow_network)
    status = "ok"
    heavy_mol = None
    if tmpl is not None:
        hm = Chem.RWMol(raw)
        for idx in sorted(range(n_heavy, raw.GetNumAtoms()), reverse=True):
            hm.RemoveAtom(idx)
        hm = hm.GetMol()
        try:
            heavy_mol = AllChem.AssignBondOrdersFromTemplate(tmpl, hm)
            notes.append(f"Bond orders from template: {tnote}")
            if source == "rcsb":
                status = "needs_confirmation"
                notes.append("Protonation/charge state taken from the CCD default; confirm if it matters for H-bond or salt-bridge calls.")
        except Exception as ex:
            notes.append(f"Template did not match the PDB heavy-atom graph ({ex}); falling back to inference.")
            heavy_mol = None
            source = "none"

    if heavy_mol is None:
        # constrained inference: needs explicit H and a declared net charge
        source = "inferred"
        status = "needs_confirmation"
        try:
            m2 = Chem.Mol(raw)
            rdDetermineBonds.DetermineBonds(m2, charge=net_charge if net_charge is not None else 0)
            heavy_mol = m2
            notes.append(
                "No template: bond orders inferred from 3D geometry"
                + (f" with declared net charge {net_charge}." if net_charge is not None else " assuming net charge 0 (NOT confirmed).")
            )
        except Exception:
            m2 = Chem.Mol(raw)
            Chem.SanitizeMol(m2, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES, catchErrors=True)
            heavy_mol = m2
            source = "unverified"
            notes.append("Bond orders could not be established. Provide an SDF/SMILES template; chemistry-dependent interactions are unreliable.")

    # Transfer heavy-atom bond orders / charges onto the full atom set (keeping the pose's own H).
    full = Chem.RWMol(raw)
    if heavy_mol.GetNumAtoms() == raw.GetNumAtoms():
        full = Chem.RWMol(heavy_mol)
    else:
        for b in heavy_mol.GetBonds():
            fb = full.GetBondBetweenAtoms(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
            if fb is None:
                full.AddBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondType())
            else:
                fb.SetBondType(b.GetBondType())
        for a in heavy_mol.GetAtoms():
            fa = full.GetAtomWithIdx(a.GetIdx())
            fa.SetFormalCharge(a.GetFormalCharge())
            fa.SetIsAromatic(a.GetIsAromatic())
            fa.SetNoImplicit(False)
        for b in full.GetBonds():
            if b.GetBeginAtom().GetIsAromatic() and b.GetEndAtom().GetIsAromatic() and \
               heavy_mol.GetBondBetweenAtoms(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) is not None:
                hb = heavy_mol.GetBondBetweenAtoms(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
                b.SetIsAromatic(hb.GetIsAromatic())
                b.SetBondType(hb.GetBondType())
    full = full.GetMol()
    try:
        Chem.SanitizeMol(full)
    except Exception as ex:
        notes.append(f"Sanitisation warning: {ex}")
        Chem.SanitizeMol(full, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES, catchErrors=True)
        status = "needs_confirmation"

    infos: list[dict] = [
        {"uid": atom_uid(a), "serial": a.serial, "name": a.name.strip(),
         "element": a.element, "source": "input"} for a in order
    ]

    # Add hydrogens only if the pose has none on this ligand.
    if not hyd:
        withH = Chem.AddHs(full, addCoords=True)
        for idx in range(full.GetNumAtoms(), withH.GetNumAtoms()):
            at = withH.GetAtomWithIdx(idx)
            parent = at.GetNeighbors()[0].GetIdx()
            hn = sum(1 for i in infos if i["source"] == "added_h" and i.get("parent") == parent) + 1
            pname = infos[parent]["name"]
            hname = f"H{pname}" if hn == 1 else f"H{pname}{hn}"
            infos.append({
                "uid": f"{infos[parent]['uid']}:H{hn}", "serial": None,
                "name": hname, "element": "H",
                "source": "added_h", "parent": parent,
            })
            # Keep the hydrogen inside the ligand residue (ProLIF groups by residue info).
            pri = withH.GetAtomWithIdx(parent).GetPDBResidueInfo()
            if pri is not None:
                ri = Chem.AtomPDBResidueInfo()
                ri.SetName(f" {hname:<3s}"[:4])
                ri.SetResidueName(pri.GetResidueName())
                ri.SetResidueNumber(pri.GetResidueNumber())
                ri.SetChainId(pri.GetChainId())
                ri.SetInsertionCode(pri.GetInsertionCode())
                ri.SetIsHeteroAtom(True)
                at.SetMonomerInfo(ri)
        full = withH
        notes.append("Ligand hydrogens were absent and were added (riding positions); heavy-atom coordinates untouched.")
    if net_charge is not None and Chem.GetFormalCharge(full) != net_charge:
        notes.append(f"Declared net charge {net_charge} differs from template charge {Chem.GetFormalCharge(full)}.")
        status = "needs_confirmation"
    return LigandChem(
        mol=full, atoms=infos, chem_source=source, status=status,
        formal_charge=Chem.GetFormalCharge(full),
        smiles=Chem.MolToSmiles(Chem.RemoveHs(full)), notes=notes,
    )
