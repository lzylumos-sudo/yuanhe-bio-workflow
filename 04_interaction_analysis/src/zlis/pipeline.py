"""End-to-end analysis for one complex PDB."""
from __future__ import annotations

import json
import re
from pathlib import Path

from . import __version__
from .chemistry import build_ligand
from .interactions import index_protein, run_analysis
from .prep import prepare_protein
from .protein_mol import build_protein_mol
from .structure import (find_ligand_candidates, read_pdb, resolve_altlocs, classify_components)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
LIBRARY = PROJECT_ROOT / "configs" / "ligand_templates.json"
CACHE = PROJECT_ROOT / "configs" / "ccd_cache.json"


def list_candidates(pdb_path: str) -> list[dict]:
    st = read_pdb(pdb_path)
    out = []
    for model, atoms in st.models.items():
        atoms, _ = resolve_altlocs(atoms)
        for c in find_ligand_candidates(atoms):
            d = c.as_dict()
            d["model"] = model
            out.append(d)
    return out


def _slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_") or "x"


def analyze(pdb_path: str, out_root: str, ligand: str | None = None, template: str | None = None,
            net_charge: int | None = None, model: int | None = None, allow_network: bool = True,
            ph: float = 7.0) -> list[dict]:
    """Analyse the chosen ligand (or every plausible ligand). Returns one result dict per ligand."""
    pdb_path = str(pdb_path)
    st = read_pdb(pdb_path)
    stem = _slug(Path(pdb_path).stem)
    results = []
    for mdl, atoms in st.models.items():
        if model is not None and mdl != model:
            continue
        atoms, notes = resolve_altlocs(atoms)
        polymer, waters, ions, _hets = classify_components(atoms)
        cands = find_ligand_candidates(atoms)
        if ligand:
            want = ligand.strip()
            cands = [c for c in cands if want in (c.label, c.group.key[3]) or
                     want == f"{c.group.key[0]}:{c.group.key[3]}{c.group.key[1]}"]
            if not cands:
                raise ValueError(f"Ligand '{ligand}' not found. Candidates: " +
                                 ", ".join(c["label"] for c in list_candidates(pdb_path)))
        elif len(cands) > 1 and sum(c.kind == "ligand" for c in cands) > 1:
            notes.append("Multiple ligand candidates; each analysed separately: " + ", ".join(c.label for c in cands))
        if not cands:
            raise ValueError("No ligand candidate found (non-polymer, >=4 heavy atoms, within 8 A of protein).")
        polymer_atoms = [a for g in polymer for a in g.atoms]
        for c in cands:
            lab = f"{c.label}" + (f"_m{mdl}" if len(st.models) > 1 else "")
            out_dir = Path(out_root) / f"{stem}__{_slug(lab)}"
            out_dir.mkdir(parents=True, exist_ok=True)
            work = out_dir / "_work"
            lig = build_ligand(c.group, template=template, library_path=LIBRARY, cache_path=CACHE,
                               allow_network=allow_network, net_charge=net_charge)
            prot_pdb, pinfo = prepare_protein(polymer_atoms, work, ph=ph)
            pmol, pstats = build_protein_mol(prot_pdb)
            pidx = index_protein(pmol, polymer_atoms, mdl)
            events = run_analysis(lig, pidx, polymer_atoms, waters, ions, mdl, c.label)
            res = {
                "zlis_version": __version__,
                "input_pdb": pdb_path, "model": mdl, "ligand": c.as_dict(),
                "chemistry": {"source": lig.chem_source, "status": lig.status, "smiles": lig.smiles,
                              "formal_charge": lig.formal_charge, "notes": lig.notes},
                "protein_prep": pinfo, "protein_charges": pstats,
                "warnings": st.warnings + notes,
                "n_waters": len(waters), "n_ions": len(ions),
                "events": events,
                "ligand_atoms": lig.atoms,
                "out_dir": str(out_dir),
            }
            (out_dir / "interactions.json").write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")
            # keep objects for downstream renderers
            res["_lig"], res["_pmol_path"], res["_polymer_atoms"], res["_group"] = lig, str(prot_pdb), polymer_atoms, c.group
            res["_waters"] = [a for g in waters for a in g.atoms if a.element.upper() == "O"]
            results.append(res)
    return results
