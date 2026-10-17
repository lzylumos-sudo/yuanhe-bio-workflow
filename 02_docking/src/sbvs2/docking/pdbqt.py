"""PDBQT preparation from the *approved* prepared structures (plan §14, §16). Meeko 0.8.0."""
from __future__ import annotations
import math
from pathlib import Path
from rdkit import Chem


def receptor_to_pdbqt(prepared_pdb: str | Path, out_stem: str | Path) -> Path:
    from meeko import PDBQTWriterLegacy  # noqa: F401  (import check)
    import subprocess, sys, shutil
    exe = shutil.which("mk_prepare_receptor.py") or str(Path(sys.executable).parent / "mk_prepare_receptor.py")
    r = subprocess.run([sys.executable, exe, "--read_pdb", str(prepared_pdb), "-o", str(out_stem), "-p", "--allow_bad_res"],
                       capture_output=True, text=True)
    out = Path(str(out_stem) + ".pdbqt")
    if r.returncode != 0 or not out.is_file() or out.stat().st_size == 0:
        raise RuntimeError(f"RECEPTOR_PDBQT_FAILED: {(r.stdout + r.stderr)[-600:]}")
    return out


def ligand_to_pdbqt(mol_h, conf_id: int, charges: list | None = None) -> tuple[str, dict]:
    """Return PDBQT text for one conformer. Charges (Gasteiger) injected from the *reviewed* state so that the
    approved charge model is what is docked; atom count must match (plan §14: complete, finite, mapped)."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    m = Chem.Mol(mol_h, confId=conf_id) if False else Chem.Mol(mol_h)
    keep = Chem.Conformer(m.GetConformer(conf_id)); m.RemoveAllConformers(); m.AddConformer(keep, assignId=True)
    if charges is not None:
        if len(charges) != m.GetNumAtoms() or not all(math.isfinite(c) for c in charges):
            raise ValueError("CHARGE_MAP_MISMATCH")
        for a, q in zip(m.GetAtoms(), charges):
            a.SetDoubleProp("_GasteigerCharge", float(q))
    prep = MoleculePreparation()
    setups = prep.prepare(m)
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise RuntimeError(f"LIGAND_PDBQT_FAILED: {err}")
    return pdbqt, {"n_atoms": m.GetNumAtoms(), "n_torsions": len(setups[0].flexibility_model.get("rigid_body_connectivity", {})) if hasattr(setups[0], "flexibility_model") else None}
