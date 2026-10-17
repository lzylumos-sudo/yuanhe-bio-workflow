"""Protein preparation. Heavy-atom coordinates are never changed."""
from __future__ import annotations

import io
from pathlib import Path

import numpy as np

from .structure import Atom, AMINO_ACIDS


def _pdb_line(i: int, a: Atom) -> str:
    el = a.element.upper() if len(a.element) == 1 else a.element.capitalize()
    name = a.name if len(a.name) == 4 else (" " + a.name.strip()).ljust(4)
    return (f"ATOM  {i % 100000:5d} {name:<4s} {a.resname:>3s} {a.chain[:1]:1s}{a.resseq:4d}{a.icode[:1] or ' ':1s}   "
            f"{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}{a.occupancy:6.2f}{a.bfactor:6.2f}          {el:>2s}")


def write_atoms_pdb(atoms: list[Atom], path: Path) -> None:
    with open(path, "w") as fh:
        for i, a in enumerate(atoms, 1):
            fh.write(_pdb_line(i, a) + "\n")
        fh.write("END\n")


def _temp_oxt(atoms: list[Atom]) -> tuple[list[Atom], list[tuple]]:
    """Add a placeholder OXT to every residue whose C is not peptide-bonded to a following N.

    Crystal structures often omit OXT. OpenMM needs it to match a terminal template.
    The placeholder is used only to place hydrogens and is deleted again afterwards,
    so exported files contain real atoms only.
    """
    from dataclasses import replace
    order, by = [], {}
    for a in atoms:
        if a.res_key not in by:
            by[a.res_key] = []
            order.append(a.res_key)
        by[a.res_key].append(a)
    added, out = [], []
    for i, key in enumerate(order):
        res = by[key]
        out.extend(res)
        names = {a.name.strip(): a for a in res}
        if "C" not in names or "OXT" in names or "CA" not in names or "O" not in names:
            continue
        nxt_n = None
        if i + 1 < len(order):
            nxt = {a.name.strip(): a for a in by[order[i + 1]]}
            nxt_n = nxt.get("N")
        c = np.array(names["C"].xyz)
        if nxt_n is not None and np.linalg.norm(c - np.array(nxt_n.xyz)) < 2.0:
            continue
        u1 = np.array(names["CA"].xyz) - c
        u2 = np.array(names["O"].xyz) - c
        d = -(u1 / np.linalg.norm(u1) + u2 / np.linalg.norm(u2))
        d /= np.linalg.norm(d)
        pos = c + 1.25 * d
        ox = replace(names["C"], name="OXT", element="O", x=float(pos[0]), y=float(pos[1]), z=float(pos[2]))
        out.append(ox)
        added.append((ox.chain, ox.resseq, ox.icode))
    return out, added


def prepare_protein(polymer_atoms: list[Atom], workdir: Path, ph: float = 7.0):
    """Return (prepared_pdb_path, info dict).

    * If the polymer already carries hydrogens they are kept as is.
    * Otherwise OpenMM Modeller adds them; every heavy atom is then verified
      to sit exactly where it was in the input.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    src = workdir / "protein_input_original.pdb"
    write_atoms_pdb(list(polymer_atoms), src)          # exactly what the user supplied
    atoms, temp_oxt = _temp_oxt(list(polymer_atoms))
    tmp_in = workdir / "protein_input_for_openmm.pdb"
    write_atoms_pdb(atoms, tmp_in)
    n_h = sum(1 for a in atoms if a.is_h)
    n_heavy = len(atoms) - n_h
    info = {"input_heavy_atoms": n_heavy, "input_h_atoms": n_h, "hydrogens_added": False,
            "ph": ph, "notes": [], "status": "ok"}
    out = workdir / "protein_prepared.pdb"
    if n_h >= 0.5 * n_heavy:
        out.write_text(src.read_text())
        info["notes"].append("Protein already contains hydrogens; kept as provided.")
        return out, info

    from openmm import app
    from openmm.app import Modeller, ForceField, PDBFile

    pdb = PDBFile(str(tmp_in))
    modeller = Modeller(pdb.topology, pdb.positions)
    # Remove any stray partial hydrogens so templates are consistent.
    if n_h:
        modeller.delete([r for r in modeller.topology.atoms() if r.element is not None and r.element.symbol == "H"])
    ff = ForceField("amber14-all.xml")
    try:
        modeller.addHydrogens(ff, pH=ph)
        if temp_oxt:
            gone = [at for at in modeller.topology.atoms() if at.name == "OXT" and
                    (at.residue.chain.id, int(at.residue.id), (at.residue.insertionCode or "").strip()) in
                    {(c, n, ic.strip()) for c, n, ic in temp_oxt}]
            modeller.delete(gone)
    except Exception as ex:
        info["status"] = "failed"
        info["notes"].append(f"Hydrogen addition failed ({ex}). Missing heavy atoms or non-standard residues?")
        raise RuntimeError(info["notes"][-1]) from ex
    with open(out, "w") as fh:
        PDBFile.writeFile(modeller.topology, modeller.positions, fh, keepIds=True)
    info["hydrogens_added"] = True
    info["temporary_OXT_placeholders_removed"] = len(temp_oxt)
    if temp_oxt:
        info["notes"].append(f"{len(temp_oxt)} chain end(s) lack OXT in the input; a placeholder was used for hydrogen "
                             "placement and removed again. Terminal carboxylate hydrogens/charges there are not modelled.")
    info["notes"].append(
        f"Hydrogens added with OpenMM Modeller (amber14, pH {ph}). Histidine/Asp/Glu/Lys states are pH-model "
        "assignments, not experimentally verified; H-bond calls involving them are provisional."
    )

    # Verify heavy atoms did not move.
    import MDAnalysis as mda
    import warnings
    warnings.filterwarnings("ignore")
    u0 = mda.Universe(str(src))
    u1 = mda.Universe(str(out))
    h0 = u0.select_atoms("not element H").positions
    h1 = u1.select_atoms("not element H").positions
    if h0.shape != h1.shape:
        info["status"] = "check"
        info["notes"].append(f"Heavy atom count changed {h0.shape[0]} -> {h1.shape[0]}.")
    else:
        d = float(np.abs(h0 - h1).max())
        info["max_heavy_atom_shift_A"] = round(d, 4)
        if d > 0.01:
            info["status"] = "check"
            info["notes"].append(f"Heavy atoms moved by up to {d:.3f} A during preparation.")
    return out, info
