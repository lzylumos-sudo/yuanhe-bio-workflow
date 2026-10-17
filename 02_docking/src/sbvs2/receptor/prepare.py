"""Receptor preparation (plan §08-§11): complete atoms -> remove chosen waters/hetero -> protonate at pH ->
restrained minimisation -> drift QC. Every decision is explicit; nothing is deleted silently.

Produces a *versioned* receptor revision. Any setting change -> new revision_hash -> approvals invalid.
"""
from __future__ import annotations
import io, math, time
from dataclasses import dataclass, field, asdict
from pathlib import Path

from ..hashing import hash_obj, sha256_bytes
from .ingest import WATER, METALS, AA3


@dataclass
class ReceptorPolicy:
    ph: float = 7.4
    ph_confirmed: bool = False
    chains: tuple = ()                 # empty = all protein chains
    water: str = "delete_all"          # delete_all | keep_listed
    keep_waters: tuple = ()            # ((chain,resseq), ...) only used for keep_listed
    remove_hetero: tuple = ()          # resnames explicitly approved for removal (e.g. buffers)
    keep_hetero: tuple = ()            # resnames kept as cofactors (treated as unsupported by minimiser -> blocked)
    altloc: str = "A"                  # whole-residue consistent altloc
    minimize: bool = True
    heavy_restraint_k: float = 1000.0  # kJ/mol/nm^2 on heavy atoms [P]
    max_iterations: int = 500
    seed: int = 20261012

    def hash(self):
        return hash_obj(asdict(self))


class PolicyError(Exception):
    pass


def _check_policy(pol: ReceptorPolicy, found_hetero: set):
    if not pol.ph_confirmed:
        raise PolicyError("PH_NOT_CONFIRMED: user must confirm pH before preparation")
    unresolved = {h for h in found_hetero if h not in WATER} - set(pol.remove_hetero) - set(pol.keep_hetero)
    if unresolved:
        raise PolicyError("HETERO_UNRESOLVED: explicit keep/remove needed for " + ",".join(sorted(unresolved)))
    if pol.keep_hetero:
        raise PolicyError("HETERO_KEEP_UNSUPPORTED: kept cofactors/metals need parameterisation (plan §24, P2); not guessed")
    if pol.water == "keep_listed" and not pol.keep_waters:
        raise PolicyError("WATER_LIST_EMPTY")


def _filter_pdb_text(text: str, pol: ReceptorPolicy) -> tuple[str, dict]:
    out, removed = [], {"water": 0, "hetero": {}, "altloc_dropped": 0, "other_chain": 0, "hydrogens_dropped": 0}
    keepw = set(pol.keep_waters)
    for ln in text.splitlines():
        rec = ln[:6].strip()
        if rec in ("ATOM", "HETATM"):
            res, ch, seq = ln[17:20].strip(), ln[21].strip() or "A", int(ln[22:26])
            alt = ln[16].strip()
            if pol.chains and ch not in pol.chains:
                removed["other_chain"] += 1; continue
            if alt and alt != pol.altloc:
                removed["altloc_dropped"] += 1; continue
            if rec == "HETATM" and res in WATER:
                if pol.water == "keep_listed" and (ch, seq) in keepw:
                    out.append(ln[:16] + " " + ln[17:]); continue
                removed["water"] += 1; continue
            if rec == "HETATM" and res in pol.remove_hetero:
                removed["hetero"][res] = removed["hetero"].get(res, 0) + 1; continue
            if rec == "HETATM":
                raise PolicyError(f"UNHANDLED_HETERO: {res}")
            el = (ln[76:78].strip() if len(ln) >= 78 else "")
            if el == "H":
                removed["hydrogens_dropped"] += 1; continue          # re-added consistently at the chosen pH
            out.append(ln[:16] + " " + ln[17:])
        elif rec in ("TER", "END", "ENDMDL"):
            continue
    return "\n".join(out) + "\nEND\n", removed


def prepare_receptor(pdb_path: str | Path, pol: ReceptorPolicy, out_dir: str | Path) -> dict:
    """Run the pipeline. Returns manifest dict; writes receptor_prepared.pdb + manifest files."""
    from openmm import app, unit, Platform, LangevinMiddleIntegrator, CustomExternalForce
    import openmm as mm
    from pdbfixer import PDBFixer
    import numpy as np

    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    raw = Path(pdb_path).read_bytes()
    text = raw.decode("utf-8", errors="replace")
    found_het = {ln[17:20].strip() for ln in text.splitlines() if ln.startswith("HETATM")}
    _check_policy(pol, found_het)
    filt, removed = _filter_pdb_text(text, pol)

    fx = PDBFixer(pdbfile=io.StringIO(filt))
    fx.findMissingResidues()
    # never model long/unknown loops here: terminal + internal missing residues are *reported* and skipped
    missing_res = {f"{k[0]}:{k[1]}": v for k, v in fx.missingResidues.items()}
    fx.missingResidues = {}
    fx.findNonstandardResidues()
    nonstd = [(r.name, r.chain.id, r.id) for r, _ in fx.nonstandardResidues]
    if nonstd:
        raise PolicyError(f"NONSTANDARD_RESIDUES: {nonstd} need an explicit mapping")
    fx.findMissingAtoms()
    missing_atoms = {f"{r.chain.id}:{r.name}{r.id}": list(a.name for a in atoms) for r, atoms in fx.missingAtoms.items()}
    missing_term = {f"{r.chain.id}:{r.name}{r.id}": list(names) for r, names in fx.missingTerminals.items()}
    fx.addMissingAtoms(seed=pol.seed)
    fx.addMissingHydrogens(pol.ph)

    n_heavy_before = sum(1 for a in fx.topology.atoms() if a.element.symbol != "H")
    ff = app.ForceField("amber14-all.xml")
    modeller = app.Modeller(fx.topology, fx.positions)
    system = ff.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
    pos0 = np.array(modeller.positions.value_in_unit(unit.angstrom))
    heavy_idx = [a.index for a in modeller.topology.atoms() if a.element.symbol != "H"]
    info = {"minimized": False}
    final_pos = modeller.positions
    if pol.minimize:
        ext = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        ext.addGlobalParameter("k", pol.heavy_restraint_k * unit.kilojoule_per_mole / unit.nanometer ** 2)
        for p in ("x0", "y0", "z0"):
            ext.addPerParticleParameter(p)
        P = modeller.positions.value_in_unit(unit.nanometer)
        for i in heavy_idx:
            ext.addParticle(i, [P[i][0], P[i][1], P[i][2]])
        system.addForce(ext)
        integ = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
        sim = app.Simulation(modeller.topology, system, integ, Platform.getPlatformByName("CPU"))
        sim.context.setPositions(modeller.positions)
        e0 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        sim.minimizeEnergy(maxIterations=pol.max_iterations)
        st = sim.context.getState(getEnergy=True, getPositions=True)
        e1 = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        final_pos = st.getPositions()
        info = {"minimized": True, "energy_before": e0, "energy_after": e1, "platform": "CPU",
                "restraint_k": pol.heavy_restraint_k, "max_iterations": pol.max_iterations}
    pos1 = np.array(final_pos.value_in_unit(unit.angstrom))

    # drift QC (P-class thresholds: backbone 0.5 A, heavy 1.0 A) - warns, never claims physical truth
    names = [a.name for a in modeller.topology.atoms()]
    bb = [i for i in heavy_idx if names[i] in ("N", "CA", "C", "O")]
    rmsd = lambda idx: float(np.sqrt(np.mean(np.sum((pos1[idx] - pos0[idx]) ** 2, axis=1))))
    drift = {"backbone_rmsd_A": rmsd(bb), "heavy_rmsd_A": rmsd(heavy_idx), "max_heavy_shift_A": float(np.max(np.linalg.norm(pos1[heavy_idx] - pos0[heavy_idx], axis=1)))}
    drift["warn_backbone"] = drift["backbone_rmsd_A"] > 0.5
    drift["warn_heavy"] = drift["heavy_rmsd_A"] > 1.0

    buf = io.StringIO()
    app.PDBFile.writeFile(modeller.topology, final_pos, buf, keepIds=True)
    pdb_out = buf.getvalue()
    # finite coordinate / H presence checks (hard)
    if not np.isfinite(pos1).all():
        raise PolicyError("NON_FINITE_COORDINATES after minimisation")
    n_h = sum(1 for a in modeller.topology.atoms() if a.element.symbol == "H")
    if n_h == 0:
        raise PolicyError("NO_HYDROGENS_ADDED")
    pdb_path_out = out_dir / "receptor_prepared.pdb"
    pdb_path_out.write_text(pdb_out)
    manifest = {
        "schema_version": 2, "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "input_sha256": sha256_bytes(raw), "policy": asdict(pol), "policy_hash": pol.hash(),
        "removed": removed, "missing_residues_reported_not_modelled": missing_res,
        "missing_atoms_completed": missing_atoms, "missing_terminals": missing_term,
        "n_heavy_atoms": n_heavy_before, "n_hydrogens": n_h, "forcefield": "amber14-all (OpenMM %s)" % mm.__version__,
        "minimization": info, "drift": drift, "prepared_sha256": sha256_bytes(pdb_out.encode()),
        "protonation_method": "PDBFixer.addMissingHydrogens (no pKa prediction; His/Asp/Glu/Lys tautomers need review)",
    }
    manifest["revision_hash"] = hash_obj({k: manifest[k] for k in ("input_sha256", "policy_hash", "prepared_sha256")})
    from ..hashing import atomic_write_json
    atomic_write_json(out_dir / "receptor_manifest.json", manifest)
    return manifest
