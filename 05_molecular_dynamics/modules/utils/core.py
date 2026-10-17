"""Project-relative configuration, validation and atomic provenance files."""
from __future__ import annotations
import copy
import hashlib
import json
import math
import os
import re
from pathlib import Path
from datetime import datetime, timezone

PROJECT_ROOT = Path(__file__).resolve().parents[2]

class WorkflowError(RuntimeError):
    pass

def relative_path(value, *, must_exist=False):
    if not isinstance(value, str) or not value or "\\" in value or ":" in value:
        raise WorkflowError("Paths must be nonempty project-relative POSIX paths.")
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        raise WorkflowError("Absolute paths and parent traversal are not allowed.")
    resolved = (PROJECT_ROOT / path).resolve()
    if not resolved.is_relative_to(PROJECT_ROOT):
        raise WorkflowError("Path escapes the project root.")
    if must_exist and not resolved.is_file():
        raise WorkflowError(f"Input missing: {value}")
    return resolved

def rel(path):
    return Path(path).resolve().relative_to(PROJECT_ROOT).as_posix()

def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)

def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))

def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def now():
    return datetime.now(timezone.utc).isoformat()

def duration_ps(value):
    match = re.fullmatch(r"\s*([0-9]+(?:\.[0-9]+)?)\s*(fs|ps|ns)\s*", str(value))
    if not match:
        raise WorkflowError(f"Invalid duration: {value!r}; use fs, ps or ns.")
    result = float(match[1]) * {"fs": .001, "ps": 1., "ns": 1000.}[match[2]]
    if not math.isfinite(result) or result <= 0:
        raise WorkflowError("Durations must be finite and positive.")
    return result

def steps(value, timestep):
    quotient = duration_ps(value) / duration_ps(timestep)
    if not math.isclose(quotient, round(quotient), abs_tol=1e-6, rel_tol=1e-10):
        raise WorkflowError("Every duration must be an integer number of timesteps.")
    return int(round(quotient))

def phase_steps(config):
    md = config["MD"]
    eq = duration_ps(md["equilibration"])
    return {
        "nvt": steps(md.get("nvt") or f"{eq/2}ps", md["timestep"]),
        "npt": steps(md.get("npt") or f"{eq/2}ps", md["timestep"]),
        "production": steps(md["production"], md["timestep"]),
    }

def positive(value, label, *, zero=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise WorkflowError(f"{label} must be a finite number.")
    if value < 0 or (value == 0 and not zero):
        raise WorkflowError(f"{label} has an invalid value.")

def input_files(config):
    files = [config["complex"]["file"]]
    if config.get("ligand") is not None:
        files.append(config["ligand"]["template_sdf"])

    return list(dict.fromkeys(files))

def forcefield_family(md, has_ligand):
    """Only explicit, compatible protein/water/small-molecule combinations."""
    protein = md["protein_forcefield"]
    water = md["water_forcefield"]
    if md["water_model"] != "tip3p":
        raise WorkflowError("Phase 1 supports three-site TIP3P water models only.")
    if protein in ("amber14-all.xml", "amber14/protein.ff14SB.xml"):
        family = "AMBER_ff14SB"
        if water != "amber14/tip3p.xml":
            raise WorkflowError("AMBER ff14SB requires amber14/tip3p.xml water and ions.")
    elif protein == "charmm36.xml":
        family = "CHARMM36"
        if water != "charmm36/water.xml":
            raise WorkflowError("CHARMM36 requires charmm36/water.xml water and ions.")
        if has_ligand:
            raise WorkflowError("CHARMM36 ligand parameterization is not implemented; mixing CHARMM with GAFF/OpenFF is not supported.")
    else:
        raise WorkflowError("Supported protein force fields: amber14-all.xml, amber14/protein.ff14SB.xml, charmm36.xml.")
    ligand = md.get("ligand_forcefield")
    if has_ligand and ligand not in ("openff-2.2.1", "gaff-2.11"):
        raise WorkflowError("Supported ligand force fields: openff-2.2.1 and gaff-2.11.")
    return family

def load_config(filename="config/default.yaml", *, inputs=True):
    import yaml
    path=relative_path(filename,must_exist=True)
    config=yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config,dict):
        raise WorkflowError("Configuration must be a YAML mapping.")
    if "docking" in config or config.get("workflow")!="md_from_complex_v1":
        raise WorkflowError("Use the MD-only configuration: workflow: md_from_complex_v1 and complex.file; legacy docking configurations are not resumed.")
    for key in ("project_name","complex","protein","ligand","MD","analysis","replicates"):
        if key not in config:
            raise WorkflowError(f"Missing configuration key: {key}")
    selector = config["ligand"]
    if selector is not None and not isinstance(selector, dict):
        raise WorkflowError("ligand must be a mapping or null for a protein-only complex.")
    system_type = config.get("system_type", "protein_protein" if selector is None else "protein_ligand")
    if system_type not in ("protein_ligand", "protein_protein", "ternary"):
        raise WorkflowError("system_type must be protein_ligand, protein_protein or ternary.")
    if (system_type == "protein_protein") != (selector is None):
        raise WorkflowError("protein_protein requires ligand: null; protein_ligand and ternary require ligand chemistry.")
    if selector is not None:
        if not selector.get("template_sdf"):
            raise WorkflowError("Provide ligand.template_sdf with defined chemical information; PDB-only bond orders/charges are not guessed.")
        if not isinstance(selector.get("resname"),str) or not 1<=len(selector["resname"])<=3:
            raise WorkflowError("ligand.resname must identify the ligand residue in the PDB.")
        for key in ("chain","resid","insertion_code"):
            if selector.get(key) is not None and not isinstance(selector[key],(str,int)):
                raise WorkflowError("Ligand residue selectors must be text, integers or null.")
        mapping=selector.get("atom_mapping")
        if mapping is not None and (not isinstance(mapping,dict) or not mapping or
           any(not isinstance(k,str) or type(v) is not int or v<1 for k,v in mapping.items()) or len(set(mapping.values()))!=len(mapping)):
            raise WorkflowError("atom_mapping must map unique PDB atom names to distinct 1-based SDF atom numbers.")
    for filename in input_files(config):
        relative_path(filename,must_exist=inputs)
    if Path(config["complex"]["file"]).suffix.lower()!=".pdb":
        raise WorkflowError("The current complex coordinate input is a PDB file.")
    md,a=config["MD"],config["analysis"]
    if md["engine"]!="openmm":
        raise WorkflowError("The current MD backend is OpenMM; GROMACS is reserved for a later backend.")
    if type(config["replicates"]) is not int or config["replicates"]<1:
        raise WorkflowError("replicates must be a positive integer.")
    for key in ("temperature","pressure","padding_nm","nonbonded_cutoff_nm","friction_per_ps"):
        positive(md[key],"MD."+key)
    for key in ("ionic_strength_molar","restraint_k_kj_mol_nm2"):
        positive(md[key],"MD."+key,zero=True)
    forcefield_family(md, selector is not None)
    if md.get("pressure_unit", "atm") not in ("atm", "bar"):
        raise WorkflowError("MD.pressure_unit must be atm or bar.")
    if md["precision"] not in ("single","mixed","double"):
        raise WorkflowError("Invalid MD precision.")
    if duration_ps(md["timestep"])>.002:
        raise WorkflowError("Timestep must be <= 2fs without hydrogen mass repartitioning.")
    phase_steps(config)
    interval=steps(md["report_interval"],md["timestep"])
    checkpoint=steps(md["checkpoint_interval"],md["timestep"])
    if checkpoint<interval or checkpoint%interval:
        raise WorkflowError("checkpoint_interval must be a multiple of report_interval.")
    production=phase_steps(config)["production"]
    if production<interval or production%interval:
        raise WorkflowError("Production must contain a whole number of reporting intervals.")
    if not 0<md["switch_distance_nm"]<md["nonbonded_cutoff_nm"]<md["padding_nm"]:
        raise WorkflowError("Require 0 < switch distance < cutoff < solvent padding.")
    if md["platform"] not in ("auto","CUDA","OpenCL","CPU","Reference"):
        raise WorkflowError("Unknown OpenMM platform.")
    for section,keys in ((md,("minimization",)),(a,("rmsd","rmsf","hbond","contacts","distance")),
                         (config["protein"],("remove_water","remove_heterogens","allow_missing_residues"))):
        for key in keys:
            if type(section[key]) is not bool:
                raise WorkflowError(f"{key} must be true or false.")
    for section,keys in ((md,("seed","cpu_threads","minimization_max_iterations")),(a,("stride",))):
        for key in keys:
            if type(section[key]) is not int or section[key]<1:
                raise WorkflowError(f"{key} must be a positive integer.")
    if not 0<=config["protein"]["ph"]<=14:
        raise WorkflowError("protein.ph must be between 0 and 14.")
    positive(a["discard_ns"],"analysis.discard_ns",zero=True)
    if a["discard_ns"]>=duration_ps(md["production"])/1000:
        raise WorkflowError("analysis.discard_ns must be smaller than production duration.")
    for key in ("contact_cutoff_angstrom","hbond_distance_angstrom","hbond_angle_degrees"):
        positive(a[key],key)
    if not 0<=a["key_residue_occupancy"]<=1 or not 0<a["hbond_angle_degrees"]<=180:
        raise WorkflowError("Invalid occupancy or hydrogen-bond angle threshold.")
    return config

def config_hash(config):
    return hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()

def run_paths(run_id):
    if not re.fullmatch(r"[A-Za-z0-9_-]+", run_id):
        raise WorkflowError("run-id permits only letters, digits, hyphens and underscores.")
    return {key: PROJECT_ROOT / "results" / run_id / folder for key, folder in
            (("md", "md"), ("analysis","analysis"), ("reports","reports"))}

def stage_context(run_id):
    paths = run_paths(run_id)
    manifest = read_json(paths["md"] / "manifest.json")
    return manifest["config"], paths, manifest

def write_pdb(topology, positions, path):
    from openmm.app import PDBFile
    with Path(path).open("w") as handle:
        PDBFile.writeFile(topology, positions, handle, keepIds=True)
