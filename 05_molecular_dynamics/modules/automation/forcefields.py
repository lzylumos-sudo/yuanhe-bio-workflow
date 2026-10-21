"""Explicitly supported combinations and early v2 configuration checks.

The historical executor remains unchanged.  This registry deliberately exposes
only combinations for which this project has archived end-to-end evidence.
"""
from __future__ import annotations
import math
from ..utils.core import WorkflowError, duration_ps, phase_steps, steps

REGISTRY = {
    "amber14-gaff211-tip3p": {
        "protein_forcefields": ["amber14-all.xml", "amber14/protein.ff14SB.xml"],
        "ligand_forcefield": "gaff-2.11", "water_forcefield": "amber14/tip3p.xml",
        "water_model": "tip3p", "ions": ["Na+", "Cl-"],
        "charge_method": "HF/6-31G* ESP and two-stage RESP; no fallback",
        "backend": "OpenMM PME; HBonds constraints; rigid water",
        "nonbonded_rules": {"combining": "Lorentz-Berthelot", "coulomb_1_4": 0.8333333333333334, "lj_1_4": 0.5,
                            "runtime_evidence": "dependencies.json compatibility_payload.nonbonded_xml_settings"},
        "version_policy": "Exact native package and included parameter-file hashes recorded in dependencies.json",
        "source": "OpenMM amber14 XML and openmmforcefields GAFF 2.11 XML/dat; restored Linux archive",
        "status": "ENABLED_REQUIRES_LOCAL_ENV_AND_V2_NATIVE_VALIDATION",
        "capabilities": ["standard_protein", "protein_protein", "one_noncovalent_ligand", "ternary_one_ligand"],
        "evidence": "Historical x86-64 Linux short-run artifacts; v2 must be validated separately on Ubuntu.",
        "limitations": ["No metal coordination parameterization", "No covalent ligand", "No multiple ligand residues",
                        "No membrane or nonaqueous solvent", "No protonation-state inference"],
    }
}


def select_forcefield(config):
    md = config["MD"]
    entry = REGISTRY["amber14-gaff211-tip3p"]
    if (md.get("protein_forcefield") not in entry["protein_forcefields"] or
        md.get("water_forcefield") != entry["water_forcefield"] or
        md.get("water_model") != entry["water_model"] or
        (config.get("ligand") is not None and md.get("ligand_forcefield") != entry["ligand_forcefield"])):
        raise WorkflowError("Unsupported force-field combination in v2. Select Amber ff14SB + GAFF 2.11 (if ligand) + amber14/TIP3P. Additional combinations require a separately validated adapter.")
    return "amber14-gaff211-tip3p", entry


def _number(value, label, *, minimum=None, maximum=None, positive=False):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise WorkflowError(label + " must be a finite number.")
    if (positive and value <= 0) or (minimum is not None and value < minimum) or (maximum is not None and value > maximum):
        raise WorkflowError(label + " is outside the supported range.")


def _integer(value, label, minimum=1):
    if type(value) is not int or value < minimum:
        raise WorkflowError(label + f" must be an integer >= {minimum}.")


def _boolean(value, label):
    if type(value) is not bool:
        raise WorkflowError(label + " must be true or false.")


def validate_config(config):
    """Validate merged configuration before preparation or simulation begins."""
    allowed_top = {"MD", "protein", "analysis", "mmgbsa", "mmpbsa", "complex", "ligand", "system_type", "replicates", "workflow", "project_name", "configuration_version"}
    unknown_top = set(config) - allowed_top
    if unknown_top:
        raise WorkflowError("Unimplemented or misspelled top-level option: " + ", ".join(sorted(unknown_top)))
    for section in ("MD", "protein", "analysis", "mmgbsa"):
        if not isinstance(config.get(section), dict):
            raise WorkflowError(section + " must be a mapping.")
    md, protein, analysis, mmgbsa = (config[k] for k in ("MD", "protein", "analysis", "mmgbsa"))
    if md.get("engine") != "openmm":
        raise WorkflowError("Only the implemented OpenMM backend is supported.")
    allowed = {
        "MD": {"engine", "temperature", "pressure", "timestep", "minimization", "minimization_max_iterations", "equilibration", "nvt", "npt", "production", "protein_forcefield", "water_forcefield", "ligand_forcefield", "water_model", "padding_nm", "ionic_strength_molar", "nonbonded_cutoff_nm", "switch_distance_nm", "friction_per_ps", "restraint_k_kj_mol_nm2", "platform", "precision", "cpu_threads", "report_interval", "checkpoint_interval", "seed", "pressure_unit", "charge_method", "qm_threads", "qm_memory", "multiplicity", "net_charge", "require_gpu", "device_index", "segment_length", "heating_duration", "initial_temperature"},
        "protein": {"file", "rebuild_hydrogens", "ph", "remove_water", "remove_heterogens", "allow_missing_residues"},
        "analysis": {"rmsd", "rmsf", "hbond", "contacts", "distance", "stride", "discard_ns", "contact_cutoff_angstrom", "hbond_distance_angstrom", "hbond_angle_degrees", "key_residue_occupancy", "distance_pairs", "groups", "sasa", "sasa_stride", "sasa_points", "sasa_probe_angstrom", "salt_bridge_cutoff_angstrom", "glue_orientation_atom_indices", "formal_window_ns", "alignment_selection", "block_ns"},
        "mmgbsa": {"enabled", "igb", "salt_molar", "max_frames", "start_ns", "entropy"},
    }
    for section, keys in allowed.items():
        unknown = set(config[section]) - keys
        if unknown:
            raise WorkflowError("Unimplemented or misspelled " + section + " option: " + ", ".join(sorted(unknown)))
    select_forcefield(config)
    _integer(config.get("replicates"), "replicates")
    for key in ("temperature", "pressure", "padding_nm", "nonbonded_cutoff_nm", "switch_distance_nm", "friction_per_ps"):
        _number(md.get(key), "MD." + key, positive=True)
    for key in ("ionic_strength_molar", "restraint_k_kj_mol_nm2"):
        _number(md.get(key), "MD." + key, minimum=0)
    for key in ("seed", "cpu_threads", "minimization_max_iterations"):
        _integer(md.get(key), "MD." + key)
    if md["seed"] + config["replicates"] + 100000 >= 2**31:
        raise WorkflowError("MD.seed plus replica/barostat offsets must fit a signed 32-bit integer.")
    _boolean(md.get("minimization"), "MD.minimization")
    if md.get("pressure_unit") not in ("atm", "bar"):
        raise WorkflowError("MD.pressure_unit must be atm or bar.")
    if md.get("platform") not in ("CPU", "CUDA", "OpenCL", "Reference", "auto"):
        raise WorkflowError("Invalid MD.platform.")
    if md.get("precision") not in ("single", "mixed", "double"):
        raise WorkflowError("Invalid MD.precision.")
    if duration_ps(md.get("timestep")) > .002:
        raise WorkflowError("Time step exceeds the supported 2 fs without hydrogen mass repartitioning.")
    for name in ("equilibration", "nvt", "npt", "production", "report_interval", "checkpoint_interval"):
        if not md.get(name):
            raise WorkflowError("Set MD." + name + " explicitly.")
        steps(md[name], md["timestep"])
    if not math.isclose(duration_ps(md["equilibration"]), duration_ps(md["nvt"]) + duration_ps(md["npt"]), abs_tol=1e-8):
        raise WorkflowError("MD.equilibration must equal MD.nvt + MD.npt; no setting is silently ignored.")
    phases = phase_steps(config)
    report = steps(md["report_interval"], md["timestep"])
    checkpoint = steps(md["checkpoint_interval"], md["timestep"])
    segment=steps(md.get("segment_length",md["checkpoint_interval"]),md["timestep"])
    if segment < report or segment % report:
        raise WorkflowError("Segment length must contain whole output intervals.")
    if checkpoint < report or checkpoint % report or phases["production"] % report:
        raise WorkflowError("Production and checkpoint intervals must contain whole reporting intervals.")
    if not 0 < md["switch_distance_nm"] < md["nonbonded_cutoff_nm"] < md["padding_nm"]:
        raise WorkflowError("Require 0 < switching distance < cutoff < solvent padding.")
    for key in ("rebuild_hydrogens", "remove_water", "remove_heterogens", "allow_missing_residues"):
        _boolean(protein.get(key), "protein." + key)
    _number(protein.get("ph"), "protein.ph", minimum=0, maximum=14)
    if not protein["remove_water"]:
        raise WorkflowError("Retained crystallographic water is not implemented in the v2 Amber/MMGBSA atom partition; set remove_water=true.")
    if protein["remove_heterogens"]:
        raise WorkflowError("Do not silently remove heterogens. Make component exclusion decisions in the input plan.")
    if protein["allow_missing_residues"]:
        raise WorkflowError("Automatic loop/sequence rebuilding is not validated. Supply a reviewed complete construct.")
    if not protein["rebuild_hydrogens"]:
        raise WorkflowError("v2 requires hydrogen rebuilding with explicit preservation of supplied histidine variants.")
    for key in ("rmsd", "rmsf", "hbond", "contacts", "distance", "sasa"):
        _boolean(analysis.get(key), "analysis." + key)
    for key in ("stride", "sasa_stride"):
        _integer(analysis.get(key), "analysis." + key)
    _integer(analysis.get("sasa_points"), "analysis.sasa_points", 20)
    _number(analysis.get("sasa_probe_angstrom"), "analysis.sasa_probe_angstrom", positive=True, maximum=9.999)
    for key in ("contact_cutoff_angstrom", "hbond_distance_angstrom"):
        _number(analysis.get(key), "analysis." + key, positive=True)
    _number(analysis.get("salt_bridge_cutoff_angstrom", 4.0), "analysis.salt_bridge_cutoff_angstrom", positive=True)
    _number(analysis.get("hbond_angle_degrees"), "analysis.hbond_angle_degrees", positive=True, maximum=180)
    _number(analysis.get("key_residue_occupancy"), "analysis.key_residue_occupancy", minimum=0, maximum=1)
    _number(analysis.get("discard_ns"), "analysis.discard_ns", minimum=0)
    total_frames = phases["production"] // report
    interval_ns = duration_ps(md["report_interval"]) / 1000
    first = max(0, math.ceil(analysis["discard_ns"] / interval_ns - 1 - 1e-10))
    first += (-first) % analysis["stride"]
    eligible = max(0, (total_frames - 1 - first) // analysis["stride"] + 1)
    if eligible < 2:
        raise WorkflowError("Analysis settings leave fewer than two production frames.")
    groups = analysis.get("groups", {})
    if not isinstance(groups, dict):
        raise WorkflowError("analysis.groups must be a mapping.")
    used_chains = set()
    for role, definition in groups.items():
        if role in ("ligand", "glue"):
            continue
        chains = definition.get("chains") if isinstance(definition, dict) else None
        if not isinstance(role, str) or not isinstance(chains, list) or not chains or any(not isinstance(c, str) for c in chains):
            raise WorkflowError("Each protein analysis group requires a named nonempty chains list.")
        if set(chains) & used_chains or len(chains) != len(set(chains)):
            raise WorkflowError("Analysis protein chain groups overlap.")
        used_chains.update(chains)
    pairs = analysis.get("distance_pairs", [])
    if not isinstance(pairs, list):
        raise WorkflowError("analysis.distance_pairs must be a list.")
    names = set()
    for pair in pairs:
        if not isinstance(pair, dict) or any(not isinstance(pair.get(k), str) or not pair[k] for k in ("name", "selection1", "selection2")):
            raise WorkflowError("Distance entries require name, selection1 and selection2 strings.")
        if pair["name"] in names:
            raise WorkflowError("Distance pair names must be unique.")
        names.add(pair["name"])
    orient = analysis.get("glue_orientation_atom_indices")
    if orient is not None and (not isinstance(orient, list) or len(orient) != 2 or any(type(i) is not int or i < 0 for i in orient) or len(set(orient)) != 2):
        raise WorkflowError("Glue orientation needs two distinct nonnegative ligand heavy-atom indices.")
    if orient is not None and config.get("ligand") is None:
        raise WorkflowError("Glue orientation requires a small molecule.")
    _boolean(mmgbsa.get("enabled"), "mmgbsa.enabled")
    _boolean(mmgbsa.get("entropy"), "mmgbsa.entropy")
    if mmgbsa["entropy"]:
        raise WorkflowError("Entropy calculation is not implemented; entropy=true is refused.")
    if type(mmgbsa.get("igb")) is not int or mmgbsa["igb"] != 5:
        raise WorkflowError("The validated MM/GBSA protocol requires igb=5 and mbondi2 radii.")
    _number(mmgbsa.get("salt_molar"), "mmgbsa.salt_molar", minimum=0)
    _integer(mmgbsa.get("max_frames"), "mmgbsa.max_frames", 2)
    _number(mmgbsa.get("start_ns"), "mmgbsa.start_ns", minimum=0)
    if mmgbsa["enabled"] and config.get("ligand") is not None:
        count = max(0, total_frames - max(0, math.ceil(mmgbsa["start_ns"] / interval_ns - 1 - 1e-10)))
        if count < 2:
            raise WorkflowError("MM/GBSA start time leaves fewer than two production frames.")
    if md.get("charge_method","resp") != "resp":
        raise WorkflowError("Formal ligand charges must use RESP; no AM1-BCC fallback.")
    if md.get("heating_duration"):
        if steps(md["heating_duration"],md["timestep"])>phases["nvt"]:raise WorkflowError("Heating must fit inside NVT.")
        _number(md.get("initial_temperature",50),"MD.initial_temperature",positive=True,maximum=md["temperature"])
    pb=config.get("mmpbsa",{})
    if pb.get("enabled"):
        from ..analysis.mmpbsa import pb_input
        pb_input(pb,2,":1",":2")
        if not 0<=float(pb.get("start_ns",400))<float(pb.get("end_ns",500)):raise WorkflowError("Invalid PB window")
        _integer(pb.get("stride",10),"mmpbsa.stride")
    return config


def list_presets():
    """Return reviewable registry entries; no unimplemented menu choices."""
    return [{"id": name, **entry} for name, entry in REGISTRY.items()] + [
        {"id": "ff19sb-opc", "status": "NOT_VALIDATED_DISABLED", "reason": "No v2 complete compatible ligand/water/ion adapter validation."},
        {"id": "charmm36m-cgenff", "status": "NOT_VALIDATED_DISABLED", "reason": "No verified CGenFF parameter assignment and compatible water/ion route."}]
