"""OpenMM explicit-solvent MD with independent seeds and transactional segments."""
from __future__ import annotations
import math
import os
from pathlib import Path
from ..utils.core import WorkflowError, atomic_json, read_json, write_pdb, duration_ps, steps, phase_steps, digest, now, forcefield_family

def pressure_quantity(settings):
    from openmm import unit
    label = settings.get("pressure_unit", "atm")
    if label not in ("atm", "bar"):
        raise WorkflowError("MD.pressure_unit must be atm or bar.")
    return settings["pressure"] * (unit.atmosphere if label == "atm" else unit.bar)

def ligand_template_generator(ligand, settings, cache):
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator, GAFFTemplateGenerator
    name = settings["ligand_forcefield"]
    if name == "openff-2.2.1":
        generator = SMIRNOFFTemplateGenerator
    elif name in ("gaff-2.11", "gaff-1.81"):
        generator = GAFFTemplateGenerator
    else:
        raise WorkflowError("Unsupported ligand force field.")
    if name not in generator.INSTALLED_FORCEFIELDS:
        raise WorkflowError(f"The installed openmmforcefields package does not supply {name}.")
    from ..parameterization.hf_resp import compute
    from openff.units import unit as offunit
    import numpy as np
    q, provenance = compute(ligand.to_rdkit(), Path(cache).parent/'resp', settings)
    ligand.partial_charges = q * offunit.elementary_charge
    atomic_json(Path(cache).parent/'charge_provenance.json', {"method":"HF/6-31G* two-stage RESP", "partial_charges_e":q.tolist(), "formal_charge_e":float(sum(q)), "calculation":provenance})
    result = generator(molecules=[ligand], forcefield=name, cache=str(cache))
    def assigned(molecule):
        if molecule.partial_charges is None or not np.allclose(molecule.partial_charges.m_as(offunit.elementary_charge),q,atol=1e-8,rtol=0):
            raise WorkflowError("GAFF attempted to replace the RESP charge vector")
        return True
    result._molecule_has_user_charges = assigned
    return result

def build_system(config, paths):
    import numpy as np
    import openmm as mm
    from openmm import app, unit
    settings = config["MD"]
    has_ligand = config.get("ligand") is not None
    family = forcefield_family(settings, has_ligand)
    folder = paths["md"] / "system"
    folder.mkdir(parents=True, exist_ok=True)
    pdb = app.PDBFile(str(paths["md"] / "prepared/protein.pdb"))
    model = app.Modeller(pdb.topology, pdb.positions)
    protein_count = model.topology.getNumAtoms()
    ligand_indices = []
    ligand, rd = None, None
    forcefield = app.ForceField(settings["protein_forcefield"], settings["water_forcefield"])
    if has_ligand:
        from openff.toolkit import Molecule
        from rdkit import Chem
        mols = list(Chem.SDMolSupplier(str(paths["md"] / "prepared/ligand.sdf"), removeHs=False))
        if len(mols) != 1 or mols[0] is None:
            raise WorkflowError("Prepared complex ligand must contain exactly one valid molecule.")
        rd = Chem.AddHs(mols[0], addCoords=True)
        ligand = Molecule.from_rdkit(rd, allow_undefined_stereo=False, hydrogens_are_explicit=True)
        for index, atom in enumerate(ligand.atoms):
            atom.name = f"{atom.symbol}{index+1}"
        ligand_top = ligand.to_topology().to_openmm()
        protein_chain_ids = {chain.id for chain in model.topology.chains()}
        available = [chain for chain in "ZYXWVUTSRQPONMLKJIHGFEDCBA9876543210" if chain not in protein_chain_ids]
        if not available:
            raise WorkflowError("No unused one-character PDB chain identifier for the ligand.")
        for chain in ligand_top.chains():
            chain.id = available[0]
        for residue in ligand_top.residues():
            residue.name = "LIG"
        ligand_indices = list(range(protein_count, protein_count + ligand.n_atoms))
        model.add(ligand_top, ligand.conformers[0].to_openmm())
        template = ligand_template_generator(ligand,settings,folder / "ligand_parameters.json")
        forcefield.registerTemplateGenerator(template.generator)
    from .amber_export import export_amber
    if has_ligand and (config.get("mmgbsa", {}).get("enabled", True) or config.get("mmpbsa", {}).get("enabled", False)):
        export_amber(forcefield, model, protein_count, len(ligand_indices), folder)
    else:
        atomic_json(folder / "amber_export.json", {
            "status": "NOT_APPLICABLE" if not has_ligand else "DISABLED",
            "reason": "No small molecule thermodynamic partition" if not has_ligand else "MMGBSA explicitly disabled"
        })
    model.addSolvent(forcefield, model=settings["water_model"],
                     padding=settings["padding_nm"]*unit.nanometer,
                     ionicStrength=settings["ionic_strength_molar"]*unit.molar,
                     positiveIon="Na+", negativeIon="Cl-", neutralize=True)
    system = forcefield.createSystem(model.topology, nonbondedMethod=app.PME,
              nonbondedCutoff=settings["nonbonded_cutoff_nm"]*unit.nanometer,
              switchDistance=settings["switch_distance_nm"]*unit.nanometer,
              constraints=app.HBonds, rigidWater=True, removeCMMotion=True)
    if config.get('mmpbsa',{}).get('enabled',False):
        import parmed
        export_system = forcefield.createSystem(model.topology, nonbondedMethod=app.PME, nonbondedCutoff=settings['nonbonded_cutoff_nm']*unit.nanometer, constraints=None, rigidWater=False)
        full = parmed.openmm.load_topology(model.topology,export_system,xyz=model.positions)
        full.save(str(folder/'system.prmtop'),overwrite=False)
    for index in range(system.getNumParticles()):
        if system.isVirtualSite(index):
            raise WorkflowError("Stage 1 analysis does not support virtual-site models.")
    atoms = list(model.topology.atoms())
    bonds = [(a.index,b.index) for a,b in model.topology.bonds()]
    adjacency = {i: [] for i in range(len(atoms))}
    for a,b in bonds:
        adjacency[a].append(b)
        adjacency[b].append(a)
    elements = [a.element.symbol if a.element else "" for a in atoms]
    charges = None
    for force in system.getForces():
        if isinstance(force, mm.NonbondedForce):
            charges = [float(force.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge))
                       for i in range(len(atoms))]
    if charges is None:
        raise WorkflowError("No particle charges found.")
    if has_ligand:
        expected = read_json(folder/'charge_provenance.json')['partial_charges_e']
        if not np.allclose(np.asarray(charges)[ligand_indices],expected,atol=1e-6,rtol=0):
            raise WorkflowError('Final OpenMM topology does not retain RESP atom charges')
    protein_indices = [i for i in range(protein_count) if atoms[i].residue.name not in ("HOH","WAT")]
    heavy_protein = [i for i in protein_indices if elements[i] != "H"]
    heavy_ligand = [i for i in ligand_indices if elements[i] != "H"]
    ca_indices = [i for i in protein_indices if atoms[i].name == "CA"]
    if len(ca_indices) < 3:
        raise WorkflowError("Need at least three protein C-alpha atoms for alignment.")
    # Explicit hydrogen-bond chemistry: RDKit features for ligand; standard protein atom types.
    ligand_acceptors = set()
    ligand_donors = set()
    if rd is not None:
        from rdkit.Chem import ChemicalFeatures
        from rdkit import RDConfig
        features = ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef")).GetFeaturesForMol(rd)
        for feature in features:
            if feature.GetFamily() == "Acceptor":
                ligand_acceptors.update(protein_count+i for i in feature.GetAtomIds())
            elif feature.GetFamily() == "Donor":
                ligand_donors.update(protein_count+i for i in feature.GetAtomIds())
    acceptor_names = {"ASP": {"OD1","OD2"}, "GLU":{"OE1","OE2"}, "ASN":{"OD1"}, "GLN":{"OE1"},
                     "SER":{"OG"}, "THR":{"OG1"}, "TYR":{"OH"}, "CYS":{"SG"}, "MET":{"SD"},
                     "HIS":{"ND1","NE2"}}
    acceptors = set(ligand_acceptors)
    donors = []
    for i in protein_indices + ligand_indices:
        atom = atoms[i]
        bonded_h = [j for j in adjacency[i] if elements[j] == "H"]
        if i < protein_count:
            accepted = atom.name in ("O","OXT") or atom.name in acceptor_names.get(atom.residue.name, set())
            if atom.residue.name == "HIS" and atom.name in ("ND1","NE2") and bonded_h:
                accepted = False
            if accepted:
                acceptors.add(i)
            is_donor = elements[i] in ("N","O","S")
        else:
            is_donor = i in ligand_donors
        if is_donor:
            donors.extend((i,j) for j in bonded_h)
    restraint = mm.CustomExternalForce("0.5*k_restraint*periodicdistance(x,y,z,x0,y0,z0)^2")
    restraint.addGlobalParameter("k_restraint", settings["restraint_k_kj_mol_nm2"])
    for key in ("x0","y0","z0"):
        restraint.addPerParticleParameter(key)
    xyz_nm = model.positions.value_in_unit(unit.nanometer)
    for i in heavy_protein + heavy_ligand:
        restraint.addParticle(i, list(xyz_nm[i]))
    system.addForce(restraint)
    barostat = mm.MonteCarloBarostat(pressure_quantity(settings), settings["temperature"]*unit.kelvin, 25)
    barostat.setFrequency(0)  # NVT first
    system.addForce(barostat)
    write_pdb(model.topology, model.positions, folder / "topology.pdb")
    (folder / "system.xml").write_text(mm.XmlSerializer.serialize(system))
    np.save(folder / "positions_nm.npy", np.asarray(xyz_nm))
    atomic_json(folder / "topology.json", {
        "atoms": len(atoms), "protein_indices": protein_indices, "ligand_indices": ligand_indices,
        "chain_ids": [a.residue.chain.id for a in atoms],
        "protein_chains": {chain: [i for i in protein_indices if atoms[i].residue.chain.id == chain]
                           for chain in sorted({atoms[i].residue.chain.id for i in protein_indices})},
        "protein_heavy": heavy_protein, "ligand_heavy": heavy_ligand, "ca_indices": ca_indices,
        "bonds": bonds, "elements": elements, "charges": charges,
        "donor_hydrogen_pairs": donors, "acceptors": sorted(acceptors),
        "residues": [f"{a.residue.chain.id}:{a.residue.name}{a.residue.id}{a.residue.insertionCode}" for a in atoms],
        "forcefields": {"protein_forcefield": settings["protein_forcefield"], "water_forcefield": settings["water_forcefield"],
                        "ligand_forcefield": settings["ligand_forcefield"] if has_ligand else None},
        "forcefield_family": family,
        "pressure_bar": pressure_quantity(settings).value_in_unit(unit.bar),
        "net_charge_e": sum(charges)
    })
    atomic_json(folder / "build.json", {"created": now(), "atoms": len(atoms), "protein_atoms": protein_count,
                 "ligand_atoms": ligand.n_atoms if ligand is not None else 0,
                 "system_type": config.get("system_type", "protein_ligand" if has_ligand else "protein_protein"),
                 "pressure": settings["pressure"], "pressure_unit": settings.get("pressure_unit", "atm"),
                 "pressure_bar": pressure_quantity(settings).value_in_unit(unit.bar),
                 "system_sha256": digest(folder / "system.xml")})

def numerical_preflight(simulation, temperature, *, check_velocities=True):
    """Reject nonfinite/pathological states; this is not an equilibration criterion."""
    import numpy as np
    import openmm as mm
    from openmm import unit
    state = simulation.context.getState(getEnergy=True, getPositions=True, getForces=True,
                                         getVelocities=check_velocities)
    positions = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
    energy = float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    if not math.isfinite(energy) or not np.isfinite(positions).all() or not np.isfinite(forces).all():
        raise WorkflowError("Numerical preflight failed: nonfinite energy, positions or forces; inspect/minimize a new run.")
    force_norms = np.linalg.norm(forces, axis=1)
    maximum_force = float(force_norms.max(initial=0))
    if maximum_force > 1.0e6:
        raise WorkflowError("Numerical preflight failed: extreme force (>1e6 kJ/mol/nm); review clashes and minimization before MD.")
    report = {"potential_energy_kj_mol": energy, "maximum_force_kj_mol_nm": maximum_force,
              "force_rms_kj_mol_nm": float(np.sqrt(np.mean(forces*forces))),
              "note": "Numerical guards only; these values do not demonstrate equilibration or convergence."}
    if check_velocities:
        velocities = state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond)
        kinetic = float(state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole))
        if not np.isfinite(velocities).all() or not math.isfinite(kinetic):
            raise WorkflowError("Numerical preflight failed: nonfinite velocities/kinetic energy.")
        system = simulation.system
        dof = 3*sum(system.getParticleMass(i)>0*unit.dalton for i in range(system.getNumParticles()))-system.getNumConstraints()
        if any(isinstance(force, mm.CMMotionRemover) for force in system.getForces()):
            dof -= 3
        if dof <= 0:
            raise WorkflowError("Numerical preflight failed: no positive kinetic degrees of freedom.")
        instantaneous_temperature = 2*kinetic/(dof*unit.MOLAR_GAS_CONSTANT_R.value_in_unit(unit.kilojoule_per_mole/unit.kelvin))
        if instantaneous_temperature < 0 or instantaneous_temperature > 20*temperature:
            raise WorkflowError("Numerical preflight failed: kinetic temperature exceeds 20 times target. Preserve legacy checkpoints and start a fresh run; velocities must be initialized after minimization.")
        report.update(kinetic_energy_kj_mol=kinetic, temperature_kelvin=instantaneous_temperature,
                      maximum_velocity_nm_ps=float(np.linalg.norm(velocities,axis=1).max(initial=0)),
                      kinetic_degrees_of_freedom=dof)
    return report

def initialize_new_state(simulation, settings, seed):
    """Only for a new run: minimize coordinates, then initialize deterministic velocities."""
    from openmm import unit
    if settings["minimization"]:
        simulation.minimizeEnergy(maxIterations=settings["minimization_max_iterations"])
    geometry = numerical_preflight(simulation, settings["temperature"], check_velocities=False)
    # Langevin integrators can offset sampled velocities using the current force.
    # Sampling before minimization retains an enormous velocity offset after
    # coordinate relaxation (demonstrated by the CHARMM36 public regression).
    simulation.context.setVelocitiesToTemperature(settings.get("initial_temperature",settings["temperature"])*unit.kelvin, seed)
    simulation.context.applyVelocityConstraints(1e-6)
    initialized = numerical_preflight(simulation, settings["temperature"])
    return {"initialization_order":"minimize_then_velocities", "seed":seed,
            "after_minimization":geometry, "initialized":initialized}

def platform_for(settings):
    import openmm as mm
    names = [settings["platform"]] if settings["platform"] != "auto" else (["CUDA","OpenCL"] if settings.get("require_gpu",True) else ["CUDA","OpenCL","CPU"])
    failures = []
    if settings.get("require_gpu",True) and settings["platform"] not in ("auto","CUDA","OpenCL"):
        raise WorkflowError("Research mode requires an actual GPU; use require_gpu=false only for explicit short tests.")
    for name in names:
        try:
            p = mm.Platform.getPlatformByName(name)
            props = {}
            if name in ("CUDA","OpenCL"):
                props["Precision"] = settings["precision"]
                props["DeviceIndex"] = str(os.environ.get("YUANHE_GPU_INDEX",settings.get("device_index","0")))
            if name == "CPU":
                props["Threads"] = str(settings["cpu_threads"])
            test = mm.System()
            test.addParticle(1)
            integrator = mm.VerletIntegrator(.001)
            context = mm.Context(test, integrator, p, props)
            context.setPositions([[0,0,0]])
            context.getState(getEnergy=True)
            del context, integrator
            return p, props
        except Exception as exc:
            failures.append(f"{name}: {type(exc).__name__}")
    raise WorkflowError("No usable OpenMM compute platform: " + ", ".join(failures))

def run_md(config, paths):
    for replica in range(1, config["replicates"]+1):
        run_replica(config, paths, replica)

def run_replica(config, paths, replica):
    import numpy as np
    import openmm as mm
    from openmm import app, unit
    md = config["MD"]
    target = paths["md"] / f"replica_{replica:02d}"
    target.mkdir(parents=True, exist_ok=True)
    manifest_file = target / "progress.json"
    base = paths["md"] / "system"
    system_hash = digest(base / "system.xml")
    state = read_json(manifest_file) if manifest_file.exists() else {
        "phase": "new", "production_steps": 0, "segments": [], "seed": md["seed"]+replica,
        "system_sha256": system_hash
    }
    if state["system_sha256"] != system_hash:
        raise WorkflowError("System changed; create a new run instead of resuming.")
    # Validate all committed history before any resumed integration, not only
    # when the trajectory eventually reaches analysis.
    committed_steps = sum(segment["steps"] for segment in state["segments"])
    if committed_steps != state["production_steps"]:
        raise WorkflowError("Committed segment steps differ from production progress.")
    for segment in state["segments"]:
        for key in ("dcd", "csv"):
            if digest(target / segment[key]) != segment[key+"_sha256"]:
                raise WorkflowError("Committed trajectory was changed or corrupted; resume refused.")
    if state["phase"] in ("production", "complete") and state.get("current_step") != committed_steps:
        raise WorkflowError("Checkpoint step metadata differs from committed production steps.")
    if state["phase"] != "new" and digest(target/state["checkpoint"]) != state["checkpoint_sha256"]:
        raise WorkflowError("Checkpoint checksum mismatch.")
    if state["phase"] == "complete":
        return
    seed = state["seed"]
    system = mm.XmlSerializer.deserialize((base / "system.xml").read_text())
    barostat = next(f for f in system.getForces() if isinstance(f, mm.MonteCarloBarostat))
    barostat.setRandomNumberSeed(seed+100000)
    if state["phase"] in ("npt", "production"):
        barostat.setFrequency(25)
    integ = mm.LangevinMiddleIntegrator(md["temperature"]*unit.kelvin,
                     md["friction_per_ps"]/unit.picosecond, duration_ps(md["timestep"])*unit.picoseconds)
    integ.setRandomNumberSeed(seed)
    integ.setConstraintTolerance(1e-6)
    platform, properties = platform_for(md)
    topology = app.PDBFile(str(base / "topology.pdb")).topology
    sim = app.Simulation(topology, system, integ, platform, properties)
    state["platform"] = platform.getName()
    state["platform_properties"] = {key: platform.getPropertyValue(sim.context,key) for key in platform.getPropertyNames()}
    if state["phase"] != "new":
        checkpoint = target / state["checkpoint"]
        if digest(checkpoint) != state["checkpoint_sha256"]:
            raise WorkflowError("Checkpoint checksum mismatch.")
        try:
            sim.loadCheckpoint(str(checkpoint))
        except Exception as exc:
            raise WorkflowError("Checkpoint incompatible with this hardware/OpenMM. Start a new run; XML states are preserved for manual recovery.") from exc
    else:
        sim.context.setPositions(np.load(base / "positions_nm.npy")*unit.nanometer)
    if state["phase"] != "new" and sim.currentStep != state["current_step"]:
        raise WorkflowError("Loaded checkpoint step differs from its committed manifest.")
    if state["phase"] == "production":
        loaded_ps=sim.context.getState().getTime().value_in_unit(unit.picosecond)
        expected_ps=state["production_steps"]*duration_ps(md["timestep"])
        if not math.isclose(loaded_ps,expected_ps,abs_tol=1e-5):
            raise WorkflowError("Loaded checkpoint time differs from committed production progress.")
    stages = phase_steps(config)
    def checkpoint(phase, label):
        # Manifest always points to a fully written checkpoint; previous file is retained.
        ck = target / (label + ".chk")
        tmp = ck.with_suffix(".chk.tmp")
        sim.saveCheckpoint(str(tmp))
        os.replace(tmp, ck)
        sim.saveState(str(target / (label + ".xml")))
        state.update(phase=phase, checkpoint=ck.name, checkpoint_sha256=digest(ck), updated=now(),
                     current_step=sim.currentStep)
        atomic_json(manifest_file, state)
    def phase_run(name, nsteps):
        handle = (target / (name+".csv")).open("w")
        reporter = app.StateDataReporter(handle, min(steps(md["report_interval"],md["timestep"]), nsteps),
                    step=True, time=True, potentialEnergy=True, temperature=True, volume=True, density=True, separator=",")
        sim.reporters = [reporter]
        try:
            if name=="nvt" and md.get("heating_duration"):
                heat=steps(md["heating_duration"],md["timestep"]);done=0
                for boundary in range(1,21):
                    end=round(heat*boundary/20)
                    if end<=done:continue
                    integ.setTemperature((md.get("initial_temperature",50)+(md["temperature"]-md.get("initial_temperature",50))*boundary/20)*unit.kelvin)
                    sim.step(end-done);done=end
                integ.setTemperature(md["temperature"]*unit.kelvin);sim.step(nsteps-heat)
            else:sim.step(nsteps)
        finally:
            sim.reporters = []
            handle.close()
        energy = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not math.isfinite(energy):
            raise WorkflowError(f"Non-finite energy during {name}.")
    if state["phase"] == "new":
        numerical = initialize_new_state(sim, md, seed)
        atomic_json(target / "initialization.json", numerical)
        state["initialization_order"] = "minimize_then_velocities"
        checkpoint("minimized", "minimized")
    if state["phase"] == "minimized":
        # Checkpoint velocities are preserved, never silently reassigned on resume.
        numerical_preflight(sim, md["temperature"])
        phase_run("nvt", stages["nvt"])
        checkpoint("nvt", "nvt")
    if state["phase"] == "nvt":
        barostat.setFrequency(25)
        sim.context.reinitialize(preserveState=True)
        phase_run("npt", stages["npt"])
        checkpoint("npt", "npt")
    if state["phase"] == "npt":
        sim.context.setParameter("k_restraint", 0.0)
        sim.currentStep = 0
        sim.context.setTime(0*unit.picosecond)
        checkpoint("production", "production_start")
    interval = steps(md["report_interval"], md["timestep"])
    chunk = steps(md.get("segment_length", md["checkpoint_interval"]), md["timestep"])
    total = stages["production"]
    while state["production_steps"] < total:
        index = len(state["segments"])
        label = f"part_{index:05d}"
        count = min(chunk, total - state["production_steps"])
        dcd = target / (label + ".dcd")
        csv = target / (label + ".csv")
        # An uncommitted interrupted part can be regenerated from the committed checkpoint.
        with csv.open("w") as log:
            sim.reporters = [
                app.DCDReporter(str(dcd), interval, enforcePeriodicBox=False),
                app.CheckpointReporter(str(target / "inflight.chk"), steps(md["checkpoint_interval"],md["timestep"])),
                app.StateDataReporter(log, interval, step=True, time=True, potentialEnergy=True,
                                     kineticEnergy=True, temperature=True, volume=True, density=True, separator=",")
            ]
            try:
                sim.step(count)
            finally:
                for reporter in sim.reporters:
                    if isinstance(reporter, app.DCDReporter):
                        reporter._out.flush()
                        reporter._out.close()
                sim.reporters = []
        energy = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not math.isfinite(energy):
            raise WorkflowError("Non-finite production energy.")
        state["production_steps"] += count
        state["segments"].append({"dcd":dcd.name, "csv":csv.name, "steps":count,
                                 "dcd_sha256":digest(dcd), "csv_sha256":digest(csv), "checkpoint_policy":"atomic segment commit; interrupted segment replayed from previous committed checkpoint"})
        checkpoint("production", label)
        print(f"replica {replica}: {state['production_steps']}/{total} production steps", flush=True)
    write_pdb(topology, sim.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(),
              target / "final.pdb")
    state["phase"] = "complete"
    state["completed"] = now()
    state["production_ns"] = state["production_steps"]*duration_ps(md["timestep"])/1000
    atomic_json(manifest_file, state)
    del sim, integ
