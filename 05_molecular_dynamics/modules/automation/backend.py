"""Audited v2 adapters around the unchanged historical OpenMM executor.

All native imports are deferred: input inspection remains usable on a Mac,
while preparation, simulation and checkpoint audits require the Linux runtime.
"""
from __future__ import annotations
import ast
import csv
import importlib
import importlib.metadata
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import math
import platform
import shutil
import sys
from pathlib import Path
from ..utils.core import WorkflowError, atomic_json, read_json, digest, duration_ps, steps, phase_steps, write_pdb, now
from .forcefields import validate_config, select_forcefield

ROOT = Path(__file__).resolve().parents[2]


def _require(condition, message):
    if not condition:
        raise WorkflowError(message)


def _path(value):
    p = Path(value)
    return p if p.is_absolute() else ROOT / p


def dependency_report(config=None):
    """Probe native imports and actual external executables, without starting MD."""
    has_ligand = config is None or config.get("ligand") is not None
    mmgbsa = has_ligand and (config is None or config.get("mmgbsa", {}).get("enabled", True))
    packages = ["yaml", "numpy", "pandas", "matplotlib.pyplot", "MDAnalysis", "MDAnalysis.lib.distances", "Bio.PDB.SASA", "openmm.app", "pdbfixer"]
    if has_ligand:
        packages += ["rdkit.Chem", "rdkit.Chem.AllChem", "rdkit.Chem.rdFMCS", "rdkit.Chem.ChemicalFeatures", "openff.toolkit", "openmmforcefields.generators"]
    if mmgbsa:
        packages += ["parmed.tools.actions"]
    imports = []
    for name in packages:
        try:
            module = importlib.import_module(name)
            distributions = {"yaml": "PyYAML", "Bio": "biopython", "openff": "openff-toolkit", "parmed": "ParmEd"}
            distribution = distributions.get(name.split(".")[0], name.split(".")[0])
            try:
                version = importlib.metadata.version(distribution)
            except importlib.metadata.PackageNotFoundError:
                parent = importlib.import_module(name.split(".")[0])
                version = str(getattr(parent, "__version__", "not exposed"))
            imports.append({"name": name, "available": True, "version": version})
        except Exception as exc:
            imports.append({"name": name, "available": False, "error": f"{type(exc).__name__}: {exc}"})
    executables = ["antechamber", "parmchk2", "sqm", "tleap"] if has_ligand else []
    if mmgbsa:
        executables += ["MMPBSA.py", "cpptraj", "sander"]
    commands = [{"name": name, "available": shutil.which(name) is not None, "path": shutil.which(name)} for name in executables]
    # Include an inventory derived from the current executor, rather than relying
    # solely on an old requirements file. These are source observations, not imports.
    observed = set()
    for source in (ROOT / "modules").rglob("*.py"):
        for node in ast.walk(ast.parse(source.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Import):
                observed.update(n.name.split(".")[0] for n in node.names)
            elif isinstance(node, ast.ImportFrom) and not node.level and node.module:
                observed.add(node.module.split(".")[0])
    stdlib = getattr(sys, "stdlib_module_names", set())
    external = sorted(observed - set(stdlib) - {"modules"})
    available = all(x["available"] for x in imports + commands)
    host_ok = platform.system() == "Linux" and platform.machine() in ("x86_64", "amd64")
    details = {"observed_import_roots": external, "python": sys.version, "python_executable": sys.executable,
               "host": {"system": platform.system(), "machine": platform.machine()}, "target": "Ubuntu/Linux x86-64",
               "target_host": host_ok, "imports": imports, "executables": commands,
               "dependencies_available": available, "ready_for_native_validation": available and host_ok,
               "scope": "Dependency presence/import probe only; no parameterization, checkpoint loading or MD was performed."}
    if available:
        try:
            from openmm import app
            from .forcefields import REGISTRY
            chosen = config["MD"] if config is not None else {"protein_forcefield": "amber14-all.xml", "water_forcefield": "amber14/tip3p.xml"}
            app.ForceField(chosen["protein_forcefield"], chosen["water_forcefield"])
            details["forcefield_xml_loadable"] = True
            if has_ligand:
                from openmmforcefields.generators import GAFFTemplateGenerator
                details["gaff211_installed"] = "gaff-2.11" in GAFFTemplateGenerator.INSTALLED_FORCEFIELDS
                details["ready_for_native_validation"] &= details["gaff211_installed"]
        except Exception as exc:
            details["forcefield_xml_loadable"] = False
            details["forcefield_error"] = f"{type(exc).__name__}: {exc}"
            details["ready_for_native_validation"] = False
    compatibility = environment_compatibility(imports, commands, config)
    details.update(compatibility)
    if any(value == "MISSING" for value in compatibility["compatibility_payload"]["forcefield_sha256"].values()):
        details["ready_for_native_validation"] = False
    details["ready"] = details["ready_for_native_validation"]
    details["status"] = "PASS" if details["ready"] else "NOT_READY"
    return details



def normalized_executable_sha256(path, prefix=None):
    """Ignore relocation edits while retaining executable content and padding.

    conda-pack rewrites shebangs and fixed-length binary prefix slots. Binary
    slots are normalized in place to their original total length, so arbitrary
    binary zero runs and other content are not discarded.
    """
    content = Path(path).read_bytes()
    known_prefix = str(Path(prefix or sys.prefix)).encode()
    if b"\x00" in content[:4096] or content.startswith(b"\x7fELF"):
        pattern = re.compile(re.escape(known_prefix) + rb"([^\x00\r\n]*)(\x00+)")
        def replace_slot(match):
            replacement = b"<ENV_PREFIX>" + match.group(1).replace(known_prefix, b"<ENV_PREFIX>")
            if len(replacement) >= len(match.group(0)):
                return match.group(0)
            return replacement + b"\x00" * (len(match.group(0)) - len(replacement))
        content = pattern.sub(replace_slot, content)
    else:
        content = content.replace(known_prefix, b"<ENV_PREFIX>")
    return hashlib.sha256(content).hexdigest()


def environment_compatibility(imports, commands, config=None):
    payload = {"python_major_minor": list(sys.version_info[:2]),
               "imports": {item["name"]: item.get("version", "UNAVAILABLE") for item in imports},
               "executable_normalized_sha256": {}, "native_extension_normalized_sha256": {}, "forcefield_sha256": {},
               "nonbonded_xml_settings": {}}
    actual_hashes = {}
    for entry in commands:
        if entry.get("available") and entry.get("path"):
            actual_hashes[entry["name"]] = digest(entry["path"])
            payload["executable_normalized_sha256"][entry["name"]] = normalized_executable_sha256(entry["path"])
        else:
            payload["executable_normalized_sha256"][entry["name"]] = "UNAVAILABLE"
    try:
        from openmm import app
        native = importlib.import_module("openmm._openmm")
        payload["native_extension_normalized_sha256"]["openmm._openmm"] = normalized_executable_sha256(native.__file__)
        actual_hashes["openmm._openmm"] = digest(native.__file__)
        data = Path(app.__file__).parent / "data"
        md = config["MD"] if config is not None else {"protein_forcefield": "amber14-all.xml", "water_forcefield": "amber14/tip3p.xml"}
        pending = [md["protein_forcefield"], md["water_forcefield"]]
        visited = set()
        while pending:
            name = pending.pop()
            if name in visited:
                continue
            visited.add(name)
            filename = data / name
            if not filename.is_file():
                payload["forcefield_sha256"][name] = "MISSING"
                continue
            payload["forcefield_sha256"][name] = digest(filename)
            xmlroot=ET.parse(filename).getroot()
            payload['nonbonded_xml_settings'][name]=[e.attrib for e in xmlroot.findall('NonbondedForce')]
            pending.extend(e.attrib["file"] for e in xmlroot.findall("Include"))
        if config is None or config.get("ligand") is not None:
            import openmmforcefields
            base = Path(openmmforcefields.__file__).parent / "ffxml/amber/gaff"
            for relative in ("dat/gaff-2.11.dat", "ffxml/gaff-2.11.xml"):
                filename = base / relative
                payload["forcefield_sha256"]["gaff/" + relative] = digest(filename) if filename.is_file() else "MISSING"
                if filename.is_file() and filename.suffix=='.xml':
                    payload['nonbonded_xml_settings']['gaff/'+relative]=[e.attrib for e in ET.parse(filename).getroot().findall('NonbondedForce')]
    except (ImportError, OSError, AttributeError) as exc:
        payload["forcefield_probe"] = type(exc).__name__
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return {"compatibility_payload": payload, "compatibility_fingerprint": hashlib.sha256(encoded).hexdigest(),
            "executable_actual_sha256": actual_hashes,
            "compatibility_note": "Fingerprint excludes paths/hostnames and normalizes conda prefix slots in executables. Actual unmodified executable hashes are retained separately. Binary checkpoint compatibility still requires a real load on the selected platform."}

def _atom_key(atom):
    r = atom.residue
    return ((r.chain.id or "").strip(), str(r.id).strip(), (r.insertionCode or "").strip(), atom.name)


def _heavy_coordinates(topology, positions):
    from openmm import unit
    xyz = positions.value_in_unit(unit.angstrom)
    result = {}
    for atom in topology.atoms():
        if atom.element is not None and atom.element.symbol != "H":
            key = _atom_key(atom)
            _require(key not in result, "Duplicate protein heavy-atom identity after normalization: " + str(key))
            result[key] = tuple(float(v) for v in xyz[atom.index])
    return result


def _residue_id_map(plan):
    residues = plan.get("residues", [])
    result = {}
    for r in residues:
        if all(k in r for k in ("chain", "resid")):
            result[(str(r["chain"]), str(r["resid"]), str(r.get("insertion_code", "")))] = r.get("id", r.get("residue_id"))
    return result


def prepare(config, paths, plan):
    """Build prepared inputs while retaining supplied histidine protonation.

    plan: protein_path; optional pose-mapped ligand_sdf; variants list containing
    chain/resid/insertion_code/variant; decisions.construct_confirmed; optional
    decisions.termini='standard_charged' or a mapping of original residue IDs to
    that value. HXT removal is an explicit chemical change, never a rename.
    """
    validate_config(config)
    decisions = plan.get("decisions", {})
    _require(decisions.get("construct_confirmed") is True or plan.get("sequence_verified") is True,
             "Confirm the intended protein construct/sequence before preparation (decisions.construct_confirmed=true). Continuous numbering alone does not prove completeness.")
    protein_path = _path(plan["protein_path"])
    _require(protein_path.is_file(), "Normalized protein coordinate file is missing.")
    _require(not plan.get("covalent_links"), "Covalent component links require a dedicated parameterization adapter.")
    if config.get("ligand") is not None:
        _require(bool(plan.get("ligand_sdf")), "One mapped, chemically resolved ligand SDF is required.")
        _require(plan.get("ligand_component_count", 1) == 1, "This backend supports exactly one noncovalent ligand component.")
    # Resolve HXT before OpenMM/PDBFixer parses names. A parser's aliases must
    # never turn this explicitly reviewed chemical edit into a silent rename.
    target = paths["md"] / "prepared"
    target.mkdir(parents=True, exist_ok=True)
    lines = protein_path.read_text().splitlines()
    residue_ids = _residue_id_map(plan)
    termini = decisions.get("termini", {})
    last_by_fragment, fragment_by_line = {}, {}
    fragment = 0
    for line_number, line in enumerate(lines):
        if line.startswith("TER"):
            fragment += 1
        elif line.startswith(("ATOM  ", "HETATM")):
            fragment_key = (line[21:22].strip(), fragment)
            fragment_by_line[line_number] = fragment_key
            last_by_fragment[fragment_key] = (line[21:22].strip(), line[22:26].strip(), line[26:27].strip())
    original_heavy, explicit_removals, cleaned = {}, [], []
    for line_number, line in enumerate(lines):
        if not line.startswith(("ATOM  ", "HETATM")):
            cleaned.append(line)
            continue
        residue_key = (line[21:22].strip(), line[22:26].strip(), line[26:27].strip())
        name = line[12:16].strip()
        element = line[76:78].strip().upper() or name.lstrip("0123456789")[:1].upper()
        key = residue_key + (name,)
        if element not in ("H", "D"):
            _require(key not in original_heavy, "Duplicate normalized protein heavy-atom identity: " + str(key))
            original_heavy[key] = tuple(float(line[i:i+8]) for i in (30, 38, 46))
        if name == "HXT":
            rid = residue_ids.get(residue_key)
            choice = termini if isinstance(termini, str) else termini.get(rid, termini.get(":".join(residue_key)))
            _require(choice == "standard_charged", "HXT terminal anomaly needs an explicit standard_charged terminus decision; it is never renamed to OXT.")
            _require(residue_key == last_by_fragment[fragment_by_line[line_number]], "HXT is not in a chain-terminal residue; manual structure review is required.")
            explicit_removals.append({"atom_identity": list(key), "element": element, "decision": "standard_charged",
                                      "action": "remove HXT record before parsing; standard terminal template builds OXT with new coordinates",
                                      "original_coordinate_angstrom": [float(line[i:i+8]) for i in (30, 38, 46)]})
        else:
            cleaned.append(line)
    preparation_input = target / "preparation_input.pdb"
    preparation_input.write_text("\n".join(cleaned) + "\n")
    from pdbfixer import PDBFixer
    from openmm import Platform, app, unit
    fixer = PDBFixer(filename=str(preparation_input), platform=Platform.getPlatformByName("CPU"))
    allowed = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
    _require(all(r.name in allowed for r in fixer.topology.residues()), "Normalized protein contains an unsupported residue or retained heterogen; an explicit parameter adapter is required.")
    parsed = _heavy_coordinates(fixer.topology, fixer.positions)
    explicit_removed_keys = {tuple(record["atom_identity"]) for record in explicit_removals if record["element"] not in ("H", "D")}
    _require(set(parsed) == set(original_heavy) - explicit_removed_keys, "PDB parser changed heavy-atom identities; review name aliases explicitly before preparation.")
    remove = [atom for atom in fixer.topology.atoms() if atom.element is not None and atom.element.symbol in ("H", "D")]
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.delete(remove)
    fixer.topology, fixer.positions = modeller.topology, modeller.positions
    fixer.findNonstandardResidues()
    _require(not fixer.nonstandardResidues, "A nonstandard protein residue remains after normalization.")
    fixer.findMissingResidues()
    gaps = {str(k): list(v) for k, v in fixer.missingResidues.items() if v}
    _require(not gaps, "Sequence gaps require a reviewed structure; automatic gap rebuilding is not performed: " + str(gaps))
    fixer.findMissingAtoms()
    missing = [{"chain": r.chain.id, "resid": r.id, "atom": a.name} for r, atoms in fixer.missingAtoms.items() for a in atoms]
    terminal_additions = [{"chain": r.chain.id, "resid": r.id, "atom": name} for r, names in fixer.missingTerminals.items() for name in names]
    fixer.addMissingAtoms(seed=config["MD"]["seed"])
    modeller = app.Modeller(fixer.topology, fixer.positions)
    variant_map = {}
    for record in plan.get("variants", []):
        variant = record["variant"]
        _require(variant in ("HID", "HIE", "HIP"), "Only explicit HID/HIE/HIP variants are supported by this adapter.")
        key = (str(record.get("chain", "")), str(record["resid"]), str(record.get("insertion_code", "")))
        _require(key not in variant_map, "Duplicate histidine variant decision.")
        variant_map[key] = variant
    variants = []
    seen = set()
    for residue in modeller.topology.residues():
        key = ((residue.chain.id or "").strip(), str(residue.id).strip(), (residue.insertionCode or "").strip())
        variant = variant_map.get(key)
        if variant is not None:
            _require(residue.name == "HIS", "A histidine variant targets a non-HIS residue.")
            seen.add(key)
        variants.append(variant)
    _require(seen == set(variant_map), "A histidine variant does not match the normalized protein.")
    ff = app.ForceField(config["MD"]["protein_forcefield"], config["MD"]["water_forcefield"])
    actual_variants = modeller.addHydrogens(ff, pH=config["protein"]["ph"], variants=variants,
                                           platform=Platform.getPlatformByName("CPU"))
    after = _heavy_coordinates(modeller.topology, modeller.positions)
    removed_heavy = {tuple(record["atom_identity"]) for record in explicit_removals if record["element"] not in ("H", "D")}
    maximum_displacement = 0.0
    for key, xyz in original_heavy.items():
        if key in removed_heavy:
            _require(key not in after, "Explicitly removed HXT still exists in prepared topology.")
            continue
        _require(key in after, "An existing heavy atom was lost or renamed: " + str(key))
        displacement = math.dist(xyz, after[key])
        maximum_displacement = max(maximum_displacement, displacement)
        _require(displacement <= 1e-4, "Protein preparation moved an existing heavy atom: " + str(key))
    target = paths["md"] / "prepared"
    target.mkdir(parents=True, exist_ok=True)
    write_pdb(modeller.topology, modeller.positions, target / "protein.pdb")
    if config.get("ligand") is not None:
        ligand = _path(plan["ligand_sdf"])
        _require(ligand.is_file(), "Mapped ligand SDF is missing.")
        shutil.copyfile(ligand, target / "ligand.sdf")
    report = {"schema_version": 2, "created": now(), "protein_input_sha256": digest(protein_path),
              "protein_output_sha256": digest(target / "protein.pdb"), "existing_heavy_atoms_preserved": len(original_heavy) - len(removed_heavy),
              "maximum_existing_heavy_displacement_angstrom": maximum_displacement,
              "explicit_chemical_removals": explicit_removals, "added_missing_heavy_atoms": missing,
              "added_terminal_atoms": terminal_additions, "added_heavy_atom_identities": [list(k) for k in after if k not in original_heavy],
              "added_heavy_atoms": [{"identity": list(k), "coordinate_angstrom": list(v)} for k, v in after.items() if k not in original_heavy],
              "histidine_variants_requested": plan.get("variants", []), "hydrogen_variants_applied": list(actual_variants),
              "ph": config["protein"]["ph"], "decisions": decisions,
              "note": "pH is a default hydrogen-placement input, not a site-specific pKa prediction. Existing heavy coordinates are retained except explicitly approved HXT removal."}
    atomic_json(target / "protein_preparation.json", report)
    return report


def build(config, paths):
    validate_config(config)
    from ..md import openmm_engine as engine
    charge_provenance = None
    original_factory = engine.ligand_template_generator
    if config.get("ligand") is not None:
        def explicit_charge_factory(ligand, settings, cache):
            nonlocal charge_provenance
            import numpy as np
            from openff.toolkit.utils import AmberToolsToolkitWrapper
            from openff.units import unit as offunit
            ligand.assign_partial_charges(partial_charge_method="am1bcc", toolkit_registry=AmberToolsToolkitWrapper(),
                                          normalize_partial_charges=True)
            q = np.asarray(ligand.partial_charges.m_as(offunit.elementary_charge), dtype=float)
            formal = float(ligand.total_charge.m_as(offunit.elementary_charge))
            _require(np.isfinite(q).all() and abs(float(q.sum()) - formal) < 1e-5, "AM1-BCC charges are nonfinite or fail formal-charge normalization.")
            charge_provenance = {"method": "AM1-BCC", "backend": "OpenFF AmberToolsToolkitWrapper / antechamber / sqm",
                                 "normalize_partial_charges": True, "formal_charge_e": formal, "partial_charges_e": q.tolist()}
            generator = original_factory(ligand, settings, cache)
            def validated_assigned_charges(molecule):
                _require(molecule.partial_charges is not None, "GAFF attempted parameterization without the explicitly assigned AmberTools charges.")
                values = np.asarray(molecule.partial_charges.m_as(offunit.elementary_charge), dtype=float)
                _require(np.isfinite(values).all() and abs(float(values.sum()) - formal) < 1e-5, "GAFF lost the verified charge vector.")
                # Unlike the legacy helper, an exactly-zero but explicitly
                # calculated charge vector must not trigger another toolkit.
                return True
            generator._molecule_has_user_charges = validated_assigned_charges
            return generator
        engine.ligand_template_generator = explicit_charge_factory
    try:
        engine.build_system(config, paths)
    finally:
        engine.ligand_template_generator = original_factory
    build_audit(config, paths, charge_provenance)
    # Validate selections against the actual built topology before any dynamics.
    import MDAnalysis as mda
    from ..analysis.trajectory import resolve_groups
    base = paths["md"] / "system"
    meta = read_json(base / "topology.json")
    resolve_groups(config, meta)
    u = mda.Universe(str(base / "topology.pdb"))
    for entry in config["analysis"].get("distance_pairs", []):
        for key in ("selection1", "selection2"):
            _require(len(u.select_atoms(entry[key])) > 0, "Empty distance selection before MD: " + entry["name"])
    orientation = config["analysis"].get("glue_orientation_atom_indices")
    if orientation is not None:
        _require(all(i < len(meta["ligand_heavy"]) for i in orientation), "Glue orientation index exceeds ligand heavy atom count.")
    return read_json(base / "build.json")



def build_audit(config, paths, charge_provenance=None):
    """Audit the assembled force-field output before minimization can begin."""
    import numpy as np
    import openmm as mm
    from openmm import app, unit
    base = paths["md"] / "system"
    meta = read_json(base / "topology.json")
    topology = app.PDBFile(str(base / "topology.pdb")).topology
    system = mm.XmlSerializer.deserialize((base / "system.xml").read_text())
    positions = np.load(base / "positions_nm.npy", allow_pickle=False)
    count = system.getNumParticles()
    _require(count == meta["atoms"] == topology.getNumAtoms() and positions.shape == (count, 3), "System/topology/coordinate atom counts differ.")
    _require(np.isfinite(positions).all(), "Built coordinates are nonfinite.")
    masses = np.asarray([system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(count)])
    _require(np.isfinite(masses).all() and (masses > 0).all(), "Invalid particle mass or unsupported massless site.")
    force = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    charges = np.asarray([force.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(count)])
    _require(np.isfinite(charges).all() and np.allclose(charges, meta["charges"], atol=1e-10), "Built partial charges are nonfinite or differ from topology metadata.")
    _require(abs(float(charges.sum())) < 1e-3, "Neutralized solvated system has a nonzero total charge.")
    protein_indices = meta["protein_indices"]
    ligand_indices = meta["ligand_indices"]
    combined = protein_indices + ligand_indices
    _require(len(combined) == len(set(combined)) and all(type(i) is int and 0 <= i < count for i in combined), "Protein/ligand atom partitions overlap or are invalid.")
    report = {"passed": True, "atoms": count, "finite_positions_and_masses": True, "net_system_charge_e": float(charges.sum()),
              "minimum_mass_dalton": float(masses.min()), "forcefield_preset": select_forcefield(config)[0]}
    if config.get("ligand") is not None:
        from rdkit import Chem
        molecules = list(Chem.SDMolSupplier(str(paths["md"] / "prepared/ligand.sdf"), removeHs=False))
        _require(len(molecules) == 1 and molecules[0] is not None, "Prepared SDF does not define one valid ligand.")
        formal = sum(atom.GetFormalCharge() for atom in molecules[0].GetAtoms())
        ligand_q = charges[ligand_indices]
        _require(abs(float(ligand_q.sum()) - formal) < 1e-4, "Parametrized ligand charge differs from the prepared SDF formal charge.")
        _require(charge_provenance is not None, "Explicit AmberTools charge provenance is missing.")
        expected = np.asarray(charge_provenance["partial_charges_e"])
        _require(expected.shape == ligand_q.shape and np.allclose(expected, ligand_q, atol=1e-5, rtol=1e-6),
                 "Applied ligand atom charges differ from the explicitly calculated AM1-BCC vector; inspect template mapping or stale parameter cache.")
        report.update(ligand_formal_charge_e=formal, ligand_partial_charge_sum_e=float(ligand_q.sum()), charge_provenance=charge_provenance)
    atomic_json(base / "build_audit_v2.json", report)
    return report

def run(config, paths):
    validate_config(config)
    from ..md.openmm_engine import run_md
    return run_md(config, paths)


def analyze(config, paths):
    validate_config(config)
    from ..analysis.trajectory import analyze as geometric_analysis
    from ..analysis.mmgbsa import analyze_mmgbsa
    geometric_analysis(config, paths)
    analyze_mmgbsa(config, paths)
    # The unchanged old analyzer emits some empty tables for disabled features;
    # retain them as data, but explicitly record status and never require them.
    for replica in range(1, config["replicates"] + 1):
        out = paths["analysis"] / f"replica_{replica:02d}"
        atomic_json(out / "analysis_status_v2.json", {name: ("ENABLED" if value else "DISABLED") for name, value in config["analysis"].items() if type(value) is bool})


def required_analysis_files(config):
    files = ["timeseries.csv", "summary.json", "quality.json", "optional_analysis.json", "MMGBSA.csv"]
    a = config["analysis"]
    if a["rmsd"]:
        files += ["RMSD.png"]
        if config.get("ligand") is not None:
            files += ["Ligand_RMSD.png"]
    if a["rmsf"]:
        files += ["rmsf.csv", "RMSF.png"]
    if a["hbond"]:
        files += ["Hbond.csv", "Hbond.png"]
    if a["contacts"]:
        files += ["contact_map.csv", "Contact_map.png", "contacts.png"]
    if a["distance"]:
        files += ["distance.csv", "distance.png"]
    if a["sasa"]:
        files += ["SASA.csv"]
    if a.get("glue_orientation_atom_indices") is not None:
        files += ["glue_orientation.png"]
    return files


def _close(actual, expected, message, tolerance=1e-4):
    _require(math.isclose(float(actual), float(expected), rel_tol=1e-7, abs_tol=tolerance), message)


def audit_run(config, paths, *, reload_checkpoint=True):
    """Check every segment and physically reload final checkpoints, without stepping.

    A caller may request a file-only audit, but its returned passed flag is false
    and its verification level explicitly says FILE_ONLY_NOT_ACCEPTED.
    """
    validate_config(config)
    import numpy as np
    from MDAnalysis.coordinates.DCD import DCDReader
    total = phase_steps(config)["production"]
    timestep = duration_ps(config["MD"]["timestep"])
    interval_ps = duration_ps(config["MD"]["report_interval"])
    interval_steps = steps(config["MD"]["report_interval"], config["MD"]["timestep"])
    base = paths["md"] / "system"
    meta = read_json(base / "topology.json")
    system_hash = digest(base / "system.xml")
    replicas = []
    artifact_hashes = {}
    for replica in range(1, config["replicates"] + 1):
        source = paths["md"] / f"replica_{replica:02d}"
        out = paths["analysis"] / f"replica_{replica:02d}"
        progress = read_json(source / "progress.json")
        _require(progress["phase"] == "complete" and progress["production_steps"] == total and progress["current_step"] == total, "Final production step audit failed.")
        _require(progress["system_sha256"] == system_hash, "Replica system fingerprint changed.")
        _require(sum(s["steps"] for s in progress["segments"]) == total, "Committed segment total differs from configured duration.")
        _require(digest(source / progress["checkpoint"]) == progress["checkpoint_sha256"], "Checkpoint hash mismatch.")
        completed = frames = 0
        details = []
        for segment in progress["segments"]:
            count = segment["steps"]
            _require(type(count) is int and count > 0 and count % interval_steps == 0, "Segment is not a whole number of reporting intervals.")
            expected_frames = count // interval_steps
            for kind in ("dcd", "csv"):
                filename = segment[kind]
                _require(Path(filename).name == filename, "Segment path escapes its replica directory.")
                _require(digest(source / filename) == segment[kind + "_sha256"], "Committed segment hash mismatch: " + filename)
                artifact_hashes[f"md/replica_{replica:02d}/{filename}"] = segment[kind + "_sha256"]
            with (source / segment["csv"]).open() as handle:
                reader = csv.reader(handle)
                header = next(reader)
                rows = list(reader)
            _require(len(rows) == expected_frames, "CSV frame count differs from segment steps.")
            values = np.asarray(rows, dtype=float)
            _require(values.ndim == 2 and np.isfinite(values).all(), "Nonfinite or malformed CSV data.")
            _require(header[0].strip('#\"') == "Step" and header[1] == "Time (ps)", "Unexpected state CSV columns.")
            expected_steps = completed + np.arange(1, expected_frames + 1) * interval_steps
            _require(np.array_equal(values[:, 0], expected_steps), "CSV production steps are discontinuous.")
            _require(np.allclose(values[:, 1], expected_steps * timestep, rtol=1e-7, atol=1e-4), "CSV time does not match its step.")
            with DCDReader(str(source / segment["dcd"])) as trajectory:
                _require(len(trajectory) == expected_frames and trajectory.n_atoms == meta["atoms"], "DCD atom/frame count mismatch.")
                _close(trajectory.dt, interval_ps, "DCD reporting interval mismatch.")
                first, last = trajectory[0].time, trajectory[-1].time
                _close(first, interval_ps, "DCD first segment-relative time mismatch.")
                _close(last, count * timestep, "DCD last segment-relative time mismatch.")
                for frame in trajectory:
                    _require(np.isfinite(frame.positions).all(), "Nonfinite DCD coordinates.")
                    _require(frame.dimensions is not None and np.isfinite(frame.dimensions).all() and (frame.dimensions[:3] > 0).all(), "Invalid periodic DCD box.")
            details.append({"dcd": segment["dcd"], "frames": expected_frames, "preceding_steps": completed,
                            "segment_relative_last_ps": float(last), "absolute_last_ps": completed * timestep + float(last)})
            completed += count
            frames += expected_frames
        _require(bool(details), "No committed production frames.")
        _close(details[-1]["absolute_last_ps"], total * timestep, "DCD final absolute time mismatch.")
        checkpoint_ps = None
        if reload_checkpoint:
            checkpoint_ps = _reload_final_checkpoint(config, base, source, progress, total)
        quality = read_json(out / "quality.json")
        _require(quality["finite_observables"] is True and quality["frames_expected"] == quality["frames_observed"] == frames,
                 "Analysis frame/finite-value audit failed.")
        for name in required_analysis_files(config):
            path = out / name
            _require(path.is_file() and path.stat().st_size > 0, "Enabled analysis output missing: " + name)
            artifact_hashes[f"analysis/replica_{replica:02d}/{name}"] = digest(path)
        with (out / "MMGBSA.csv").open() as handle:
            result = next(csv.DictReader(handle))
        expected_status = "NOT_APPLICABLE" if config.get("ligand") is None else ("COMPUTED" if config["mmgbsa"]["enabled"] else "DISABLED")
        _require(result.get("status") == expected_status, "MM/GBSA result status differs from requested capability.")
        if expected_status == "COMPUTED":
            for key in ("mean_kcal_mol", "std_kcal_mol", "sem_kcal_mol"):
                _require(math.isfinite(float(result[key])), "Nonfinite MM/GBSA estimate.")
            _require(int(result["frames"]) >= 2, "Too few MM/GBSA frames.")
        for key in ("progress.json", progress["checkpoint"]):
            artifact_hashes[f"md/replica_{replica:02d}/{key}"] = digest(source / key)
        replicas.append({"replica": replica, "steps": total, "frames": frames, "binary_checkpoint_reloaded": reload_checkpoint,
                         "checkpoint_time_ps": checkpoint_ps, "segments": details, "mmgbsa": expected_status})
    return {"schema_version": 2, "audited": now(), "passed": bool(reload_checkpoint),
            "verification_level": "NATIVE_CHECKPOINT_AND_TRAJECTORY" if reload_checkpoint else "FILE_ONLY_NOT_ACCEPTED",
            "production_ns": total * timestep / 1000, "replicas": replicas, "artifact_sha256": artifact_hashes,
            "system_sha256": system_hash,
            "dcd_time_definition": "Absolute time equals each segment-relative DCD time plus verified preceding committed steps times dt; CSV is global.",
            "claims": "Numerical workflow completion only; no convergence, affinity, pose, recruitment, or degradation inference."}


def _reload_final_checkpoint(config, base, source, progress, expected_steps):
    import numpy as np
    import openmm as mm
    from openmm import app, unit
    md = config["MD"]
    actual_platform = progress["platform"]
    _require(actual_platform in ("CPU", "CUDA", "OpenCL", "Reference"), "Checkpoint platform is not recognized.")
    system = mm.XmlSerializer.deserialize((base / "system.xml").read_text())
    barostat = next(f for f in system.getForces() if isinstance(f, mm.MonteCarloBarostat))
    barostat.setRandomNumberSeed(progress["seed"] + 100000)
    barostat.setFrequency(25)
    dt_ps = duration_ps(md["timestep"])
    integrator = mm.LangevinMiddleIntegrator(md["temperature"] * unit.kelvin, md["friction_per_ps"] / unit.picosecond, dt_ps * unit.picosecond)
    integrator.setRandomNumberSeed(progress["seed"])
    integrator.setConstraintTolerance(1e-6)
    properties = {"Threads": str(md["cpu_threads"])} if actual_platform == "CPU" else ({"Precision": md["precision"]} if actual_platform in ("CUDA", "OpenCL") else {})
    topology = app.PDBFile(str(base / "topology.pdb")).topology
    simulation = app.Simulation(topology, system, integrator, mm.Platform.getPlatformByName(actual_platform), properties)
    try:
        simulation.loadCheckpoint(str(source / progress["checkpoint"]))
        loaded = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
        _require(simulation.currentStep == expected_steps, "Actual binary checkpoint step mismatch.")
        ps = loaded.getTime().value_in_unit(unit.picosecond)
        _close(ps, expected_steps * dt_ps, "Actual binary checkpoint time mismatch.")
        _require(np.isfinite(loaded.getPositions(asNumpy=True).value_in_unit(unit.nanometer)).all(), "Nonfinite checkpoint positions.")
        _require(np.isfinite(loaded.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond)).all(), "Nonfinite checkpoint velocities.")
        _require(math.isfinite(loaded.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)), "Nonfinite checkpoint potential energy.")
        return float(ps)
    except Exception as exc:
        if isinstance(exc, WorkflowError):
            raise
        raise WorkflowError("Binary checkpoint cannot be reloaded in this runtime. Do not claim acceptance or silently regenerate velocities; use the original compatible environment or an explicitly audited state migration.") from exc
    finally:
        del simulation, integrator, system
