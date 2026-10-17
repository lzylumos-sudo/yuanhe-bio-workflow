"""Export complete unconstrained AMBER parameters before restraints are added."""
from pathlib import Path
import numpy as np
from ..utils.core import WorkflowError, atomic_json, digest

def export_amber(forcefield,model,protein_count,ligand_count,folder):
    import openmm as mm
    from openmm import app,unit
    import parmed as pmd
    from parmed.tools.actions import changeRadii
    dry=app.Modeller(model.topology,model.positions)
    # A separate unconstrained system retains all bond/angle parameters needed by Amber.
    system=forcefield.createSystem(dry.topology,nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False,removeCMMotion=False)
    if any(not isinstance(f,(mm.HarmonicBondForce,mm.HarmonicAngleForce,mm.PeriodicTorsionForce,mm.NonbondedForce)) for f in system.getForces()):
        raise WorkflowError('Unsupported AMBER export force; MM/GBSA cannot silently drop custom forces.')
    structure=pmd.openmm.load_topology(dry.topology,system,xyz=dry.positions)
    structure.box=None
    for atom in structure.atoms:
        if atom.atomic_number<=0: raise WorkflowError('Virtual particles are not supported by this AMBER export.')
    amber=pmd.amber.AmberParm.from_structure(structure)
    changeRadii(amber,'mbondi2').execute()
    amber.save(str(folder/'solute.prmtop'),overwrite=True)
    amber.save(str(folder/'solute.inpcrd'),overwrite=True)
    if ligand_count:
        receptor=pmd.amber.AmberParm.from_structure(structure[:protein_count])
        ligand=pmd.amber.AmberParm.from_structure(structure[protein_count:protein_count+ligand_count])
        for name,obj in [('receptor',receptor),('ligand',ligand)]:
            changeRadii(obj,'mbondi2').execute(); obj.save(str(folder/(name+'.prmtop')),overwrite=True)
    atomic_json(folder/'amber_export.json',{'status':'complete','solute_atoms':len(structure.atoms),
        'protein_atoms':protein_count,'ligand_atoms':ligand_count,'constraints_for_export':False,
        'restraints_in_export':False,'radii':'mbondi2','prmtop_sha256':digest(folder/'solute.prmtop'),
        'partition':'all protein atoms vs small molecule' if ligand_count else 'not configured; no small molecule',
        'note':'Same force-field templates as MD. Export uses an unconstrained isolated solute to preserve bonded terms.'})
