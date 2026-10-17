"""Public-test energy roundtrip for the actual OpenMM -> Amber parameter export."""
import json,sys
from pathlib import Path
import numpy as np
import openmm as mm
from openmm import app,unit
from openff.toolkit import Molecule
from rdkit import Chem
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root))
from modules.utils.core import read_json,atomic_json
from modules.md.openmm_engine import ligand_template_generator
run=root/'results'/sys.argv[1]
import yaml
config=yaml.safe_load((run/'config.yaml').read_text());base=run/'md/system'
protein=app.PDBFile(str(run/'md/prepared/protein.pdb'))
model=app.Modeller(protein.topology,protein.positions)
rd=next(iter(Chem.SDMolSupplier(str(run/'md/prepared/ligand.sdf'),removeHs=False)))
ligand=Molecule.from_rdkit(rd,hydrogens_are_explicit=True,allow_undefined_stereo=False)
top=ligand.to_topology().to_openmm()
for residue in top.residues():residue.name='LIG'
model.add(top,ligand.conformers[0].to_openmm())
ff=app.ForceField(config['MD']['protein_forcefield'],config['MD']['water_forcefield'])
gen=ligand_template_generator(ligand,config['MD'],base/'ligand_parameters.json')
ff.registerTemplateGenerator(gen.generator)
original=ff.createSystem(model.topology,nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False,removeCMMotion=False)
amber=app.AmberPrmtopFile(str(base/'solute.prmtop'))
converted=amber.createSystem(nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False,removeCMMotion=False)
energies=[]
for system in [original,converted]:
    integ=mm.VerletIntegrator(0.001)
    ctx=mm.Context(system,integ,mm.Platform.getPlatformByName('CPU'),{'Threads':'4'})
    ctx.setPositions(model.positions)
    energies.append(float(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)))
    del ctx,integ
delta=abs(energies[0]-energies[1])
passed=bool(np.isfinite(energies).all() and delta<0.01)
result={'passed':passed,'platform':'CPU','method':'Identical unconstrained isolated-solute coordinates; original template-generated system versus exported Amber prmtop','original_kj_mol':energies[0],'amber_kj_mol':energies[1],'absolute_difference_kj_mol':delta,'tolerance_kj_mol':0.01}
atomic_json(root/'installation_evidence/amber_export_roundtrip.json',result)
print(json.dumps(result))
if not passed:raise SystemExit(1)
