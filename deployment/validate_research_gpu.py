#!/usr/bin/env python3
"""Small engineering GPU/continuation/PB check using a REAL HF-RESP ligand topology.
No 500 ns or target-system inference. Inputs and outputs stay outside the source tree.
"""
import argparse,json,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'05_molecular_dynamics'))

def main():
 p=argparse.ArgumentParser();p.add_argument('--precision',choices=['single','mixed','double'],default='mixed');p.add_argument('--parameter-dir',required=True,type=Path);p.add_argument('--output',required=True,type=Path);args=p.parse_args()
 import numpy as np,parmed,openmm as mm
 from openmm import app,unit
 from modules.md.openmm_engine import platform_for
 from modules.analysis.research import analyze
 from modules.analysis.mmpbsa import run as pb_run
 out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
 parameter=args.parameter_dir.resolve()
 for name in ['ligand.mol2','ligand.frcmod']:
  import shutil
  shutil.copyfile(parameter/name,out/name)
 ligand=parmed.load_file(str(out/'ligand.mol2'),structure=True);ligand.coordinates+=np.array([8.,0.,0.]);ligand.save(str(out/'positioned.mol2'),overwrite=False)
 (out/'tleap.in').write_text('source leaprc.protein.ff14SB\nsource leaprc.gaff2\nsource leaprc.water.tip3p\nset default PBradii mbondi2\nloadamberparams ligand.frcmod\nlig = loadmol2 positioned.mol2\nprot = sequence { ACE ALA NME }\ncomplex = combine { prot lig }\nsolvatebox complex TIP3PBOX 12.0\nsaveamberparm complex system.prmtop system.inpcrd\nsavepdb complex system.pdb\nquit\n')
 with (out/'tleap.log').open('w') as log:subprocess.run(['tleap','-f','tleap.in'],cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True)
 top=app.AmberPrmtopFile(str(out/'system.prmtop'));xyz=app.AmberInpcrdFile(str(out/'system.inpcrd'))
 system=top.createSystem(nonbondedMethod=app.PME,nonbondedCutoff=.9*unit.nanometer,constraints=app.HBonds)
 settings={'platform':'auto','require_gpu':True,'precision':args.precision,'device_index':'0','cpu_threads':2}
 platform,props=platform_for(settings)
 def new():
  integrator=mm.LangevinMiddleIntegrator(300*unit.kelvin,1/unit.picosecond,.002*unit.picosecond);integrator.setRandomNumberSeed(17)
  return app.Simulation(top.topology,system,integrator,platform,props)
 sim=new();sim.context.setPositions(xyz.positions);sim.context.setPeriodicBoxVectors(*xyz.boxVectors);sim.minimizeEnergy(maxIterations=1000);sim.context.setVelocitiesToTemperature(300*unit.kelvin,17);sim.step(500)
 sim.currentStep=0;sim.context.setTime(0*unit.picosecond)
 sim.reporters=[app.DCDReporter(str(out/'part1.dcd'),100,enforcePeriodicBox=False)]
 sim.step(1000);sim.saveCheckpoint(str(out/'checkpoint.chk'));checkpoint_time=sim.context.getState().getTime().value_in_unit(unit.picosecond)
 sim.reporters[0]._out.close();sim.reporters=[];sim.step(100)
 expected=np.array(sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer));del sim
 resumed=new();resumed.loadCheckpoint(str(out/'checkpoint.chk'));loaded_time=resumed.context.getState().getTime().value_in_unit(unit.picosecond)
 if abs(loaded_time-checkpoint_time)>1e-10:raise ValueError('Checkpoint time reset')
 resumed.step(100);actual=np.array(resumed.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer))
 error=float(np.max(np.abs(actual-expected)))
 if error>1e-4:raise ValueError('Checkpoint trajectory differs from uninterrupted continuation')
 # Return to the original checkpoint before writing the continuation segment.
 resumed.loadCheckpoint(str(out/'checkpoint.chk'));resumed.reporters=[app.DCDReporter(str(out/'part2.dcd'),100,enforcePeriodicBox=False)];resumed.step(1000);resumed.reporters[0]._out.close();resumed.reporters=[]
 devices={k:platform.getPropertyValue(resumed.context,k) for k in platform.getPropertyNames()}
 result={'scope':'engineering only: solvated alanine dipeptide + HF-RESP ethanol','platform':platform.getName(),'properties':devices,'checkpoint_time_ps':checkpoint_time,'loaded_time_ps':loaded_time,'final_time_ps':resumed.context.getState().getTime().value_in_unit(unit.picosecond),'continuation_max_error_nm':error,'atoms':system.getNumParticles()}
 (out/'gpu_acceptance.json').write_text(json.dumps(result,indent=2)+'\n');del resumed
 times=np.arange(1,21)*.0002
 analyze(out/'system.prmtop',[out/'part1.dcd',out/'part2.dcd'],out/'analysis',alignment='resname ACE ALA NME and name CA C N',protein='resname ACE ALA NME',times_ns=times)
 pb=pb_run(out/'system.prmtop',[out/'part1.dcd',out/'part2.dcd'],out/'pb',{'formal':False,'start_ns':.0002,'end_ns':.004,'stride':5,'entropy':False,'membrane':False},times_ns=times,protein='resname ACE ALA NME')
 print(json.dumps(result,indent=2));print('PB:',pb['components'])
if __name__=='__main__':main()
