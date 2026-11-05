"""Isolated workers started only by the local UI task manager."""
from __future__ import annotations
import argparse,json,math,os,subprocess,sys,time,traceback
from pathlib import Path
from webui import workflow as w
from webui.studio import within

def receptor(d,c,p,out):
 from pdbfixer import PDBFixer
 from openmm import Platform,unit,Context,VerletIntegrator,LocalEnergyMinimizer
 from openmm.app import PDBFile,ForceField,NoCutoff
 w.require(d,'prepare');src=d/'prepare/receptor_full.pdb'
 ph=float(p.get('pH',7.4))
 if not math.isfinite(ph) or not 0<=ph<=14:raise ValueError('pH Invalid')
 # Use the same double-precision placement backend as the membrane workflow.
 fixer=PDBFixer(filename=str(src),platform=Platform.getPlatformByName('Reference'))
 fixer.findMissingResidues();fixer.findNonstandardResidues();fixer.findMissingAtoms()
 report={'pH':ph,'input_sha256':w.digest(src),'missing_residues':str(fixer.missingResidues),'nonstandard_residues':str(fixer.nonstandardResidues),'missing_heavy_atoms':{r.id:[x.name for x in a] for r,a in fixer.missingAtoms.items()},'missing_terminal_atoms':{r.id:a for r,a in fixer.missingTerminals.items()},'mutations':'preserved','scientific_review_required':True}
 w.atomic(out/'result.json',report)
 if fixer.missingResidues or fixer.nonstandardResidues:raise ValueError('Missing or nonstandard residues require review; no automatic rebuilding or deletion')
 if fixer.missingAtoms and not p.get('repair_sidechains'):raise ValueError('Missing side-chain heavy atoms; review the report and explicitly allow completion')
 before={(a.residue.chain.id,a.residue.id,a.name):v.value_in_unit(unit.nanometer) for a,v in zip(fixer.topology.atoms(),fixer.positions) if a.element and a.element.symbol!='H'}
 fixer.addMissingAtoms(seed=c['seed'])
 forcefield=ForceField('amber14-all.xml')
 fixer.addMissingHydrogens(ph,forcefield=forcefield)
 # PDBFixer's soft placement can leave intramolecular side-chain contacts
 # that distance-based receptor readers interpret as spurious covalent bonds.
 # Relax only added atoms with the same physical protein force field used by MD.
 system=forcefield.createSystem(fixer.topology,nonbondedMethod=NoCutoff,constraints=None)
 movable=0
 for atom in fixer.topology.atoms():
  key=(atom.residue.chain.id,atom.residue.id,atom.name)
  if key in before:system.setParticleMass(atom.index,0)
  else:movable+=1
 integrator=VerletIntegrator(.001*unit.picoseconds)
 context=Context(system,integrator,Platform.getPlatformByName('Reference'))
 context.setPositions(fixer.positions)
 LocalEnergyMinimizer.minimize(context,tolerance=10*unit.kilojoule_per_mole/unit.nanometer,maxIterations=2000)
 state=context.getState(getPositions=True,getEnergy=True)
 if not math.isfinite(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)):
  raise ValueError('Nonfinite receptor completion energy; inspect the retained preparation report')
 fixer.positions=state.getPositions()
 del context,integrator
 report['added_atom_relaxation']={'forcefield':'amber14-all.xml','fixed_atoms':'all existing heavy atoms','movable_atoms':movable,'max_iterations':2000,'tolerance_kJ_mol_nm':10,'purpose':'repair new side-chain and hydrogen geometry; not whole-receptor equilibration'}
 pos=fixer.positions.value_in_unit(unit.nanometer)
 for a in fixer.topology.atoms():
  key=(a.residue.chain.id,a.residue.id,a.name)
  if key in before:pos[a.index]=before[key]
 target=out/'receptor.pdb'
 with target.open('w') as f:PDBFile.writeFile(fixer.topology,pos*unit.nanometer,f,keepIds=True)
 report.update(receptor=str(target.relative_to(d)),sha256=w.digest(target),placement_platform='Reference',note='Construct and existing heavy-atom coordinates preserved; new side chains and termini require review; no WT reversion to wild type')
 return report

def environment(d,c,p,out):
 from webui.readiness import collect
 result=collect(c,w.ROOT,publish=lambda value:w.atomic(out/'result.json',value))
 if not result['all_runtime_checks_pass']:raise RuntimeError('Deployment readiness blocked; inspect result.json')
 return result


def main():
 parser=argparse.ArgumentParser();parser.add_argument('--run-id',required=True);parser.add_argument('--operation',required=True);a=parser.parse_args()
 if not a.operation or any(x in a.operation for x in ('/','\\','..')):raise ValueError('Invalid operation identifier')
 d=w.run_dir(a.run_id);out=d/'operations'/a.operation;req=w.read(out/'request.json');rec=w.read(out/'operation.json')
 if rec['status']!='running':return 1
 c=req['config'];p=req['payload'];kind=req['kind']
 os.environ['PATH']=str(Path(sys.executable).parent)+os.pathsep+os.environ.get('PATH','')
 lease_context=None
 try:
  if c.get('research') and kind=='md_build':
   from research.resources import detect,gpu_lease,disk_guard
   disk_guard(d,c['research'])
   # Desktop redraws can briefly cross the idle threshold during the UI click.
   # Recheck for a bounded interval without relaxing memory/utilization limits
   # or preempting another process. A busy device is never silently selected.
   deadline=time.monotonic()+10
   while True:
    snapshot=detect(d)
    available=[g for g in snapshot['gpus'] if g['allocatable']]
    selected=None
    for gpu in available:
     try:
      lease_context=gpu_lease(gpu);lease_context.__enter__()
      selected=gpu
      break
     except BlockingIOError:
      lease_context=None
    if selected is not None:break
    if time.monotonic()>=deadline:
     detail=snapshot.get('gpu_detection_error') or [{k:g[k] for k in ('index','free_mib','total_mib','utilization')} for g in snapshot['gpus']]
     raise ValueError(f'No idle GPU after 10-second availability check; task was not queued. Retry from the UI when resources are free. Observed: {detail}')
    print('MD_BUILD: checking GPU availability; existing workloads are preserved',flush=True)
    time.sleep(1)
   os.environ['CUDA_VISIBLE_DEVICES']=selected['uuid'];os.environ['YUANHE_GPU_INDEX']='0'
   p={**p,'gpu_id':'0'}
  if kind=='receptor':result=receptor(d,c,p,out)
  elif kind=='environment':result=environment(d,c,p,out)
  elif kind=='md_validation':
   from scientific.validation import run
   result=run(d,c,p,out,w.ROOT)
  elif kind=='install_membrane':
   subprocess.run(['bash',str(w.ROOT/'deployment/install_membrane.sh')],check=True,timeout=7200)
   result={'status':'installed','md_python':'/opt/yuanhe/envs/membrane/bin/python','scientific_validation':'not_run'}
  else:
   from scientific import membrane
   result={'md_build':membrane.build,'md_run':membrane.run,'md_analyze':membrane.analyze}[kind](d,c,p,out)
  w.atomic(out/'result.json',result)
  if w.read(out/'operation.json')['status']=='running':w.atomic(out/'operation.json',{**rec,'status':'completed','ended':time.time()})
  return 0
 except Exception as e:
  traceback.print_exc()
  if kind=='md_build':
   ident=w.read(out/'result.json',{}).get('system_id')
   if ident:
    folder=d/'md_systems'/ident;meta=w.read(folder/'system.json',{})
    w.atomic(folder/'system.json',{**meta,'status':'failed','error':str(e)})
  if w.read(out/'operation.json')['status']=='running':w.atomic(out/'operation.json',{**rec,'status':'failed','error':str(e),'ended':time.time()})
  return 1
 finally:
  if lease_context:lease_context.__exit__(None,None,None)
if __name__=='__main__':raise SystemExit(main())
