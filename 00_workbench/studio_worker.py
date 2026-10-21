"""Isolated workers started only by the local UI task manager."""
from __future__ import annotations
import argparse,json,math,os,subprocess,sys,time,traceback
from pathlib import Path
from webui import workflow as w
from webui.studio import within

def receptor(d,c,p,out):
 from pdbfixer import PDBFixer
 from openmm import Platform,unit
 from openmm.app import PDBFile,Modeller
 w.require(d,'prepare');src=d/'prepare/receptor_full.pdb'
 ph=float(p.get('pH',7.4))
 if not math.isfinite(ph) or not 0<=ph<=14:raise ValueError('pH Invalid')
 fixer=PDBFixer(filename=str(src),platform=Platform.getPlatformByName('CPU'))
 fixer.findMissingResidues();fixer.findNonstandardResidues();fixer.findMissingAtoms()
 report={'pH':ph,'input_sha256':w.digest(src),'missing_residues':str(fixer.missingResidues),'nonstandard_residues':str(fixer.nonstandardResidues),'missing_heavy_atoms':{r.id:[x.name for x in a] for r,a in fixer.missingAtoms.items()},'missing_terminal_atoms':{r.id:a for r,a in fixer.missingTerminals.items()},'mutations':'preserved','scientific_review_required':True}
 w.atomic(out/'result.json',report)
 if fixer.missingResidues or fixer.nonstandardResidues:raise ValueError('Missing or nonstandard residues require review; no automatic rebuilding or deletion')
 if fixer.missingAtoms and not p.get('repair_sidechains'):raise ValueError('Missing side-chain heavy atoms; review the report and explicitly allow completion')
 before={(a.residue.chain.id,a.residue.id,a.name):v.value_in_unit(unit.nanometer) for a,v in zip(fixer.topology.atoms(),fixer.positions) if a.element and a.element.symbol!='H'}
 fixer.addMissingAtoms(seed=c['seed']);fixer.addMissingHydrogens(ph)
 pos=fixer.positions.value_in_unit(unit.nanometer)
 for a in fixer.topology.atoms():
  key=(a.residue.chain.id,a.residue.id,a.name)
  if key in before:pos[a.index]=before[key]
 target=out/'receptor.pdb'
 with target.open('w') as f:PDBFile.writeFile(fixer.topology,pos*unit.nanometer,f,keepIds=True)
 report.update(receptor=str(target.relative_to(d)),sha256=w.digest(target),note='Construct and existing heavy-atom coordinates preserved; new side chains and termini require review; no WT reversion to wild type')
 return report

def environment(d,c,p,out):
 probes={
 'chemistry':[c['chem_python'],'-c','import rdkit,yaml,gemmi,numpy; print(rdkit.__version__)'],
 'screening':[c.get('screen_python') or c['chem_python'],'-c','import MDAnalysis,prolif,posebusters; print("OK")'],
 'peptide_cuda':[c['peptide_python'],'-c','import torch; from torch_scatter import scatter_add; from torch_cluster import radius_graph; assert torch.cuda.is_available(); x=torch.tensor([[0.,0.,0.],[1.,0.,0.]],device="cuda"); print(radius_graph(x,2)); print(scatter_add(torch.ones(2,device="cuda"),torch.tensor([0,0],device="cuda"))); torch.cuda.synchronize(); print(torch.cuda.get_device_name())'],
 'reinvent':[c['reinvent_executable'],'--help'],
 'receptor':[c['prep_python'],'-c','import pdbfixer,openmm; print(openmm.__version__)'],
 'vina':[c['vina'],'--version'],
 'gromacs':[c['gmx'],'--version'],
 'membrane':[c['md_python'],'-c','import openmm,parmed,MDAnalysis,psi4,resp; from openff.toolkit import Molecule; from openmmforcefields.generators import GAFFTemplateGenerator; import shutil; assert shutil.which("antechamber") and shutil.which("parmchk2"); print(openmm.__version__)']}
 results={}
 for key,argv in probes.items():
  try:
   env=os.environ.copy()
   env['PATH']=str(Path(argv[0]).parent)+os.pathsep+env.get('PATH','')
   r=subprocess.run(argv,capture_output=True,text=True,timeout=90,env=env);results[key]={'passed':r.returncode==0,'detail':(r.stdout+r.stderr)[-6000:]}
  except Exception as e:results[key]={'passed':False,'detail':str(e)}
  w.atomic(out/'result.json',{'checks':results});print(key,results[key],flush=True)
 return {'checks':results,'runtime_passed':all(x['passed'] for x in results.values()),'scientific_validation':'not_run'}

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--run-id',required=True);parser.add_argument('--operation',required=True);a=parser.parse_args()
 if not a.operation or any(x in a.operation for x in ('/','\\','..')):raise ValueError('Invalid operation identifier')
 d=w.run_dir(a.run_id);out=d/'operations'/a.operation;req=w.read(out/'request.json');rec=w.read(out/'operation.json')
 if rec['status']!='running':return 1
 c=req['config'];p=req['payload'];kind=req['kind']
 os.environ['PATH']=str(Path(sys.executable).parent)+os.pathsep+os.environ.get('PATH','')
 lease_context=None
 try:
  if c.get('research') and kind in {'md_build','environment','md_validation'}:
   from research.resources import detect,gpu_lease,disk_guard
   disk_guard(d,c['research'])
   available=[g for g in detect(d)['gpus'] if g['allocatable']]
   if not available:raise ValueError('No allocatable GPU for operation')
   lease_context=gpu_lease(available[0]);lease_context.__enter__()
   os.environ['CUDA_VISIBLE_DEVICES']=available[0]['uuid'];os.environ['YUANHE_GPU_INDEX']='0'
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
