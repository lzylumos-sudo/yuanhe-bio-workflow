"""All-atom membrane building, GROMACS execution/restart, and trajectory analysis.

OpenMM builds the solvated membrane; GAFF2/HF-RESP parameterises the ligand;
ParmEd exports the same System to GROMACS. No surrogate energies or trajectories.
"""
from __future__ import annotations
import json,math,os,re,subprocess,time,uuid,hashlib
from pathlib import Path
from . import artifacts as w
from .artifacts import within

AMINO={'ALA','ARG','ASN','ASP','ASH','CYS','CYX','GLN','GLU','GLH','GLY','HIS','HID','HIE','HIP','ILE','LEU','LYS','LYN','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
def number(p,key,default,lo,hi):
 v=float(p.get(key,default))
 if not math.isfinite(v) or not lo<=v<=hi:raise ValueError(f'{key} must be within {lo}–{hi}')
 return v

def rotation(normal):
 import numpy as np
 n=np.asarray(normal,dtype=float)
 if n.shape!=(3,) or not np.isfinite(n).all() or np.linalg.norm(n)<1e-8:raise ValueError('Invalid membrane normal')
 n=n/np.linalg.norm(n);z=np.array([0.,0.,1.]);v=np.cross(n,z);c=float(n@z)
 if c>1-1e-10:return np.eye(3)
 if c<-1+1e-10:return np.diag([1.,-1.,-1.])
 k=np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
 return np.eye(3)+k+k@k/(1+c)

def system_dir(d,ident):
 if not re.fullmatch(r'[A-Za-z0-9_-]{1,100}',str(ident)):raise ValueError('Invalid system identifier')
 return d/'md_systems'/ident

def execute(argv,cwd,timeout=604800,input_text=None):
 print(json.dumps({'argv':list(map(str,argv)),'cwd':str(cwd)},ensure_ascii=False),flush=True)
 p=subprocess.run(list(map(str,argv)),cwd=cwd,input=input_text,text=True,timeout=timeout,pass_fds=w.lease_fds())
 if p.returncode:raise RuntimeError(f'External program exited {p.returncode}; inspect the log')

def add_restraints(top):
 """Restrain solute heavy atoms using molecule-local indices, never global indices."""
 lines=top.read_text().splitlines();result=[];section='';rest=[]
 def flush():
  if rest:
   result.extend(['#ifdef POSRES','[ position_restraints ]','; atom type fx fy fz']+[f'{i} 1 POSRES_FC POSRES_FC POSRES_FC' for i in rest]+['#endif',''])
   rest.clear()
 for line in lines:
  match=re.match(r'\s*\[\s*([^]]+)\s*\]',line)
  if match:
   nxt=match.group(1).strip()
   if nxt in ('moleculetype','system','molecules'):flush()
   section=nxt
  if section=='atoms' and line.strip() and not line.lstrip().startswith((';','#','[')):
   fields=line.split(';')[0].split()
   if len(fields)>=8 and fields[3] in AMINO|{'LIG'} and float(fields[7])>2:rest.append(int(fields[0]))
  result.append(line)
 flush();top.write_text('\n'.join(result)+'\n')

def normalize_lipid_rounding(topology,system,lipid):
 """Correct only tiny residue charge rounding in nominally neutral lipids."""
 from openmm import unit
 force=next(f for f in system.getForces() if type(f).__name__=='NonbondedForce')
 original_charges=[float(force.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)) for i in range(force.getNumParticles())]
 edits=[];lipid_count=0
 for residue in topology.residues():
  if residue.name not in {lipid,lipid[:3]}:continue
  lipid_count+=1
  atoms=list(residue.atoms());charge=sum(float(force.getParticleParameters(a.index)[0].value_in_unit(unit.elementary_charge)) for a in atoms)
  if abs(charge)>1e-5:raise ValueError('Lipid charge differs materially from neutral; automatic rounding correction refused')
  if abs(charge)<1e-12:continue
  atom=next((a for a in atoms if a.name=='P'),None)
  if atom is None:raise ValueError('Lipid phosphate atom missing')
  q,sigma,epsilon=force.getParticleParameters(atom.index)
  force.setParticleParameters(atom.index,q-charge*unit.elementary_charge,sigma,epsilon)
  edits.append({'residue_index':residue.index,'residue_name':residue.name,'atom_index':atom.index,'atom_name':atom.name,'original_residue_charge_e':charge,'delta_charge_e':-charge})
 if not lipid_count:raise ValueError('No lipid residues matched the selected membrane model')
 # Coulomb exceptions store charge products, not references to particle
 # charges. Update those products to preserve their original 1-4 scaling.
 updated={e['atom_index'] for e in edits};exception_count=0
 for i in range(force.getNumExceptions()):
  a,b,qprod,sigma,epsilon=force.getExceptionParameters(i)
  if a not in updated and b not in updated:continue
  old=original_charges[a]*original_charges[b]
  new=float(force.getParticleParameters(a)[0].value_in_unit(unit.elementary_charge))*float(force.getParticleParameters(b)[0].value_in_unit(unit.elementary_charge))
  if old==0:
   if new!=0:raise ValueError('Cannot infer zero-charge exception scaling')
   continue
  force.setExceptionParameters(i,a,b,qprod*(new/old),sigma,epsilon);exception_count+=1
 return {'updated_coulomb_exceptions':exception_count,'lipid_residue_count':lipid_count,'policy':'Neutral-lipid rounding only, absolute residue charge <=1e-5 e; residual assigned to phosphorus; other parameters unchanged','edits':edits,'total_delta_e':sum(e['delta_charge_e'] for e in edits)}

def build(d,c,p,out):
 """Thin UI adapter: verify stage artifacts, then use the standalone engine."""
 from webui import workflow
 workflow.require(d,'screening');workflow.require(d,'docking')
 from .validation import effective_candidates
 rows=effective_candidates(d)
 row=next((r for r in rows if r['compound_id']==p.get('compound_id')),None)
 if not row or not row.get('sdf'):raise ValueError('Candidate has no traceable screening pose')
 row=dict(row)
 matches=[x for x in w.read(d/'docking/input_molecules.json',[]) if x.get('smiles')==row.get('smiles')]
 if len(matches)!=1:raise ValueError('Screening candidate must map to exactly one upstream molecule')
 row['molecule_id']=matches[0]['id']
 from .md_cli import export_bundle
 receptor=d/'docking/receptor.pdb'
 if c.get('research'):
  proof=row.get('research_receptor',{});receptor=Path(proof.get('path',''))
  if not receptor.is_file() or w.digest(receptor)!=proof.get('sha256'):raise ValueError('Selected pose receptor provenance changed')
 bundle=export_bundle(d/'md_handoffs'/('H_'+uuid.uuid4().hex[:12]),row,receptor,within(d,row['sdf']),p,{'gmx':c['gmx'],'seed':c['seed'],'timeout_seconds':c.get('timeout_seconds',86400)}, {'run_id':d.name,'source_revision':w.read(d/'active_source_revision.json'),'screening_stage_sha256':w.digest(d/'screening/stage.json'),'docking_stage_sha256':w.digest(d/'docking/stage.json')})
 return build_prepared(d,c,{**p,'handoff_bundle':str(bundle.relative_to(d))},out,row,within(d,row['sdf']),receptor,w.read(d/'active_source_revision.json'))

def build_prepared(d,c,p,out,row,ligand_path,receptor,source_revision=None):
 """Build from explicit screened pose/receptor files; no UI service required."""
 import numpy as np
 from rdkit import Chem
 from openff.toolkit import Molecule
 from openff.units import unit as offunit
 from openmm import unit,Platform,XmlSerializer
 from openmm.app import PDBFile,Modeller,ForceField,PME,HBonds
 from openmmforcefields.generators import GAFFTemplateGenerator
 import parmed
 if row.get('qualification')!='formal' and not p.get('exploratory_ack'):raise ValueError('Exploratory acknowledgement required')
 if not p.get('orientation_reviewed') or not p.get('orientation_source','').strip():raise ValueError('Review membrane orientation and provide its source or selection rationale')
 center=np.asarray(p.get('center'),dtype=float)
 if center.shape!=(3,) or not np.isfinite(center).all():raise ValueError('Invalid membrane center')
 rot=rotation(p.get('normal'))
 lipid=p.get('lipid','POPE')
 if lipid not in {'POPE','POPC','DPPC','DOPC'}:raise ValueError('Unsupported lipid')
 salt=number(p,'salt_M',.15,0,1);padding=number(p,'padding_nm',1.2,1,3)
 temp=number(p,'temperature_K',310,273,350)
 production=number(p,'production_ns',500,.01,1000);equil=number(p,'equilibration_ns',2,.1,100)
 output_interval=number(p,'output_interval_ps',10,.002,1000);checkpoint_minutes=number(p,'checkpoint_minutes',2,.1,60)
 qm_threads=p.get('qm_threads',4);velocity_seed=p.get('replica_seed',17)
 if type(qm_threads)!=int or qm_threads<1:raise ValueError('qm_threads must be a positive integer')
 if type(velocity_seed)!=int or not 1<=velocity_seed<=2147483647:raise ValueError('Initial velocity seed must be a positive integer')
 gpu_id=str(p.get('gpu_id','0'))
 if not gpu_id.isdigit():raise ValueError('Choose one nonnegative GPU device index')
 if not str(p.get('alignment_selection','protein and backbone')).strip():raise ValueError('Stable protein alignment selection is required')
 p={**p,'output_interval_ps':output_interval,'checkpoint_minutes':checkpoint_minutes,'qm_threads':qm_threads,'replica_seed':velocity_seed,'gpu_id':gpu_id}
 ident='MD_'+uuid.uuid4().hex[:12];folder=system_dir(d,ident);folder.mkdir(parents=True)
 import shutil
 shutil.copyfile(ligand_path,folder/'ligand.sdf');shutil.copyfile(receptor,folder/'receptor.pdb')
 meta={'source_revision':source_revision,'molecule_id':row.get('molecule_id'),'id':ident,'compound_id':row['compound_id'],'state_id':row.get('state_id'),'pose_id':row.get('pose_id'),'qualification':row.get('qualification'),'qualification_provenance':row.get('qualification_provenance'),'status':'building','source_hashes':{'receptor':w.digest(receptor),'ligand':w.digest(ligand_path)},'parameters':{**p,'temperature_K':temp,'production_ns':production,'equilibration_ns':equil},'force_fields':['amber14-all.xml','amber14/tip3p.xml','gaff-2.11'],'membrane_model':lipid+' single-lipid idealised bilayer','orientation_rotation':rot.tolist(),'orientation_center_A':center.tolist(),'scientific_validation':'requires_review'}
 w.atomic(folder/'system.json',meta);w.atomic(out/'result.json',{'system_id':ident})
 mol=Chem.SDMolSupplier(str(ligand_path),removeHs=False)[0]
 if mol is None:raise ValueError('Ligand SDF cannot be read')
 canonical=Chem.MolToSmiles(Chem.RemoveHs(mol),isomericSmiles=True)
 expected=Chem.MolToSmiles(Chem.MolFromSmiles(row['smiles']),isomericSmiles=True)
 if canonical!=expected:raise ValueError('Screening identity and selected SDF chemistry differ')
 mol=Chem.AddHs(mol,addCoords=True)
 off=Molecule.from_rdkit(mol,allow_undefined_stereo=False)
 off.name='LIG'
 for atom in off.atoms:atom.metadata['residue_name']='LIG'
 import sys
 shared=next((parent/'05_molecular_dynamics' for parent in Path(__file__).resolve().parents if (parent/'05_molecular_dynamics/modules').is_dir()),None)
 if shared is None:raise RuntimeError('Shared 05_molecular_dynamics source is required')
 sys.path.insert(0,str(shared))
 from modules.parameterization.hf_resp import compute
 if p.get('mmpbsa',{}).get('enabled'):
  from modules.analysis.mmpbsa import pb_input
  pb=p['mmpbsa']
  if not pb.get('membrane'):raise ValueError('MurJ requires membrane-aware PB settings')
  start=float(pb.get('start_ns',400));end=float(pb.get('end_ns',500));stride=pb.get('stride',10)
  if not 0<=start<end<=production or type(stride)!=int or stride<1:raise ValueError('PB sampling window must lie within production; stride must be a positive integer')
  pb_input(pb,2,':1',':2')
 print('MD_BUILD: HF/6-31G* ESP and two-stage RESP',flush=True)
 q,resp_record=compute(mol,folder/'resp',p)
 off.partial_charges=q*offunit.elementary_charge
 charges=off.partial_charges.m_as(offunit.elementary_charge)
 if not np.isfinite(charges).all() or abs(float(sum(charges))-Chem.GetFormalCharge(mol))>1e-4:raise ValueError('Ligand partial charges do not match its formal charge')
 meta['ligand_chemistry']={'isomeric_smiles':Chem.MolToSmiles(Chem.RemoveHs(mol),isomericSmiles=True),'formal_charge':Chem.GetFormalCharge(mol),'partial_charge_sum':float(sum(charges)),'charge_method':'HF/6-31G* two-stage RESP','resp':resp_record,'protonation':'as screened; no pH state enumeration','exploratory':row.get('qualification')!='formal'}
 w.atomic(folder/'system.json',meta)
 off.to_file(str(folder/'parameterized_ligand.sdf'),file_format='sdf')
 gaff=GAFFTemplateGenerator(molecules=off,forcefield='gaff-2.11',cache=str(folder/'gaff-cache.json'))
 def checked_user_charges(molecule):
  if molecule.partial_charges is None or not np.allclose(molecule.partial_charges.m_as(offunit.elementary_charge),q,atol=1e-8,rtol=0):raise ValueError('GAFF must retain RESP charges')
  return True
 gaff._molecule_has_user_charges=checked_user_charges
 ff=ForceField('amber14-all.xml','amber14/tip3p.xml');ff.registerTemplateGenerator(gaff.generator)
 pdb=PDBFile(str(receptor));model=Modeller(pdb.topology,pdb.positions)
 coords=(np.asarray(pdb.positions.value_in_unit(unit.angstrom))-center)@rot.T
 model.positions=coords*unit.angstrom
 # Receptor has already been explicitly prepared and reviewed by the UI.
 lt=off.to_topology().to_openmm()
 for res in lt.residues():res.name='LIG'
 lcoords=(off.conformers[0].m_as(offunit.angstrom)-center)@rot.T
 model.add(lt,lcoords*unit.angstrom)
 import random
 random.seed(c['seed']);np.random.seed(c['seed'])
 print('MD_BUILD: adding lipid bilayer, solvent and ions',flush=True)
 # OpenMM's insertion relaxation uses a short minimization and an unseeded
 # 1 fs integrator. Use a reproducible, more conservative 0.5 fs relaxation
 # with the same elapsed relaxation time, and retain this policy in metadata.
 import openmm.app.modeller as modeller_module
 old_integrator=modeller_module.LangevinIntegrator;old_minimizer=modeller_module.LocalEnergyMinimizer
 class InsertionIntegrator(old_integrator):
  def __init__(self,temperature,friction,step):
   super().__init__(temperature,friction,step/2);self.setRandomNumberSeed(int(c['seed']))
  def step(self,n):return super().step(n*2)
 class InsertionMinimizer:
  @staticmethod
  def minimize(context,tolerance,maxIterations):return old_minimizer.minimize(context,tolerance,max(300,maxIterations))
 from modules.md.openmm_engine import platform_for
 platform,properties=platform_for({'platform':p.get('platform','auto'),'require_gpu':True,'precision':p.get('precision','mixed'),'device_index':p.get('gpu_id','0'),'cpu_threads':4})
 for key,value in properties.items():platform.setPropertyDefaultValue(key,value)
 meta['membrane_builder']={'platform':platform.getName(),'properties':properties,'seed':c['seed'],'insertion_timestep_fs':0.5,'initial_minimization_max_iterations':300,'insertion_elapsed_time':'same as OpenMM 8.5.1 default; twice as many integration steps'}
 w.atomic(folder/'system.json',meta)
 try:
  modeller_module.LangevinIntegrator=InsertionIntegrator;modeller_module.LocalEnergyMinimizer=InsertionMinimizer
  model.addMembrane(ff,lipidType=lipid,minimumPadding=padding*unit.nanometer,ionicStrength=salt*unit.molar,membraneCenterZ=0*unit.nanometer,platform=platform)
 finally:
  modeller_module.LangevinIntegrator=old_integrator;modeller_module.LocalEnergyMinimizer=old_minimizer
 print('MD_BUILD: exporting parameterized system',flush=True)
 # Preserve bond/angle terms for ParmEd. GROMACS applies H-bond constraints
 # and SETTLE from the exported topology at runtime.
 system=ff.createSystem(model.topology,nonbondedMethod=PME,nonbondedCutoff=1*unit.nanometer,constraints=None,rigidWater=False)
 force=next(f for f in system.getForces() if type(f).__name__=='NonbondedForce')
 ligand_ids=[a.index for a in model.topology.atoms() if a.residue.name=='LIG']
 actual=[float(force.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)) for i in ligand_ids]
 if not np.allclose(actual,q,atol=1e-6,rtol=0):raise ValueError('Final membrane topology changed RESP charges')
 meta['lipid_charge_rounding']=normalize_lipid_rounding(model.topology,system,lipid)
 w.atomic(folder/'system.json',meta)
 structure=parmed.openmm.load_topology(model.topology,system,xyz=model.positions)
 # Save the complete pre-export model so converter failures retain evidence.
 (folder/'system.xml').write_text(XmlSerializer.serialize(system))
 with (folder/'system.pdb').open('w') as f:PDBFile.writeFile(model.topology,model.positions,f,keepIds=True)
 structure.save(str(folder/'system.prmtop'),overwrite=False)
 structure.save(str(folder/'system.top'),overwrite=False);structure.save(str(folder/'system.gro'),overwrite=False)
 converted=parmed.load_file(str(folder/'system.top'),xyz=str(folder/'system.gro'))
 if not np.allclose([a.charge for a in converted.atoms],[a.charge for a in structure.atoms],atol=1e-6,rtol=0):raise ValueError('Topology conversion changed atom charges')
 if converted.defaults.comb_rule != 2 or abs(converted.defaults.fudgeQQ-1/1.2)>1e-6 or abs(converted.defaults.fudgeLJ-.5)>1e-6:raise ValueError('Incompatible GROMACS nonbonded rules')
 add_restraints(folder/'system.top')
 (folder/'system.xml').write_text(XmlSerializer.serialize(system))
 with (folder/'system.pdb').open('w') as f:PDBFile.writeFile(model.topology,model.positions,f,keepIds=True)
 groups={'System':[],'Protein':[],'Backbone':[],'Ligand':[],'Solute':[]}
 for atom in model.topology.atoms():
  i=atom.index+1;groups['System'].append(i)
  if atom.residue.name in AMINO:
   groups['Protein'].append(i);groups['Solute'].append(i)
   if atom.name in ('N','CA','C'):groups['Backbone'].append(i)
  elif atom.residue.name=='LIG':groups['Ligand'].append(i);groups['Solute'].append(i)
 if not groups['Protein'] or not groups['Ligand']:raise ValueError('Exported protein or ligand atom group is empty')
 with (folder/'index.ndx').open('w') as f:
  for key,values in groups.items():f.write('[ '+key+' ]\n'+'\n'.join(' '.join(map(str,values[i:i+15])) for i in range(0,len(values),15))+'\n')
 import importlib.metadata as im
 versions={}
 for package in ['openmm','parmed','openff-toolkit','openmmforcefields','rdkit','MDAnalysis']:
  try:versions[package]=im.version(package)
  except im.PackageNotFoundError:versions[package]='unknown'
 meta['software_versions']=versions
 from collections import Counter
 meta['residue_counts']=dict(Counter(r.name for r in model.topology.residues()))
 from yuanhe.protocols import membrane as protocol_identity
 meta['protocol']=protocol_identity(meta['parameters'])
 meta['software_identity']=__import__('yuanhe.version',fromlist=['identity']).identity()
 meta['forcefield_parameter_files']={str(x):w.digest(x) for x in folder.rglob('*') if x.is_file() and x.suffix in {'.xml','.itp','.top','.prmtop','.frcmod'}}
 meta['net_charge_e']=sum(float(system.getForce(i).getParticleParameters(j)[0].value_in_unit(unit.elementary_charge)) for i in range(system.getNumForces()) if system.getForce(i).__class__.__name__=='NonbondedForce' for j in range(system.getNumParticles()))
 if abs(meta['net_charge_e'])>1e-3:raise ValueError('Built system is not neutral')
 if meta['residue_counts'].get('LIG')!=1:raise ValueError('Exactly one ligand required')
 meta.update(status='built',atom_count=system.getNumParticles(),groups=groups,files={name:w.digest(folder/name) for name in ['system.top','system.gro','system.prmtop','system.xml','system.pdb','index.ndx','ligand.sdf','receptor.pdb']})
 w.atomic(folder/'system.json',meta)
 return {'system_id':ident,'atom_count':system.getNumParticles(),'note':'Membrane and topology built; review the system before minimization/Equilibration'}

def mdp(phase,p,seed):
 temp=float(p['temperature_K']);ns=float(p['production_ns'] if phase=='production' else p['equilibration_ns'])
 em=phase=='em';nvt=phase=='nvt';restrained=phase!='production'
 steps=5000 if em else max(1,round((.1 if nvt else ns if phase=='production' else ns/2)*500000))
 lines=['integrator = '+('steep' if em else 'md'),f'nsteps = {steps}','dt = 0.002','emtol = 1000','emstep = 0.01','cutoff-scheme = Verlet','nstlist = 20','rlist = 1.0','coulombtype = PME','rcoulomb = 1.0','vdwtype = Cut-off','rvdw = 1.0','DispCorr = EnerPres','pbc = xyz','constraints = h-bonds','constraint-algorithm = lincs','nstenergy = 500','nstlog = 500','nstxout-compressed = '+str(round(float(p.get('output_interval_ps',10))/.002)),'compressed-x-grps = System','nstxout = 0','nstvout = 0','nstfout = 0']
 if restrained:lines+=['define = -DPOSRES -DPOSRES_FC='+('100' if phase=='npt_soft' else '1000'),'refcoord-scaling = com']
 if not em:
  lines+=[f'ld-seed = {seed}','tcoupl = V-rescale','tc-grps = System','tau-t = 1.0',f'ref-t = {temp}','continuation = '+('no' if nvt else 'yes'),'gen-vel = '+('yes' if nvt else 'no')]
  if nvt:lines +=['gen-temp = 50',f'gen-seed = {seed}','annealing = single','annealing-npoints = 2','annealing-time = 0 100',f'annealing-temp = 50 {temp}']
  lines+=['pcoupl = '+('no' if nvt else 'Parrinello-Rahman' if phase=='production' else 'C-rescale')]
  if not nvt:lines+=['pcoupltype = semiisotropic','tau-p = 5.0','ref-p = 1.0 1.0','compressibility = 4.5e-5 4.5e-5']
 return '\n'.join(lines)+'\n'

def check_completion(folder,phase,recipe):
 log=(folder/(phase+'.log')).read_text(errors='replace')
 if phase=='em':
  if not re.search(r'converged to Fmax\s*<',log,re.I):raise RuntimeError('Minimization did not meet Fmax convergence criteria; equilibration is blocked')
 else:
  expected=int(re.search(r'^nsteps\s*=\s*(\d+)',recipe,re.M).group(1))
  steps=re.findall(r'^\s*Step\s+Time\s*\n\s*(\d+)\s+[\d.eE+-]+',log,re.M)
  if not steps or int(steps[-1])<expected:raise RuntimeError('Simulation ended before the planned step count; checkpoint retained for resume')
 if not (folder/(phase+'.gro')).is_file():raise RuntimeError('Final structure missing')
 if phase=='production' and not (folder/'production.xtc').is_file():raise RuntimeError('Production trajectory missing')


def checkpoint_probe(folder,c,meta):
 """Independent 20 ps test after equilibration; production is not interrupted."""
 import shutil
 probe=folder/'checkpoint_acceptance';probe.mkdir(exist_ok=True)
 evidence=probe/'acceptance.json'
 if evidence.exists():
  rec=w.read(evidence)
  if rec.get('status')=='passed' and all(w.digest(probe/n)==v for n,v in rec['files'].items()):return rec
  raise ValueError('Existing checkpoint acceptance evidence is inconsistent')
 gmx=c['gmx'];recipe=mdp('production',meta['parameters'],c['seed'])
 recipe=re.sub(r'^nsteps = \d+', 'nsteps = 10000',recipe,flags=re.M)
 (probe/'probe.mdp').write_text(recipe)
 if not (probe/'probe.tpr').exists():
  execute([gmx,'grompp','-f','probe.mdp','-c','../npt_soft.gro','-t','../npt_soft.cpt','-p','../system.top','-n','../index.ndx','-o','probe.tpr'],probe)
 args=[gmx,'mdrun','-deffnm','probe','-ntmpi','1','-ntomp',str(int(os.environ.get('YUANHE_MD_THREADS','12'))),'-nb','gpu','-pme','gpu','-bonded','gpu','-cpt','1']
 if not (probe/'first_checkpoint.json').exists():
  execute(args+['-nsteps','5000'],probe)
  if not (probe/'probe.cpt').exists():raise ValueError('Checkpoint not produced')
  dump=subprocess.run([gmx,'dump','-cp','probe.cpt'],cwd=probe,capture_output=True,text=True,check=True)
  (probe/'first_checkpoint.txt').write_text(dump.stdout+dump.stderr)
  if not re.search(r'step\s*=\s*5000\b',dump.stdout):raise ValueError('Unexpected first checkpoint step')
  shutil.copyfile(probe/'probe.cpt',probe/'step5000.cpt')
  w.atomic(probe/'first_checkpoint.json',{'step':5000,'sha256':w.digest(probe/'step5000.cpt'),'tpr_sha256':w.digest(probe/'probe.tpr')})
 first=w.read(probe/'first_checkpoint.json')
 if w.digest(probe/'probe.tpr')!=first['tpr_sha256']:raise ValueError('Probe TPR changed')
 execute(args+['-cpi','probe.cpt','-append'],probe)
 check_completion(probe,'probe',recipe)
 log=(probe/'probe.log').read_text(errors='replace')
 if re.search(r'LINCS WARNING|Fatal error',log):raise ValueError('Checkpoint test has constraint warnings')
 rec={'status':'passed','scope':'independent 20 ps checkpoint continuity test; not production or biological validation','first_step':5000,'planned_final_step':10000,'production_trajectory_interrupted':False,'files':{n:w.digest(probe/n) for n in ['probe.tpr','probe.cpt','probe.log','probe.gro','probe.edr','probe.xtc','first_checkpoint.json','step5000.cpt']}}
 w.atomic(evidence,rec);return rec

def run(d,c,p,out):
 folder=system_dir(d,p.get('system_id'));meta=w.read(folder/'system.json')
 if not meta or meta.get('status') not in {'built','completed'}:raise ValueError('System construction has not succeeded')
 if not p.get('system_reviewed'):raise ValueError('Review the membrane system and ligand correspondence')
 for name,sha in meta['files'].items():
  if w.digest(folder/name)!=sha:raise ValueError('System inputs changed; resume blocked')
 seed=meta['parameters'].get('replica_seed',c['seed'])
 if type(seed)!=int or not 1<=seed<=2147483647:raise ValueError('Invalid MD velocity seed')
 c={**c,'seed':seed}
 import shutil
 gmx=shutil.which(c['gmx']) or str(Path(c['gmx']).resolve());progress=w.read(folder/'progress.json',{'phases':{}})
 fingerprint=hashlib.sha256(json.dumps({'files':meta['files'],'parameters':meta['parameters'],'gmx':gmx,'gmx_sha256':w.digest(gmx)},sort_keys=True).encode()).hexdigest()
 if progress.get('fingerprint',fingerprint)!=fingerprint:raise ValueError('Simulation protocol changed; resume blocked')
 progress['fingerprint']=fingerprint
 progress['execution_source_revision']=w.read(d/'active_source_revision.json',meta.get('source_revision'))
 progress['engine_sha256']=w.digest(gmx)
 execute([gmx,'--version'],folder,60)
 phases=['em','nvt','npt','npt_soft','production'];previous='system.gro'
 for phase in phases:
  record=progress['phases'].get(phase,{})
  if record.get('status')=='completed':
   for name,sha in record['files'].items():
    if w.digest(folder/name)!=sha:raise ValueError('Completed simulation artifact validation failed')
   previous=phase+'.gro';continue
  if phase=='production':
   print('MD: validating independent checkpoint continuation before production',flush=True)
   checkpoint_probe(folder,c,meta)
  print('MD_PHASE: '+phase,flush=True)
  recipe=mdp(phase,meta['parameters'],c['seed']);recipe_path=folder/(phase+'.mdp')
  if recipe_path.exists() and recipe_path.read_text()!=recipe:raise ValueError('MDP changed; cannot resume')
  recipe_path.write_text(recipe)
  if not (folder/(phase+'.tpr')).exists():
   argv=[gmx,'grompp','-f',phase+'.mdp','-c',previous,'-r','system.gro','-p','system.top','-n','index.ndx','-o',phase+'.tpr']
   prior=phases[phases.index(phase)-1] if phase!='em' else None
   if phase not in ('em','nvt') and (folder/(prior+'.cpt')).exists():argv+=['-t',prior+'.cpt']
   execute(argv,folder)
  tpr_hash=w.digest(folder/(phase+'.tpr'))
  if record.get('tpr_sha256',tpr_hash)!=tpr_hash:raise ValueError('TPR validation failed')
  progress['phases'][phase]={'status':'running','tpr_sha256':tpr_hash};w.atomic(folder/'progress.json',progress)
  argv=[gmx,'mdrun','-deffnm',phase,'-ntmpi','1','-ntomp',str(int(os.environ.get('YUANHE_MD_THREADS','12'))),'-cpt',str(meta['parameters'].get('checkpoint_minutes',2)),'-nb','gpu','-gpu_id',str(os.environ.get('YUANHE_GPU_INDEX',meta['parameters'].get('gpu_id','0')))]
  if phase!='em':argv+=['-pme','gpu','-bonded','gpu']
  if (folder/(phase+'.cpt')).exists():argv+=['-cpi',phase+'.cpt','-append']
  elif record.get('status')=='running' and (folder/(phase+'.log')).exists():raise ValueError('No checkpoint after interruption; create a new system to preserve partial output')
  execute(argv,folder,int(c.get('timeout_seconds',86400)))
  check_completion(folder,phase,recipe)
  progress['phases'][phase]={'status':'completed','tpr_sha256':tpr_hash,'files':{phase+s:w.digest(folder/(phase+s)) for s in ('.gro','.tpr','.mdp','.log','.cpt','.xtc','.edr') if (folder/(phase+s)).is_file()}}
  w.atomic(folder/'progress.json',progress);previous=phase+'.gro'
 if meta.get('protocol'):
  from yuanhe.configuration import integration_steps
  expected=integration_steps(meta['parameters']['production_ns'],2)
  actual=int(re.search(r'^nsteps\s*=\s*(\d+)',mdp('production',meta['parameters'],c['seed']),re.M).group(1))
  if actual!=expected or meta['protocol']['production_steps']!=expected:raise ValueError('Production steps differ from frozen protocol')
 meta['status']='completed';w.atomic(folder/'system.json',meta)
 return {'system_id':meta['id'],'status':'completed','scientific_validation':'requires_trajectory_review'}

def energy_series(path):
 legends={};rows=[]
 for line in path.read_text().splitlines():
  match=re.match(r'@\s+s(\d+)\s+legend\s+"([^"]+)"',line)
  if match:legends[int(match.group(1))]=match.group(2)
  if not line.strip() or line.lstrip().startswith(('#','@')):continue
  values=[float(v) for v in line.split()]
  if not all(math.isfinite(v) for v in values):raise ValueError('Energy file contains nonfinite values')
  if len(values)!=len(legends)+1:raise ValueError('Energy columns do not match labels')
  rows.append({'time_ns':values[0]/1000,**{legends[i]:values[i+1] for i in sorted(legends)}})
 if not rows or not {'Temperature','Pressure','Potential'}<=set(legends.values()):raise ValueError('Temperature, pressure or potential energy record missing')
 return rows


def analyze(d,c,p,out):
 import numpy as np
 import MDAnalysis as mda
 from MDAnalysis.lib.distances import distance_array,minimize_vectors
 folder=system_dir(d,p.get('system_id'));meta=w.read(folder/'system.json',{})
 if meta.get('status')!='completed':raise ValueError('Production simulation incomplete')
 for name,sha in meta['files'].items():
  if w.digest(folder/name)!=sha:raise ValueError('Analysis input validation failed')
 progress=w.read(folder/'progress.json',{})
 for record in progress.get('phases',{}).values():
  for name,sha in record.get('files',{}).items():
   if w.digest(folder/name)!=sha:raise ValueError('Simulation output changed; analysis refused')
 dest=folder/'analysis';dest.mkdir(exist_ok=True)
 import sys
 shared=next((parent/'05_molecular_dynamics' for parent in Path(__file__).resolve().parents if (parent/'05_molecular_dynamics/modules').is_dir()),None)
 if shared is None:raise RuntimeError('Shared 05_molecular_dynamics source is required')
 sys.path.insert(0,str(shared))
 from modules.analysis.research import analyze as research_analyze
 formal=research_analyze(folder/'system.prmtop',[folder/'production.xtc'],dest/'research',alignment=meta['parameters'].get('alignment_selection','protein and backbone'))
 pb_result=None
 if meta['parameters'].get('mmpbsa',{}).get('enabled'):
  from modules.analysis.mmpbsa import run as pb_run
  pb_settings=meta['parameters']['mmpbsa']
  if not pb_settings.get('membrane'):raise ValueError('Explicit membrane system requires membrane-aware PB parameters')
  pb_result=pb_run(folder/'system.prmtop',[folder/'production.xtc'],dest/'mmpbsa',pb_settings)
 execute([c['gmx'],'trjconv','-s','production.tpr','-f','production.xtc','-o','analysis/whole.xtc','-n','index.ndx','-pbc','mol','-center'],folder,input_text='Protein\nSystem\n')
 # Preserved PDB avoids TPR-version coupling and retains construct residue IDs.
 u=mda.Universe(str(folder/'system.pdb'),str(dest/'whole.xtc'))
 if len(u.atoms)!=meta['atom_count']:raise ValueError('Trajectory/topology atom count mismatch')
 protein=u.atoms[np.array(meta['groups']['Protein'])-1];back=u.atoms[np.array(meta['groups']['Backbone'])-1];ligand=u.atoms[np.array(meta['groups']['Ligand'])-1]
 ligheavy=ligand.select_atoms('not name H*');ca=protein.select_atoms('name CA')
 if len(back)<3 or not len(ligheavy):raise ValueError('Analysis selection is empty')
 u.trajectory[0];reference=back.positions.copy();ref_center=reference.mean(0);q=reference-ref_center
 ref_lig=ligheavy.positions.copy();delta=ref_lig.mean(0)-ref_center;ref_lig+=minimize_vectors(delta[None,:],u.dimensions)[0]-delta
 stride=max(1,len(u.trajectory)//200);series=[];frames=[];aligned=[];contact={};lipid_p=u.select_atoms('resname POPE POPC DPPC DOPC POP DPP DOP and name P')
 frame_ids=sorted(set(range(0,len(u.trajectory),stride))|{len(u.trajectory)-1})
 for frame_id in frame_ids:
  ts=u.trajectory[frame_id]
  b=back.positions.copy();center=b.mean(0);P=b-center;U,S,V=np.linalg.svd(P.T@q);R=U@V
  if np.linalg.det(R)<0:U[:,-1]*=-1;R=U@V
  bx=P@R+ref_center;l=ligheavy.positions.copy();delta=l.mean(0)-center;l+=minimize_vectors(delta[None,:],ts.dimensions)[0]-delta;lx=(l-center)@R+ref_center
  cap=(ca.positions-center)@R+ref_center;aligned.append(cap)
  dist=distance_array(protein.positions,ligheavy.positions,box=ts.dimensions);touched=set(str(a.resid) for a,x in zip(protein,dist.min(axis=1)) if x<4.5)
  for r in touched:contact[r]=contact.get(r,0)+1
  thickness=None
  if len(lipid_p)>4:
   z=lipid_p.positions[:,2];mid=np.median(z);lo=z[z<mid];hi=z[z>=mid]
   if len(lo) and len(hi):thickness=float(np.median(hi)-np.median(lo))
  series.append({'time_ns':float(ts.time/1000),'backbone_rmsd_A':float(np.sqrt(np.mean(np.sum((bx-reference)**2,axis=1)))),'ligand_rmsd_A':float(np.sqrt(np.mean(np.sum((lx-ref_lig)**2,axis=1)))),'minimum_contact_A':float(dist.min()),'membrane_thickness_proxy_A':thickness})
  frames.append({'time_ns':float(ts.time/1000),'positions':np.concatenate([cap,lx]).round(3).tolist()})
 arr=np.asarray(aligned);rmsf=np.sqrt(np.mean(np.sum((arr-arr.mean(0))**2,axis=2),axis=0))
 summary={'molecule_id':meta.get('molecule_id'),'pose_id':meta.get('pose_id'),'qualification':meta['qualification'],'formal_qualification_changed':False,'production_span_ns':series[-1]['time_ns']-series[0]['time_ns'],'production_end_ns':series[-1]['time_ns'],'first_saved_frame_ns':series[0]['time_ns'],'system_id':meta['id'],'compound_id':meta['compound_id'],'sampled_frames':len(series),'trajectory_frames':len(u.trajectory),'alignment':'protein backbone least squares; ligand not independently fitted','rmsf':[{'resid':int(a.resid),'rmsf_A':float(v)} for a,v in zip(ca,rmsf)],'contact_occupancy_4_5A':{r:n/len(series) for r,n in contact.items()},'limitations':['Contact occupancy uses a 4.5 Å distance cutoff, not a hydrogen-bond definition','Membrane thickness is a phosphorus-layer separation proxy','Analysis references the first production frame; compare replicas separately'],'source_hashes':{x:w.digest(folder/x) for x in ['production.tpr','production.xtc','system.pdb']}}
 execute([c['gmx'],'energy','-f','production.edr','-o','analysis/thermodynamics.xvg'],folder,input_text='Temperature\nPressure\nPotential\nTotal-Energy\nDensity\nBox-X\nBox-Y\nBox-Z\n0\n')
 thermo=energy_series(dest/'thermodynamics.xvg')
 summary['thermodynamics']={key:{'mean':float(np.mean([r[key] for r in thermo])),'std':float(np.std([r[key] for r in thermo])),'min':min(r[key] for r in thermo),'max':max(r[key] for r in thermo)} for key in thermo[0] if key!='time_ns'}
 summary['limitations'].append('Review temperature, pressure and energy time series; finite values and means do not establish equilibration or convergence')
 summary['units']={'time':'ns','RMSD_RMSF_distance':'angstrom','Temperature':'K','Pressure':'bar','Potential':'kJ/mol','Density':'kg/m3','Box':'nm'}
 summary['checkpoint_acceptance']=w.read(folder/'checkpoint_acceptance/acceptance.json',{})
 summary['quality_observations']={'finite_trajectory':bool(np.isfinite(np.asarray([f['positions'] for f in frames])).all()),'expected_duration_ns':meta['parameters']['production_ns'],'actual_duration_ns':series[-1]['time_ns'],'constraint_warnings':bool(re.search(r'LINCS WARNING|Fatal error',(folder/'production.log').read_text(errors='replace'))),'last_backbone_rmsd_A':series[-1]['backbone_rmsd_A'],'last_ligand_rmsd_A':series[-1]['ligand_rmsd_A'],'interpretation':'Descriptive engineering checks; no activity or binding validation is inferred'}
 import csv
 for name,rows in [('trajectory_metrics',series),('thermodynamics',thermo),('residue_rmsf',summary['rmsf'])]:
  with (dest/(name+'.csv')).open('w',newline='') as f:
   writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
 for name,index in [('first',0),('last',len(u.trajectory)-1)]:
  u.trajectory[index]
  with mda.Writer(str(dest/(name+'_production.pdb')),n_atoms=len(u.atoms)) as writer:writer.write(u.atoms)
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,axes=plt.subplots(3,2,figsize=(11,10),layout='constrained')
 for key,label in [('backbone_rmsd_A','Backbone'),('ligand_rmsd_A','Ligand')]:axes[0,0].plot([r['time_ns'] for r in series],[r[key] for r in series],label=label)
 axes[0,0].set(xlabel='Time (ns)',ylabel='RMSD (angstrom)');axes[0,0].legend()
 axes[0,1].plot([r['resid'] for r in summary['rmsf']],[r['rmsf_A'] for r in summary['rmsf']]);axes[0,1].set(xlabel='Residue',ylabel='CA RMSF (angstrom)')
 for ax,key,label in [(axes[1,0],'Temperature','Temperature (K)'),(axes[1,1],'Pressure','Pressure (bar)'),(axes[2,0],'Potential','Potential energy (kJ/mol)')]:
  ax.plot([r['time_ns'] for r in thermo],[r[key] for r in thermo],lw=.8);ax.set(xlabel='Time (ns)',ylabel=label)
 axes[2,1].plot([r['time_ns'] for r in series],[r['minimum_contact_A'] for r in series]);axes[2,1].set(xlabel='Time (ns)',ylabel='Minimum protein-ligand distance (angstrom)')
 fig.suptitle(meta['id']+' | Exploratory MD; qualification unchanged')
 fig.savefig(dest/'analysis_panels.png',dpi=180);fig.savefig(dest/'analysis_panels.pdf');plt.close(fig)
 w.atomic(dest/'thermodynamics.json',thermo)
 summary['last100ns']=formal
 summary['mmpbsa']=pb_result
 summary['preview_alignment']='Backbone fit for interactive preview; core statistics use last100ns.alignment_selection'
 w.atomic(dest/'summary.json',summary);w.atomic(dest/'series.json',series);w.atomic(dest/'frames.json',{'atoms':[{'chain':'REC','name':'CA','resid':int(a.resid)} for a in ca]+[{'chain':'LIG','name':a.name} for a in ligheavy],'frames':frames})
 return summary
