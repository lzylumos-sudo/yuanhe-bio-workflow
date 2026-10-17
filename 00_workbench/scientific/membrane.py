"""All-atom membrane building, GROMACS execution/restart, and trajectory analysis.

OpenMM builds the solvated membrane; GAFF2/AM1-BCC parameterises the ligand;
ParmEd exports the same System to GROMACS. No surrogate energies or trajectories.
"""
from __future__ import annotations
import json,math,os,re,subprocess,time,uuid,hashlib
from pathlib import Path
from webui import workflow as w
from webui.studio import within

AMINO={'ALA','ARG','ASN','ASP','ASH','CYS','CYX','GLN','GLU','GLH','GLY','HIS','HID','HIE','HIP','ILE','LEU','LYS','LYN','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
def number(p,key,default,lo,hi):
 v=float(p.get(key,default))
 if not math.isfinite(v) or not lo<=v<=hi:raise ValueError(f'{key} Must be within {lo}–{hi} and')
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
 if p.returncode:raise RuntimeError(f'External program exited {p.returncode}; See log')

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

def build(d,c,p,out):
 import numpy as np
 from rdkit import Chem
 from openff.toolkit import Molecule
 from openff.units import unit as offunit
 from openmm import unit,Platform,XmlSerializer
 from openmm.app import PDBFile,Modeller,ForceField,PME,HBonds
 from openmmforcefields.generators import GAFFTemplateGenerator
 import parmed
 w.require(d,'screening')
 rows=w.read(d/'screening/candidates.json',[])
 row=next((r for r in rows if r['compound_id']==p.get('compound_id')),None)
 w.require(d,'docking')
 if not row or not row.get('sdf'):raise ValueError('Candidate has no traceable screening pose')
 if row.get('qualification')!='formal' and not p.get('exploratory_ack'):raise ValueError('This is a preliminary candidate; explicitly enable exploratory simulation')
 if not p.get('orientation_reviewed') or not p.get('orientation_source','').strip():raise ValueError('Review membrane orientation and provide its source or selection rationale')
 center=np.asarray(p.get('center'),dtype=float)
 if center.shape!=(3,) or not np.isfinite(center).all():raise ValueError('Invalid membrane center')
 rot=rotation(p.get('normal'))
 lipid=p.get('lipid','POPE')
 if lipid not in {'POPE','POPC','DPPC','DOPC'}:raise ValueError('Unsupported lipid')
 salt=number(p,'salt_M',.15,0,1);padding=number(p,'padding_nm',1.2,1,3)
 temp=number(p,'temperature_K',310,273,350)
 production=number(p,'production_ns',1,.01,1000);equil=number(p,'equilibration_ns',1,.1,100)
 ident='MD_'+uuid.uuid4().hex[:12];folder=system_dir(d,ident);folder.mkdir(parents=True)
 ligand_path=within(d,row['sdf']);receptor=d/'docking/receptor.pdb'
 import shutil
 shutil.copyfile(ligand_path,folder/'ligand.sdf');shutil.copyfile(receptor,folder/'receptor.pdb')
 meta={'id':ident,'compound_id':row['compound_id'],'pose_id':row.get('pose_id'),'qualification':row.get('qualification'),'status':'building','source_hashes':{'receptor':w.digest(receptor),'ligand':w.digest(ligand_path)},'parameters':{**p,'temperature_K':temp,'production_ns':production,'equilibration_ns':equil},'force_fields':['amber14-all.xml','amber14/tip3p.xml','gaff-2.11'],'membrane_model':lipid+' single-lipid idealised bilayer','orientation_rotation':rot.tolist(),'orientation_center_A':center.tolist(),'scientific_validation':'requires_review'}
 w.atomic(folder/'system.json',meta);w.atomic(out/'result.json',{'system_id':ident})
 mol=Chem.SDMolSupplier(str(ligand_path),removeHs=False)[0]
 if mol is None:raise ValueError('Ligand SDF cannot be read')
 mol=Chem.AddHs(mol,addCoords=True)
 off=Molecule.from_rdkit(mol,allow_undefined_stereo=False)
 off.name='LIG'
 for atom in off.atoms:atom.metadata['residue_name']='LIG'
 off.assign_partial_charges('am1bcc')
 off.to_file(str(folder/'parameterized_ligand.sdf'),file_format='sdf')
 gaff=GAFFTemplateGenerator(molecules=off,forcefield='gaff-2.11',cache=str(folder/'gaff-cache.json'))
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
 model.addMembrane(ff,lipidType=lipid,minimumPadding=padding*unit.nanometer,ionicStrength=salt*unit.molar,membraneCenterZ=0*unit.nanometer,platform=Platform.getPlatformByName('CPU'))
 system=ff.createSystem(model.topology,nonbondedMethod=PME,nonbondedCutoff=1*unit.nanometer,constraints=HBonds,rigidWater=True)
 structure=parmed.openmm.load_topology(model.topology,system,xyz=model.positions)
 structure.save(str(folder/'system.top'),overwrite=False);structure.save(str(folder/'system.gro'),overwrite=False)
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
 meta.update(status='built',atom_count=system.getNumParticles(),groups=groups,files={name:w.digest(folder/name) for name in ['system.top','system.gro','system.xml','system.pdb','index.ndx','ligand.sdf','receptor.pdb']})
 w.atomic(folder/'system.json',meta)
 return {'system_id':ident,'atom_count':system.getNumParticles(),'note':'Membrane and topology built; review the system before minimization/Equilibration'}

def mdp(phase,p,seed):
 temp=float(p['temperature_K']);ns=float(p['production_ns'] if phase=='production' else p['equilibration_ns'])
 em=phase=='em';nvt=phase=='nvt';restrained=phase!='production'
 steps=5000 if em else max(1,round((.1 if nvt else ns if phase=='production' else ns/2)*500000))
 lines=['integrator = '+('steep' if em else 'md'),f'nsteps = {steps}','dt = 0.002','emtol = 1000','emstep = 0.01','cutoff-scheme = Verlet','nstlist = 20','rlist = 1.0','coulombtype = PME','rcoulomb = 1.0','vdwtype = Cut-off','rvdw = 1.0','DispCorr = EnerPres','pbc = xyz','constraints = h-bonds','constraint-algorithm = lincs','nstenergy = 500','nstlog = 500','nstxout-compressed = 500','compressed-x-grps = System','nstxout = 0','nstvout = 0','nstfout = 0']
 if restrained:lines+=['define = -DPOSRES -DPOSRES_FC='+('100' if phase=='npt_soft' else '1000'),'refcoord-scaling = com']
 if not em:
  lines+=[f'ld-seed = {seed}','tcoupl = V-rescale','tc-grps = System','tau-t = 1.0',f'ref-t = {temp}','continuation = '+('no' if nvt else 'yes'),'gen-vel = '+('yes' if nvt else 'no')]
  if nvt:lines +=[f'gen-temp = {temp}',f'gen-seed = {seed}']
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


def run(d,c,p,out):
 folder=system_dir(d,p.get('system_id'));meta=w.read(folder/'system.json')
 if not meta or meta.get('status') not in {'built','completed'}:raise ValueError('System construction has not succeeded')
 if not p.get('system_reviewed'):raise ValueError('Review the membrane system and ligand correspondence')
 for name,sha in meta['files'].items():
  if w.digest(folder/name)!=sha:raise ValueError('System inputs changed; resume blocked')
 gmx=c['gmx'];progress=w.read(folder/'progress.json',{'phases':{}})
 fingerprint=hashlib.sha256(json.dumps({'files':meta['files'],'parameters':meta['parameters'],'gmx':gmx},sort_keys=True).encode()).hexdigest()
 if progress.get('fingerprint',fingerprint)!=fingerprint:raise ValueError('Simulation protocol changed; resume blocked')
 progress['fingerprint']=fingerprint
 execute([gmx,'--version'],folder,60)
 phases=['em','nvt','npt','npt_soft','production'];previous='system.gro'
 for phase in phases:
  record=progress['phases'].get(phase,{})
  if record.get('status')=='completed':
   for name,sha in record['files'].items():
    if w.digest(folder/name)!=sha:raise ValueError('Completed simulation artifact validation failed')
   previous=phase+'.gro';continue
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
  argv=[gmx,'mdrun','-deffnm',phase,'-ntmpi','1','-ntomp','4']
  if (folder/(phase+'.cpt')).exists():argv+=['-cpi',phase+'.cpt','-append']
  elif record.get('status')=='running' and (folder/(phase+'.log')).exists():raise ValueError('No checkpoint after interruption; create a new system to preserve partial output')
  execute(argv,folder,int(c.get('timeout_seconds',86400)))
  check_completion(folder,phase,recipe)
  progress['phases'][phase]={'status':'completed','tpr_sha256':tpr_hash,'files':{phase+s:w.digest(folder/(phase+s)) for s in ('.gro','.tpr','.mdp','.log','.cpt','.xtc','.edr') if (folder/(phase+s)).is_file()}}
  w.atomic(folder/'progress.json',progress);previous=phase+'.gro'
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
 progress=w.read(folder/'progress.json',{})
 for record in progress.get('phases',{}).values():
  for name,sha in record.get('files',{}).items():
   if w.digest(folder/name)!=sha:raise ValueError('Simulation output changed; analysis refused')
 dest=folder/'analysis';dest.mkdir(exist_ok=True)
 execute([c['gmx'],'trjconv','-s','production.tpr','-f','production.xtc','-o','analysis/whole.xtc','-n','index.ndx','-pbc','mol','-center'],folder,input_text='Protein\nSystem\n')
 u=mda.Universe(str(folder/'production.tpr'),str(dest/'whole.xtc'))
 protein=u.atoms[np.array(meta['groups']['Protein'])-1];back=u.atoms[np.array(meta['groups']['Backbone'])-1];ligand=u.atoms[np.array(meta['groups']['Ligand'])-1]
 ligheavy=ligand.select_atoms('not name H*');ca=protein.select_atoms('name CA')
 if len(back)<3 or not len(ligheavy):raise ValueError('Analysis selection is empty')
 u.trajectory[0];reference=back.positions.copy();ref_center=reference.mean(0);q=reference-ref_center
 ref_lig=ligheavy.positions.copy();ref_lig=ref_center+minimize_vectors(ref_lig-ref_center,u.dimensions)
 stride=max(1,len(u.trajectory)//200);series=[];frames=[];aligned=[];contact={};lipid_p=u.select_atoms('resname POPE POPC DPPC DOPC and name P')
 for ts in u.trajectory[::stride]:
  b=back.positions.copy();center=b.mean(0);P=b-center;U,S,V=np.linalg.svd(P.T@q);R=U@V
  if np.linalg.det(R)<0:U[:,-1]*=-1;R=U@V
  bx=P@R+ref_center;l=center+minimize_vectors(ligheavy.positions-center,ts.dimensions);lx=(l-center)@R+ref_center
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
 summary={'system_id':meta['id'],'compound_id':meta['compound_id'],'sampled_frames':len(series),'trajectory_frames':len(u.trajectory),'alignment':'protein backbone least squares; ligand not independently fitted','rmsf':[{'resid':int(a.resid),'rmsf_A':float(v)} for a,v in zip(ca,rmsf)],'contact_occupancy_4_5A':{r:n/len(series) for r,n in contact.items()},'limitations':['Contact occupancy uses a 4.5 Å distance cutoff, not a hydrogen-bond definition','Membrane thickness is a phosphorus-layer separation proxy','Analysis references the first production frame; compare replicas separately'],'source_hashes':{x:w.digest(folder/x) for x in ['production.tpr','production.xtc']}}
 execute([c['gmx'],'energy','-f','production.edr','-o','analysis/thermodynamics.xvg'],folder,input_text='Temperature\nPressure\nPotential\nTotal-Energy\nDensity\nBox-X\nBox-Y\nBox-Z\n0\n')
 thermo=energy_series(dest/'thermodynamics.xvg')
 summary['thermodynamics']={key:{'mean':float(np.mean([r[key] for r in thermo])),'std':float(np.std([r[key] for r in thermo])),'min':min(r[key] for r in thermo),'max':max(r[key] for r in thermo)} for key in thermo[0] if key!='time_ns'}
 summary['limitations'].append('Review temperature, pressure and energy time series; finite values and means do not establish equilibration or convergence')
 w.atomic(dest/'thermodynamics.json',thermo)
 w.atomic(dest/'summary.json',summary);w.atomic(dest/'series.json',series);w.atomic(dest/'frames.json',{'atoms':[{'chain':'REC','name':'CA','resid':int(a.resid)} for a in ca]+[{'chain':'LIG','name':a.name} for a in ligheavy],'frames':frames})
 return summary
