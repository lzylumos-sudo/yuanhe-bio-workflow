"""Local UI services. All expensive scientific work runs in isolated workers."""
from __future__ import annotations
import base64,hashlib,json,math,os,shutil,subprocess,sys,threading,time,uuid
from pathlib import Path
from . import workflow as w

def config(rid):
 d=w.run_dir(rid);c=w.read(d/'config.json')
 if c is None:raise ValueError('Run not found')
 if w.digest(d/'config.json')!=w.read(d/'run.json')['config_sha256']:raise ValueError('Configuration validation failed')
 for stage,_,_ in w.STAGES:c.update(w.read(d/'reviews'/f'{stage}.json',{}))
 return d,c

def within(d,p):
 p=(d/p).resolve()
 if not p.is_relative_to(d.resolve()) or not p.is_file():raise ValueError('Artifact missing or path outside permitted scope')
 return p

def upload(payload):
 name=Path(payload.get('name','')).name
 if not name.lower().endswith('.pdb'):raise ValueError('Structure uploads currently accept PDB; Please export PDB format')
 content=payload.get('text','')
 if not isinstance(content,str) or len(content.encode())>8*1024*1024:raise ValueError('Maximum structure file size 8 MB')
 if not any(l.startswith('ATOM  ') for l in content.splitlines()):raise ValueError('File contains no PDB atom records')
 dest=w.ROOT/'inputs'/f'{uuid.uuid4().hex}.pdb';dest.parent.mkdir(exist_ok=True)
 dest.write_text(content,encoding='utf-8')
 return {'structure':str(dest.relative_to(w.ROOT)),'name':name,'sha256':w.digest(dest)}

def core(rid,seed_id):
 from rdkit import Chem
 from rdkit.Chem.Draw import rdMolDraw2D
 d,c=config(rid);w.require(d,'peptide2mol')
 row=next((r for r in w.read(d/'peptide2mol/seeds.json',[]) if r['id']==seed_id),None)
 if not row:raise ValueError('Seed not found')
 m=Chem.SDMolSupplier(str(within(d,row['raw_sdf'])),removeHs=True)[0]
 if m is None:raise ValueError('Unreadable structure')
 draw=rdMolDraw2D.MolDraw2DSVG(800,500);draw.drawOptions().addAtomIndices=False
 for a in m.GetAtoms():a.SetProp('atomNote',str(a.GetIdx()+1))
 rdMolDraw2D.PrepareAndDrawMolecule(draw,m)
 points=[{'index':i+1,'x':draw.GetDrawCoords(i).x,'y':draw.GetDrawCoords(i).y,'element':m.GetAtomWithIdx(i).GetSymbol()} for i in range(m.GetNumAtoms())]
 draw.FinishDrawing()
 return {'svg':draw.GetDrawingText(),'atoms':points,'seed_id':seed_id,'smiles':row['smiles']}

def save_core(rid,payload):
 d,c=config(rid);seed=payload.get('seed_id');info=core(rid,seed)
 ids=payload.get('indices',[])
 if not ids or any(type(i)!=int or i<1 or i>len(info['atoms']) for i in ids):raise ValueError('Select valid core atoms')
 if not payload.get('reason','').strip():raise ValueError('Record the chemical rationale for preserving this core')
 rows=[r for r in c.get('core_selections',[]) if r['seed_id']!=seed]
 rows.append({'seed_id':seed,'key_atom_indices':sorted(set(ids)),'review_reason':payload['reason'][:2000]})
 w.review(rid,'reinvent4',{'reinvent_executable':c['reinvent_executable'],'core_selections':rows})
 return {'selections':rows}

def structure(path):
 atoms=[]
 for l in path.read_text().splitlines():
  if l.startswith(('ATOM  ','HETATM')):
   atoms.append({'p':[float(l[30:38]),float(l[38:46]),float(l[46:54])],'chain':l[21:22],'resid':int(l[22:26]) if l[22:26].strip().lstrip('-').isdigit() else l[22:26].strip(),'resname':l[17:20].strip(),'name':l[12:16].strip(),'element':l[76:78].strip()})
 return atoms

def box(rid,payload):
 d,c=config(rid);w.require(d,'prepare')
 atoms=structure(d/'prepare/receptor_full.pdb')
 ids=payload.get('residues',[50,51,53,54,237,240,241,244])
 if not isinstance(ids,list) or not ids or any(type(i)!=int for i in ids):raise ValueError('Invalid residue numbers')
 region=[a for a in atoms if a['chain']==c['receptor_chain'] and a['resid'] in ids and a['element']!='H']
 if set(ids)-{a['resid'] for a in region}:raise ValueError('Some residues have no coordinates')
 padding=float(payload.get('padding',5))
 if not math.isfinite(padding) or not 1<=padding<=15:raise ValueError('Padding must be 1–15 Å')
 lo=[min(a['p'][i] for a in region) for i in range(3)];hi=[max(a['p'][i] for a in region) for i in range(3)]
 b={'center':[(lo[i]+hi[i])/2 for i in range(3)],'size':[max(14,hi[i]-lo[i]+2*padding) for i in range(3)]}
 return {'box':b,'atoms':atoms+[{**a,'role':'peptide'} for a in structure(d/'prepare/peptide_reference.pdb')],'review_required':True}

def operations(rid):
 d,c=config(rid);out=[]
 for p in sorted((d/'operations').glob('*/operation.json'),reverse=True):
  x=w.read(p);x['id']=p.parent.name
  if x['status']=='running' and rid not in w._active and not w.live_record(x,rid):x={**x,'status':'interrupted','error':'Service restarted; confirm the old process has stopped before resuming'}
  log=p.parent/'process.log'
  if log.exists():
   with log.open('rb') as f:f.seek(max(0,log.stat().st_size-20000));x['log']=f.read().decode('utf-8','replace')
  x['result']=w.read(p.parent/'result.json',{})
  out.append(x)
 return out[:30]

def launch(rid,kind,payload):
 if kind not in {'receptor','environment','install_membrane','md_build','md_run','md_analyze','md_validation'}:raise ValueError('Unknown operation')
 d,c=config(rid);w.verify_source(d)
 if c.get('research') and kind=='md_run':raise ValueError('Use the research MD stage queue to run/resume selected systems')
 if c.get('research') and kind=='md_build':
  if float(payload.get('production_ns',0))!=c['research']['md_production_ns']:raise ValueError('MD length differs from frozen research configuration')
  existing=[w.read(p,{}) for p in (d/'md_systems').glob('*/system.json')]
  if any(x.get('compound_id')==payload.get('compound_id') for x in existing):raise ValueError('Existing system for this compound; resume it instead of creating a repeat')
  if len(existing)>=c['research']['md_limit']:raise ValueError('MD candidate budget exhausted')
 with w._lock:
  if w._active:raise ValueError('A task is running; finish or stop it first')
  if kind=='receptor':w.require(d,'prepare')
  if kind=='md_build':w.require(d,'screening')
  if kind.startswith('md_') and not c.get('md_python'):raise ValueError('Select the membrane simulation environment in settings Python')
  lease=w.acquire_compute()
  ident=str(time.time_ns())+'_'+uuid.uuid4().hex[:8];out=d/'operations'/ident;out.mkdir(parents=True)
  w.atomic(out/'request.json',{'kind':kind,'payload':payload,'config':c,'run_id':rid,'source_revision':w.read(d/'active_source_revision.json')})
  w.atomic(out/'operation.json',{'kind':kind,'status':'running','started':time.time()})
  py=(c.get('md_python') if kind.startswith('md_') else c.get('prep_python') if kind=='receptor' else c['chem_python']) or sys.executable
  log=(out/'process.log').open('w')
  try:proc=subprocess.Popen([py,str(w.UI/'studio_worker.py'),'--run-id',rid,'--operation',ident],stdout=log,stderr=subprocess.STDOUT,cwd=w.UI,start_new_session=True,pass_fds=(lease.fileno(),),env={**os.environ,'PATH':str(Path(py).parent)+os.pathsep+os.environ.get('PATH',''),'OPENMM_CPU_THREADS':'12','OMP_NUM_THREADS':'12','YUANHE_COMPUTE_FD':str(lease.fileno())})
  except Exception as e:
   log.close();lease.close();w.atomic(out/'operation.json',{'kind':kind,'status':'failed','error':str(e)});raise
  w._active[rid]=proc
  w.atomic(out/'operation.json',{**w.read(out/'operation.json'),**w.process_identity(proc.pid)})
  def wait():
   code=proc.wait();log.close();lease.close()
   with w._lock:
    w._active.pop(rid,None);r=w.read(out/'operation.json')
    if r['status']=='running':w.atomic(out/'operation.json',{**r,'status':'failed','error':f'Process exited {code}, without complete results','ended':time.time()})
  threading.Thread(target=wait,daemon=True).start()
 return {'id':ident,'kind':kind}

def apply_receptor(rid,payload):
 d,c=config(rid);op=payload.get('operation','')
 if not op or '/' in op or '\\' in op:raise ValueError('Invalid operation identifier')
 folder=d/'operations'/op
 if w.read(folder/'operation.json',{}).get('status')!='completed':raise ValueError('Receptor preparation incomplete')
 result=w.read(folder/'result.json',{})
 path=within(d,result.get('receptor',''))
 if w.digest(path)!=result['sha256']:raise ValueError('Receptor changed')
 if not payload.get('reviewed'):raise ValueError('Review protonation, residues and termini')
 values={} if c.get('research') else {k:c[k] for k in w.REVIEW_FIELDS['docking'] if k in c}
 values['docking_receptor']=str(path.relative_to(d));values['docking_box']=payload['box']
 values['receptor_review']={'sha256':result['sha256'],'pH':result['pH'],'reviewed':True,'operation':op}
 w.review(rid,'docking',values)
 return {'saved':True}

def system_summary(meta):
 meta=dict(meta);groups=meta.pop('groups',{})
 if groups:meta['group_atom_counts']={k:len(v) for k,v in groups.items()}
 if 'lipid_charge_rounding' in meta:
  rounding=meta['lipid_charge_rounding'];edits=rounding.get('edits',[])
  meta['lipid_charge_rounding']={k:v for k,v in rounding.items() if k!='edits'}
  meta['lipid_charge_rounding']['corrected_lipids']=len(edits)
  meta['lipid_charge_rounding']['max_abs_correction_e']=max((abs(e['delta_charge_e']) for e in edits),default=0)
 return meta

def md_list(rid):
 d,c=config(rid)
 rows=[{'id':p.parent.name,**system_summary(w.read(p)),'progress':w.read(p.parent/'progress.json',{})} for p in sorted((d/'md_systems').glob('*/system.json'))]
 for ident,record in w.read(d/'inherited_md.json',{}).items():
  source=w.run_dir(record['source_run']);rows.append({'id':ident,**system_summary(w.read(source/'md_systems'/ident/'system.json')),'reused_from':source.name})
 return sorted(rows,key=lambda x:(x.get('status') not in {'built','completed'},x['id']))

def membrane_preview(rid,payload):
 import numpy as np
 from rdkit import Chem
 from scientific.membrane import rotation
 d,c=config(rid);w.require(d,'docking')
 center=np.asarray(payload.get('center'),dtype=float)
 if center.shape!=(3,) or not np.isfinite(center).all():raise ValueError('Invalid membrane center')
 rot=rotation(payload.get('normal'));atoms=structure(d/'docking/receptor.pdb')
 row=next((x for x in w.read(d/'screening/candidates.json',[]) if x['compound_id']==payload.get('compound_id')),None)
 if row and row.get('sdf'):
  mol=Chem.SDMolSupplier(str(within(d,row['sdf'])),removeHs=True)[0]
  if mol:
   atoms.extend({'p':list(mol.GetConformer().GetAtomPosition(a.GetIdx())),'chain':'LIG','name':'CA','element':a.GetSymbol()} for a in mol.GetAtoms())
 for atom in atoms:atom['p']=((np.array(atom['p'])-center)@rot.T).tolist()
 return {'atoms':atoms,'box':{'center':[0,0,0],'size':[100,100,30]},'caption':'Membrane orientation preview: wireframes represent layers around the center ±15 Å ; actual thickness depends on the selected lipid'}

def md_data(rid,system_id):
 d,c=config(rid)
 if not system_id or any(x in system_id for x in ('/','\\','..')):raise ValueError('Invalid system identifier')
 folder=d/'md_systems'/system_id
 thermo=w.read(folder/'analysis/thermodynamics.json',[])
 thermo=thermo[::max(1,math.ceil(len(thermo)/2000))]
 return {'system':system_summary(w.read(folder/'system.json',{})),'analysis':w.read(folder/'analysis/summary.json',{}),'trajectory':w.read(folder/'analysis/frames.json',{}),'series':w.read(folder/'analysis/series.json',{}),'thermodynamics':thermo}

def md_structure(rid,system_id):
 from scientific.membrane import system_dir
 d,c=config(rid);folder=system_dir(d,system_id);meta=w.read(folder/'system.json',{})
 if meta.get('status') not in {'built','completed'}:raise ValueError('System has not been built successfully')
 if w.digest(folder/'system.pdb')!=meta['files']['system.pdb']:raise ValueError('System structure changed')
 protein=set(meta['groups']['Protein']);ligand=set(meta['groups']['Ligand']);atoms=[]
 for i,atom in enumerate(structure(folder/'system.pdb'),1):
  if i in ligand and atom['element']!='H':atom['chain']='LIG';atoms.append(atom)
  elif i in protein and atom['name']=='CA':atom['role']='protein';atoms.append(atom)
  elif i not in protein and i not in ligand and atom['name']=='P':atom['role']='lipid';atoms.append(atom)
  elif atom['resname'] in {'NA','CL','Na','Cl'}:atom['role']='ion';atoms.append(atom)
 return {'atoms':atoms,'caption':'Actual built coordinates: protein CA (green), ligand heavy atoms (gold), lipid phosphorus (purple), ions (blue). Water and remaining atoms are omitted only from this preview. Full coordinates remain in system.pdb.','system_id':system_id,'source_sha256':meta['files']['system.pdb']}

def environment_save(payload):
 allowed={'chem_python','screen_python','prep_python','md_python','peptide_python','reinvent_executable','vina','gmx'}
 if set(payload)-allowed:raise ValueError('Unknown environment setting')
 if any(not isinstance(v,str) or len(v)>2000 for v in payload.values()):raise ValueError('Invalid executable path')
 cfg=w.read(w.ROOT/'deployment/runtime.json',{});cfg.update(payload);w.atomic(w.ROOT/'deployment/runtime.json',cfg)
 return {'saved':True,'note':'Used for new runs; existing runs retain frozen settings'}

# Environment diagnosis must be available even when no research task exists.
_diagnostics_lock = threading.Lock()
_diagnostics_state = {'status':'idle','checks':{}}

def diagnostics_status():
 import json
 with _diagnostics_lock:return json.loads(json.dumps(_diagnostics_state))

def diagnostics_start(payload):
 import json,threading
 from .readiness import collect
 global _diagnostics_state
 if set(payload)-{'run_id'}:raise ValueError('Unknown diagnostics setting')
 cfg=config(payload['run_id'])[1] if payload.get('run_id') else w.defaults()
 with w._lock:
  if w._active:raise ValueError('Finish or stop the active computation before diagnostics')
  with _diagnostics_lock:
   if _diagnostics_state['status']=='running':raise ValueError('Diagnostics are already running')
   lease=w.acquire_compute()
   ident=str(time.time_ns())
   report=w.ROOT/'installation_evidence'/('readiness_'+ident+'.json')
   _diagnostics_state={'id':ident,'status':'running','checks':{},'report':str(report)}
 def publish(value):
  global _diagnostics_state
  snapshot=json.loads(json.dumps({**value,'id':ident,'report':str(report)}))
  w.atomic(report,snapshot)
  with _diagnostics_lock:_diagnostics_state=snapshot
 def worker():
  try:collect(cfg,w.ROOT,ui_python=sys.executable,publish=publish)
  except Exception as exc:publish({'status':'failed','checks':diagnostics_status().get('checks',{}),'error':str(exc),'all_runtime_checks_pass':False})
  finally:lease.close()
 threading.Thread(target=worker,name='deployment-readiness',daemon=True).start()
 return {'id':ident,'status':'running','report':str(report)}

def export_run(rid,payload):
 import zipfile
 d,c=config(rid)
 with w._lock:
  if w._active:raise ValueError('Wait for the active task to finish before exporting a consistent snapshot')
  lease=w.acquire_compute()
  target=d/'exports'/('run_'+str(time.time_ns())+'.zip');target.parent.mkdir(exist_ok=True)
  full=bool(payload.get('full'));omitted=[];manifest={}
  try:
   paths=[]
   for path in sorted(d.rglob('*')):
    if not path.is_file() or path.is_symlink() or target.parent in path.parents:continue
    rel=str(path.relative_to(d))
    if not full and (path.suffix.lower() in {'.xtc','.trr','.cpt','.tpr','.edr'} or path.stat().st_size>20*1024*1024):
     omitted.append({'path':rel,'sha256':w.digest(path),'bytes':path.stat().st_size});continue
    paths.append(path)
   with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=3,allowZip64=True) as archive:
    for path in paths:
     rel=str(path.relative_to(d));manifest[rel]=w.digest(path);archive.write(path,rel)
    release=w.read(w.ROOT/'deployment/release.json',{})
    archive.writestr('EXPORT_MANIFEST.json',json.dumps({'run_id':rid,'created':time.time(),'mode':'full' if full else 'review','files':manifest,'omitted':omitted,'release':release,'note':'Large files omitted from review archives remain in the original run directory; full archives still exclude weights and installed environments'},ensure_ascii=False,indent=2))
   return {'url':f'/artifacts/{rid}/exports/{target.name}?download=1','files':len(manifest),'omitted':len(omitted),'sha256':w.digest(target)}
  except Exception:
   target.unlink(missing_ok=True);raise
  finally:lease.close()
