"""Actual AmberTools PB solver via MMPBSA.py; one processed trajectory for all terms."""
from __future__ import annotations
import csv,hashlib,io,json,os,re,shutil,subprocess
import time as clock_time
from pathlib import Path
import numpy as np
from .research import statistics

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def pb_input(settings,frames,receptor_mask,ligand_mask):
    values={'indi':float(settings.get('interior_dielectric',1)), 'exdi':float(settings.get('exterior_dielectric',80)),
            'istrng':float(settings.get('salt_molar',.15)), 'ipb':2,'inp':2,'radiopt':0}
    if not all(np.isfinite(v) for v in values.values()) or values['indi']<=0 or values['exdi']<=0 or values['istrng']<0:raise ValueError('Invalid PB parameters')
    if settings.get('entropy',False):raise ValueError('Entropy is not implemented; explicitly use entropy=false')
    if settings.get('membrane'):
        m=settings.get('membrane_model',{})
        if not m.get('reviewed') or not m.get('source') or m.get('normal')!=[0,0,1]:raise ValueError('Membrane PB requires reviewed z-aligned slab geometry and source')
        for key in ['thickness_A','dielectric','poretype']:
            if key not in m:raise ValueError('Missing membrane PB parameter: '+key)
        if not 5<float(m['thickness_A'])<80 or not 1<=float(m['dielectric'])<=20:raise ValueError('Invalid membrane slab thickness/dielectric')
        values.update(memopt=1,emem=float(m['dielectric']),mthick=float(m['thickness_A']),mctrdz=0.,poretype=int(m['poretype']),bcopt=10)
    text=f'&general\n startframe=1, endframe={frames}, interval=1, verbose=2, keep_files=2,\n receptor_mask="{receptor_mask}", ligand_mask="{ligand_mask}",\n/\n&pb\n'
    return text+',\n'.join(f' {k}={v}' for k,v in values.items())+',\n/\n',values

def parse_frames(path):
    lines=Path(path).read_text().splitlines();section=False;rows=[];header=None
    for line in lines:
        if line.strip().startswith('DELTA Energy Terms'):section=True;header=None;continue
        if not section:continue
        fields=next(csv.reader([line]))
        if fields and fields[0].strip().lower().startswith('frame'):
            header=[x.strip() for x in fields];continue
        if header and len(fields)==len(header):
            try:rows.append({k:float(v) for k,v in zip(header,fields)})
            except ValueError:break
    if not rows:raise ValueError('No DELTA frame energies parsed from real MMPBSA output')
    return rows

def run(topology,trajectories,output,settings,*,times_ns=None,protein='protein',ligand='resname LIG'):
    import MDAnalysis as mda,parmed
    from MDAnalysis.lib.distances import minimize_vectors
    from parmed.tools.actions import changeRadii
    exe=shutil.which('MMPBSA.py')
    if not exe:raise RuntimeError('AmberTools MMPBSA.py and PBSA required')
    out=Path(output).resolve()
    hashes={'topology':digest(topology),**{str(t):digest(t) for t in trajectories}}
    request={'hashes':hashes,'settings':settings,'protein':protein,'ligand':ligand,
             'times_sha256':hashlib.sha256(np.asarray(times_ns,dtype=float).tobytes()).hexdigest() if times_ns is not None else None}
    if out.exists():
        prior=json.loads((out/'request.json').read_text()) if (out/'request.json').is_file() else None
        if prior!=request:raise ValueError('PB output exists for different or unknown inputs; choose a new output directory')
        if (out/'results.json').is_file():
            result=json.loads((out/'results.json').read_text())
            for name,expected in result['output_hashes'].items():
                if digest(out/name)!=expected:raise ValueError('PB cached artifact checksum mismatch')
            return result
        out.rename(out.with_name(out.name+'.failed_'+str(clock_time.time_ns())))
    out.mkdir(parents=True,exist_ok=False)
    (out/'request.json').write_text(json.dumps(request,indent=2)+'\n')
    u=mda.Universe(str(topology),list(map(str,trajectories)))
    prot=u.select_atoms(protein);lig=u.select_atoms(ligand)
    if not len(prot) or not len(lig):raise ValueError('Explicit protein/ligand partition missing')
    start=float(settings.get('start_ns',400));end=float(settings.get('end_ns',500));stride=int(settings.get('stride',10))
    if stride<1 or not 0<=start<end:raise ValueError('Invalid PB sampling window')
    time=np.array(times_ns if times_ns is not None else [ts.time/1000 for ts in u.trajectory])
    if len(time)!=len(u.trajectory) or np.any(np.diff(time)<=0):raise ValueError('PB time axis invalid')
    if time[-1]<end-1e-6 or (settings.get('formal',True) and time[-1]<500-1e-6):raise ValueError('Full requested production window not available')
    indices=np.flatnonzero((time>=start-1e-9)&(time<=end+1e-9))[::stride]
    if len(indices)<2:raise ValueError('At least two PB frames required')
    ids=np.concatenate([prot.indices,lig.indices]);group=u.atoms[ids]
    structure=parmed.load_file(str(topology));structure.box=None
    # Boolean selection preserves MD atom order; receptor must precede ligand in dry topology.
    if np.any(np.diff(ids)<=0):raise ValueError('PB export requires receptor then ligand atom order')
    for name,selection in [('solute',ids),('receptor',prot.indices),('ligand',lig.indices)]:
        mask=np.zeros(len(structure.atoms),dtype=bool);mask[selection]=True
        dry=structure[mask.tolist()];dry.box=None
        amber=parmed.amber.AmberParm.from_structure(dry.copy(parmed.Structure));changeRadii(amber,'mbondi2').execute()
        if not np.allclose([a.charge for a in amber.atoms],[structure.atoms[i].charge for i in selection],atol=1e-8,rtol=0):raise ValueError('PB topology changed MD charges')
        amber.save(str(out/(name+'.prmtop')),overwrite=False)
        checked=parmed.load_file(str(out/(name+'.prmtop')))
        if len(checked.atoms)!=len(selection):raise ValueError('PB topology round trip failed')
    membrane=bool(settings.get('membrane'));lipids=None
    if membrane:
        lipids=u.select_atoms(settings.get('membrane_model',{}).get('phosphate_selection','resname POPE POPC DPPC DOPC POP DPP DOP and name P'))
        if len(lipids)<10:raise ValueError('Membrane PB needs explicit bilayer phosphate coordinates')
    with mda.Writer(str(out/'solute.dcd'),n_atoms=len(ids)) as writer:
        for i in indices:
            ts=u.trajectory[int(i)]
            if ts.dimensions is None:raise ValueError('PB PBC box missing')
            group.unwrap(compound='fragments',reference=None,inplace=True)
            anchor=prot.fragments[0].center_of_geometry()
            for frag in prot.fragments[1:]:
                delta=frag.center_of_geometry()-anchor;frag.positions+=minimize_vectors(delta[None,:],ts.dimensions)[0]-delta
            delta=lig.center_of_geometry()-prot.center_of_geometry();lig.positions+=minimize_vectors(delta[None,:],ts.dimensions)[0]-delta
            if membrane:
                # Keep membrane normal on z. Do not use the protein-fitted RMSD coordinates.
                z=lipids.positions[:,2];center=float(np.median(z));group.positions-=np.array([0,0,center])
            writer.write(group)
    nres=len(prot.residues);total=nres+len(lig.residues)
    text,parameters=pb_input(settings,len(indices),f':1-{nres}',f':{nres+1}-{total}')
    (out/'mmpbsa.in').write_text(text)
    command=[exe,'-O','-i','mmpbsa.in','-cp','solute.prmtop','-rp','receptor.prmtop','-lp','ligand.prmtop','-y','solute.dcd','-o','FINAL_RESULTS.dat','-eo','frame_energies.csv']
    env=os.environ.copy()
    amber_home=Path(exe).resolve().parents[1]
    if not env.get('AMBERHOME') and (amber_home/'dat'/'mmpbsa').is_dir():env['AMBERHOME']=str(amber_home)
    with (out/'backend.log').open('w') as log:subprocess.run(command,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    rows=parse_frames(out/'frame_energies.csv');
    if len(rows)!=len(indices):raise ValueError('PB energy frame count differs from selected trajectory')
    total_key=next((k for k in rows[0] if k.upper() in {'DELTA TOTAL','TOTAL'}),None)
    if not total_key or not any(k.upper() in {'EPB','PB'} for k in rows[0]):raise ValueError('PB energy component absent; cannot label calculation MM/PBSA')
    result={'method':'AmberTools MMPBSA.py / PBSA','unit':'kcal/mol','entropy_included':False,'parameters':parameters,'settings':settings,
            'selected_frame_indices':indices.tolist(),'times_ns':time[indices].tolist(),'components':{k:statistics([r[k] for r in rows]) for k in rows[0] if not k.lower().startswith('frame')},
            'uncertainty':'sample_sd measures frame variation; correlation_adjusted_sem estimates mean error under stationarity, not convergence',
            'interpretation':'Auxiliary end-state comparison; not experimental binding free energy',
            'coordinate_policy':'whole complex; membrane slab z preserved and centered from bilayer phosphates' if membrane else 'whole complex nearest images',
            'input_hashes':hashes,
            'output_hashes':{name:digest(out/name) for name in ['mmpbsa.in','solute.dcd','FINAL_RESULTS.dat','frame_energies.csv','solute.prmtop','receptor.prmtop','ligand.prmtop']},
            'topology_hashes':{n:digest(out/(n+'.prmtop')) for n in ['solute','receptor','ligand']}}
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n');u.trajectory.close();return result

def analyze_mmpbsa(config,paths):
    settings=config.get('mmpbsa',{})
    if not settings.get('enabled',False) or config.get('ligand') is None:return
    from ..utils.core import read_json,duration_ps
    for replica in range(1,config['replicates']+1):
        source=paths['md']/f'replica_{replica:02d}';p=read_json(source/'progress.json')
        if p['phase']!='complete':raise ValueError('PB analysis requires committed complete trajectory')
        trajectories=[source/s['dcd'] for s in p['segments']]
        for segment in p['segments']:
            if digest(source/segment['dcd'])!=segment['dcd_sha256']:raise ValueError('Trajectory checksum mismatch')
        import MDAnalysis as mda
        base=paths['md']/'system';u=mda.Universe(str(base/'topology.pdb'),list(map(str,trajectories)))
        times=(np.arange(len(u.trajectory))+1)*duration_ps(config['MD']['report_interval'])/1000;u.trajectory.close()
        # Complete solvated AMBER topology is exported from the same MD System.
        run(base/'system.prmtop',trajectories,paths['analysis']/f'replica_{replica:02d}'/'mmpbsa',settings,times_ns=times)
