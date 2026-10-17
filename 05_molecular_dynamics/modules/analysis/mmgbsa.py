"""Real AmberTools MMPBSA.py single-trajectory MM/GBSA; entropy excluded."""
import csv,json,math,os,re,shutil,subprocess
from pathlib import Path
import numpy as np
from ..utils.core import WorkflowError,atomic_json,read_json,duration_ps
from .trajectory import write_csv

def analyze_mmgbsa(config,paths):
    from openmm import app
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import minimize_vectors
    settings=config.get('mmgbsa',{})
    base=paths['md']/'system'; meta=read_json(base/'topology.json')
    nsolute=len(meta['protein_indices'])+len(meta['ligand_indices'])
    has_ligand=bool(meta['ligand_indices'])
    for replica in range(1,config['replicates']+1):
        out=paths['analysis']/f'replica_{replica:02d}'; out.mkdir(parents=True,exist_ok=True)
        if not has_ligand:
            status={'status':'NOT_APPLICABLE','reason':'No small molecule. A protein-protein thermodynamic partition is not inferred.','entropy_included':False}
            write_csv(out/'MMGBSA.csv',[status])
            optional=read_json(out/'optional_analysis.json'); optional['MMGBSA']=status; atomic_json(out/'optional_analysis.json',optional)
            continue
        if not settings.get('enabled',True):
            status={'status':'DISABLED','reason':'Explicit configuration disabled MMGBSA.','entropy_included':False}
            write_csv(out/'MMGBSA.csv',[status])
            optional=read_json(out/'optional_analysis.json'); optional['MMGBSA']=status; atomic_json(out/'optional_analysis.json',optional)
            continue
        executable=shutil.which('MMPBSA.py')
        if not executable: raise WorkflowError('AmberTools MMPBSA.py is required for real MM/GBSA.')
        calc=out/'mmgbsa'; calc.mkdir(exist_ok=True)
        for name in ['solute.prmtop','receptor.prmtop','ligand.prmtop']:
            if not (base/name).is_file(): raise WorkflowError('AMBER topology export is missing: '+name)
            shutil.copyfile(base/name,calc/name)
        source=paths['md']/f'replica_{replica:02d}'; progress=read_json(source/'progress.json')
        u=mda.Universe(str(base/'topology.pdb'),[str(source/s['dcd']) for s in progress['segments']])
        if hasattr(u.atoms,'bonds'): u.delete_bonds(u.bonds)
        u.add_TopologyAttr('bonds',meta['bonds'])
        group=u.atoms[:nsolute]; protein=u.atoms[meta['protein_indices']]; ligand=u.atoms[meta['ligand_indices']]
        start=float(settings.get('start_ns',config['analysis'].get('discard_ns',0)))
        interval=duration_ps(config['MD']['report_interval'])/1000
        eligible=[i for i in range(len(u.trajectory)) if (i+1)*interval>=start]
        max_frames=int(settings.get('max_frames',200))
        if max_frames<2: raise WorkflowError('mmgbsa.max_frames must be at least 2.')
        chosen=np.unique(np.linspace(0,len(eligible)-1,min(len(eligible),max_frames),dtype=int)) if eligible else []
        indices=[eligible[i] for i in chosen]
        if len(indices)<2: raise WorkflowError('At least two MM/GBSA frames are required.')
        times=[]
        with mda.Writer(str(calc/'solute.dcd'),n_atoms=nsolute) as writer:
            for i in indices:
                ts=u.trajectory[i]; box=ts.dimensions.copy()
                group.unwrap(compound='fragments',reference=None,inplace=True)
                anchor=protein.fragments[0].center_of_geometry()
                for frag in protein.fragments[1:]:
                    delta=frag.center_of_geometry()-anchor; frag.positions+=minimize_vectors(delta[None,:],box)[0]-delta
                delta=ligand.center_of_geometry()-protein.center_of_geometry()
                ligand.positions+=minimize_vectors(delta[None,:],box)[0]-delta
                writer.write(group); times.append((i+1)*interval)
        u.trajectory.close()
        prmtop=app.AmberPrmtopFile(str(calc/'solute.prmtop'))
        protein_residues=len({a.residue for a in list(prmtop.topology.atoms())[:len(meta['protein_indices'])]})
        total_residues=prmtop.topology.getNumResidues()
        receptor_mask=f':1-{protein_residues}'
        ligand_mask=f':{protein_residues+1}-{total_residues}'
        igb=int(settings.get('igb',5)); salt=float(settings.get('salt_molar',0.15))
        if igb!=5: raise WorkflowError('Validated MM/GBSA protocol currently requires igb=5 and mbondi2 radii.')
        (calc/'mmgbsa.in').write_text(f'&general\n startframe=1, endframe={len(times)}, interval=1, verbose=1, keep_files=0,\n receptor_mask="{receptor_mask}", ligand_mask="{ligand_mask}",\n/\n&gb\n igb={igb}, saltcon={salt},\n/\n')
        command=[executable,'-O','-i','mmgbsa.in','-cp','solute.prmtop','-rp','receptor.prmtop','-lp','ligand.prmtop','-y','solute.dcd','-o','FINAL_RESULTS.dat','-eo','frame_energies.csv']
        with (calc/'backend.log').open('w') as log:
            result=subprocess.run(command,cwd=calc,stdout=log,stderr=subprocess.STDOUT,env=os.environ.copy())
        if result.returncode: raise WorkflowError('MMPBSA.py failed; local backend.log retains the diagnostic.')
        text=(calc/'FINAL_RESULTS.dat').read_text()
        matches=re.findall(r'^DELTA TOTAL\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)',text,re.M)
        if len(matches)!=1: raise WorkflowError('Cannot uniquely parse Amber MM/GBSA DELTA TOTAL.')
        mean,sd,sem=map(float,matches[0])
        if not all(map(math.isfinite,[mean,sd,sem])): raise WorkflowError('Non-finite MM/GBSA results.')
        row={'status':'COMPUTED','method':'AmberTools MMPBSA.py single-trajectory MMGBSA','mean_kcal_mol':mean,'std_kcal_mol':sd,'sem_kcal_mol':sem,
             'frames':len(times),'start_ns':times[0],'end_ns':times[-1],'igb':igb,'salt_molar':salt,'radii':'mbondi2','entropy_included':False,
             'partition':'all proteins (receptor) versus small molecule (ligand)','interpretation':'end-state energy estimate; not experimental affinity or ternary cooperativity'}
        write_csv(out/'MMGBSA.csv',[row]); atomic_json(calc/'protocol.json',{**row,'source_frame_indices':indices,'times_ns':times,'command':[Path(command[0]).name]+command[1:]})
        optional=read_json(out/'optional_analysis.json') if (out/'optional_analysis.json').is_file() else {}
        optional['MMGBSA']=row; atomic_json(out/'optional_analysis.json',optional)
    overall=read_json(paths['analysis']/'analysis.json')
    overall['MMGBSA']='COMPUTED' if has_ligand and settings.get('enabled',True) else ('NOT_APPLICABLE' if not has_ligand else 'DISABLED')
    atomic_json(paths['analysis']/'analysis.json',overall)
