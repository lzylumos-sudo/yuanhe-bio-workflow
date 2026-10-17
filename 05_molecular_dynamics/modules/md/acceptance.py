"""Check actual production steps/time/checkpoint and required numerical outputs."""
import csv,math
import numpy as np
from ..utils.core import read_json,digest,duration_ps,phase_steps,WorkflowError

def verify_run(config,paths):
    import openmm as mm
    from openmm import unit
    required=['RMSD.png','RMSF.png','Ligand_RMSD.png','Hbond.csv','Contact_map.png','SASA.csv','MMGBSA.csv']
    replicas=[]; total=phase_steps(config)['production']
    target_ns=total*duration_ps(config['MD']['timestep'])/1000
    for r in range(1,config['replicates']+1):
        source=paths['md']/f'replica_{r:02d}'; out=paths['analysis']/f'replica_{r:02d}'
        p=read_json(source/'progress.json')
        if p['phase']!='complete' or p['production_steps']!=total or p['current_step']!=total:
            raise WorkflowError('Production step audit failed.')
        if digest(source/p['checkpoint'])!=p['checkpoint_sha256']: raise WorkflowError('Final checkpoint hash failed.')
        state=mm.XmlSerializer.deserialize((source/(p['checkpoint'].removesuffix('.chk')+'.xml')).read_text())
        time_ns=state.getTime().value_in_unit(unit.nanosecond)
        if not math.isclose(time_ns,target_ns,abs_tol=1e-6): raise WorkflowError('Final checkpoint XML time differs from configured production.')
        if not np.isfinite(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)).all(): raise WorkflowError('Final positions are non-finite.')
        quality=read_json(out/'quality.json')
        if not quality['finite_observables'] or quality['frames_expected']!=quality['frames_observed']: raise WorkflowError('Trajectory/thermodynamic audit failed.')
        for name in required:
            if not (out/name).is_file() or (out/name).stat().st_size==0: raise WorkflowError('Required output missing: '+name)
        with (out/'MMGBSA.csv').open() as f: mmgbsa=next(csv.DictReader(f))
        if config.get('ligand') and config.get('mmgbsa',{}).get('enabled',True) and mmgbsa['status']!='COMPUTED': raise WorkflowError('Real MM/GBSA was not computed.')
        replicas.append({'replica':r,'production_steps':total,'production_ns':time_ns,'frames':quality['frames_observed'],'checkpoint_sha256':p['checkpoint_sha256'],'platform':p['platform'],'mmgbsa':mmgbsa['status']})
    return {'passed':True,'production_ns':target_ns,'replicates':replicas,'mmgbsa_status':replicas[0]['mmgbsa'],
            'claims':'Numerical workflow completion only; no inference of convergence, affinity, pose correctness, recruitment, or degradation.'}
