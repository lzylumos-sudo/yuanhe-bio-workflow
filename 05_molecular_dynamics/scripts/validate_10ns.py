#!/usr/bin/env python3
"""Independent completed-trajectory audit; never advances or edits a simulation."""
from __future__ import annotations
import argparse,csv,json,math,sys
from pathlib import Path

PROJECT_ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PROJECT_ROOT))
from modules.utils.core import (WorkflowError,atomic_json,read_json,digest,load_config,
    run_paths,relative_path,input_files,config_hash,duration_ps,phase_steps,now)
from modules.md.acceptance import verify_run
from run_md import implementation,verify_files
from state_transfer_audit import verify_state_transfers


def require(condition,message):
    if not condition: raise WorkflowError(message)


def close(actual,expected,message,tol=1e-4):
    require(math.isclose(float(actual),float(expected),rel_tol=1e-7,abs_tol=tol),message)


def audit_completed(run_id):
    import numpy as np
    import openmm as mm
    from openmm import app,unit
    from MDAnalysis.coordinates.DCD import DCDReader
    paths=run_paths(run_id); folder=paths['md'].parent
    require(not (folder/'.running').exists(),'Run is still locked/running; audit a completed run.')
    config=load_config((folder/'config.yaml').relative_to(PROJECT_ROOT).as_posix())
    manifest=read_json(folder/'manifest.json')
    require(manifest.get('status')=='complete','Run manifest is not complete.')
    require(manifest.get('accepted') is True,'Run was not accepted by its normal workflow.')
    require(config_hash(config)==manifest['config_sha256'],'Scientific configuration hash changed.')
    current={name:digest(relative_path(name,must_exist=True)) for name in input_files(config)}
    require(current==manifest['input_sha256'],'Original staged input hash changed.')
    require(implementation()==manifest['md_implementation_sha256'],'MD implementation differs from the run fingerprint.')
    state_transfer_audit=verify_state_transfers(folder,config,manifest)
    require(not state_transfer_audit.get('candidate_only',False),'Final acceptance cannot use candidate-only transfer evidence.')
    for name in ['prepare','build','md','analysis']:
        record=manifest['stages'].get(name,{})
        require(record.get('status')=='complete',name+' stage is incomplete.')
        verify_files(record)
    numerical=verify_run(config,paths)
    timestep=duration_ps(config['MD']['timestep'])
    interval_ps=duration_ps(config['MD']['report_interval'])
    interval_steps=round(interval_ps/timestep)
    total=phase_steps(config)['production']
    require(total%interval_steps==0,'Production end does not coincide with a saved trajectory frame.')
    source_hashes={'config.yaml':digest(folder/'config.yaml'),'manifest.json':digest(folder/'manifest.json'),
      'system.xml':digest(paths['md']/'system/system.xml'),**current}
    source_hashes.update(state_transfer_audit['source_sha256'])
    replicas=[]
    for replica in range(1,config['replicates']+1):
        base=paths['md']/'system'; target=paths['md']/f'replica_{replica:02d}'
        p=read_json(target/'progress.json')
        require(p['phase']=='complete' and p['production_steps']==total and p['current_step']==total,'Incorrect final production step.')
        require(p['system_sha256']==digest(base/'system.xml'),'Progress system hash differs.')
        require(sum(s['steps'] for s in p['segments'])==total,'Segment steps do not sum to final production.')
        completed=0; frames=0; details=[]
        atoms=read_json(base/'topology.json')['atoms']
        for segment in p['segments']:
            nsteps=segment['steps']; expected_frames=nsteps//interval_steps
            require(nsteps>0 and nsteps%interval_steps==0,'Segment has incomplete report interval.')
            for kind in ['dcd','csv']:
                require(digest(target/segment[kind])==segment[kind+'_sha256'],'Committed '+kind+' hash mismatch.')
            with (target/segment['csv']).open() as handle:
                reader=csv.reader(handle); header=next(reader)
                rows=list(reader)
            require(len(rows)==expected_frames,'CSV frame count differs from committed steps.')
            values=np.asarray(rows,dtype=float)
            require(bool(np.isfinite(values).all()),'CSV contains non-finite thermodynamics.')
            require(header[0].strip('#"')=='Step' and header[1]=='Time (ps)','Unexpected state CSV columns.')
            expected_steps=completed+np.arange(1,expected_frames+1)*interval_steps
            require(bool(np.array_equal(values[:,0],expected_steps)),'CSV production steps are discontinuous.')
            require(bool(np.allclose(values[:,1],expected_steps*timestep,rtol=1e-7,atol=1e-4)),'CSV times disagree with steps.')
            with DCDReader(str(target/segment['dcd'])) as trajectory:
                require(len(trajectory)==expected_frames and trajectory.n_atoms==atoms,'DCD frame/atom count mismatch.')
                close(trajectory.dt,interval_ps,'DCD time interval mismatch.')
                first=trajectory[0].time; last=trajectory[-1].time
                close(first,interval_ps,'DCD first segment-relative time mismatch.')
                close(last,nsteps*timestep,'DCD last segment-relative time mismatch.')
                for frame in trajectory:
                    require(bool(np.isfinite(frame.positions).all()),'DCD coordinates are not finite.')
                    require(frame.dimensions is not None and bool(np.isfinite(frame.dimensions).all()) and bool((frame.dimensions[:3]>0).all()),'DCD periodic box is invalid.')
            absolute_end_ps=completed*timestep+last
            close(absolute_end_ps,(completed+nsteps)*timestep,'Mapped DCD absolute time differs from committed end.')
            details.append({'dcd':segment['dcd'],'frames':expected_frames,'committed_start_step':completed,
                'raw_dcd_first_ps':first,'raw_dcd_last_ps':last,'mapped_absolute_last_ps':absolute_end_ps,
                'dcd_sha256':segment['dcd_sha256'],'csv_sha256':segment['csv_sha256']})
            completed+=nsteps;frames+=expected_frames
        require(bool(details),'No committed trajectory segments.')
        close(details[-1]['mapped_absolute_last_ps'],total*timestep,'DCD production end time mismatch.')
        # Load, inspect and destroy a separate CPU Context. Never call step().
        require(p['platform']=='CPU','This audit only reloads CPU checkpoints; GPU execution is not silently selected.')
        system=mm.XmlSerializer.deserialize((base/'system.xml').read_text())
        barostat=next(f for f in system.getForces() if isinstance(f,mm.MonteCarloBarostat))
        barostat.setRandomNumberSeed(p['seed']+100000);barostat.setFrequency(25)
        md=config['MD']
        integrator=mm.LangevinMiddleIntegrator(md['temperature']*unit.kelvin,md['friction_per_ps']/unit.picosecond,timestep*unit.picosecond)
        integrator.setRandomNumberSeed(p['seed']);integrator.setConstraintTolerance(1e-6)
        topology=app.PDBFile(str(base/'topology.pdb')).topology
        simulation=app.Simulation(topology,system,integrator,mm.Platform.getPlatformByName('CPU'),{'Threads':str(md['cpu_threads'])})
        simulation.loadCheckpoint(str(target/p['checkpoint']))
        loaded=simulation.context.getState(getPositions=True,getVelocities=True)
        require(simulation.currentStep==total,'Actual binary checkpoint step count mismatch.')
        checkpoint_ps=loaded.getTime().value_in_unit(unit.picosecond)
        close(checkpoint_ps,total*timestep,'Actual binary checkpoint time mismatch.')
        require(bool(np.isfinite(loaded.getPositions(asNumpy=True).value_in_unit(unit.nanometer)).all()),'Checkpoint positions are non-finite.')
        require(bool(np.isfinite(loaded.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond)).all()),'Checkpoint velocities are non-finite.')
        del simulation,integrator,system
        source_hashes[f'replica_{replica:02d}/progress.json']=digest(target/'progress.json')
        source_hashes[f'replica_{replica:02d}/'+p['checkpoint']]=digest(target/p['checkpoint'])
        replicas.append({'replica':replica,'steps':total,'frames':frames,'checkpoint_time_ns':checkpoint_ps/1000,
            'mapped_dcd_end_ns':details[-1]['mapped_absolute_last_ps']/1000,'segments':details})
    return {'run_id':run_id,'audited_utc':now(),'passed':True,'timestep_ps':timestep,
        'production_ns':total*timestep/1000,'production_steps':total,'replicas':replicas,'numerical_acceptance':numerical,
        'state_transfer_audit':state_transfer_audit,
        'checks':{'completed_stages_and_config':True,'input_and_stage_hashes':True,'segment_dcd_csv_hashes':True,
         'continuous_csv_steps_and_times':True,'actual_dcd_frames_atoms_times_and_coordinates':True,
         'actual_binary_checkpoint_loaded_on_CPU':True,'checkpoint_steps_times_finite_positions_velocities':True,
         'required_analysis_and_real_mmgbsa':True,'documented_state_transfer_and_rng_epoch':True},'source_sha256':source_hashes,
        'dcd_time_definition':'Each OpenMM segment DCD starts at one report interval. Absolute production time = independently verified committed preceding-step offset + DCD time; checked against global state CSV.',
        'scope':'10ns numerical workflow completion only; no convergence, affinity or correct-pose claim.'}


def require_10ns(record):
    close(record['timestep_ps'],0.002,'10 ns demo acceptance requires the configured 2 fs timestep.',tol=1e-12)
    require(record['production_steps']>=5_000_000 and record['production_ns']>=10-1e-8,
        f"10 ns acceptance refused: only {record['production_steps']} production steps / {record['production_ns']} ns completed.")


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run-id',required=True)
    args=parser.parse_args()
    record=audit_completed(args.run_id);require_10ns(record)
    record['checks']['at_least_5000000_steps_at_2fs']=True
    record['validator_sha256']=digest(Path(__file__))
    out=PROJECT_ROOT/'installation_evidence/10ns_final_audit.json'
    atomic_json(out,record)
    folder=run_paths(args.run_id)['md'].parent
    report='\n'.join(['# MD test PASS','',f'Run ID: {args.run_id}',
        f'Production: {record["production_steps"]:,} steps per replica at 2 fs = {record["production_ns"]:g} ns.',
        f'Replicas: {len(record["replicas"])}.',
        'Completed input/build/analysis hashes, actual CPU binary checkpoint step/time, every DCD/CSV hash, finite coordinates, continuous global CSV times, and offset-mapped DCD final times passed.',
        'MM/GBSA is real AmberTools single-trajectory endpoint analysis; entropy excluded.',
        'Scope: numerical workflow completion only; no convergence, binding-affinity or correct-pose claim.',
        'Detailed evidence: installation_evidence/10ns_final_audit.json',''])
    if record['state_transfer_audit']['verified_transfers']:
        report+='\nA documented CPU10-to-CPU16 XML physical-state transfer preserved coordinates, velocities, box and committed time/segments. Langevin/barostat RNG streams and barostat adaptation were restarted; this is not bitwise binary continuation.\n'
    (folder/'MD_TEST_PASS.md').write_text(report)
    (PROJECT_ROOT/'MD_TEST_PASS.md').write_text(report)
    atomic_json(PROJECT_ROOT/'validation_status.json',{'status':'PASS','scope':'10ns numerical workflow completion only',
        'run_id':args.run_id,'production_ns':record['production_ns'],'production_steps':record['production_steps'],
        'replicates':len(record['replicas']),'audited_utc':record['audited_utc'],
        'evidence':'installation_evidence/10ns_final_audit.json','evidence_sha256':digest(out)})
    print(json.dumps({'status':'PASS','run_id':args.run_id,'evidence':str(out.relative_to(PROJECT_ROOT))}))


if __name__=='__main__':
    try:main()
    except (WorkflowError,FileNotFoundError,KeyError,ValueError) as exc:
        print('NOT ACCEPTED: '+str(exc),file=sys.stderr);sys.exit(2)
