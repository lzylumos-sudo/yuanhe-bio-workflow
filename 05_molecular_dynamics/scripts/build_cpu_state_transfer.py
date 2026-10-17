"""Build an audited zero-step CPU state bridge with an explicit new RNG epoch.

The caller freezes the chosen run before invoking prepare_transfer. This module
never signals a process or writes a run file. The target directory must be empty.
Original CPU binary checkpoints are loaded only with their original RNG count.
"""
from __future__ import annotations
from datetime import datetime,timezone
from pathlib import Path
import gc,hashlib,json,math,os,re,shutil,struct,time,xml.etree.ElementTree as ET


def _sha(path):
    value=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):value.update(block)
    return value.hexdigest()


def _require(condition,message):
    if not condition:raise RuntimeError(message)


def _now():return datetime.now(timezone.utc).isoformat()


def _rng_streams(blob,seed):
    # Read-only diagnostic pinned to the inspected OpenMM 8.6.1 CpuRandom format.
    # Checkpoint bytes are never patched or used as a migration implementation.
    prefix=struct.pack('<ii',1,seed);matches=[];start=0
    while True:
        offset=blob.find(prefix,start)
        if offset<0:break
        if offset+12<=len(blob):
            count=struct.unpack_from('<i',blob,offset+8)[0]
            if 1<=count<=1024:matches.append({'offset':offset,'streams':count})
        start=offset+1
    _require(len(matches)==1,'Ambiguous version-specific RNG stream diagnostic; refusing bridge.')
    return matches[0]


def prepare_transfer(root:Path,run:Path,out:Path,target_threads:int=16,new_seed:int=20262609):
    """Return a successful probe record, or write failed probe.json and raise.

    Preconditions: caller has frozen run; root is its module; out is an empty
    evidence directory outside run or a new nested run/state_transfers directory.
    All restoration/testing occurs on copies.
    out/native16.chk and native16.xml remain at the original production step.
    A separate Context executes 25 test steps and saves only validation artifacts.
    No old binary checkpoint is loaded into a differently sized CPU Context.
    """
    import numpy as np
    import openmm as mm
    from openmm import app,unit
    import yaml
    root=Path(root).resolve();run=Path(run).resolve();out=Path(out).resolve()
    _require(run.is_relative_to(root/'results'),'Run must belong to the explicitly supplied module.')
    embedded_audit=out.is_relative_to(run/'state_transfers') and len(out.relative_to(run/'state_transfers').parts)>=2
    _require(not out.is_relative_to(run) or embedded_audit,
        'Output inside a run is limited to a new nested state_transfers evidence directory.')
    _require(target_threads==16,'This reviewed builder is restricted to the authorized CPU16 bridge.')
    _require(isinstance(new_seed,int) and 0<new_seed<2**31-100000,'Invalid reproducible runtime seed.')
    _require(not out.exists() or not any(out.iterdir()),'Bridge evidence output must be empty.')
    out.mkdir(parents=True,exist_ok=True)
    report={'passed':False,'status':'RUNNING','started':_now(),'run_id':run.name,
        'openmm_version':mm.__version__,'new_threads':target_threads,
        'integrator_seed':new_seed,'barostat_seed':new_seed+100000,
        'caller_must_freeze_run':True,'cpu_affinity':sorted(os.sched_getaffinity(0)),
        'scope':'Physical state transfer with explicit new thermostat/barostat RNG epoch; not full stochastic-state continuity.'}
    def publish():
        temporary=out/'probe.json.tmp';temporary.write_text(json.dumps(report,indent=2)+'\n');temporary.replace(out/'probe.json')
    publish()
    try:
        _require(mm.__version__=='8.6.1','The RNG diagnostic was reviewed for OpenMM 8.6.1 only.')
        source_progress=run/'md/replica_01/progress.json'
        progress_bytes=source_progress.read_bytes();progress=json.loads(progress_bytes)
        config_bytes=(run/'config.yaml').read_bytes();config=yaml.safe_load(config_bytes);md=config['MD']
        _require(md['platform']=='CPU','Source platform is not CPU.')
        old_threads=int(md['cpu_threads']);old_seed=int(progress['seed'])
        _require(old_seed!=new_seed,'A new explicit RNG epoch requires a distinct seed.')
        _require(progress['phase']=='production','Source must be a committed production checkpoint.')
        _require(progress['production_steps']==progress['current_step']==sum(s['steps'] for s in progress['segments']),
            'Committed steps are inconsistent.')
        checkpoint=progress['checkpoint'];_require(Path(checkpoint).name==checkpoint,'Unsafe checkpoint filename.')
        paths={'system.xml':run/'md/system/system.xml','topology.pdb':run/'md/system/topology.pdb',
            'source.chk':run/'md/replica_01'/checkpoint,
            'source.xml':run/'md/replica_01'/(Path(checkpoint).stem+'.xml')}
        snapshots=out/'source';snapshots.mkdir()
        source_sha={name:_sha(path) for name,path in paths.items()}
        _require(source_sha['source.chk']==progress['checkpoint_sha256'],'Source checkpoint checksum mismatch.')
        _require(source_sha['system.xml']==progress['system_sha256'],'Source System checksum mismatch.')
        for name,path in paths.items():
            shutil.copyfile(path,snapshots/name)
            _require(_sha(snapshots/name)==source_sha[name],'Source changed during snapshot.')
        (snapshots/'config.yaml').write_bytes(config_bytes)
        (snapshots/'progress.json').write_bytes(progress_bytes)
        report.update(old_threads=old_threads,old_seed=old_seed,old_barostat_seed=old_seed+100000,
            source_checkpoint=checkpoint,source_checkpoint_sha256=source_sha['source.chk'],
            source_xml_sha256=source_sha['source.xml'],source_system_sha256=source_sha['system.xml'],
            source_config_sha256=hashlib.sha256(config_bytes).hexdigest(),
            source_progress_sha256=hashlib.sha256(progress_bytes).hexdigest(),
            step=progress['current_step'],source_sha256=source_sha)
        source_rng=_rng_streams((snapshots/'source.chk').read_bytes(),old_seed)
        _require(source_rng['streams']==old_threads,'Source binary RNG count differs from CPU thread configuration.')
        report['source_rng_diagnostic']=source_rng
        value=re.fullmatch(r'\s*([0-9.]+)\s*(fs|ps|ns)\s*',str(md['timestep']))
        _require(value is not None,'Unsupported timestep representation.')
        dt_ps=float(value[1])*{'fs':0.001,'ps':1.0,'ns':1000.0}[value[2]]
        platform=mm.Platform.getPlatformByName('CPU')
        def simulation(threads,seed):
            system=mm.XmlSerializer.deserialize((snapshots/'system.xml').read_text())
            barostat=next(f for f in system.getForces() if isinstance(f,mm.MonteCarloBarostat))
            barostat.setRandomNumberSeed(seed+100000);barostat.setFrequency(25)
            integrator=mm.LangevinMiddleIntegrator(md['temperature']*unit.kelvin,
                md['friction_per_ps']/unit.picosecond,dt_ps*unit.picosecond)
            integrator.setRandomNumberSeed(seed);integrator.setConstraintTolerance(1e-6)
            sim=app.Simulation(app.PDBFile(str(snapshots/'topology.pdb')).topology,system,integrator,platform,{'Threads':str(threads)})
            return sim,barostat
        def snapshot(state):
            return {'positions':np.array(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer),copy=True),
                'velocities':np.array(state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond),copy=True),
                'box':np.array(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer),copy=True),
                'step':state.getStepCount(),'time_ps':state.getTime().value_in_unit(unit.picosecond),
                'parameters':dict(state.getParameters())}
        def context_snapshot(context):
            return snapshot(context.getState(getPositions=True,getVelocities=True,getParameters=True,getIntegratorParameters=True))
        def compare(a,b):
            result={key+'_max_abs_delta':float(np.max(np.abs(a[key]-b[key]))) for key in ['positions','velocities','box']}
            result.update(step_delta=a['step']-b['step'],time_ps_delta=a['time_ps']-b['time_ps'],parameters_exact=a['parameters']==b['parameters'])
            for key in ['positions','velocities','box']:result[key+'_exact']=bool(np.array_equal(a[key],b[key]))
            result['exactly_equal']=all(result[k] for k in ['positions_exact','velocities_exact','box_exact','parameters_exact']) and result['step_delta']==0 and result['time_ps_delta']==0
            return result
        def effective(sim,barostat):
            integrator=sim.integrator
            return {'actual_threads':int(platform.getPropertyValue(sim.context,'Threads')),
                'integrator_seed':integrator.getRandomNumberSeed(),'barostat_seed':barostat.getRandomNumberSeed(),
                'temperature_K':integrator.getTemperature().value_in_unit(unit.kelvin),
                'friction_per_ps':integrator.getFriction().value_in_unit(unit.picosecond**-1),
                'timestep_ps':integrator.getStepSize().value_in_unit(unit.picosecond),
                'constraint_tolerance':integrator.getConstraintTolerance(),
                'barostat_pressure_atm':barostat.getDefaultPressure().value_in_unit(unit.atmosphere),
                'barostat_temperature_K':barostat.getDefaultTemperature().value_in_unit(unit.kelvin),
                'barostat_frequency':barostat.getFrequency()}
        old,old_baro=simulation(old_threads,old_seed)
        old.loadCheckpoint(str(snapshots/'source.chk'))
        original=context_snapshot(old.context);old_settings=effective(old,old_baro)
        _require(original['step']==progress['current_step'],'Actual source checkpoint step mismatch.')
        _require(math.isclose(original['time_ps'],original['step']*dt_ps,abs_tol=1e-5),'Actual source checkpoint time mismatch.')
        source_state=mm.XmlSerializer.deserialize((snapshots/'source.xml').read_text())
        xml_snapshot=snapshot(source_state)
        xml_comparison=compare(xml_snapshot,original)
        _require(xml_comparison['exactly_equal'],'Source XML is not identical to the committed binary physical state.')
        xml_root=ET.parse(snapshots/'source.xml').getroot();integrator_node=xml_root.find('IntegratorParameters')
        report.update(source_xml_vs_binary=xml_comparison,old_effective_settings=old_settings,time_ps=original['time_ps'],
            source_xml_metadata={'root_attributes':dict(xml_root.attrib),'state_data_flags':source_state.getDataTypes(),
                'integrator_parameters_present':integrator_node is not None,
                'integrator_parameters_xml':ET.tostring(integrator_node,encoding='unicode') if integrator_node is not None else None})
        del old,old_baro;gc.collect()
        new,new_baro=simulation(target_threads,new_seed)
        new.loadState(str(snapshots/'source.xml'))
        restored=context_snapshot(new.context);comparison=compare(restored,original);new_settings=effective(new,new_baro)
        _require(comparison['exactly_equal'],'Physical state changed during XML transfer.')
        _require(new_settings['actual_threads']==target_threads,'Incorrect target Threads property.')
        _require(new_settings['integrator_seed']==new_seed and new_settings['barostat_seed']==new_seed+100000,
            'XML load changed the explicitly assigned new RNG seeds.')
        unchanged=['temperature_K','friction_per_ps','timestep_ps','constraint_tolerance','barostat_pressure_atm','barostat_temperature_K','barostat_frequency']
        _require(all(new_settings[key]==old_settings[key] for key in unchanged),'Public physical integration settings changed.')
        for key in ['positions_exact','velocities_exact','box_exact','parameters_exact']:report[key]=comparison[key]
        report.update(state_comparison=comparison,new_effective_settings=new_settings,actual_threads=target_threads,
            new_seed_confirmed_after_XML_load=True,public_integrator_parameters_equal=True,
            old_binary_never_loaded_into_new_thread_count=True)
        new.saveCheckpoint(str(out/'native16.chk'));new.saveState(str(out/'native16.xml'))
        native_rng=_rng_streams((out/'native16.chk').read_bytes(),new_seed)
        _require(native_rng['streams']==target_threads,'New checkpoint does not contain native16 RNG streams.')
        report.update(native16_checkpoint_sha256=_sha(out/'native16.chk'),native16_xml_sha256=_sha(out/'native16.xml'),
            native16_rng_diagnostic=native_rng,native16_step=new.currentStep,native16_advanced_steps=0)
        _require(new.currentStep==report['step'],'Zero-step bridge was accidentally advanced.')
        del new,new_baro;gc.collect();publish()
        verified,verified_baro=simulation(target_threads,new_seed)
        verified.loadCheckpoint(str(out/'native16.chk'))
        reload_comparison=compare(context_snapshot(verified.context),original)
        roundtrip=verified.context.createCheckpoint()==(out/'native16.chk').read_bytes()
        _require(reload_comparison['exactly_equal'] and roundtrip,'Native16 binary roundtrip verification failed.')
        report.update(native16_roundtrip_verified=True,native16_reload_comparison=reload_comparison,
            native16_roundtrip_exact=roundtrip,native16_reload_effective_settings=effective(verified,verified_baro))
        start=time.monotonic();verified.step(25);continued=context_snapshot(verified.context)
        energy=verified.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        finite=bool(np.isfinite(continued['positions']).all() and np.isfinite(continued['velocities']).all() and math.isfinite(energy))
        _require(continued['step']==original['step']+25 and math.isclose(continued['time_ps'],original['time_ps']+25*dt_ps,abs_tol=1e-7) and finite,
            'Independent native16 short continuation failed.')
        verified.saveCheckpoint(str(out/'validation_25steps.chk'));verified.saveState(str(out/'validation_25steps.xml'))
        _require(_rng_streams((out/'validation_25steps.chk').read_bytes(),new_seed)['streams']==target_threads,'Continued checkpoint lost native16 RNG streams.')
        report['continuation_test']={'steps':25,'actual_threads':target_threads,'final_step':continued['step'],
            'final_time_ps':continued['time_ps'],'potential_energy_kj_mol':energy,'finite':finite,'wall_seconds':time.monotonic()-start,
            'separate_context_and_checkpoint':True}
        del verified,verified_baro;gc.collect()
        _require(_sha(out/'native16.chk')==report['native16_checkpoint_sha256'],'Pristine bridge checkpoint changed.')
        _require(_sha(out/'native16.xml')==report['native16_xml_sha256'],'Pristine bridge XML changed.')
        unchanged_sources={name:_sha(path)==source_sha[name] for name,path in paths.items()}
        unchanged_sources['config.yaml']=(run/'config.yaml').read_bytes()==config_bytes
        unchanged_sources['progress.json']=source_progress.read_bytes()==progress_bytes
        _require(all(unchanged_sources.values()),'Source run changed during transfer; caller must freeze it and retry in a new directory.')
        report.update(passed=True,status='PASS_PHYSICAL_STATE_MIGRATION_WITH_NEW_RNG_EPOCH',completed=_now(),
            source_files_unchanged=unchanged_sources,no_minimization_equilibration_or_velocity_reassignment=True,
            reset_internal_state=['CPU Langevin RNG streams','barostat RNG','barostat proposal adaptation and internal counters'],
            unchanged_science=['force-field System source','positions','velocities','periodic box','time','step','Context parameters',*unchanged],
            warning='XML restores physical state and public parameters; RNG/internal adaptation restarts. No claim of exact full stochastic state or future trajectory equivalence.')
        publish();return report
    except BaseException as exc:
        report.update(passed=False,status='FAILED',completed=_now(),error_type=type(exc).__name__,error=str(exc));publish();raise
