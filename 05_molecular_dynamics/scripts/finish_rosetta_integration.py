#!/usr/bin/env python3
"""Continue one authorized public Rosetta integration after its independent build."""
from __future__ import annotations
import argparse,csv,fcntl,json,os,shutil,subprocess,sys,time
from pathlib import Path

PROJECT_ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PROJECT_ROOT))
from modules.utils.core import atomic_json,read_json,digest,load_config,config_hash,run_paths,now,WorkflowError
from run_md import verify_files
from validate_10ns import audit_completed

RUN_ID='integration_rosetta_single_file_v2'
EXPECTED_INPUT_SHA='0349394d0cf319df0cf86a6751f5b1e65cc4ccc633453134bd46437008b2fbbc'


def require(value,message):
    if not value:raise WorkflowError(message)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--poll-seconds',type=float,default=30)
    parser.add_argument('--timeout-hours',type=float,default=24)
    args=parser.parse_args()
    require(args.poll_seconds>=5 and args.timeout_hours>0,'Invalid observer interval/timeout.')
    paths=run_paths(RUN_ID);folder=paths['md'].parent
    evidence=PROJECT_ROOT/'installation_evidence'/('integration_observer_'+RUN_ID)
    evidence.mkdir(parents=True,exist_ok=True)
    lock=(evidence/'observer.lock').open('a')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise WorkflowError('An integration observer is already active.')
    config=load_config((folder/'config.yaml').relative_to(PROJECT_ROOT).as_posix())
    expected_config=config_hash(config)
    require(config['MD']['platform']=='CPU' and config['MD']['cpu_threads']==2,'Only the authorized CPU two-thread integration is allowed.')
    require(config['MD']['production']=='4ps' and config['replicates']==1,'Only the authorized one-replica 4 ps integration is allowed.')
    require(digest(folder/'inputs/complex.pdb')==EXPECTED_INPUT_SHA,'The authorized public single-PDB input changed.')
    record={'run_id':RUN_ID,'status':'WAITING_FOR_BUILD','started':now(),'pid':os.getpid(),
        'input_sha256':EXPECTED_INPUT_SHA,'config_sha256':expected_config,
        'observer_sha256':digest(Path(__file__)),'platform':'CPU','threads':2,
        'scope':'One existing public integration: preserve saved config, continue once after build, never certify 10 ns.'}
    status_path=evidence/'status.json'
    def publish():
        record['updated']=now();atomic_json(status_path,record)
    publish();started=time.monotonic()
    try:
        while True:
            require(time.monotonic()-started<args.timeout_hours*3600,'Observer timed out; no workflow was stopped or restarted.')
            manifest=read_json(folder/'manifest.json')
            require(manifest['config_sha256']==expected_config,'Committed scientific configuration changed.')
            require(config_hash(load_config((folder/'config.yaml').relative_to(PROJECT_ROOT).as_posix()))==expected_config,'Saved configuration changed.')
            require(manifest.get('status') not in ['failed','interrupted','error'],'Build/workflow did not complete; inspect its original error, no automatic retry.')
            record['manifest_status']=manifest.get('status');record['build_status']=manifest['stages'].get('build',{}).get('status')
            if manifest['stages'].get('build',{}).get('status')=='complete' and not (folder/'.running').exists():break
            publish();time.sleep(args.poll_seconds)
        for stage in ['prepare','build']:verify_files(manifest['stages'][stage])
        import openmm as mm
        import parmed as pmd
        system=mm.XmlSerializer.deserialize((paths['md']/'system/system.xml').read_text())
        topology=read_json(paths['md']/'system/topology.json')
        ligand=pmd.load_file(str(paths['md']/'system/ligand.prmtop'))
        charge=sum(a.charge for a in ligand.atoms)
        require(system.getNumParticles()==topology['atoms'],'Built System/topology atom count mismatch.')
        require(len(ligand.atoms)==len(topology['ligand_indices'])==129 and abs(charge)<1e-4,'Parameterized MZ1 atom count or net charge differs.')
        build={'run_id':RUN_ID,'status':'PASS','stage':'build','audited_utc':now(),
            'system_atoms':topology['atoms'],'ligand_atoms':129,'ligand_amber_net_charge_e':charge,
            'ligand_forcefield':'GAFF 2.11 / AM1-BCC','own_environment_and_no_upstream_parameter_cache':True,
            'system_sha256':digest(paths['md']/'system/system.xml'),
            'ligand_prmtop_sha256':digest(paths['md']/'system/ligand.prmtop'),
            'prepared_and_built_artifact_hashes_verified':True}
        atomic_json(PROJECT_ROOT/'installation_evidence/rosetta_md_build_validation.json',build)
        if manifest.get('status')!='complete':
            record['status']='RUNNING_SAME_CONFIG_COMPLETION';publish()
            command=[sys.executable,str(PROJECT_ROOT/'scripts/run_md.py'),'--resume','--run-id',RUN_ID]
            record['command']=['python','scripts/run_md.py','--resume','--run-id',RUN_ID]
            with (evidence/'continuation.log').open('w') as output:
                result=subprocess.run(command,cwd=PROJECT_ROOT,stdout=output,stderr=subprocess.STDOUT)
            record['continuation_exit_code']=result.returncode
            require(result.returncode==0,'Integration continuation failed; inspect continuation.log. No automatic retry.')
        record['status']='INDEPENDENT_SHORT_RUN_AUDIT';publish()
        audit=audit_completed(RUN_ID)
        require(audit['production_steps']==2000 and abs(audit['production_ns']-0.004)<1e-10,'Actual production does not equal the authorized 4 ps.')
        audit['scope']='4 ps numerical interoperability only; not 10 ns completion, convergence, binding affinity or native-system stability.'
        audit['build_validation']=build
        atomic_json(evidence/'completion_audit.json',audit)
        relative_audit=(evidence/'completion_audit.json').relative_to(PROJECT_ROOT).as_posix()
        metadata_path=PROJECT_ROOT/'installation_evidence/input_contract_tests.json'
        metadata_lock=(evidence/'input_contract_update.lock').open('a')
        fcntl.flock(metadata_lock,fcntl.LOCK_EX)
        data=read_json(metadata_path)
        for case in data['cases']:
            if case['run_id']==RUN_ID:
                case.update(status='PASS',stage='acceptance',build_status='PASS',simulation_platform='CPU',
                    system_atoms=topology['atoms'],ligand_atoms=129,ligand_forcefield=build['ligand_forcefield'],
                    ligand_amber_net_charge_e=charge,production_steps=2000,production_ns=0.004,
                    frames=audit['replicas'][0]['frames'],checkpoint_sha256=audit['numerical_acceptance']['replicates'][0]['checkpoint_sha256'],
                    mmgbsa=audit['numerical_acceptance']['mmgbsa_status'],mmgbsa_entropy_included=False,
                    full_integration_evidence=relative_audit)
                break
        else:raise WorkflowError('The original preparation case is absent; refusing to invent its provenance.')
        data['scope']='All three actual upstream single-PDB cases completed independent preparation, force-field build, minimization, NVT, NPT, 4 ps production and applicable analysis. These are numerical interoperability tests, not 10 ns completion or scientific binding/convergence validation.'
        atomic_json(metadata_path,data);metadata_lock.close()
        report=PROJECT_ROOT/'INSTALLATION_REPORT.md';marker='## Actual Rosetta single-PDB full 4 ps integration'
        if marker not in report.read_text():
            with report.open('a') as handle:
                handle.write('\n'+marker+'\n\n')
                handle.write(f'`results/{RUN_ID}`: PASS. Independent GAFF 2.11 / AM1-BCC in the MD runtime produced a {topology["atoms"]:,}-atom solvated system with 129 ligand atoms and net ligand Amber charge {charge:.8f} e. The existing saved CPU two-thread configuration completed minimization, 1 ps NVT, 1 ps NPT, and 2,000 production steps = 4 ps, followed by real AmberTools MM/GBSA and all applicable analysis. Independent verification reloaded the CPU binary checkpoint and checked each DCD/CSV hash, frame/time continuity and analysis artifact hash.\n\n')
                handle.write(f'Evidence: `{relative_audit}`. Only the upstream complex.pdb was copied. This is a short interoperability test; the module 10 ns validation status and MD_TEST_PASS.md are not modified by this observer.\n')
        record.update(status='PASS',completed=now(),production_steps=2000,production_ns=0.004,
            completion_audit=relative_audit,completion_audit_sha256=digest(evidence/'completion_audit.json'))
        publish();print(json.dumps(record),flush=True);return 0
    except BaseException as exc:
        record.update(status='FAILED',completed=now(),error_type=type(exc).__name__,error=str(exc))
        candidate=folder/'SHORT_RUN_TEST.md'
        if candidate.exists() and f'Run ID: {RUN_ID}' in candidate.read_text().splitlines():
            archived=evidence/('unverified_short_run_report_'+str(time.time_ns())+'.md')
            shutil.copyfile(candidate,archived)
            if digest(candidate)==digest(archived):candidate.unlink();record['unverified_report_archived']=archived.name
        publish();print(json.dumps(record),file=sys.stderr,flush=True);return 2


if __name__=='__main__':
    raise SystemExit(main())
