"""Frozen, resumable CPU Vina execution with a wall-clock cutoff."""
import concurrent.futures
import datetime
import json
import os
import statistics
import threading
import time
from pathlib import Path
from filelock import FileLock
from scientific import artifacts as a
from .common import ROOT, now, sha, command


def protocol():
    from sbvs2.docking.protocol import Protocol
    return Protocol.preset('standard',seeds=(20261023,),exhaustiveness=8,num_modes=5,energy_range=5.)


def box_for(root):
    from sbvs2.docking.box import Box
    b=a.read(Path(root)/'target/box.json')
    return Box(tuple(b['center']),tuple(b['size']),'1YCR peptide contact box')


def freeze(root,development=False):
    root=Path(root);base=root/'development' if development else root/'inputs'
    dest=base/'preregistered_protocol.json'
    if dest.exists():return dest
    rows=a.read(base/'prepared.json')
    if not rows:raise ValueError('No prepared anonymous library')
    env=a.read(root/'environment.json');p=protocol();b=box_for(root)
    if not development:
        from .scoring import filter_library
        filter_library(root)
    manifest={'schema':'yuanhe-frozen-benchmark/1','run_id':root.name,
        'purpose':'throughput_development_only' if development else 'retrospective_benchmark',
        'target_type':'soluble_PPI','frozen_at':now(),'receptor':str(root/'target/receptor.pdbqt'),
        'receptor_sha256':a.digest(root/'target/receptor.pdbqt'),'engine':env['vina_path'],'engine_sha256':a.digest(env['vina_path']),
        'protocol_hash':p.hash(),'box_hash':b.hash(),'seeds':[20261023],'exhaustiveness':8,'num_modes':5,'energy_range':5,
        'cpu_per_worker':env['cpu_per_worker'],'workers':env['workers'],'top_k':[10,20,50],
        'pose_selection':'minimum Vina affinity per unique parent; same pose used by all arms',
        'A0':'Vina affinity ascending, anonymous parent ID ascending',
        'A1':'existing peptide residue-contact recall descending, Vina affinity ascending, parent ID ascending',
        'A2':'same rule using frozen rigidly perturbed peptide contact atoms',
        'peptide_filter':'none; contact-guided reranking only',
        'acquisition':'product gate mandatory; research scientific ranking annotated unknown; deliverable requires actual evidence',
        'state_rule':__import__('scientific.mdm2.library',fromlist=['STATE_RULE']).STATE_RULE,
        'timeout_seconds':300,'max_attempts':2,
        'ligands':[{'parent_id':r['parent_id'],'state_id':r['state_id'],'path':str(base/'ligands'/r['sdf']),
                    'sha256':r['sha256']} for r in rows if r['status']=='ok'],
        'input_parent_count':len(rows),'preparation_failures':[r for r in rows if r['status']!='ok'],
        'files':{str(p.relative_to(root)):a.digest(p) for p in [root/'target/box.json',root/'target/peptide_features.json',root/'target/perturbed_features.json',root/'target/receptor_original.pdb',base/'prepared.json']+([] if development else [base/'library.json',base/'blinded_library.sdf',base/'filters.json',base/'product_filter_config.json',base/'filter_profiles.json',root/'target/chemical_features.json',root/'target_audit.json',root/'environment.json'])},
        'git_commit':command(['git','-C',str(ROOT),'rev-parse','HEAD'])['stdout'],
        'git_status':command(['git','-C',str(ROOT),'status','--porcelain'])['stdout'],
        'label_seal_sha256':None if development else a.digest(root/'evaluation/labels.json'),
        'source_files':{str(p.relative_to(ROOT)):a.digest(p) for pattern in ['00_workbench/scientific/mdm2/*.py','00_workbench/scientific/peptide_contacts.py','02_docking/src/sbvs2/docking/*.py','02_docking/src/yuanhe_preflight/*.py','02_docking/src/yuanhe_preflight/configs/*.json','01_molecule_generation/src/chemistry.py'] for p in ROOT.glob(pattern)},
        'pretrained_model_overlap':'not_applicable_no_trained_scoring_model; retrospective public structures and compounds'}
    a.atomic(dest,manifest);dest.with_suffix('.sha256').write_text(a.digest(dest)+'\n')
    return dest


def monitor(root,stop,output):
    rows=[]
    while not stop.is_set():
        def read(p):
            try:return Path(p).read_text().strip()
            except OSError:return None
        rows.append({'time':now(),'monotonic':time.monotonic(),'cpu_stat':read('/sys/fs/cgroup/cpu.stat'),
                     'memory_current':read('/sys/fs/cgroup/memory.current'),
                     'gpu':command(['nvidia-smi','--query-gpu=utilization.gpu,memory.used','--format=csv,noheader,nounits'])})
        a.atomic(output,rows);stop.wait(5)


def run(root,development=False):
    from sbvs2.docking.vina_official import run_benchmark_job
    root=Path(root);base=root/'development' if development else root/'inputs'
    frozen_path=base/'preregistered_protocol.json';frozen=a.read(frozen_path)
    if not frozen:raise ValueError('Freeze before compute')
    expected=frozen_path.with_suffix('.sha256').read_text().strip()
    if a.digest(frozen_path)!=expected:raise ValueError('Frozen protocol changed')
    for rel,digest in frozen['files'].items():
        if a.digest(root/rel)!=digest:raise ValueError('Frozen input changed: '+rel)
    # Development source is allowed to evolve before the formal freeze. Formal
    # computation refuses any edits to its actual execution/scoring sources.
    if not development:
        for rel,digest in frozen['source_files'].items():
            if a.digest(ROOT/rel)!=digest:raise ValueError('Frozen code changed: '+rel)
    manifest=a.read(root/'run_manifest.json')
    cutoff=datetime.datetime.fromisoformat(manifest['new_docking_cutoff']).timestamp()
    out=base/'jobs';out.mkdir(exist_ok=True);started=time.time();records=[]
    stop=threading.Event();thread=threading.Thread(target=monitor,args=(root,stop,base/'resource_samples.json'),daemon=True)
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
    def task(lig):
        if time.time()>=cutoff:return {'parent_id':lig['parent_id'],'status':'not_started_budget'}
        identity={k:frozen[k] for k in ['receptor_sha256','engine_sha256','protocol_hash','box_hash','cpu_per_worker']}
        identity.update(ligand_sha256=lig['sha256'],seed=20261023,preparation=frozen['state_rule'])
        cache_key=sha(identity);directory=out/lig['parent_id']/cache_key
        directory.mkdir(parents=True,exist_ok=True)
        with FileLock(str(directory/'job.lock'),timeout=1):
            for attempt in range(frozen['max_attempts']):
                jd=directory/f'attempt_{attempt+1}';rp=jd/'job_record.json'
                if rp.exists():
                    rec=a.read(rp)
                    if rec['status']=='ok' and a.digest(jd/'poses.sdf')==rec['output_sha256']:
                        return {**rec,'parent_id':lig['parent_id'],'state_id':lig['state_id'],'job_dir':str(jd),'cache_identity':cache_key,'cache_hit':True,'reuse_evidence':a.read(jd/'reuse_evidence.json')}
                    continue
                if jd.exists():continue # Interrupted attempt retained; never overwrite.
                if time.time()>=cutoff:return {'parent_id':lig['parent_id'],'status':'not_started_budget'}
                rec=run_benchmark_job(frozen['engine'],frozen['receptor'],lig['path'],box_for(root),protocol(),20261023,jd,
                    manifest=frozen_path,manifest_sha256=expected,cpu=frozen['cpu_per_worker'],
                    timeout=min(frozen['timeout_seconds'],max(1,cutoff-time.time())),meta={'parent_id':lig['parent_id'],'state_id':lig['state_id']})
                if rec['status']=='ok':break
            else:
                rec={'status':'failed','error':'all attempts failed or interrupted'}
            return {**rec,'parent_id':lig['parent_id'],'state_id':lig['state_id'],'job_dir':str(jd),'cache_identity':cache_key}
    with FileLock(str(base/'run.lock'),timeout=1):
        thread.start()
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=frozen['workers']) as pool:
                futures={pool.submit(task,lig):lig for lig in frozen['ligands']}
                for future in concurrent.futures.as_completed(futures):
                    try:record=future.result()
                    except Exception as exc:record={'parent_id':futures[future]['parent_id'],'status':'failed','error':str(exc)}
                    records.append(record);a.atomic(base/'job_records.json',records)
                    progress={'phase':'development_throughput' if development else 'formal_docking','completed':len(records),
                              'success':sum(r['status']=='ok' for r in records),'planned':len(frozen['ligands']),'updated_at':now()}
                    a.atomic(root/'progress.json',progress);print(json.dumps(progress),flush=True)
        finally:stop.set();thread.join(timeout=10)
    elapsed=time.time()-started;ok=[r for r in records if r['status']=='ok']
    times=[r['audit']['elapsed_s'] for r in records if r.get('audit')]
    summary={'phase':frozen['purpose'],'started_at':datetime.datetime.fromtimestamp(started,datetime.timezone.utc).isoformat(),
             'finished_at':now(),'elapsed_seconds':elapsed,'input_parents':frozen['input_parent_count'],
             'prepared':len(frozen['ligands']),'attempted':len(records),'success':len(ok),
             'failure_rate_including_preparation':1-len(ok)/frozen['input_parent_count'],
             'mean_job_seconds':statistics.mean(times) if times else None,'median_job_seconds':statistics.median(times) if times else None,
             'projected_550_seconds':elapsed/max(1,len(records))*550,
             'cache_hits':sum(r.get('cache_hit',False) for r in records),'cross_run_reused':sum(bool(r.get('reuse_evidence')) for r in records),
             'workers':frozen['workers'],'cpu_per_worker':frozen['cpu_per_worker'],'resource_samples':str(base/'resource_samples.json')}
    a.atomic(root/('throughput.json' if development else 'docking_summary.json'),summary)
    return summary


def dispatch(root,action):
    if action=='throughput':
        from .library import prepare_development
        prepare_development(root);freeze(root,True);return run(root,True)
    if action=='freeze':return {'manifest':str(freeze(root))}
    return run(root)
