"""Two independent membrane jobs. Existing scientific engine and checkpoints are unchanged."""
import concurrent.futures,os,queue,shutil,subprocess,sys,threading
from pathlib import Path
from .resources import detect,disk_guard,gpu_lease

def run(d,c,out):
    from webui import workflow as w
    rc=c['research'];ids=c.get('md_system_ids',[])
    if not c.get('gromacs_system_reviewed') or not ids:raise ValueError('Reviewed membrane systems required')
    if len(ids)>rc['md_limit'] or len(ids)!=len(set(ids)):raise ValueError('Invalid MD candidate count or duplicate system')
    compounds=set()
    for ident in ids:
        if Path(ident).name!=ident:raise ValueError('Invalid system ID')
        meta=w.read(d/'md_systems'/ident/'system.json',{})
        if not meta or meta['compound_id'] in compounds:raise ValueError('Only one trajectory per compound')
        compounds.add(meta['compound_id'])
        if meta['parameters']['production_ns']!=rc['md_production_ns']:raise ValueError('Prepared production length differs from frozen research configuration')
    disk_guard(d,rc)
    gpus=[g for g in detect(d)['gpus'] if g['allocatable']][:rc['md_concurrency']]
    if not gpus:raise ValueError('No authorized allocatable GPU; MD remains queued')
    tasks=queue.Queue()
    for ident in ids:tasks.put(ident)
    state={ident:{'status':'queued','target_ns':rc['md_production_ns']} for ident in ids}
    records=[];lock=threading.Lock()
    package=out/'standalone_md';package.mkdir(exist_ok=True)
    for name in ('__init__.py','artifacts.py','membrane.py','md_cli.py'):
        dest=package/'scientific'/name;dest.parent.mkdir(exist_ok=True);shutil.copyfile(w.UI/'scientific'/name,dest)
    driver=package/'run.py'
    driver.write_text('import json,sys\nfrom pathlib import Path\nsys.path.insert(0,str(Path(__file__).resolve().parent))\nfrom scientific import membrane,artifacts\nr=json.loads(Path(sys.argv[1]).read_text());d=Path(r["run_dir"]);out=Path(r["response"]).parent\nmembrane.run(d,r["config"],{"system_id":r["system_id"],"system_reviewed":True},out)\na=membrane.analyze(d,r["config"],{"system_id":r["system_id"]},out)\nartifacts.atomic(Path(r["response"]),a)\n')
    def persist():w.atomic(d/'md_queue.json',{'tasks':state,'requested_concurrency':rc['md_concurrency'],'available_concurrency':len(gpus)})
    persist()
    def worker(gpu):
        with gpu_lease(gpu) as lease:
            while True:
                try:ident=tasks.get_nowait()
                except queue.Empty:return
                folder=d/'md_systems'/ident
                response=package/(ident+'.result.json');request=package/(ident+'.request.json')
                try:
                    estimate=rc.get('estimated_md_gb_per_candidate')
                    disk_guard(d,rc,0 if estimate is None else float(estimate)*1024**3)
                    with lock:state[ident].update(status='running',gpu=gpu['index']);persist()
                    w.atomic(request,{'run_dir':str(d),'config':c,'system_id':ident,'response':str(response)})
                    # UUID exposes exactly one device; both GROMACS run and checkpoint probe see it as device 0.
                    with (folder/'queue.log').open('a') as log:
                        subprocess.run([sys.executable,'-I',str(driver),str(request)],check=True,
                            stdout=log,stderr=subprocess.STDOUT,pass_fds=(*w.lease_fds(),lease.fileno()),
                            env={**os.environ,'CUDA_VISIBLE_DEVICES':gpu['uuid'],'YUANHE_GPU_INDEX':'0','PYTHONPATH':''})
                    analysis=w.read(response);quality=analysis.get('quality_observations',{})
                    complete=(analysis.get('production_end_ns',0)>=rc['md_production_ns'] and quality.get('finite_trajectory') is True and not quality.get('constraint_warnings',True))
                    if not complete:raise ValueError('Production duration or trajectory quality check failed')
                    rec={'system_id':ident,'compound_id':analysis['compound_id'],'analysis':analysis,'md_completed':True,
                         'md_analysis_completed':True,'files':{str(x.relative_to(d)):w.digest(x) for x in folder.rglob('*') if x.is_file()}}
                    with lock:
                        records.append(rec);state[ident].update(status='md_completed',completed_ns=analysis['production_end_ns']);persist();w.atomic(out/'systems.json',records)
                except Exception as exc:
                    with lock:state[ident].update(status='failed',error=str(exc));persist()
                finally:tasks.task_done()
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(gpus)) as pool:
        futures=[pool.submit(worker,gpu) for gpu in gpus]
        for future in futures:future.result()
    if len(records)!=len(ids):raise RuntimeError('Some MD jobs remain incomplete; checkpoints retained')
    return 0
