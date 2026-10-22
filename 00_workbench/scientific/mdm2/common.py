import csv
import datetime
import hashlib
import json
import os
import subprocess
from pathlib import Path
from scientific import artifacts as a

ROOT=Path(__file__).resolve().parents[3]


def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()

def sha(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

def write_csv(path,rows,fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fields=fields or list(dict.fromkeys(k for r in rows for k in r))
    tmp=path.with_suffix('.tmp')
    with tmp.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        writer.writerows({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in rows)
    tmp.replace(path)

def command(argv):
    p=subprocess.run(argv,capture_output=True,text=True,timeout=30)
    return {'returncode':p.returncode,'stdout':p.stdout.strip(),'stderr':p.stderr.strip()}

def preflight(root,vina,started_at=None):
    import importlib.metadata as im
    import platform,shutil,sys
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    def read(path):
        p=Path(path);return p.read_text().strip() if p.is_file() else None
    affinity=len(os.sched_getaffinity(0));max_cpu=read('/sys/fs/cgroup/cpu.max')
    quota=affinity
    if max_cpu and max_cpu.split()[0]!='max':quota=min(quota,int(max_cpu.split()[0])//int(max_cpu.split()[1]))
    cpu=min(4,max(1,quota-1));workers=min(5,max(1,(quota-2)//cpu))
    from sbvs2.docking.vina_official import vina_version
    environment={'checked_at':now(),'os':platform.platform(),'python':sys.version,'executable':sys.executable,
        'affinity_threads':affinity,'cpu_max':max_cpu,'effective_cpu_quota':quota,'workers':workers,'cpu_per_worker':cpu,
        'memory_max_bytes':read('/sys/fs/cgroup/memory.max'),'memory_current_bytes':read('/sys/fs/cgroup/memory.current'),
        'disk':dict(zip(('total','used','free'),shutil.disk_usage(root))),
        'gpu':command(['nvidia-smi','--query-gpu=name,memory.total,driver_version','--format=csv,noheader']),
        'packages':{p:im.version(p) for p in ['rdkit','numpy','scipy','pandas','biopython','gemmi','meeko','openmm','pdbfixer','molscrub','psutil']},
        'vina_version':vina_version(vina),'vina_sha256':a.digest(vina),'vina_path':str(Path(vina).resolve()),
        'git_commit':command(['git','-C',str(ROOT),'rev-parse','HEAD']),
        'git_status':command(['git','-C',str(ROOT),'status','--short','--branch']),
        'generation_models':'not_required_disabled','md':'not_requested','gpu_scoring':'not_requested'}
    a.atomic(root/'environment.json',environment)
    start=datetime.datetime.fromisoformat(started_at) if started_at else datetime.datetime.now(datetime.timezone.utc)
    if start.tzinfo is None:raise ValueError('started_at requires timezone')
    manifest={'run_id':root.name,'target_type':'soluble_PPI','started_at':start.isoformat(),'hard_budget_minutes':180,
        'report_reserve_minutes':15,'hard_deadline':(start+datetime.timedelta(minutes=180)).isoformat(),
        'new_docking_cutoff':(start+datetime.timedelta(minutes=165)).isoformat(),'status':'preparing',
        'scope':'retrospective quick benchmark; no generation or MD',
        'git_commit':environment['git_commit']['stdout'],'phase_times':{'preflight_completed':now()}}
    if not (root/'run_manifest.json').exists():a.atomic(root/'run_manifest.json',manifest)
    return environment
