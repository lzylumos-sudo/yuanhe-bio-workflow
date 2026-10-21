"""Conservative runtime detection and cross-process GPU leases; no preemption."""
import contextlib, csv, fcntl, os, shutil, subprocess, tempfile
from pathlib import Path

def detect(path):
    affinity=len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else (os.cpu_count() or 1)
    quota=None
    try:
        q,period=Path('/sys/fs/cgroup/cpu.max').read_text().split()
        if q!='max': quota=max(1,int(q)//int(period));affinity=min(affinity,quota)
    except (OSError,ValueError): pass
    memory={}
    try:
        memory={line.split(':')[0]:int(line.split()[1])*1024 for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith(('MemTotal:','MemAvailable:'))}
    except OSError: pass
    for name in ('memory.max','memory.current'):
        try: memory[name]=Path('/sys/fs/cgroup/'+name).read_text().strip()
        except OSError: pass
    gpus=[];error=None
    try:
        p=subprocess.run(['nvidia-smi','--query-gpu=index,uuid,name,memory.total,memory.free,utilization.gpu,driver_version','--format=csv,noheader,nounits'],capture_output=True,text=True,timeout=10,check=True)
        allowed=os.environ.get('CUDA_VISIBLE_DEVICES',os.environ.get('NVIDIA_VISIBLE_DEVICES','all'))
        for row in csv.reader(p.stdout.splitlines()):
            index,uid,name,total,free,util,driver=map(str.strip,row)
            if allowed not in ('all',) and index not in allowed.split(',') and uid not in allowed.split(','):continue
            gpus.append({'index':index,'uuid':uid,'name':name,'total_mib':int(total),'free_mib':int(free),'utilization':int(util),'driver':driver,'allocatable':int(util)<10 and int(free)>int(total)*.8})
    except (OSError,ValueError,subprocess.SubprocessError) as exc:error=str(exc)
    physical=None
    try:
        lines=subprocess.run(['lscpu','-p=CORE,SOCKET'],capture_output=True,text=True,timeout=5,check=True).stdout.splitlines()
        physical=len({x for x in lines if x and not x.startswith('#')})
    except (OSError,subprocess.SubprocessError):pass
    return {'logical_cpus':os.cpu_count(),'physical_cores':physical,'cpu_quota':quota,'effective_cpus':affinity,
            'memory':memory,'disk_free_bytes':shutil.disk_usage(path).free,'gpus':gpus,'gpu_detection_error':error,
            'cuda_compatibility':'requires backend runtime validation','container':Path('/.dockerenv').exists(),
            'scheduler_gpu_allocation':{k:os.environ.get(k) for k in ('CUDA_VISIBLE_DEVICES','NVIDIA_VISIBLE_DEVICES','SLURM_JOB_GPUS')}}

def disk_guard(path,cfg,additional_bytes=0):
    if shutil.disk_usage(path).free < cfg['disk_reserve_gb']*1024**3+additional_bytes:
        raise RuntimeError('Insufficient free disk space; existing trajectories preserved')

@contextlib.contextmanager
def gpu_lease(gpu):
    # UUID prevents ambiguity between CUDA visibility masks in separate processes.
    path=Path(tempfile.gettempdir())/('yuanhe-gpu-'+gpu['uuid']+'.lock')
    handle=path.open('a+')
    try:
        fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        yield handle
    finally:handle.close()
