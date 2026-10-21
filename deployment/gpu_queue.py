#!/usr/bin/env python3
"""One process per CUDA device, with cross-process leases and explicit resume argv.
The queue executes only when invoked by the operator. Never launches on import.
"""
import argparse,concurrent.futures,fcntl,hashlib,json,os,subprocess,threading,time
from pathlib import Path

def main():
 p=argparse.ArgumentParser();p.add_argument('manifest',type=Path);p.add_argument('--gpus',required=True,help='CUDA device indices or UUIDs, comma separated');p.add_argument('--output',required=True,type=Path);a=p.parse_args()
 devices=a.gpus.split(',')
 if len(set(devices))!=len(devices) or any(not d.strip() for d in devices):raise ValueError('Distinct explicit GPU identifiers required')
 manifest=json.loads(a.manifest.read_text());jobs=manifest['jobs'];ids=[j['id'] for j in jobs]
 if len(set(ids))!=len(ids):raise ValueError('Duplicate job IDs')
 for job in jobs:
  if not job['id'].replace('_','').replace('-','').isalnum():raise ValueError('Invalid job ID')
  for key in ['argv','resume_argv']:
   if not isinstance(job.get(key),list) or not job[key] or not all(isinstance(v,str) for v in job[key]):raise ValueError('Each job requires explicit argv and resume_argv arrays')
 a.output.mkdir(parents=True,exist_ok=True);pending=iter(jobs);mutex=threading.Lock()
 lease_root=Path(os.environ.get('YUANHE_GPU_LOCK_DIR','/tmp/yuanhe-gpu-leases'));lease_root.mkdir(parents=True,exist_ok=True)
 def worker(device):
  failures=[]
  # CUDA_VISIBLE_DEVICES maps the selected physical device to process-local index 0.
  env={**os.environ,'CUDA_VISIBLE_DEVICES':device,'YUANHE_GPU_INDEX':'0'}
  with (lease_root/(hashlib.sha256(device.encode()).hexdigest()+'.lock')).open('a+') as lease:
   fcntl.flock(lease,fcntl.LOCK_EX)
   while True:
    with mutex:job=next(pending,None)
    if job is None:break
    folder=a.output/job['id'];folder.mkdir(exist_ok=True);record=folder/'job.json'
    key=hashlib.sha256(json.dumps(job,sort_keys=True).encode()).hexdigest()
    if record.exists() and json.loads(record.read_text())['request_sha256']!=key:raise ValueError('Queue request changed; use a new output directory')
    command=job['resume_argv'] if record.exists() else job['argv']
    state={'request_sha256':key,'gpu':device,'process_local_gpu':'0','argv':command,'cwd':job['cwd'],'started':time.time(),'status':'running'}
    def save():
     tmp=record.with_suffix('.tmp');tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(record)
    save()
    with (folder/'execution.log').open('ab') as log:
     process=subprocess.run(command,cwd=job['cwd'],env=env,stdout=log,stderr=subprocess.STDOUT,check=False)
    state.update(status='complete' if process.returncode==0 else 'failed',returncode=process.returncode,finished=time.time());save()
    if process.returncode:failures.append(job['id'])
  return failures
 with concurrent.futures.ThreadPoolExecutor(max_workers=len(devices)) as pool:
  failures=sum(list(pool.map(worker,devices)),[])
 if failures:raise SystemExit('Failed queue jobs (logs retained): '+', '.join(failures))
if __name__=='__main__':main()
