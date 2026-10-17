import os,signal,subprocess,time
from pathlib import Path
from ..paths import write_json,sha
class Cancelled(Exception):pass
def execute(argv,folder,timeout,cancel=lambda:False):
 folder=Path(folder);start=time.time();env=os.environ.copy();env.update(OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
 libs=[str(p) for p in Path('/opt/sbvs/env/lib/python3.11/site-packages/nvidia').glob('*/lib')]
 env['LD_LIBRARY_PATH']=':'.join(libs+['/opt/sbvs/env/lib','/usr/lib/wsl/lib',env.get('LD_LIBRARY_PATH','')])
 with (folder/'stdout.log').open('w') as stdout,(folder/'stderr.log').open('w') as stderr:
  p=subprocess.Popen([str(a) for a in argv],stdout=stdout,stderr=stderr,env=env,start_new_session=True,shell=False)
  audit={'argv':[str(a) for a in argv],'pid':p.pid,'started':start}
  write_json(folder/'execution.json',audit)
  reason=None
  gpu_samples=[];last_gpu_sample=0
  while p.poll() is None:
   if str(argv[0]).endswith('/gnina') and '--no_gpu' not in argv and time.time()-last_gpu_sample>1:
    last_gpu_sample=time.time()
    try:
     query=subprocess.run(['/usr/lib/wsl/lib/nvidia-smi','--query-compute-apps=pid,process_name,used_gpu_memory','--format=csv,noheader'],capture_output=True,text=True,timeout=3)
     selected=[line for line in query.stdout.splitlines() if line.split(',')[0].strip()==str(p.pid) or '/opt/sbvs/vendor/gnina' in line]
     gpu_samples.append({'time':time.time(),'process_matches':selected,'query_returncode':query.returncode})
    except (OSError,subprocess.TimeoutExpired) as e:gpu_samples.append({'error':str(e)})
   if cancel():reason='cancelled'
   elif time.time()-start>timeout:reason='timeout'
   if reason:
    os.killpg(p.pid,signal.SIGTERM)
    try:p.wait(5)
    except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
    break
   time.sleep(.15)
  audit.update(returncode=p.returncode,finished=time.time(),elapsed_seconds=time.time()-start,termination_reason=reason)
  write_json(folder/'execution.json',audit)
  if gpu_samples:write_json(folder/'gpu_samples.json',gpu_samples)
 if reason=='cancelled':raise Cancelled('Cancelled by user')
 if reason=='timeout':raise TimeoutError('Per-ligand timeout')
 if p.returncode:raise RuntimeError(f'Engine exited {p.returncode}: '+(folder/'stderr.log').read_text()[-1800:])
 return audit
