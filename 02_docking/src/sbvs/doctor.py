import os,platform,subprocess,importlib.metadata
from pathlib import Path
from .paths import ROOT,project,write_json,read_json
def doctor():
 packages={}
 for p in ('rdkit','meeko','biopython','numpy','scipy','flask','playwright'):
  try:packages[p]=importlib.metadata.version(p)
  except importlib.metadata.PackageNotFoundError:packages[p]=None
 engines={}
 for name in ('vina','gnina'):
  exe=Path('/opt/sbvs/vendor')/name
  if exe.exists():
   env=os.environ.copy();env['LD_LIBRARY_PATH']=':'.join([str(p) for p in Path('/opt/sbvs/env/lib/python3.11/site-packages/nvidia').glob('*/lib')]+['/opt/sbvs/env/lib','/usr/lib/wsl/lib'])
   r=subprocess.run([str(exe),'--version'],capture_output=True,text=True,timeout=45,env=env);engines[name]={'version_output':r.stdout+r.stderr,'returncode':r.returncode}
  else:engines[name]={'installed':False}
 try:
  r=subprocess.run(['/usr/lib/wsl/lib/nvidia-smi','--query-gpu=name,driver_version,memory.total','--format=csv'],capture_output=True,text=True,timeout=15);gpu={'returncode':r.returncode,'output':r.stdout+r.stderr}
 except OSError as e:gpu={'error':str(e)}
 result={'project_uuid':project()['uuid'],'python':platform.python_version(),'packages':packages,'cpu_count':os.cpu_count(),'cpu_available':len(os.sched_getaffinity(0)),'engines':engines,'gpu_detection':gpu,'note':'GPU detection is not GNINA CNN validation','receptor_ready':(ROOT/'data/prepared/demo_1drf/receptor.pdbqt').exists()}
 previous=ROOT/'reports/doctor.json'
 if previous.exists():
  old=read_json(previous)
  if old.get('gpu_detection')!=gpu:
   caps=ROOT/'reports/gnina_capabilities.json'
   if caps.exists():
    data=read_json(caps);data.update(gpu_passed=False,cpu_passed=False,revalidation_required='Hardware detection changed; run both GNINA smoke tests');write_json(caps,data)
 write_json(ROOT/'reports/doctor.json',result);return result
