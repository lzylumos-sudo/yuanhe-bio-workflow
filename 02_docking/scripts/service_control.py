import argparse,os,sys,subprocess,json,time,signal
from pathlib import Path
ROOT=Path(os.environ['SBVS_PROJECT_ROOT']);sys.path.insert(0,str(ROOT/'src'))
from sbvs.paths import read_json,project,write_json
import psutil,requests
from filelock import FileLock,Timeout
def ensure_no_compute():
 lock=FileLock('/opt/sbvs/state/compute.lock')
 try:lock.acquire(timeout=0)
 except Timeout:raise RuntimeError('A standalone project computation is active; finish it or explicitly cancel it before migration')
 finally:
  if lock.is_locked:lock.release()
p=argparse.ArgumentParser();p.add_argument('action',choices=['start','stop']);p.add_argument('--port',type=int,default=5902);p.add_argument('--cancel-active',action='store_true');args=p.parse_args();state=Path('/opt/sbvs/state/server.json')
if args.action=='start':
 with (ROOT/'logs/web.log').open('a') as log:
  proc=subprocess.Popen(['bash',str(ROOT/'run.sh'),'serve','--port',str(args.port)],stdout=log,stderr=log,start_new_session=True,env=os.environ.copy())
 print('Started project service PID',proc.pid)
else:
 if not state.exists():ensure_no_compute();print('No service state');sys.exit(0)
 info=read_json(state)
 if info['project_uuid']!=project()['uuid']:raise RuntimeError('Wrong project service identity')
 if not psutil.pid_exists(info['pid']):ensure_no_compute();print('Service already stopped');sys.exit(0)
 proc=psutil.Process(info['pid']);cmd=proc.cmdline()
 if 'sbvs.cli' not in cmd or 'serve' not in cmd or not proc.exe().startswith('/opt/sbvs/env/'):raise RuntimeError('PID identity mismatch; refusing signal')
 port=int(cmd[cmd.index('--port')+1]) if '--port' in cmd else 5902
 url=f'http://127.0.0.1:{port}';headers={'X-SBVS-Request':'local-ui'}
 h=requests.get(url+'/api/health',timeout=5).json()
 if h['project_uuid']!=project()['uuid']:raise RuntimeError('Service UUID mismatch')
 idle=requests.post(url+'/api/drain',headers=headers,timeout=5).json()['idle']
 if not idle and args.cancel_active:
  import sqlite3
  db=sqlite3.connect('/opt/sbvs/state/jobs.sqlite')
  ids=[r[0] for r in db.execute("SELECT id FROM jobs WHERE status IN ('queued','preparing','running')")];db.close()
  for rid in ids:requests.post(url+'/api/jobs/'+rid+'/cancel',headers=headers,timeout=5).raise_for_status()
  for _ in range(60):
   if requests.post(url+'/api/drain',headers=headers,timeout=5).json()['idle']:idle=True;break
   time.sleep(1)
 if not idle:raise RuntimeError('Tasks active; service draining. Explicit cancellation is required to stop now')
 ensure_no_compute()
 proc.terminate()
 try:proc.wait(15)
 except psutil.TimeoutExpired:raise RuntimeError('Server did not stop; project not terminated')
 state.unlink(missing_ok=True);print('Project service stopped')
