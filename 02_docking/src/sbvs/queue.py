import sqlite3,threading,time,traceback
from pathlib import Path
from filelock import FileLock
from .paths import ROOT,read_json,write_json
from . import runner
class Queue:
 def __init__(self):
  self.lock=FileLock('/opt/sbvs/state/worker.lock');self.lock.acquire(timeout=0)
  self.db=sqlite3.connect('/opt/sbvs/state/jobs.sqlite',check_same_thread=False);self.mutex=threading.RLock();self.accepting=True;self.stopping=False
  self.db.execute('CREATE TABLE IF NOT EXISTS jobs (id TEXT PRIMARY KEY,path TEXT,status TEXT,cancel INTEGER DEFAULT 0,resume INTEGER DEFAULT 0)')
  for manifest in (ROOT/'runs').glob('*/run_manifest.json'):
   m=read_json(manifest);status=m['status']
   if status in ('running','preparing'):status='interrupted'
   self.db.execute('INSERT OR IGNORE INTO jobs VALUES(?,?,?,?,?)',(m['run_id'],str(manifest.parent.relative_to(ROOT)),status,0,0))
  self.db.commit()
  # Exclusive process lock proves the prior scheduler is no longer running.
  for rid,path in self.db.execute("SELECT id,path FROM jobs WHERE status IN ('preparing','running')").fetchall():
   import psutil,os,signal
   out=ROOT/path
   for audit in (out/'entries').glob('*/execution.json'):
    info=read_json(audit)
    if 'finished' in info:continue
    try:
     proc=psutil.Process(info['pid'])
     if proc.cmdline()==info['argv'] and abs(proc.create_time()-info['started'])<5 and proc.exe().startswith('/opt/sbvs/vendor/'):
      os.killpg(proc.pid,signal.SIGTERM)
      try:proc.wait(5)
      except psutil.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL)
    except psutil.NoSuchProcess:pass
   m=read_json(out/'run_manifest.json');m['status']='interrupted';write_json(out/'run_manifest.json',m)
  self.db.execute("UPDATE jobs SET status='interrupted' WHERE status IN ('preparing','running')");self.db.commit()
  self.thread=threading.Thread(target=self.work,daemon=True);self.thread.start()
 def add(self,path,resume=False):
  if not self.accepting:raise ValueError('Service draining')
  m=read_json(path/'run_manifest.json')
  with self.mutex:
   self.db.execute('INSERT OR REPLACE INTO jobs VALUES (?,?,?,?,?)',(m['run_id'],str(path.relative_to(ROOT)),'queued',0,int(resume)));self.db.commit()
  return m['run_id']
 def get(self,rid):
  with self.mutex:r=self.db.execute('SELECT * FROM jobs WHERE id=?',(rid,)).fetchone()
  if not r:raise KeyError(rid)
  return dict(zip(['id','path','status','cancel','resume'],r))
 def all(self):
  with self.mutex:ids=[r[0] for r in self.db.execute('SELECT id FROM jobs ORDER BY rowid DESC LIMIT 200')]
  return [self.get(rid) for rid in ids]
 def update(self,rid,status):
  with self.mutex:self.db.execute('UPDATE jobs SET status=? WHERE id=?',(status,rid));self.db.commit()
 def cancel(self,rid):
  with self.mutex:self.db.execute('UPDATE jobs SET cancel=1 WHERE id=?',(rid,));self.db.commit()
 def work(self):
  while not self.stopping:
   with self.mutex:r=self.db.execute("SELECT id,path,resume,cancel FROM jobs WHERE status='queued' ORDER BY rowid LIMIT 1").fetchone()
   if not r:time.sleep(.3);continue
   rid,path,resume,cancel=r
   if cancel:
    m=read_json(ROOT/path/'run_manifest.json');m['status']='cancelled';write_json(ROOT/path/'run_manifest.json',m);self.update(rid,'cancelled');continue
   self.update(rid,'preparing')
   try:runner.run(ROOT/path,bool(resume),lambda:bool(self.get(rid)['cancel']),lambda m:self.update(rid,m['status']))
   except Exception:
    (ROOT/path/'worker_error.log').write_text(traceback.format_exc());self.update(rid,'failed')
 def drain(self):
  self.accepting=False
  return not any(x['status'] in ('queued','preparing','running') for x in self.all())
