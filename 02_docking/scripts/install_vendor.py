"""Resumable official downloads, atomically installed after validation."""
import hashlib,json,os,subprocess,sys,time,platform
from pathlib import Path
ROOT=Path(os.environ['SBVS_PROJECT_ROOT'])
if platform.system()!='Linux' or platform.machine() not in {'x86_64','amd64'}:raise SystemExit('This official download installer is Linux x86-64 only; use the pinned platform environment otherwise')
SPECS={
'vina':('1.2.7','https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_linux_x86_64','f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644'),
}

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
for name in sys.argv[1:]:
 version,url,expected=SPECS[name]; dst=Path(os.environ.get('YUANHE_VENDOR_DIR','/opt/sbvs/vendor'))/name; dst.parent.mkdir(parents=True,exist_ok=True); partial=dst.with_suffix('.part')
 if not dst.exists():
  if not (partial.exists() and expected and sha(partial)==expected):
   subprocess.run(['curl','-fL','--retry','3','--connect-timeout','20','--max-time','3600','-C','-',url,'-o',str(partial)],check=True)
  digest=sha(partial)
  if partial.stat().st_size<1000000 or (expected and digest!=expected):raise RuntimeError('Size/hash mismatch '+name)
  partial.replace(dst)
 digest=sha(dst)
 if expected and digest!=expected:raise RuntimeError('Installed hash mismatch')
 dst.chmod(0o755)
 info={'version':version,'url':url,'sha256':digest,'repository_sha256_verified':bool(expected),'path':str(dst),'installed_at':time.time()}
 (ROOT/'logs/install').mkdir(parents=True,exist_ok=True)
 for flag in ['--version','--help']:
  env=os.environ.copy()
  r=subprocess.run([str(dst),flag],capture_output=True,text=True,timeout=60,env=env)
  if r.returncode or (flag=='--version' and 'AutoDock Vina v1.2.7' not in r.stdout+r.stderr):raise RuntimeError('Official Vina version verification failed')
  (ROOT/'logs/install'/f'{name}_{flag[2:]}.txt').write_text(r.stdout+r.stderr)
  info[flag[2:]+'_returncode']=r.returncode
 lock=Path(os.environ.get('YUANHE_VENDOR_RECEIPT',str(ROOT/'env/vendor.lock.json')))
 lock.parent.mkdir(parents=True,exist_ok=True)
 data=json.loads(lock.read_text()) if lock.exists() else {};data[name]=info
 lock.write_text(json.dumps(data,indent=2));print(json.dumps(info),flush=True)
