"""Resumable official downloads, atomically installed after validation."""
import hashlib,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(os.environ['SBVS_PROJECT_ROOT'])
SPECS={
'vina':('1.2.7','https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_linux_x86_64',None),
'gnina':('1.3.3','https://github.com/gnina/gnina/releases/download/v1.3.3/gnina.cuda12.8.static','3340c1f49cd3c7c84d8699182a1c6af13c7fa2a22448d1204640446106f72172')}
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
for name in sys.argv[1:]:
 version,url,expected=SPECS[name]; dst=Path('/opt/sbvs/vendor')/name; partial=dst.with_suffix('.part')
 if not dst.exists():
  if not (partial.exists() and expected and sha(partial)==expected):
   subprocess.run(['curl','-fL','--retry','3','--connect-timeout','20','--max-time','3600','-C','-',url,'-o',str(partial)],check=True)
  digest=sha(partial)
  if partial.stat().st_size<1000000 or (expected and digest!=expected):raise RuntimeError('Size/hash mismatch '+name)
  partial.replace(dst)
 digest=sha(dst)
 if expected and digest!=expected:raise RuntimeError('Installed hash mismatch')
 dst.chmod(0o755)
 info={'version':version,'url':url,'sha256':digest,'publisher_sha256_verified':bool(expected),'path':str(dst),'installed_at':time.time()}
 for flag in ['--version','--help']:
  env=os.environ.copy();env['LD_LIBRARY_PATH']=':'.join([str(p) for p in Path('/opt/sbvs/env/lib/python3.11/site-packages/nvidia').glob('*/lib')]+['/opt/sbvs/env/lib','/usr/lib/wsl/lib'])
  r=subprocess.run([str(dst),flag],capture_output=True,text=True,timeout=60,env=env)
  (ROOT/'logs/install'/f'{name}_{flag[2:]}.txt').write_text(r.stdout+r.stderr)
  info[flag[2:]+'_returncode']=r.returncode
 lock=ROOT/'env/vendor.lock.json'; data=json.loads(lock.read_text()) if lock.exists() else {};data[name]=info
 lock.write_text(json.dumps(data,indent=2));print(json.dumps(info),flush=True)
