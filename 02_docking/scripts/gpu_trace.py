import os,subprocess,json,re
from pathlib import Path
from rdkit import Chem
from sbvs.paths import ROOT,read_json,write_json
from sbvs.engines.gnina import value
source=ROOT/read_json(ROOT/'reports/gnina_gpu_smoke.json')['run'];m=read_json(source/'run_manifest.json');entry=source/'entries'/next(iter(m['completed_entries']));audit=read_json(entry/'execution.json');argv=audit['argv'][:]
out=ROOT/'reports/gpu_trace';out.mkdir(exist_ok=True);argv[argv.index('-l')+1]=str(entry/'pose_001.sdf');argv[argv.index('-o')+1]=str(out/'scored.sdf');argv+=['--score_only']
env=os.environ.copy();env['LD_LIBRARY_PATH']=':'.join([str(p) for p in Path('/opt/sbvs/env/lib/python3.11/site-packages/nvidia').glob('*/lib')]+['/opt/sbvs/env/lib','/usr/lib/wsl/lib'])
with (out/'stdout.log').open('w') as stdout,(out/'stderr.log').open('w') as stderr:
 r=subprocess.run(['strace','-f','-yy','-e','trace=openat,ioctl','-o',str(out/'device.trace')]+argv,env=env,stdout=stdout,stderr=stderr,timeout=120)
assert r.returncode==0
mol=next(iter(Chem.SDMolSupplier(str(out/'scored.sdf'),removeHs=False)));scores={key:value(mol,[key],True) for key in ['CNNscore','CNNaffinity','minimizedAffinity']}
lines=[line for line in (out/'device.trace').read_text().splitlines() if '/dev/dxg' in line]
successful=[l for l in lines if 'ioctl(' in l and '= 0' in l]
write_json(ROOT/'reports/gpu_device_trace.json',{'status':'passed' if successful else 'not_tested','requested_device':0,'wsl_gpu_device':'/dev/dxg','successful_gpu_ioctls':len(successful),'trace':'reports/gpu_trace/device.trace','finite_scores':scores,'argv':argv});print(json.dumps({'gpu_ioctl_success':len(successful),'scores':scores}));assert successful
