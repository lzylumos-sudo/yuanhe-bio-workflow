import os,time,json,subprocess,sys
from pathlib import Path
import requests,psutil
from sbvs.paths import ROOT,read_json,write_json
BASE='http://127.0.0.1:5902';HEAD={'X-SBVS-Request':'local-ui'};results={}
def submit(text,preset='class_full'):
 r=requests.post(BASE+'/api/jobs',headers=HEAD,data={'smiles':text,'preset':preset},timeout=10);r.raise_for_status();return r.json()['id']
def get(rid):return requests.get(BASE+'/api/jobs/'+rid,timeout=10).json()
def wait(rid,limit=180):
 start=time.time()
 while time.time()-start<limit:
  d=get(rid)
  if d['status'] in ('completed','partial_failed','failed','cancelled','interrupted'):return d
  time.sleep(.2)
 raise TimeoutError(rid)
# Partial failure, duplicate names, Unicode and blank validation.
r=requests.post(BASE+'/api/jobs',headers=HEAD,data={'smiles':'  '},timeout=10);assert r.status_code==400
rid=submit('CCO 重复 名字\nINVALID bad\nCCCO 重复 名字\n[Fe] unsupported','class_quick');d=wait(rid)
assert d['status']=='partial_failed' and len(d['ranking'])==2 and len(d['failures'])==2
results['mixed_batch']={'status':'passed','run_id':rid,'successful':2,'failed':2}
# Actual cancellation, not a mocked process.
rid=submit((ROOT/'data/demo/ligands_20.smi').read_text());start=time.time()
while time.time()-start<40:
 if get(rid)['status']=='running':break
 time.sleep(.1)
requests.post(BASE+'/api/jobs/'+rid+'/cancel',headers=HEAD,timeout=10).raise_for_status();d=wait(rid);assert d['status']=='cancelled'
results['cancel']={'status':'passed','run_id':rid}
# Hard service interruption after at least one durable entry.
rid=submit((ROOT/'data/demo/ligands_20.smi').read_text());start=time.time()
while time.time()-start<60:
 d=get(rid)
 if d['counts']['successful_entries']>=1 and d['status']=='running':break
 time.sleep(.1)
else:raise RuntimeError('Unable to catch a real active task')
manifest=read_json(ROOT/'runs'/rid/'run_manifest.json');before=dict(manifest['completed_entries'])
server=read_json('/opt/sbvs/state/server.json');proc=psutil.Process(server['pid'])
assert 'sbvs.cli' in proc.cmdline() and 'serve' in proc.cmdline() and proc.exe().startswith('/opt/sbvs/env/')
proc.kill();proc.wait(10)
subprocess.run([sys.executable,str(ROOT/'scripts/service_control.py'),'start'],check=True)
for _ in range(100):
 try:
  d=get(rid)
  if d['status']=='interrupted':break
 except requests.RequestException:pass
 time.sleep(.2)
assert d['status']=='interrupted'
requests.post(BASE+'/api/jobs/'+rid+'/resume',headers=HEAD,timeout=10).raise_for_status();d=wait(rid,180);assert d['status']=='completed' and len(d['ranking'])==20
after=read_json(ROOT/'runs'/rid/'run_manifest.json')['completed_entries']
assert all(after[k]==v for k,v in before.items())
results['interruption_resume']={'status':'passed','run_id':rid,'verified_reused_entries':len(before),'final_successful':len(d['ranking'])}
write_json(ROOT/'reports/recovery_acceptance.json',results);print(json.dumps(results,indent=2))
