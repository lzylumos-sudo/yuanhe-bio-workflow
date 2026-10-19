import subprocess,sys,os
from filelock import FileLock
from sbvs.paths import ROOT,write_json
# Represents the exact lock held by real CLI/worker computations; does not alter any process.
with FileLock('/opt/sbvs/state/compute.lock'):
 r=subprocess.run([sys.executable,str(ROOT/'scripts/service_control.py'),'stop'],capture_output=True,text=True)
 assert r.returncode!=0 and 'standalone project computation' in r.stderr
write_json(ROOT/'reports/stop_guard.json',{'status':'passed','note':'Stop refuses while the shared computation lock is held; no process signalled'})
print('Active computation stop guard passed')
