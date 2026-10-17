import os,json,shutil,subprocess,sys
from pathlib import Path
from sbvs.paths import ROOT,write_json
target=ROOT/'runtime/path tests/English Project'
target.mkdir(parents=True,exist_ok=True)
for name in ('src','configs','data','env','scripts'):
 shutil.copytree(ROOT/name,target/name,dirs_exist_ok=True)
for name in ('project.json','run.sh'):shutil.copyfile(ROOT/name,target/name)
for name in ('runs','logs','reports'): (target/name).mkdir(exist_ok=True)
env=os.environ.copy();env['SBVS_PROJECT_ROOT']=str(target);env['PYTHONPATH']=str(target/'src')
r=subprocess.run(['bash',str(target/'run.sh'),'screen','--config','class_quick','--ligands','data/demo/cco.smi'],env=env,capture_output=True,text=True,timeout=180)
audit={'status':'passed' if r.returncode==0 else 'failed','path':str(target.relative_to(ROOT)),'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr,'scope':'same WSL runtime, relocated source/data with spaces and Chinese; no VHDX copy and no second machine test'}
write_json(ROOT/'reports/path_acceptance.json',audit);print(json.dumps(audit));assert r.returncode==0
