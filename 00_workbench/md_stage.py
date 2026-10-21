"""Run selected prepared membrane systems and analysis from pipeline stage 6."""
import argparse,json,shutil,subprocess,sys,os
from pathlib import Path
from webui.studio import config
from webui import workflow as w
from scientific import membrane

def main():
 p=argparse.ArgumentParser();p.add_argument('--run-id',required=True);p.add_argument('--output',required=True);a=p.parse_args()
 d,c=config(a.run_id);out=Path(a.output).resolve()
 if out!= (d/'gromacs').resolve():raise ValueError('Invalid output path')
 if c.get('research'):
  from research.md_queue import run
  return run(d,c,out)
 if not c.get('gromacs_system_reviewed'):raise ValueError('System review required')
 ids=c.get('md_system_ids',[])
 if not ids:raise ValueError('Select a system in the membrane wizard first')
 records=[]
 for ident in ids:
  # Exercise the independently exportable engine in a child interpreter.
  # It has no webui import path and uses explicit run/config/system inputs.
  package=out/'standalone_md';package.mkdir(exist_ok=True)
  for name in ['__init__.py','artifacts.py','membrane.py','md_cli.py']:
   target=package/'scientific'/name;target.parent.mkdir(exist_ok=True)
   source=w.UI/'scientific'/name
   if target.exists() and w.digest(target)!=w.digest(source):raise ValueError('Standalone MD source changed')
   if not target.exists():shutil.copy2(source,target)
  request=package/(ident+'.request.json');response=package/(ident+'.result.json')
  w.atomic(request,{'run_dir':str(d),'config':c,'system_id':ident,'response':str(response),'source_files':{str(x.relative_to(package)):w.digest(x) for x in (package/'scientific').glob('*.py')}})
  driver=package/'run.py'
  driver.write_text('import json,sys\nfrom pathlib import Path\nsys.path.insert(0,str(Path(__file__).resolve().parent))\nfrom scientific import membrane,artifacts\nr=json.loads(Path(sys.argv[1]).read_text());d=Path(r["run_dir"]);out=Path(r["response"]).parent\nmembrane.run(d,r["config"],{"system_id":r["system_id"],"system_reviewed":True},out)\nresult=membrane.analyze(d,r["config"],{"system_id":r["system_id"]},out)\nartifacts.atomic(Path(r["response"]),result)\n')
  from scientific.artifacts import lease_fds
  subprocess.run([sys.executable,'-I',str(driver),str(request)],check=True,pass_fds=lease_fds(),env={**os.environ,'PYTHONPATH':''})
  analysis=w.read(response)
  folder=membrane.system_dir(d,ident)
  records.append({'system_id':ident,'compound_id':analysis['compound_id'],'files':{str(x.relative_to(d)):w.digest(x) for x in folder.rglob('*') if x.is_file()},'analysis':analysis})
  w.atomic(out/'systems.json',records)
 return 0
if __name__=='__main__':raise SystemExit(main())
