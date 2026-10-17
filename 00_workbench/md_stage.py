"""Run selected prepared membrane systems and analysis from pipeline stage 6."""
import argparse
from pathlib import Path
from webui.studio import config
from webui import workflow as w
from scientific import membrane

def main():
 p=argparse.ArgumentParser();p.add_argument('--run-id',required=True);p.add_argument('--output',required=True);a=p.parse_args()
 d,c=config(a.run_id);out=Path(a.output).resolve()
 if out!= (d/'gromacs').resolve():raise ValueError('Invalid output path')
 if not c.get('gromacs_system_reviewed'):raise ValueError('System review required')
 ids=c.get('md_system_ids',[])
 if not ids:raise ValueError('Select a system in the membrane wizard first')
 records=[]
 for ident in ids:
  record=membrane.run(d,c,{'system_id':ident,'system_reviewed':True},out)
  analysis=membrane.analyze(d,c,{'system_id':ident},out)
  folder=membrane.system_dir(d,ident)
  records.append({'system_id':ident,'compound_id':analysis['compound_id'],'files':{str(x.relative_to(d)):w.digest(x) for x in folder.rglob('*') if x.is_file()},'analysis':analysis})
  w.atomic(out/'systems.json',records)
 return 0
if __name__=='__main__':raise SystemExit(main())
