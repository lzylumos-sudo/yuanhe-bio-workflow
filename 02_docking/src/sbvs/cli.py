import argparse,json,sys,os,subprocess
from pathlib import Path
from .paths import ROOT,rootpath,read_json,write_json,validate_runtime
from .config import load
from . import runner
def main():
 p=argparse.ArgumentParser(prog='sbvs');sub=p.add_subparsers(dest='command',required=True)
 sub.add_parser('doctor')
 s=sub.add_parser('serve');s.add_argument('--host',default='127.0.0.1');s.add_argument('--port',type=int,default=5902)
 s=sub.add_parser('demo');s.add_argument('--preset',default='class_quick')
 s=sub.add_parser('prepare-receptor');s.add_argument('--input',required=True);s.add_argument('--chains',default='A');s.add_argument('--out',required=True);s.add_argument('--policy');s.add_argument('--classroom',action='store_true')
 s=sub.add_parser('prepare-ligands');s.add_argument('--input',required=True);s.add_argument('--out',required=True)
 s=sub.add_parser('make-box');s.add_argument('--reference',required=True);s.add_argument('--padding',type=float,default=5);s.add_argument('--out',required=True)
 s=sub.add_parser('screen');s.add_argument('--config',required=True);s.add_argument('--ligands',required=True);s.add_argument('--out')
 for name in ('resume','report'):
  s=sub.add_parser(name);s.add_argument('--run',required=True)
 s=sub.add_parser('export-md');s.add_argument('--run',required=True);s.add_argument('--selection',required=True)
 s=sub.add_parser('rescore');s.add_argument('--run',required=True);s.add_argument('--engine',choices=['vina'],default='vina');s.add_argument('--operation',choices=['redock'],required=True);s.add_argument('--top-k',type=int,default=10);s.add_argument('--device',choices=['cpu'],default='cpu')
 s=sub.add_parser('redock-reference');s.add_argument('--preset',default='class_full')
 args=p.parse_args();validate_runtime()
 if args.command=='doctor':
  from .doctor import doctor
  print(json.dumps(doctor(),indent=2));return 0
 if args.command=='serve':
  from .web.app import make_app
  from waitress import serve
  from werkzeug.middleware.dispatcher import DispatcherMiddleware
  app=make_app();base=os.environ.get('APP_BASE_PATH','').rstrip('/')
  if base:
   if not base.startswith('/'):raise ValueError('APP_BASE_PATH must start /')
   app=DispatcherMiddleware(lambda e,s:(s('404 Not Found',[]),[b'Not Found'])[1],{base:app})
  write_json(Path('/opt/sbvs/state/server.json'),{'pid':os.getpid(),'project_uuid':read_json(ROOT/'project.json')['uuid'],'root':str(ROOT)})
  serve(app,host=args.host,port=args.port,threads=6);return 0
 if args.command=='prepare-receptor':
  from .preparation.receptor import prepare
  print(json.dumps(prepare(rootpath(args.input),rootpath(args.out),args.chains.split(','),read_json(rootpath(args.policy)) if args.policy else None,args.classroom),indent=2));return 0
 if args.command=='prepare-ligands':
  from .preparation.ligands import inputs,prepare
  from .ranking import csvwrite
  out=rootpath(args.out);out.mkdir(parents=True,exist_ok=True);rows=[]
  for item in inputs(rootpath(args.input)):
   try:rows.append(dict(prepare(item,out/f'P{item["index"]:06d}'),status='success'))
   except Exception as e:rows.append({'name':item['name'],'status':'failed','error':str(e)})
  csvwrite(out/'ligand_manifest.csv',rows);print(json.dumps(rows,ensure_ascii=False));return int(all(r['status']=='failed' for r in rows))
 if args.command=='make-box':
  from .preparation.box import make_box
  print(json.dumps(make_box(rootpath(args.reference),args.padding,rootpath(args.out))));return 0
 if args.command in ('screen','demo'):
  cfg=load(args.config if args.command=='screen' else args.preset)
  source=rootpath(args.ligands) if args.command=='screen' else ROOT/'data/demo/ligands_20.smi'
  out=runner.create(source,cfg,args.out if args.command=='screen' else None);m=runner.run(out);print(json.dumps({'path':str(out),'status':m['status'],'counts':m['counts']}));return 0 if m['status']=='completed' else 2
 if args.command=='resume':
  m=runner.run(rootpath(args.run),True);print(json.dumps(m['counts']));return 0 if m['status']=='completed' else 2
 if args.command=='report':
  from .export import report
  print(report(args.run));return 0
 if args.command=='export-md':
  from .export import handoff
  print(handoff(args.run,args.selection));return 0
 if args.command=='rescore':
  from .stages import vina_stage
  out=vina_stage(args.run,args.operation,args.top_k,args.device);result=runner.run(out);print(out);return 0 if result['status']=='completed' else 2
 if args.command=='redock-reference':
  from rdkit import Chem
  from .qc import fixed_frame_rmsd
  from .export import read_rows
  reference=next(iter(Chem.SDMolSupplier(str(ROOT/'data/reference/FOL_crystal.sdf'),removeHs=False)))
  source=ROOT/'data/reference/FOL_redock.smi';source.write_text(Chem.MolToSmiles(reference)+' FOL_redock\n')
  out=runner.create(source,load(args.preset));result=runner.run(out);values=[]
  for row in read_rows(out/'poses.csv'):
   mol=next(iter(Chem.SDMolSupplier(str(out/row['pose_file']),removeHs=False)));values.append({'pose_id':row['pose_id'],'rmsd_A':fixed_frame_rmsd(mol,reference)})
  audit={'run':str(out.relative_to(ROOT)),'method':'symmetry-aware CalcRMS without ligand alignment','values':values,'top1':values[0] if values else None,'best_top3':min(values[:3],key=lambda r:r['rmsd_A']) if values else None,'best_all':min(values,key=lambda r:r['rmsd_A']) if values else None,'reference_threshold_A':2,'enrichment':'not_evaluated_no_labels'};write_json(ROOT/'reports/FOL_redocking.json',audit);print(json.dumps(audit));return 0 if result['status']=='completed' else 2
if __name__=='__main__':
 try:sys.exit(main())
 except Exception as e:
  print(f'{type(e).__name__}: {e}',file=sys.stderr);sys.exit(1)
