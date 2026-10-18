"""Standalone MD module: immutable handoff -> build -> run -> analyze.

Copy the scientific package and a handoff bundle into another workflow. The
same engine is used by Research Studio; a web server is not required.
"""
from pathlib import Path
import argparse,json,shutil,sys,os,fcntl,time
from . import artifacts as a
SCHEMA='yuanhe.md-input/1'

def validate(root):
 root=Path(root).resolve();m=a.read(root/'handoff.json')
 if not m or m.get('schema')!=SCHEMA:raise ValueError('Unsupported MD input schema')
 for key in ('receptor','ligand'):
  row=m['artifacts'][key];p=a.within(root,row['path'])
  if a.digest(p)!=row['sha256']:raise ValueError('Input checksum mismatch: '+key)
 for key in ('molecule_id','compound_id','pose_id','smiles','qualification'):
  if not m['candidate'].get(key):raise ValueError('Missing candidate identity: '+key)
 if m['candidate']['qualification']!='formal' and not m['parameters'].get('exploratory_ack'):raise ValueError('Exploratory scope must be explicit')
 return m

def export_bundle(target,candidate,receptor,ligand,parameters,execution,provenance):
 target=Path(target)
 if target.exists():raise ValueError('Refusing to overwrite an existing handoff')
 target.mkdir(parents=True);inputs=target/'inputs';inputs.mkdir()
 artifacts={}
 for key,source in [('receptor',receptor),('ligand',ligand)]:
  name=key+Path(source).suffix;shutil.copyfile(source,inputs/name)
  artifacts[key]={'path':'inputs/'+name,'sha256':a.digest(inputs/name)}
 m={'schema':SCHEMA,'candidate':candidate,'artifacts':artifacts,'parameters':parameters,
    'execution':execution,'provenance':provenance,'scientific_claim':'exploratory simulation input; qualification unchanged'}
 a.atomic(target/'handoff.json',m);validate(target)
 return target

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['validate','build','run','analyze']);parser.add_argument('--bundle',required=True,type=Path);parser.add_argument('--system-id');parser.add_argument('--review-system',action='store_true');args=parser.parse_args();root=args.bundle.resolve();m=validate(root)
 if args.action=='validate':print(json.dumps({'status':'valid','schema':SCHEMA,'candidate':m['candidate']}));return
 from . import membrane
 with (root/'.compute.lock').open('a+') as lock:
  try:fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
  except BlockingIOError:raise RuntimeError('This bundle already has an active computation')
  os.environ['YUANHE_COMPUTE_FD']=str(lock.fileno())
  c=m['execution'];out=root/'operations'/(args.action+'_'+str(time.time_ns()));out.mkdir(parents=True,exist_ok=True)
  if args.action=='build':
   result=membrane.build_prepared(root,c,m['parameters'],out,m['candidate'],a.within(root,m['artifacts']['ligand']['path']),a.within(root,m['artifacts']['receptor']['path']),m['provenance'].get('source_revision'))
  else:
   if not args.system_id:raise ValueError('--system-id is required')
   result=getattr(membrane,args.action)(root,c,{'system_id':args.system_id,'system_reviewed':args.review_system},out)
  a.atomic(out/'result.json',result);print(json.dumps(result))
if __name__=='__main__':main()
