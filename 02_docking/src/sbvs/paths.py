import os,json,hashlib
from pathlib import Path
ROOT=Path(os.environ.get('SBVS_PROJECT_ROOT',Path(__file__).resolve().parents[2])).resolve()
def rootpath(value):
 p=Path(value); p=(p if p.is_absolute() else ROOT/p).resolve()
 if not p.is_relative_to(ROOT):raise ValueError('Path outside project')
 return p
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def write_json(path,data):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8');temp.replace(path)
def read_json(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def project():return read_json(ROOT/'project.json')
def validate_runtime():
 if Path('/opt/sbvs/project_uuid').read_text().strip()!=project()['uuid']:raise RuntimeError('Project UUID differs from dedicated runtime')
