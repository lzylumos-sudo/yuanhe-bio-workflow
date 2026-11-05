"""Small artifact primitives shared by standalone scientific modules and the UI."""
import hashlib,json,os,uuid
from pathlib import Path
from decimal import Decimal

def integration_steps(duration_ns,timestep_fs):
 """Exact units for the standalone MD engine, without importing the workbench."""
 duration=Decimal(str(duration_ns));step=Decimal(str(timestep_fs))
 if not duration.is_finite() or not step.is_finite() or duration<=0 or step<=0:
  raise ValueError('Duration and timestep must be finite and positive (ns, fs)')
 count=duration*Decimal(1000000)/step
 if count!=count.to_integral_value():raise ValueError('Duration is not an integer number of integration steps')
 return int(count)

def read(path,default=None):
 p=Path(path)
 return json.loads(p.read_text()) if p.is_file() else default

def atomic(path,value):
 p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
 t=p.with_name(p.name+'.'+uuid.uuid4().hex+'.tmp')
 t.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False));t.replace(p)

def digest(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
 return h.hexdigest()

def within(root,rel):
 p=(Path(root)/rel).resolve()
 if not p.is_relative_to(Path(root).resolve()) or not p.is_file():raise ValueError('Missing artifact or path escapes bundle')
 return p

def lease_fds():
 # Preserve all resource and output locks through nested engine processes.
 values=os.environ.get('YUANHE_COMPUTE_FDS',os.environ.get('YUANHE_COMPUTE_FD',''))
 fds=[]
 for value in values.split(','):
  try:
   fd=int(value);os.fstat(fd)
   if fd not in fds:fds.append(fd)
  except (ValueError,OSError):continue
 return tuple(fds)
