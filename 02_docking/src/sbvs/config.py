import os,math,yaml
from .paths import ROOT,rootpath,sha
def load(value='class_quick'):
 path=ROOT/'configs'/f'{value}.yaml'
 if not path.exists():path=rootpath(value)
 cfg=yaml.safe_load(path.read_text(encoding='utf-8'))
 if str(value) in ('class_quick','class_full','research_gnina','research_refine'):cfg['cpu_per_job']=min(int(cfg['cpu_per_job']),len(os.sched_getaffinity(0)))
 validate(cfg);return cfg
def validate(c):
 if c['engine'] not in ('vina','gnina'):raise ValueError('Unknown engine')
 if c.get('ligand_states','preserve') not in ('preserve','enumerate'):raise ValueError('Unknown state policy')
 if c.get('ligand_states')=='enumerate':
  if not 0<=float(c.get('ph_low',6.4))<=float(c.get('ph_high',8.4))<=14:raise ValueError('Invalid pH range')
  if not 1<=int(c.get('max_states',8))<=32:raise ValueError('max_states must be 1..32')
 if c.get('workers',1)!=1:raise ValueError('This scheduler uses one controlled worker')
 if not 1<=int(c['cpu_per_job'])<=len(os.sched_getaffinity(0)):raise ValueError('CPU budget exceeded')
 if not c['seeds'] or any(not isinstance(s,int) or not 0<s<2147483647 for s in c['seeds']):raise ValueError('Invalid seeds')
 for key in ('exhaustiveness','num_modes','per_ligand_timeout_seconds'):
  if not isinstance(c[key],(int,float)) or not math.isfinite(c[key]) or c[key]<=0:raise ValueError('Invalid '+key)
 for key in ('receptor','receptor_pdb','box'):
  if not rootpath(c[key]).is_file():raise ValueError('Missing '+key)
 box=validate_box(__import__('json').loads(rootpath(c['box']).read_text()))
 if box.get('receptor_sha256') and box['receptor_sha256']!=sha(rootpath(c['receptor'])):raise ValueError('Box belongs to a different receptor')
def validate_box(b):
 for k in ('center','size'):
  if len(b[k])!=3 or any(not math.isfinite(float(v)) for v in b[k]):raise ValueError('Invalid box')
 if any(v<=0 or v>126 for v in b['size']):raise ValueError('Box size must be in (0,126] angstrom')
 return b
