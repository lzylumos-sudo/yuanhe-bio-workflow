"""Read-only preflight dependency/assets inventory; never starts a search or MD."""
import argparse,importlib.util,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'02_docking/src'))
sys.path.insert(0,str(ROOT/'00_workbench'))


def load_config(path=None):
 """Read the same effective settings as new workbench tasks, without writing."""
 from yuanhe_preflight.config import absolute_paths,validate
 if path is None:
  from webui import workflow
  cfg=workflow.defaults()['predocking']
  runtime=ROOT/'deployment/runtime.json'
  source=str(runtime) if runtime.is_file() else 'workbench research defaults'
 else:
  source=str(Path(path).resolve())
  cfg=json.loads(Path(path).read_text())
  # Accept a standalone preflight file, runtime/run configuration, or resolved run.
  if isinstance(cfg,dict):
   cfg=cfg.get('config',cfg)
  if isinstance(cfg,dict) and 'retro' not in cfg:
   research=cfg.get('research')
   cfg=research['predocking'] if isinstance(research,dict) and 'predocking' in research else cfg.get('predocking')
 return absolute_paths(validate(cfg),ROOT),source


def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--config',help='Preflight or workbench configuration; default: effective deployment/runtime.json settings')
 a=p.parse_args()
 import hashlib
 def file_hash(path):
  h=hashlib.sha256()
  with Path(path).open("rb") as handle:
   for chunk in iter(lambda:handle.read(1<<20),b""):h.update(chunk)
  return h.hexdigest()
 cfg,source=load_config(a.config)
 missing=[];assets=[]
 for package in ('rdkit','filelock','aizynthfinder'):
  if importlib.util.find_spec(package) is None:missing.append({'item':package,'effect':'local retrosynthesis unavailable' if package=='aizynthfinder' else 'preflight unavailable','environment':sys.executable})
 for key in ('model_path','templates_path','stock_path'):
  path=Path(cfg['retro'].get(key) or '')
  if not path.is_file():missing.append({'item':key,'effect':'retrosynthesis unavailable'})
  else:assets.append({'item':key,'path':str(path),'sha256':file_hash(path)})
 print(json.dumps({'missing':missing,'assets':assets,'config_source':source,'configured_retro_python':cfg['retro']['python'],
  'note':'Package check describes this doctor interpreter. Run it with the configured retro Python too. No model loading/search or external request performed.'},ensure_ascii=False,indent=2))
 return 1 if missing else 0
if __name__=='__main__':raise SystemExit(main())
