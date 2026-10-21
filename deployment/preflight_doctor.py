"""Read-only preflight dependency/assets inventory; never starts a search or MD."""
import argparse,importlib.util,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'02_docking/src'))
def main():
 p=argparse.ArgumentParser();p.add_argument('--config');a=p.parse_args()
 from yuanhe_preflight.config import defaults,validate
 from yuanhe_preflight.common import read,file_hash
 cfg=validate(read(a.config)) if a.config else defaults(ROOT)
 missing=[];assets=[]
 for package in ('rdkit','filelock','aizynthfinder'):
  if importlib.util.find_spec(package) is None:missing.append({'item':package,'effect':'local retrosynthesis unavailable' if package=='aizynthfinder' else 'preflight unavailable','environment':sys.executable})
 if not cfg['catalogs']:missing.append({'item':'supplier_catalogs','effect':'no purchase evidence'})
 if not cfg['stock_catalogs']:missing.append({'item':'trusted_terminal_supplier_catalogs','effect':'routes cannot qualify'})
 for key in ('model_path','templates_path','stock_path'):
  path=Path(cfg['retro'].get(key) or '')
  if not path.is_file():missing.append({'item':key,'effect':'retrosynthesis unavailable'})
  else:assets.append({'item':key,'path':str(path),'sha256':file_hash(path)})
 print(json.dumps({'missing':missing,'assets':assets,'configured_retro_python':cfg['retro']['python'],
  'note':'Package check describes this doctor interpreter. Run it with the configured retro Python too. No model loading/search or external request performed.'},ensure_ascii=False,indent=2))
 return 1 if missing else 0
if __name__=='__main__':raise SystemExit(main())
