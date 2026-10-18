#!/usr/bin/env python3
"""Package an independently installable component with exact file checksums."""
import argparse,hashlib,json,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
MODULES={'docking':'02_docking','screening':'03_virtual_screening'}
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--module',choices=MODULES,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 if a.output.exists():raise ValueError('Refusing to overwrite an existing package')
 base=ROOT/MODULES[a.module];files={}
 for name in ['src','configs','env','docs','README.md','pyproject.toml']:
  source=base/name
  for f in ([source] if source.is_file() else source.rglob('*')):
   if f.is_file() and not f.is_symlink() and '__pycache__' not in f.parts and f.suffix not in {'.pyc','.pyo'}:files['module/'+str(f.relative_to(base))]=f
 for name in ['LICENSE','LICENSE_SCOPE.txt','THIRD_PARTY_NOTICES.md']:files[name]=ROOT/name
 manifest={n:hashlib.sha256(f.read_bytes()).hexdigest() for n,f in files.items()}
 a.output.parent.mkdir(parents=True,exist_ok=True)
 with zipfile.ZipFile(a.output,'x',zipfile.ZIP_DEFLATED) as z:
  for n,f in files.items():z.write(f,n)
  z.writestr('MANIFEST.json',json.dumps({'schema':'yuanhe.component/1','module':a.module,'files':manifest},indent=2))
  z.writestr('README.txt',f'Yuanhe {a.module} independent source package.\nInstall module/ into a compatible Python environment. Read its environment and usage documentation.\nGNINA and model binaries are external dependencies. Set explicit input, output and engine paths.\nNo run data, credentials or weights are bundled. Archive validation is not a new engine execution.\n')
 print(json.dumps({'module':a.module,'files':len(files),'sha256':hashlib.sha256(a.output.read_bytes()).hexdigest()}))
if __name__=='__main__':main()
