#!/usr/bin/env python3
"""Export the shared MD core without the UI, histories, data or credentials."""
import argparse,hashlib,json,zipfile
from pathlib import Path
from package_licenses import license_files
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
 if a.output.exists():raise ValueError('Refusing to overwrite an existing package')
 a.output.parent.mkdir(parents=True,exist_ok=True)
 files={f'scientific/{n}':ROOT/'00_workbench/scientific'/n for n in ['__init__.py','artifacts.py','membrane.py','md_cli.py']}
 for n in ['deployment/membrane-environment.yml','docs/07_MODULE_CONTRACTS.md']:files[n]=ROOT/n
 for f in (ROOT/'05_molecular_dynamics/modules').rglob('*.py'):
  files[str(f.relative_to(ROOT))]=f
 files['docs/15_VINA_RESP_RESEARCH.md']=ROOT/'docs/15_VINA_RESP_RESEARCH.md'
 for f in (ROOT/'00_workbench/yuanhe').glob('*.py'):files['yuanhe/'+f.name]=f
 for n in ['VERSION','configs/module_registry.yaml','configs/protocol_registry.yaml']:files[n]=ROOT/n
 for f in (ROOT/'configs/schemas').glob('*.json'):files[str(f.relative_to(ROOT))]=f
 files.update(license_files())
 manifest={n:hashlib.sha256(p.read_bytes()).hexdigest() for n,p in files.items()}
 with zipfile.ZipFile(a.output,'x',zipfile.ZIP_DEFLATED) as z:
  for name,path in files.items():z.write(path,name)
  z.writestr('MANIFEST.json',json.dumps({'schema':'yuanhe.md-module/1','files':manifest},indent=2))
  z.writestr('README.txt','Standalone membrane MD. Same scientific core as Research Studio.\nInstall deployment/membrane-environment.yml and GROMACS with CUDA support.\nRead docs/07_MODULE_CONTRACTS.md.\npython -m scientific.md_cli validate --bundle /path/to/handoff\nThe package contains no weights, credentials, run data or implicit candidate approval.\n')
 print(json.dumps({'path':str(a.output),'sha256':hashlib.sha256(a.output.read_bytes()).hexdigest(),'files':len(files)}))
if __name__=='__main__':main()
