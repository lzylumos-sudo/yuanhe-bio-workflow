"""Export the pure peptide-contact comparison module and its licensing."""
import argparse,hashlib,json,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);args=p.parse_args()
 if args.output.exists():raise ValueError('Output already exists')
 files={f'scientific/{n}':ROOT/'00_workbench/scientific'/n for n in ('__init__.py','artifacts.py','peptide_contacts.py')}
 for n in ('LICENSE','LICENSE_SCOPE.txt','THIRD_PARTY_NOTICES.md'):files[n]=ROOT/n
 args.output.parent.mkdir(parents=True,exist_ok=True)
 with zipfile.ZipFile(args.output,'x',zipfile.ZIP_DEFLATED) as z:
  for name,src in files.items():z.write(src,name)
  z.writestr('MANIFEST.json',json.dumps({'schema':'yuanhe.peptide-evidence/1','files':{n:hashlib.sha256(p.read_bytes()).hexdigest() for n,p in files.items()}},indent=2))
  z.writestr('README.txt','Requires Python >=3.11, numpy and RDKit.\npython -m scientific.peptide_contacts --manifest input.json --output result\nManifest keys: receptor, reference_receptor (original peptide coordinate frame), peptide, candidates (compound_id,pose_id,sdf), selection (primary/preliminary/exploratory arrays with compound_id and raw.affinity). Paths are absolute or relative to manifest.\nThis distance-contact metric is not a pharmacophore or activity predictor.\n')
if __name__=='__main__':main()
