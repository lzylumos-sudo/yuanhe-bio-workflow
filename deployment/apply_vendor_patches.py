#!/usr/bin/env python3
"""Apply bounded-sampling/stereo patches only to the recorded vendor revision."""
import argparse,hashlib,json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PATCHES=Path(__file__).with_name('vendor_patches')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(vendor=None):
 rows=json.loads((PATCHES/'manifest.json').read_text())
 if vendor:rows=[row for row in rows if row['vendor']==vendor]
 if not rows:raise ValueError('No recorded patches for vendor')
 for row in rows:
  p=ROOT/'01_molecule_generation/vendor'/row['vendor']/row['file']
  accepted = [row['before_sha256'],row['after_sha256']] + [x['before_sha256'] for x in row.get('upgrades', [])]
  if not p.is_file() or sha(p) not in accepted:
   raise SystemExit(f'Vendor source revision differs: {p}; inspect before applying')
 for row in rows:
  vendor=ROOT/'01_molecule_generation/vendor'/row['vendor'];p=vendor/row['file']
  if sha(p)==row['after_sha256']:continue
  patch = next((x['patch'] for x in row.get('upgrades', []) if sha(p)==x['before_sha256']), row['patch'])
  command=['patch','--batch','-p1','-i',str(PATCHES/patch)]
  subprocess.run(command+['--dry-run'],cwd=vendor,check=True)
  subprocess.run(command,cwd=vendor,check=True)
  if sha(p)!=row['after_sha256']:raise RuntimeError(f'Patched source hash mismatch: {p}')
  print(row['vendor']+': patched and source hash verified')
if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--vendor',choices=['Peptide2Mol','REINVENT4'])
 main(parser.parse_args().vendor)
