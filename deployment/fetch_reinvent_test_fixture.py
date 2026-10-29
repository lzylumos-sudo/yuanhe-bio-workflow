"""Fetch only the pinned REINVENT4 fixture for isolated CI model regression."""
import argparse
from pathlib import Path
import shutil
import fetch_generation_models as fetch
import apply_vendor_patches as patches

ROOT = Path(__file__).resolve().parents[1]

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--destination', type=Path, required=True)
    args = ap.parse_args()
    record = fetch.plan(args.destination.resolve())
    record['models'] = [r for r in record['models'] if r['name'] == 'REINVENT4']
    fetch.execute(record)
    shutil.copytree(args.destination/'vendor/REINVENT4', ROOT/'01_molecule_generation/vendor/REINVENT4',
                    ignore=shutil.ignore_patterns('.git','__pycache__'))
    patches.main('REINVENT4')

if __name__ == '__main__':
    main()
