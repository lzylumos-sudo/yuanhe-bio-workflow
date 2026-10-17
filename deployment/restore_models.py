#!/usr/bin/env python3
"""Restore locally supplied models without installing or running them."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', required=True, type=Path)
    args = parser.parse_args()
    candidates = [args.source_root, args.source_root / '01_molecule_generation']
    candidates += list(args.source_root.glob('*/v2'))
    matches = [p for p in candidates if (p / 'vendor/Peptide2Mol').is_dir()
               and (p / 'vendor/REINVENT4').is_dir()]
    if len(matches) != 1:
        raise SystemExit('Provide a source root containing exactly one generation directory with both vendor models')
    source = matches[0]
    target = ROOT / '01_molecule_generation'
    config = json.loads((ROOT / 'deployment/model_dependencies.json').read_text())
    entries = [('peptide2mol', 'Peptide2Mol', 'checkpoint', 'checkpoint_sha256'),
               ('reinvent4', 'REINVENT4', 'prior', 'prior_sha256')]
    # Preflight both sources and destinations before the first copy.
    for key, name, path_key, hash_key in entries:
        model = config[key]
        weights = source / model[path_key]
        if not weights.is_file() or sha256(weights) != model[hash_key]:
            raise SystemExit(f'Missing or mismatched source weights: {weights}')
        if (target / 'vendor' / name).exists():
            raise SystemExit(f'Destination already exists; refusing overwrite: {name}')
        if not (source / 'vendor' / name / 'LICENSE').is_file():
            raise SystemExit(f'Third-party license missing: {name}')
    for key, name, path_key, hash_key in entries:
        destination = target / 'vendor' / name
        shutil.copytree(source / 'vendor' / name, destination,
                        ignore=shutil.ignore_patterns('.git', '.env', '.env.*', '._*', '.DS_Store',
                                                      '__pycache__', '.pytest_cache', '*.egg-info'))
        if sha256(target / config[key][path_key]) != config[key][hash_key]:
            raise SystemExit(f'Copied weight verification failed: {name}')
        print(f'{name}: copied; weight SHA-256 verified; inference not run')

if __name__ == '__main__':
    main()
