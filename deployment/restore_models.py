#!/usr/bin/env python3
"""Restore verified local weights with fresh, pinned upstream model code."""
import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

def restore_source(source, destination, model):
    """Use local committed objects if available; never copy a vendor worktree."""
    cache = None
    if (source / '.git').exists():
        check = subprocess.run(['git', '-C', str(source), 'cat-file', '-e',
                                model['repo_commit'] + '^{commit}'], capture_output=True)
        if check.returncode == 0:
            cache = str(source.resolve())
    destination.mkdir(parents=True)
    commands = [
        ['git', 'init', str(destination)],
        ['git', '-C', str(destination), 'config', 'core.autocrlf', 'false'],
        ['git', '-C', str(destination), 'remote', 'add', 'origin', model['source_url']],
        ['git', '-C', str(destination), 'fetch', '--depth=1', cache or 'origin', model['repo_commit']],
        ['git', '-C', str(destination), 'checkout', '--detach', 'FETCH_HEAD'],
    ]
    for command in commands:
        subprocess.run(command, check=True, timeout=900)
    actual = subprocess.check_output(['git', '-C', str(destination), 'rev-parse', 'HEAD'], text=True).strip()
    if actual != model['repo_commit'] or not (destination / 'LICENSE').is_file():
        raise ValueError('Restored source revision or upstream LICENSE mismatch: ' + destination.name)


def restore(source_root, root=ROOT):
    source_root, root = Path(source_root), Path(root)
    candidates = [source_root, source_root / '01_molecule_generation']
    candidates += list(source_root.glob('*/v2'))
    matches = [p for p in candidates if (p / 'vendor/Peptide2Mol').is_dir()
               and (p / 'vendor/REINVENT4').is_dir()]
    if len(matches) != 1:
        raise SystemExit('Provide a source root containing exactly one generation directory with both vendor models')
    source = matches[0]
    target = root / '01_molecule_generation'
    config = json.loads((root / 'deployment/model_dependencies.json').read_text())
    entries = [('peptide2mol', 'Peptide2Mol', 'checkpoint', 'checkpoint_sha256'),
               ('reinvent4', 'REINVENT4', 'prior', 'prior_sha256')]
    # Preflight both sources and destinations before the first copy.
    for key, name, path_key, hash_key in entries:
        model = config[key]
        weights = source / model[path_key]
        if not weights.is_file() or sha256(weights) != model[hash_key]:
            raise SystemExit(f'Missing or mismatched source weights: {weights}')
        if (target / 'vendor' / name).exists() or (target / 'vendor' / name).is_symlink():
            raise SystemExit(f'Destination already exists; refusing overwrite: {name}')
        if not (source / 'vendor' / name / 'LICENSE').is_file():
            raise SystemExit(f'Third-party license missing: {name}')
    for key, name, path_key, hash_key in entries:
        destination = target / 'vendor' / name
        model = config[key]
        restore_source(source / 'vendor' / name, destination, model)
        weights = target / model[path_key]
        weights.parent.mkdir(parents=True, exist_ok=True)
        if weights.exists() and sha256(weights) != model[hash_key]:
            raise ValueError('Pinned source already contains different weights: ' + name)
        if not weights.exists():
            shutil.copy2(source / model[path_key], weights)
        if sha256(weights) != model[hash_key]:
            raise SystemExit(f'Copied weight verification failed: {name}')
        print(f'{name}: pinned source restored; local weight SHA-256 verified; inference not run')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', required=True, type=Path,
                        help='Existing weights; pinned code is fetched if its Git objects are unavailable locally')
    args = parser.parse_args()
    restore(args.source_root)

if __name__ == '__main__':
    main()
