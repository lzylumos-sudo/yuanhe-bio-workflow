#!/usr/bin/env python3
"""Bind installed environments/assets to current research defaults, atomically."""
import argparse
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / '00_workbench'))
sys.path.insert(0, str(ROOT / '02_docking/src'))
from research.config import resolve
from yuanhe_preflight.config import validate, absolute_paths


def configured(existing, prefix, gmx, assets, root=ROOT):
    """Retain valid existing policies; only bind installation-specific assets."""
    cfg = copy.deepcopy(existing)
    try:
        research = resolve(cfg.get('research'))
    except (ValueError, KeyError, TypeError) as exc:
        raise ValueError('Existing research configuration is incompatible. Preserve it with '
                         'the old deployment; configure a fresh checkout. No runtime file changed: ' + str(exc)) from exc
    validate(assets)
    pc = research['predocking']
    if 'engine' not in assets or 'engine' not in pc:
        raise ValueError('Current synthesis-only engine configuration is required')
    for key in ('python', 'model_path', 'templates_path', 'stock_path', 'model_version',
                'templates_version', 'stock_version', 'model_license', 'templates_license', 'stock_license'):
        if not assets['retro'].get(key):
            raise ValueError('Missing synthesis installation asset: ' + key)
        pc['retro'][key] = assets['retro'][key]
    pc['engine']['toxicity']['python'] = assets['engine']['toxicity']['python']
    for key in ('scscore_path', 'scscore_sha256'):
        if not assets['engine']['synthesis'].get(key):
            raise ValueError('Missing synthesis installation asset: ' + key)
        pc['engine']['synthesis'][key] = assets['engine']['synthesis'][key]
    research['predocking'] = absolute_paths(pc, root)
    cfg['research'] = resolve(research)
    cfg['predocking'] = copy.deepcopy(cfg['research']['predocking'])
    runtime = str(prefix / 'envs/runtime/bin/python')
    cfg.update(chem_python=runtime, screen_python=runtime, prep_python=runtime, md_python=runtime,
               peptide_python=str(prefix / 'envs/peptide2mol/bin/python'),
               reinvent_executable=str(prefix / 'envs/reinvent4/bin/reinvent'),
               vina=str(prefix / 'bin/vina'), gmx=str(gmx))
    return cfg


def write_config(prefix, gmx, root=ROOT):
    dest = root / 'deployment/runtime.json'
    assets = prefix / 'predocking-assets/predocking-config.json'
    if not assets.is_file():
        raise ValueError('Run install_workbench.sh synthesis first; missing ' + str(assets))
    proc = subprocess.run([str(gmx), '--version'], capture_output=True, text=True, timeout=60)
    if proc.returncode or not any('GPU support:' in line and 'CUDA' in line
                                  for line in (proc.stdout + proc.stderr).splitlines()):
        raise ValueError('Selected GROMACS does not report a working CUDA build')
    cfg = configured(json.loads(dest.read_text()) if dest.exists() else {}, prefix, gmx,
                     json.loads(assets.read_text()), root)
    evidence = prefix / 'evidence'
    evidence.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        shutil.copy2(dest, evidence / ('runtime.before-' + str(time.time_ns()) + '.json'))
    temp = dest.with_name('runtime.' + str(time.time_ns()) + '.tmp')
    try:
        temp.write_text(json.dumps(cfg, indent=2, allow_nan=False) + '\n')
        os.replace(temp, dest)
    finally:
        temp.unlink(missing_ok=True)
    return dest


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--prefix', type=Path, required=True)
    ap.add_argument('--gmx', type=Path, required=True)
    args = ap.parse_args()
    if not args.prefix.is_absolute() or not args.gmx.is_absolute():
        ap.error('prefix and gmx must be absolute')
    print('Saved ' + str(write_config(args.prefix, args.gmx)) + '; readiness and real acceptance remain required.')


if __name__ == '__main__':
    main()
