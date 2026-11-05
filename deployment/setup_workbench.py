#!/usr/bin/env python3
"""Validate site settings and prerequisites, then optionally install the full workbench."""
import argparse
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

import deploy_workbench
from restore_models import sha256

ROOT = Path(__file__).resolve().parents[1]


def load_config(path):
    value = json.loads(Path(path).read_text(encoding='utf-8'))
    required = {'prefix', 'model_source', 'gmx', 'host', 'port'}
    metadata = {'schema_version', 'software_version'}
    if not isinstance(value, dict) or not required <= set(value) or set(value) - required - metadata:
        raise ValueError('Use the machine settings and version fields in deployment/site.example.json')
    if set(value) & metadata:
        if not metadata <= set(value) or value['schema_version'] != 'yuanhe.site/2':
            raise ValueError('Unsupported or incomplete site configuration schema')
        if not isinstance(value['software_version'], str) or not value['software_version'].strip():
            raise ValueError('software_version must name the intended release')
    for key in required - {'port'}:
        if not isinstance(value[key], str) or not value[key].strip():
            raise ValueError('Expected a nonempty string: ' + key)
    if type(value['port']) is not int:
        raise ValueError('port must be an integer')
    return value


def preflight(config, root=ROOT):
    checks = []
    def add(name, ok, detail):
        checks.append(dict(check=name, ok=bool(ok), detail=detail))
    if config.get('software_version'):
        actual = (root / 'VERSION').read_text().strip()
        add('release_version', actual == config['software_version'],
            'Configuration expects ' + config['software_version'] + '; source is ' + actual)
    add('platform', platform.system() == 'Linux' and platform.machine() in ('x86_64', 'AMD64'),
        'Linux x86_64 or Linux inside WSL2 required')
    add('python', sys.version_info >= (3, 11), 'Bootstrap Python 3.11+ required')
    try:
        verification = subprocess.run([sys.executable, '-B', '-c',
            'import sys;sys.path.insert(0,' + repr(str(root / '00_workbench')) + ');'
            'from yuanhe.protocols import validate_sources;print(validate_sources())'],
            capture_output=True, text=True, timeout=30)
        add('scientific_source_integrity', verification.returncode == 0,
            (verification.stdout + verification.stderr)[-2000:])
    except (OSError, subprocess.SubprocessError) as exc:
        add('scientific_source_integrity', False, str(exc))
    for tool in ('bash', 'git', 'curl', 'patch', 'make', 'gcc', 'g++', 'tar',
                 os.environ.get('MICROMAMBA', 'micromamba')):
        add('tool:' + tool, shutil.which(tool), 'Install prerequisite; see DEPLOY.md')
    for name, argv, predicate in (
        ('driver', ['nvidia-smi'], lambda text: 'NVIDIA' in text),
        ('gromacs', [config['gmx'], '--version'],
         lambda text: any('GPU support:' in line and 'CUDA' in line for line in text.splitlines())),
    ):
        try:
            result = subprocess.run(argv, capture_output=True, text=True, timeout=60)
            add(name, result.returncode == 0 and predicate(result.stdout + result.stderr),
                'Requires working NVIDIA driver' if name == 'driver' else 'Requires CUDA GROMACS')
        except (OSError, subprocess.SubprocessError) as exc:
            add(name, False, str(exc))
    source = Path(config['model_source'])
    candidates = [source, source / '01_molecule_generation', *source.glob('*/v2')]
    matches = [p for p in candidates if all((p / 'vendor' / n).is_dir() for n in ('Peptide2Mol', 'REINVENT4'))]
    add('model_bundle', len(matches) == 1, 'Run fetch_generation_models.py or supply a compatible existing bundle')
    if len(matches) == 1:
        pins = json.loads((root / 'deployment/model_dependencies.json').read_text())
        for key, name, field in (('peptide2mol', 'Peptide2Mol', 'checkpoint'), ('reinvent4', 'REINVENT4', 'prior')):
            weight = matches[0] / pins[key][field]
            try:
                ok = weight.is_file() and sha256(weight) == pins[key][field + '_sha256']
                ok = ok and (matches[0] / 'vendor' / name / 'LICENSE').is_file()
                add(name, ok, 'Pinned weights and upstream LICENSE required')
            except OSError as exc:
                add(name, False, str(exc))
    record = None
    try:
        record = deploy_workbench.plan(Path(config['prefix']), source, Path(config['gmx']),
                                        config['host'], config['port'], root)
        add('installation_plan', True, 'Clean source and new destination verified')
    except (ValueError, OSError, subprocess.SubprocessError) as exc:
        add('installation_plan', False, str(exc))
    return dict(ready=all(item['ok'] for item in checks), checks=checks, plan=record,
                scientific_acceptance='not_run')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    try:
        result = preflight(load_config(args.config))
        print(json.dumps(result, indent=2), flush=True)
        if not result['ready']:
            print('Fix all failed checks using DEPLOY.md, then repeat. No installation started.', file=sys.stderr)
            return 1
        if args.execute:
            deploy_workbench.execute(result['plan'])
        else:
            print('Preflight passed; add --execute to install. No files changed.')
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
