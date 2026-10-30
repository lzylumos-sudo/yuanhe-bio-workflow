#!/usr/bin/env python3
"""Plan or install a separate workbench deployment; never switch a live service."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / '00_workbench'))
from yuanhe.version import git_identity, software_version


def source_identity(root=ROOT):
    if (root / 'SOURCE_MANIFEST.json').is_file():
        from verify_release import verify
        return verify(root)['source_commit']
    identity = git_identity(root)
    if identity['commit'] == 'unknown' or identity['dirty'] is not False:
        raise ValueError('Use a clean Git checkout or a verified administrator package')
    return identity['commit']


def plan(prefix, model_source, gmx, host, port, root=ROOT):
    if not prefix.is_absolute() or prefix == Path('/') or prefix.exists():
        raise ValueError('Choose an absolute NEW installation prefix; existing directories are preserved')
    if not model_source.is_absolute() or not model_source.is_dir():
        raise ValueError('Provide the existing licensed model source directory as an absolute path')
    if not gmx.is_absolute() or not gmx.is_file() or not os.access(gmx, os.X_OK):
        raise ValueError('Provide the existing executable CUDA GROMACS path')
    if not host or any(c.isspace() for c in host) or not 1 <= port <= 65535:
        raise ValueError('Supply the approved private bind address and port')
    if (root / 'deployment/runtime.json').exists() or (root / '00_workbench/workflow_runs').exists():
        raise ValueError('Use a fresh checkout; preserve the existing runtime and run directories')
    for name in ('Peptide2Mol', 'REINVENT4'):
        if (root / '01_molecule_generation/vendor' / name).exists():
            raise ValueError('Use a fresh checkout without restored vendor directories')
    commit = source_identity(root)
    commands = [('restore_models', [sys.executable, str(root / 'deployment/restore_models.py'),
                                   '--source-root', str(model_source)])]
    for step in ('runtime', 'peptide', 'reinvent', 'synthesis', 'vina', 'configure', 'check'):
        commands.append((step, ['bash', str(root / 'deployment/install_workbench.sh'), step, str(prefix)]))
    return {'schema': 'yuanhe.deployment/1', 'software_version': software_version(),
            'source_commit': commit, 'root': str(root), 'prefix': str(prefix), 'gmx': str(gmx),
            'host': host, 'port': port, 'commands': commands,
            'service_switched': False, 'scientific_acceptance': 'not_run'}


def execute(record, root=ROOT):
    if platform.system() != 'Linux' or platform.machine() not in ('x86_64', 'AMD64'):
        raise ValueError('Installation requires Linux x86_64 / WSL2')
    for tool in ('bash', 'git', 'curl', 'patch', os.environ.get('MICROMAMBA', 'micromamba')):
        if not shutil.which(tool):
            raise ValueError('Missing administrator prerequisite: ' + tool)
    probe = subprocess.run([record['gmx'], '--version'], capture_output=True, text=True, timeout=60)
    if probe.returncode or not any('GPU support:' in line and 'CUDA' in line
                                  for line in (probe.stdout + probe.stderr).splitlines()):
        raise ValueError('GROMACS must report CUDA support before installation starts')
    prefix = Path(record['prefix'])
    prefix.mkdir(parents=True, exist_ok=False)
    evidence = prefix / 'evidence'
    evidence.mkdir()
    state = {**record, 'started': time.time(), 'status': 'installing', 'steps': []}
    def save():
        temp = evidence / 'deployment.json.tmp'
        temp.write_text(json.dumps(state, indent=2) + '\n')
        temp.replace(evidence / 'deployment.json')
    save()
    env = dict(os.environ, YUANHE_GMX=record['gmx'])
    env.pop('PYTHONPATH', None)
    env.pop('PYTHONHOME', None)
    env['PYTHONNOUSERSITE'] = '1'
    try:
        for name, argv in record['commands']:
            state['current_step'] = name
            save()
            print('Installing step: ' + name, flush=True)
            with (evidence / (name + '.log')).open('x') as log:
                proc = subprocess.run(argv, cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT)
            state['steps'].append({'name': name, 'returncode': proc.returncode})
            if proc.returncode:
                raise RuntimeError('Step failed: ' + name + '; inspect ' + str(evidence / (name + '.log')))
            save()
        readiness = json.loads((evidence / 'readiness.json').read_text())
        if not readiness.get('all_runtime_checks_pass'):
            raise RuntimeError('Readiness did not pass; service handoff refused')
        launch = evidence / 'start-workbench.sh'
        launch.write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec bash ' +
                          shlex.quote(str(root / 'start_linux.sh')) + ' --host ' +
                          shlex.quote(record['host']) + ' --port ' + str(record['port']) + '\n')
        launch.chmod(0o700)
        state.update(status='ready_for_service_switch', launch_script=str(launch),
                     runtime_sha256=hashlib.sha256((root / 'deployment/runtime.json').read_bytes()).hexdigest())
    except BaseException as exc:
        state.update(status='failed', error=str(exc))
        raise
    finally:
        state['ended'] = time.time()
        save()
    print('Readiness passed. Inspect ' + str(evidence / 'deployment.json'))
    print('After checking existing jobs, switch the site service to: bash ' + str(launch))
    print('No service was stopped or started. Real UI/model/MD acceptance is still required.')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--prefix', required=True, type=Path)
    ap.add_argument('--model-source', required=True, type=Path)
    ap.add_argument('--gmx', required=True, type=Path)
    ap.add_argument('--host', default='127.0.0.1')
    ap.add_argument('--port', type=int, default=8777)
    ap.add_argument('--execute', action='store_true', help='Install after inspecting the default read-only plan')
    args = ap.parse_args()
    try:
        record = plan(args.prefix, args.model_source, args.gmx, args.host, args.port)
        if args.execute:
            execute(record)
        else:
            print(json.dumps(record, indent=2))
            print('Plan only; no files changed. Add --execute to install into this new prefix.')
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
