"""Run an explicitly reviewed, self-contained TPR; do not infer molecular qualification."""
import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from . import artifacts as a


def capture(argv, folder, name, timeout=120):
    proc = subprocess.run([str(x) for x in argv], cwd=folder, capture_output=True, text=True,
                          errors='replace', timeout=timeout, pass_fds=a.lease_fds())
    text = proc.stdout+'\n'+proc.stderr
    (folder/name).write_text(text, encoding='utf-8')
    if proc.returncode:
        raise RuntimeError('GROMACS command failed; inspect '+name)
    return text


def value(text, name, cast=float):
    m = re.search(r'^\s*'+re.escape(name)+r'\s*=\s*(\S+)', text, re.M)
    if not m:
        raise ValueError('TPR metadata missing: '+name)
    return cast(m[1])


def run(d, cfg, params, source):
    folder = d/'prepared_md'; folder.mkdir(exist_ok=True)
    gmx = shutil.which(cfg.get('gmx', '')) or cfg.get('gmx', '')
    if not gmx or not Path(gmx).is_file():
        raise ValueError('Configure a GROMACS executable in deployment settings')
    source_hash = a.digest(source); engine_hash = a.digest(gmx)
    identity = {'input_sha256': source_hash, 'engine_sha256': engine_hash, 'parameters': params}
    old = a.read(folder/'identity.json')
    if old and old != identity:
        raise ValueError('MD input, engine or execution settings changed; resume refused')
    if not old:
        shutil.copy2(source, folder/'production.tpr'); a.atomic(folder/'identity.json', identity)
    if a.digest(folder/'production.tpr') != source_hash:
        raise ValueError('Staged TPR changed')
    capture([gmx, '--version'], folder, 'engine.txt')
    dump = capture([gmx, 'dump', '-s', 'production.tpr'], folder, 'tpr_dump.txt')
    steps = value(dump, 'nsteps', int); initial = value(dump, 'init-step', int)
    dt = value(dump, 'dt'); tinit = value(dump, 'tinit')
    if steps <= 0 or initial < 0 or not math.isfinite(dt) or dt <= 0 or not math.isfinite(tinit):
        raise ValueError('A finite, bounded dynamics TPR is required')
    if value(dump, 'integrator', str) not in {'md', 'md-vv', 'md-vv-avek', 'sd', 'bd'}:
        raise ValueError('Upload a production dynamics TPR, not an energy minimization input')
    if value(dump, 'nstxout-compressed', int) <= 0:
        raise ValueError('TPR must record a compressed trajectory')
    plan = {'initial_step': initial, 'final_step': initial+steps, 'dt_ps': dt,
            'final_time_ps': tinit+(initial+steps)*dt, 'duration_ns': steps*dt/1000}
    a.atomic(folder/'plan.json', plan)
    if (folder/'completion.json').exists():
        done = a.read(folder/'completion.json')
        if all(a.digest(folder/k) == v for k, v in done['files'].items()):
            return done
        raise ValueError('Completed MD artifacts changed')
    argv = [gmx, 'mdrun', '-s', 'production.tpr', '-deffnm', 'production',
            '-ntmpi', '1', '-ntomp', str(params['threads']), '-nb', params['device'], '-cpt', '1']
    if (folder/'production.cpt').exists():
        argv += ['-cpi', 'production.cpt', '-append']
    elif (folder/'production.log').exists():
        raise ValueError('Interrupted output has no checkpoint; preserve this job and create a new one')
    # The worker process group and inherited lease also cover this child.
    with (folder/'execution.log').open('a', encoding='utf-8') as log:
        proc = subprocess.run(argv, cwd=folder, stdout=log, stderr=subprocess.STDOUT,
                              timeout=int(cfg.get('timeout_seconds', 86400)), pass_fds=a.lease_fds(),
                              env={**os.environ, 'OMP_NUM_THREADS': str(params['threads'])})
    if proc.returncode:
        raise RuntimeError('MD execution failed; inspect execution.log and retained checkpoint')
    log = (folder/'production.log').read_text(errors='replace')
    observed = re.findall(r'^\s*Step\s+Time\s*\n\s*(\d+)\s+([\d.eE+-]+)', log, re.M)
    if 'Finished mdrun' not in log or not observed or int(observed[-1][0]) != plan['final_step']:
        raise RuntimeError('MD did not finish its planned final step')
    if abs(float(observed[-1][1])-plan['final_time_ps']) > max(dt, .001):
        raise RuntimeError('Final simulation time differs from the TPR plan')
    if re.search(r'Fatal error|LINCS WARNING|constraint failure|not finite|nan detected', log, re.I):
        raise RuntimeError('MD log contains a fatal or constraint warning; acceptance blocked')
    checkpoint = capture([gmx, 'dump', '-cp', 'production.cpt'], folder, 'checkpoint_audit.txt')
    if value(checkpoint, 'step', int) != plan['final_step']:
        raise RuntimeError('Final checkpoint does not match planned completion')
    capture([gmx, 'check', '-f', 'production.xtc'], folder, 'trajectory_check.txt')
    required = ['production.tpr', 'production.log', 'production.gro', 'production.cpt',
                'production.xtc', 'production.edr', 'checkpoint_audit.txt', 'trajectory_check.txt', 'plan.json']
    done = {'status': 'completed', 'scope': 'prepared-input MD execution acceptance',
            'candidate_qualified': False, 'scientific_validation': 'not_established',
            'plan': plan, 'files': {name: a.digest(folder/name) for name in required}}
    a.atomic(folder/'completion.json', done)
    return done
