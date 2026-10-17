# -*- coding: utf-8 -*-
"""webui.jobs -- run the pipeline entry points from the browser and stream logs.

Only whitelisted commands can be started (the four entry scripts of the
repository); arguments are validated per command, so the web console can never
be used to run an arbitrary shell command.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import threading
import time
import uuid

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

KINDS = {
    'evidence': dict(
        label='Evidence reconstruction',
        script='run_evidence.py',
        docs='Recalculate from three PDB structures and verify all 24 literature priors (approximately 30 s)',
        args={'out': dict(type='path', default='evidence')},
        suffix=[('out', '--out')],
    ),
    'pipeline': dict(
        label='Design pipeline',
        script='run_pipeline.py',
        docs='Interface → Pharmacophore → Pocket → Scaffold enumeration → Pose scoring → Rank export (approximately 50 s)',
        args={'out': dict(type='path', default='results')},
        suffix=[('out', '--out')],
    ),
    'filters': dict(
        label='Structure-derived filters',
        script='run_filters.py',
        docs='F1 Anchor specificity / F2 Like-charge repulsion / F3 Cross-subpocket / F4 Three-structure ensemble (approximately 80 s)',
        args={'out': dict(type='path', default='results'),
              'targets': dict(type='int', default=0, minimum=0, maximum=5000),
              'max_spread': dict(type='float', default=5.0, minimum=0.5, maximum=50.0)},
        suffix=[('out', '--out'), ('targets', '--targets'), ('max_spread', '--max-spread')],
    ),
    'tests': dict(
        label='Smoke tests',
        script='tests/test_smoke.py',
        docs='Structure parsing / Priors / Evidence chain / Assembly / Enumeration / Pocket / End-to-end self-check',
        args={},
        suffix=[],
    ),
}

_SAFE_PATH = re.compile(r'^[\w./\\-]{1,64}$')
MAX_LINES = 6000

_lock = threading.Lock()
_jobs: dict = {}
_order: list = []


def _now():
    return time.time()


def list_kinds():
    out = []
    for key, spec in KINDS.items():
        out.append(dict(kind=key, label=spec['label'], script=spec['script'],
                        docs=spec['docs'],
                        args={k: {kk: vv for kk, vv in v.items()} for k, v in spec['args'].items()}))
    return out


def list_jobs(limit=30):
    """Recent jobs, without log lines or the process handle (JSON-safe)."""
    with _lock:
        ids = list(_order)[-limit:]
        out = []
        for job_id in reversed(ids):
            job = _jobs[job_id]
            out.append(dict(id=job['id'], kind=job['kind'], label=job['label'],
                            status=job['status'], started=job['started'],
                            ended=job['ended'], exit=job['exit'], error=job['error'],
                            cmd=job['cmd'], params=job.get('params') or {},
                            total_lines=len(job['lines'])))
        return out


def get_job(job_id, since=0):
    with _lock:
        job = _jobs.get(job_id)
        if job is None:
            return None
        lines = [line for line in job['lines'] if line['i'] >= max(0, int(since))]
        return dict(id=job['id'], kind=job['kind'], label=job['label'],
                    status=job['status'], cmd=job['cmd'], started=job['started'],
                    ended=job['ended'], exit=job['exit'], error=job['error'],
                    lines=lines, cursor=job.get('next_line', len(job['lines'])),
                    total_lines=len(job['lines']))


def _validate(kind, payload):
    spec = KINDS[kind]
    argv, out = [], {}
    for name, rules in spec['args'].items():
        raw = (payload or {}).get(name, rules.get('default'))
        if rules['type'] == 'path':
            if not _SAFE_PATH.fullmatch(str(raw)) or '..' in str(raw) or os.path.isabs(str(raw)) or str(raw).startswith(('\\', '/')):
                raise ValueError(f'Invalid output directory: {raw!r}')
            resolved = os.path.realpath(os.path.join(ROOT, str(raw)))
            if os.path.commonpath([resolved, os.path.realpath(ROOT)]) != os.path.realpath(ROOT):
                raise ValueError('Output directory is outside the project')
            out[name] = str(raw)
        elif rules['type'] == 'int':
            val = int(float(raw))
            if not (rules['minimum'] <= val <= rules['maximum']):
                raise ValueError(f'{name} Out of range [{rules["minimum"]}, {rules["maximum"]}]')
            out[name] = val
        elif rules['type'] == 'float':
            val = float(raw)
            if not (rules['minimum'] <= val <= rules['maximum']):
                raise ValueError(f'{name} Out of range [{rules["minimum"]}, {rules["maximum"]}]')
            out[name] = val
    for key, flag in spec['suffix']:
        argv += [flag, str(out[key])]
    return argv, out


def start(kind, payload=None):
    if kind not in KINDS:
        raise ValueError(f'Unknown job type: {kind}')
    raise RuntimeError('Legacy enumeration is archived as read-only; use the workbench model workflow')
    spec = KINDS[kind]
    with _lock:
        running_same = [j for j in _jobs.values()
                        if j['kind'] == kind and j['status'] in ('queued', 'running')]
        if running_same:
            raise RuntimeError(f'This job is already running ({running_same[0]["id"]})')
        n_running = sum(1 for j in _jobs.values() if j['status'] in ('queued', 'running'))
        if n_running >= 2:
            raise RuntimeError('At most 2 concurrent jobs; wait for an active job to finish')
    argv, params = _validate(kind, payload or {})
    cmd = [sys.executable, os.path.join(ROOT, spec['script'])] + argv
    job_id = uuid.uuid4().hex[:10]
    job = dict(id=job_id, kind=kind, label=spec['label'], cmd=cmd, params=params,
               status='queued', started=None, ended=None, exit=None, error=None,
               lines=[], proc=None)
    with _lock:
        active = [j for j in _jobs.values() if j['status'] in ('queued', 'running')]
        if active:
            raise RuntimeError('A job is running; shared output directories require serial writes')
        _jobs[job_id] = job
        _order.append(job_id)
    threading.Thread(target=_run, args=(job_id,), name=f'job-{job_id}', daemon=True).start()
    return job_id


def _append(job, text, stream='stdout'):
    line = text.rstrip('\r\n')
    with _lock:
        index = job.get('next_line', 0)
        job['lines'].append(dict(i=index, t=round(_now(), 2), text=line))
        job['next_line'] = index + 1
        if len(job['lines']) > MAX_LINES:
            del job['lines'][:len(job['lines']) - MAX_LINES]


def _run(job_id):
    with _lock:
        job = _jobs[job_id]
        if job['status'] == 'cancelled':
            return
        job['status'] = 'running'
        job['started'] = _now()
    spec = KINDS[job['kind']]
    env = dict(os.environ)
    env['PYTHONIOENCODING'] = 'utf-8'
    env['PYTHONUNBUFFERED'] = '1'
    _append(job, f'$ {os.path.basename(spec["script"])} '
                 f'{" ".join(job["cmd"][2:])}'.strip())
    try:
        proc = subprocess.Popen(job['cmd'], cwd=ROOT, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, bufsize=1,
                                encoding='utf-8', errors='replace', env=env)
    except OSError as exc:
        with _lock:
            job['status'] = 'failed'
            job['error'] = f'{type(exc).__name__}: {exc}'
            job['ended'] = _now()
            job['exit'] = -1
        _append(job, f'[Startup failed] {job["error"]}', 'stderr')
        return
    with _lock:
        job['proc'] = proc
        if job['status'] == 'cancelled':
            proc.terminate()
    try:
        assert proc.stdout is not None
        for line in proc.stdout:
            _append(job, line)
        code = proc.wait()
    except Exception as exc:                                # noqa: BLE001
        code = -1
        job['error'] = f'{type(exc).__name__}: {exc}'
        try:
            proc.kill()
        except OSError:
            pass
    with _lock:
        job['exit'] = code
        job['ended'] = _now()
        job['proc'] = None
        if job['status'] != 'cancelled':
            job['status'] = 'done' if code == 0 else 'failed'
    _append(job, f'[{"Completed" if code == 0 else "Failed"}] exit={code} '
                 f'Duration {job["ended"] - job["started"]:.1f}s')


def cancel(job_id):
    with _lock:
        job = _jobs.get(job_id)
        if job is None:
            raise KeyError(job_id)
        proc = job.get('proc')
        if job['status'] not in ('queued', 'running'):
            return False
        job['status'] = 'cancelled'
    try:
        if proc is not None:
            proc.terminate()
    except OSError:
        pass
    return True
