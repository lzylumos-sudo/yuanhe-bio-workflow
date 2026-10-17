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
        label='证据链重建',
        script='run_evidence.py',
        docs='从三个 PDB 现算全部结构数字，并逐条核对 24 条文献先验（约 30 s）',
        args={'out': dict(type='path', default='evidence')},
        suffix=[('out', '--out')],
    ),
    'pipeline': dict(
        label='设计管线',
        script='run_pipeline.py',
        docs='界面 → 药效团 → 口袋 → 骨架枚举 → 位姿打分 → 排名导出（约 50 s）',
        args={'out': dict(type='path', default='results')},
        suffix=[('out', '--out')],
    ),
    'filters': dict(
        label='结构衍生过滤器',
        script='run_filters.py',
        docs='F1 锚点特异性 / F2 同电荷排斥 / F3 跨亚口袋 / F4 三结构集合（约 80 s）',
        args={'out': dict(type='path', default='results'),
              'targets': dict(type='int', default=0, minimum=0, maximum=5000),
              'max_spread': dict(type='float', default=5.0, minimum=0.5, maximum=50.0)},
        suffix=[('out', '--out'), ('targets', '--targets'), ('max_spread', '--max-spread')],
    ),
    'tests': dict(
        label='冒烟测试',
        script='tests/test_smoke.py',
        docs='结构解析 / 先验 / 证据链 / 组装 / 枚举 / 口袋 / 端到端自检',
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
                raise ValueError(f'非法的输出目录: {raw!r}')
            resolved = os.path.realpath(os.path.join(ROOT, str(raw)))
            if os.path.commonpath([resolved, os.path.realpath(ROOT)]) != os.path.realpath(ROOT):
                raise ValueError('输出目录超出项目范围')
            out[name] = str(raw)
        elif rules['type'] == 'int':
            val = int(float(raw))
            if not (rules['minimum'] <= val <= rules['maximum']):
                raise ValueError(f'{name} 超出范围 [{rules["minimum"]}, {rules["maximum"]}]')
            out[name] = val
        elif rules['type'] == 'float':
            val = float(raw)
            if not (rules['minimum'] <= val <= rules['maximum']):
                raise ValueError(f'{name} 超出范围 [{rules["minimum"]}, {rules["maximum"]}]')
            out[name] = val
    for key, flag in spec['suffix']:
        argv += [flag, str(out[key])]
    return argv, out


def start(kind, payload=None):
    if kind not in KINDS:
        raise ValueError(f'未知任务类型: {kind}')
    raise RuntimeError('旧版枚举入口已归档为只读，请使用统一工作台的真实模型流程')
    spec = KINDS[kind]
    with _lock:
        running_same = [j for j in _jobs.values()
                        if j['kind'] == kind and j['status'] in ('queued', 'running')]
        if running_same:
            raise RuntimeError(f'该任务已在运行中（{running_same[0]["id"]}）')
        n_running = sum(1 for j in _jobs.values() if j['status'] in ('queued', 'running'))
        if n_running >= 2:
            raise RuntimeError('同时最多运行 2 个任务，请等待当前任务结束')
    argv, params = _validate(kind, payload or {})
    cmd = [sys.executable, os.path.join(ROOT, spec['script'])] + argv
    job_id = uuid.uuid4().hex[:10]
    job = dict(id=job_id, kind=kind, label=spec['label'], cmd=cmd, params=params,
               status='queued', started=None, ended=None, exit=None, error=None,
               lines=[], proc=None)
    with _lock:
        active = [j for j in _jobs.values() if j['status'] in ('queued', 'running')]
        if active:
            raise RuntimeError('已有任务运行；共享输出目录只能串行写入')
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
        _append(job, f'[启动失败] {job["error"]}', 'stderr')
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
    _append(job, f'[{"完成" if code == 0 else "失败"}] exit={code} '
                 f'用时 {job["ended"] - job["started"]:.1f}s')


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
