"""Bounded local job leases for Linux/WSL workers and their descendants.

One pipeline computation and one CPU receptor preparation may overlap. The
original installation-wide lock remains an exclusive maintenance barrier, so
old servers, exports and environment changes cannot race the new workers.
"""
from __future__ import annotations

import hashlib
import os
from pathlib import Path
import tempfile


class Lease:
    def __init__(self):
        self.handles = []

    def take(self, path, shared=False, message='Resource is in use'):
        import fcntl
        handle = path.open('a+')
        try:
            fcntl.flock(handle, (fcntl.LOCK_SH if shared else fcntl.LOCK_EX) | fcntl.LOCK_NB)
        except BlockingIOError:
            handle.close()
            raise RuntimeError(message) from None
        except BaseException:
            handle.close()
            raise
        self.handles.append(handle)

    def fileno(self):
        return self.handles[0].fileno()

    def fds(self):
        return tuple(h.fileno() for h in self.handles)

    def environment(self):
        return {'YUANHE_COMPUTE_FD': str(self.fileno()),
                'YUANHE_COMPUTE_FDS': ','.join(map(str, self.fds()))}

    def close(self):
        # Do not LOCK_UN: inherited descriptors must retain the locks until
        # the last descendant exits, including after a UI service restart.
        for handle in reversed(self.handles):
            handle.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def acquire(root, run=None, lane=None, stage=None, review=None):
    key = hashlib.sha256(str(Path(root).resolve()).encode()).hexdigest()[:24]
    prefix = Path(tempfile.gettempdir()) / ('yuanhe-compute-' + key)
    lease = Lease()
    try:
        lease.take(Path(str(prefix) + '.lock'), shared=lane is not None or review is not None,
                   message='A background process or maintenance operation is active')
        if lane is None and review is None:
            return lease
        run_key = hashlib.sha256(str(Path(run).resolve()).encode()).hexdigest()[:24]
        def scoped(name):
            return Path(str(prefix) + '-' + run_key + '-' + name + '.lock')
        if review is not None:
            lease.take(scoped('review-' + review), message='This stage is still active; settings are immutable')
            return lease
        if lane not in {'compute', 'receptor'}:
            raise ValueError('Unknown resource lane')
        lease.take(scoped(lane), message='This run already has an active ' + lane + ' job')
        if stage:
            lease.take(scoped('review-' + stage), message='Stage settings are being updated')
        lease.take(Path(str(prefix) + '-' + lane + '.lock'),
                   message=lane.capitalize() + ' capacity is occupied; retry after that job finishes')
        return lease
    except BaseException:
        lease.close()
        raise


def receptor_environment():
    # Leave the GPU to model/MD workers and bound CPU library thread pools.
    cpus = len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else (os.cpu_count() or 1)
    threads = str(max(1, min(2, cpus)))
    return {'CUDA_VISIBLE_DEVICES': '', 'OPENMM_DEFAULT_PLATFORM': 'CPU',
            **{key: threads for key in ('OPENMM_CPU_THREADS', 'OMP_NUM_THREADS',
                                        'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                                        'NUMEXPR_NUM_THREADS')}}
