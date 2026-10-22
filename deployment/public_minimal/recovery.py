"""Transactional recovery for the bounded public engineering control (POSIX)."""
import contextlib
import fcntl
import hashlib
import json
import os
import shutil
import signal
import subprocess
import uuid


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_dir = path.parent / '.atomic'
    temporary_dir.mkdir(exist_ok=True)
    temporary = temporary_dir / (path.name + '.' + uuid.uuid4().hex + '.tmp')
    with temporary.open('w') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


EXCLUDED = {'.atomic', 'cache', 'scratch', '.recovery', '.run.lock', 'RUN_STATE.json',
            'FAILED.json', 'COMPLETED.json', 'artifact_hashes.json',
            'pip-check.txt', 'gpu.json', 'source_hashes.json'}


def inventory(out):
    result = {}
    for path in out.rglob('*'):
        relative = path.relative_to(out)
        if relative.parts[0] in EXCLUDED:
            continue
        if path.is_symlink():
            raise ValueError('Symlink in run artifacts: ' + str(relative))
        if path.is_file():
            result[relative.as_posix()] = digest(path)
    return result


def verify(out, expected):
    for name, sha in expected.items():
        path = out / name
        if not path.is_file() or path.is_symlink() or digest(path) != sha:
            raise ValueError('Recorded artifact changed or missing: ' + name)


@contextlib.contextmanager
def lease(out):
    with (out / '.run.lock').open('a+') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Another process still owns this run') from None
        yield handle.fileno()


@contextlib.contextmanager
def interruption_handlers():
    def stop(signum, frame):
        raise InterruptedError('Interrupted by signal ' + str(signum))
    saved = {sig: signal.signal(sig, stop) for sig in (signal.SIGINT, signal.SIGTERM)}
    try:
        yield
    finally:
        for sig, handler in saved.items():
            signal.signal(sig, handler)


def execute(command, log, lock_fd, timeout=600):
    with log.open('w') as stream:
        child = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT,
                                 start_new_session=True, pass_fds=(lock_fd,))
        try:
            code = child.wait(timeout=timeout)
            if code:
                raise subprocess.CalledProcessError(code, command)
        except BaseException:
            for sig in (signal.SIGTERM, signal.SIGKILL):
                try:
                    os.killpg(child.pid, sig)
                except ProcessLookupError:
                    pass
                if sig == signal.SIGTERM:
                    try:
                        child.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        pass
            child.wait()
            raise


class Recovery:
    def __init__(self, out, identity, resume=False):
        self.out = out
        self.path = out / 'RUN_STATE.json'
        self.directory = out / '.recovery'
        if resume:
            self.state = json.loads(self.path.read_text())
            if self.state.get('schema') != 1 or self.state['identity'] != identity:
                raise ValueError('Code, interpreter, environment or run path changed; resume refused')
        else:
            if self.path.exists() or inventory(out):
                raise ValueError('A fresh output directory is required; use --resume for recorded runs')
            self.state = {'schema': 1, 'identity': identity, 'completed': [],
                          'artifacts': {}, 'active': None, 'status': 'running'}
            self.save()
        self.directory.mkdir(exist_ok=True)
        self.restore()

    def save(self):
        write(self.path, self.state)

    def restore(self):
        active = self.state['active']
        if not active:
            verify(self.out, self.state['artifacts'])
            if inventory(self.out) != self.state['artifacts']:
                raise ValueError('Unrecorded artifacts in completed stages; resume refused')
            return
        baseline = self.directory / active['baseline']
        verify(baseline, self.state['artifacts'])
        if active['phase'] == 'dynamics':
            allowed = lambda name: name.startswith('dynamics/md/replica_') or name in {
                'dynamics.log', 'dynamics_timing.json', 'md_execution.json'}
            current = inventory(self.out)
            protected = {k: v for k, v in self.state['artifacts'].items() if not allowed(k)}
            verify(self.out, protected)
            if any(k not in protected and not allowed(k) for k in current):
                raise ValueError('Unexpected artifacts outside interrupted MD workspace')
            archive = self.directory / ('interrupted-' + uuid.uuid4().hex)
            archive.mkdir()
            for name in current:
                if allowed(name):
                    dest = archive / name
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(self.out / name, dest)
        else:
            if not active.get('restoring'):
                verify(self.out, self.state['artifacts'])
            archive = self.directory / ('interrupted-' + uuid.uuid4().hex)
            archive.mkdir()
            for name in inventory(self.out):
                dest = archive / name
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(self.out / name, dest)
            active['restoring'] = True
            self.save()
            for child in self.out.iterdir():
                if child.name not in EXCLUDED:
                    if child.is_dir():
                        shutil.rmtree(child)
                    else:
                        child.unlink()
            for name in self.state['artifacts']:
                dest = self.out / name
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(baseline / name, dest)
        failed = self.out / 'FAILED.json'
        if failed.exists():
            shutil.copy2(failed, archive / 'FAILED.json')
        # Keep active until a successful commit, including during recovery itself.

    def phase(self, name, callback):
        if name in self.state['completed']:
            print('Verified; skipping:', name, flush=True)
            return
        if self.state['active'] and self.state['active']['phase'] != name:
            raise ValueError('Stage order differs from interrupted run')
        if not self.state['active']:
            baseline_name = 'baseline-' + uuid.uuid4().hex
            baseline = self.directory / baseline_name
            baseline.mkdir()
            verify(self.out, self.state['artifacts'])
            for path in self.state['artifacts']:
                dest = baseline / path
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(self.out / path, dest)
            self.state['active'] = {'phase': name, 'baseline': baseline_name}
        self.state['active'].pop('restoring', None)
        self.state['status'] = 'running'
        self.save()
        try:
            callback()
            verify(self.out, self.state['artifacts'])
            artifacts = inventory(self.out)
        except BaseException as exc:
            self.state['status'] = 'interrupted' if isinstance(exc, (KeyboardInterrupt, InterruptedError)) else 'failed'
            self.save()
            write(self.out / 'FAILED.json', {'phase': name, 'status': self.state['status'],
                  'error': str(exc), 'type': type(exc).__name__})
            raise
        self.state['artifacts'] = artifacts
        self.state['completed'].append(name)
        self.state['active'] = None
        self.save()

    def finish(self):
        failed = self.out / 'FAILED.json'
        if failed.exists():
            failed.replace(self.directory / ('failure-' + uuid.uuid4().hex + '.json'))
        self.state['status'] = 'completed'
        self.save()
