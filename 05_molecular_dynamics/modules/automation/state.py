"""Durable job records and artifact fingerprints; never infer completion from names."""
from __future__ import annotations
import hashlib, json, os, platform, socket
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

def now():
    return datetime.now(timezone.utc).isoformat()

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''): h.update(chunk)
    return h.hexdigest()

def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False).encode()).hexdigest()

def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))

def write(path, data):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f'.{os.getpid()}.tmp')
    with tmp.open('w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False, allow_nan=False)
        f.write('\n'); f.flush(); os.fsync(f.fileno())
    os.replace(tmp, path)

def update(job, status, **fields):
    path = Path(job)/'state.json'
    data = read(path) if path.exists() else {'schema_version': 2, 'created': now(), 'stages': {}}
    data.update(status=status, updated=now(), **fields)
    write(path, data)
    with (Path(job)/'events.jsonl').open('a', encoding='utf-8') as f:
        f.write(json.dumps({'time': now(), 'status':status, **fields}, ensure_ascii=False)+'\n')
    from . import governance
    records=governance()
    if records:records.sync_soluble(job,data)
    return data

def files_under(root, folders):
    root=Path(root); result={}
    for name in folders:
        path=root/name
        candidates=path.rglob('*') if path.is_dir() else [path]
        for p in candidates:
            if p.is_file() and '__pycache__' not in p.parts and not p.name.endswith('.tmp'):
                result[p.relative_to(root).as_posix()]=digest(p)
    return result

def verify(root, files):
    root=Path(root).resolve()
    for name,expected in files.items():
        p=(root/name).resolve()
        if not p.is_relative_to(root) or not p.is_file() or digest(p)!=expected:
            raise ValueError('Committed artifacts corrupted or changed; reuse refused:'+name)

@contextmanager
def job_lock(job):
    lock=Path(job)/'.running'
    try: lock.mkdir()
    except FileExistsError:
        raise LockBusy('Task is locked. Check the host and process recorded in owner.json ; do not delete checkpoints directly.')
    write(lock/'owner.json', {'pid':os.getpid(),'hostname':socket.gethostname(),'platform':platform.platform(),'started':now()})
    try: yield
    finally:
        (lock/'owner.json').unlink(missing_ok=True)
        lock.rmdir()

class LockBusy(ValueError):
    pass
