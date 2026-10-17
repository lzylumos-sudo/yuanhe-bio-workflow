"""Record actual source bytes; a declared commit is never treated as verification."""
import hashlib
import subprocess
from pathlib import Path

def source_fingerprint(repo, expected):
    repo = Path(repo)
    commit = None
    if (repo / '.git').exists():
        p = subprocess.run(['git', '-C', str(repo), 'rev-parse', 'HEAD'], capture_output=True, text=True, timeout=15)
        if p.returncode:
            raise RuntimeError('Cannot verify vendor Git HEAD: ' + p.stderr[-300:])
        commit = p.stdout.strip()
        if commit != expected:
            raise ValueError('Vendor commit mismatch: ' + commit)
    h = hashlib.sha256(); count = 0
    for p in sorted(repo.rglob('*')):
        rel = p.relative_to(repo)
        if any(x in {'.git','__pycache__','.pytest_cache','logs','priors','ckpts'} for x in rel.parts):
            continue
        if p.is_file() and p.suffix in {'.py','.yaml','.yml','.toml','.sh'}:
            h.update(rel.as_posix().encode()); h.update(b'\0'); h.update(hashlib.sha256(p.read_bytes()).digest()); count += 1
    return {'expected_commit': expected, 'verified_commit': commit, 'source_sha256': h.hexdigest(),
            'source_files': count, 'note': 'Fingerprint includes local patches; absent Git metadata is not commit verification.'}
