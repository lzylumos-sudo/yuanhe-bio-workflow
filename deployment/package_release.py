#!/usr/bin/env python3
"""Package only a committed Git tree with source provenance and file hashes."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[1]

def build(output, ref='HEAD', root=ROOT):
    output = Path(output)
    if output.exists(): raise FileExistsError(output)
    def git(*args):
        return subprocess.check_output(['git', '-C', str(root), *args])
    commit = git('rev-parse', '--verify', '--end-of-options', ref + '^{commit}').decode().strip()
    payload = git('archive', '--format=zip', commit)
    with zipfile.ZipFile(io.BytesIO(payload)) as source:
        version = source.read('VERSION').decode().strip()
        if not re.fullmatch(r'\d+\.\d+\.\d+(?:-[0-9A-Za-z.-]+)?', version): raise ValueError('Invalid version')
        files = {n: hashlib.sha256(source.read(n)).hexdigest() for n in source.namelist() if not n.endswith('/')}
        manifest = {'schema': 'yuanhe.source-package/1', 'software_version': version, 'git_commit': commit,
                    'git_tree': git('rev-parse', commit + '^{tree}').decode().strip(), 'files': files}
        prefix = 'yuanhe-bio-platform-' + version + '/'
        output.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(output, 'x', zipfile.ZIP_DEFLATED) as target:
            for entry in source.infolist():
                data = source.read(entry.filename)
                entry.filename = prefix + entry.filename
                target.writestr(entry, data)
            target.writestr(prefix + 'SOURCE_MANIFEST.json', json.dumps(manifest, indent=2) + '\n')
    return {'software_version': version, 'git_commit': commit, 'files': len(files),
            'sha256': hashlib.sha256(output.read_bytes()).hexdigest(), 'archive': output.name}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True); p.add_argument('--ref', default='HEAD')
    a = p.parse_args(); print(json.dumps(build(a.output, a.ref)))

if __name__ == '__main__': main()
