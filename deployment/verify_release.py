#!/usr/bin/env python3
"""Verify the named administrator package before installing environments."""
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]

def verify(root=ROOT):
    root = Path(root).resolve()
    value = json.loads((root / 'SOURCE_MANIFEST.json').read_text())
    if value.get('schema') != 'yuanhe.source-package/1': raise ValueError('Unknown package schema')
    if not re.fullmatch('[0-9a-f]{40}', value.get('git_commit', '')): raise ValueError('Invalid source commit')
    if (root / 'VERSION').read_text().strip() != value['software_version']: raise ValueError('Software version mismatch')
    for name, expected in value['files'].items():
        path = (root / name).resolve()
        if not path.is_relative_to(root) or not path.is_file(): raise ValueError('Missing or unsafe file: ' + name)
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected: raise ValueError('Modified file: ' + name)
    return {'verified': True, 'software_version': value['software_version'],
            'source_commit': value['git_commit'], 'verified_files': len(value['files'])}

if __name__ == '__main__': print(json.dumps(verify(), indent=2))
