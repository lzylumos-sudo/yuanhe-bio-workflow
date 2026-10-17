#!/usr/bin/env python3
"""Check current tracked/distributable files, not historical Git objects."""
import ast
import hashlib
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
HAN = re.compile(r'[\u3400-\u9fff\uf900-\ufaff]')

def main():
    names = subprocess.check_output(['git', 'ls-files', '-z', '--cached', '--others', '--exclude-standard'], cwd=ROOT).decode().split('\0')
    errors = []
    count = 0
    for name in sorted(set(filter(None, names))):
        path = ROOT / name
        if not path.is_file():
            continue
        count += 1
        if HAN.search(name):
            errors.append(name + ': non-English filename')
        try:
            source = path.read_text(encoding='utf-8-sig')
        except UnicodeDecodeError:
            continue
        if HAN.search(source):
            errors.append(name + ': untranslated text')
        try:
            if path.suffix == '.py':
                ast.parse(source, filename=name)
            elif path.suffix == '.json':
                value = json.loads(source)
                if HAN.search(json.dumps(value, ensure_ascii=False)):
                    errors.append(name + ': escaped untranslated text')
        except (SyntaxError, ValueError) as exc:
            errors.append(name + ': ' + str(exc))
    provenance = json.loads((ROOT / 'deployment/license_provenance.json').read_text())
    if hashlib.sha256((ROOT / 'LICENSE').read_bytes()).hexdigest() != provenance['sha256']:
        errors.append('LICENSE: upstream digest mismatch')
    for error in errors:
        print(error)
    print(f'Checked {count} files; {len(errors)} errors.')
    return bool(errors)

if __name__ == '__main__':
    raise SystemExit(main())
