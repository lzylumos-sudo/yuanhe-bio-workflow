#!/usr/bin/env python3
"""Portable repository integrity checks; does not validate scientific engines."""
import ast
import json
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]

def check(root=ROOT):
    names = subprocess.check_output(
        ['git', 'ls-files', '-z', '--cached', '--others', '--exclude-standard'],
        cwd=root, text=True).split('\0')
    errors = []
    files = [root / n for n in sorted(set(filter(None, names))) if (root / n).is_file()]
    for path in files:
        name = path.relative_to(root).as_posix()
        if re.search(r'(^|/)(\.env|id_rsa|id_ed25519)$|\.(pem|key)$', name):
            errors.append(f'{name}: private configuration/key filename')
        try:
            source = path.read_text(encoding='utf-8-sig')
        except UnicodeDecodeError:
            continue
        if path.suffix == '.md':
            # Ignore code fences and inline code such as chemical SMILES.
            prose = re.sub(r'```.*?```|`[^`\n]*`', '', source, flags=re.S)
            for url in re.findall(r'\]\(([^\s)]+)\)', prose):
                if '://' in url or url.startswith(('#', 'mailto:')):
                    continue
                target = unquote(url.split('#')[0])
                if target and not (path.parent / target).exists():
                    errors.append(f'{name}: broken local link {url}')
            if re.search(r'\b(?:competition|hackathon)\b', prose, re.I):
                errors.append(f'{name}: outdated event framing')
        if name == 'deployment/validation_5070/INSTALL_AND_VALIDATE.txt':
            if re.search(r'\bcompetition\b', source, re.I):
                errors.append(f'{name}: outdated event framing')
    server = root / '00_workbench/webui/server.py'
    version = next(node.value.value for node in ast.parse(server.read_text()).body
                   if isinstance(node, ast.Assign) and any(
                       isinstance(t, ast.Name) and t.id == 'VERSION' for t in node.targets))
    if f'Release {version}' not in (root/'README.md').read_text():
        errors.append('README.md: release differs from UI/API version')
    if f'Version {version}' not in (root/'00_workbench/webui/static/index.html').read_text():
        errors.append('index.html: footer differs from UI/API version')
    for name in ('LICENSE', 'LICENSE_SCOPE.txt', 'THIRD_PARTY_NOTICES.txt'):
        if (root/name).read_bytes() != (root/'00_workbench/webui/static/licenses'/('LICENSE.txt' if name == 'LICENSE' else name)).read_bytes():
            errors.append(f'{name}: UI license copy differs')
    for item in json.loads((root/'02_docking/docs/licenses/sources.json').read_text()):
        import hashlib
        path = root/'02_docking'/item['path']
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']:
            errors.append(f"{item['path']}: missing or modified upstream license")
    for error in errors:
        print(error)
    print(json.dumps({'checked_files': len(files), 'errors': len(errors), 'version': version}))
    return bool(errors)

if __name__ == '__main__':
    sys.exit(check())
