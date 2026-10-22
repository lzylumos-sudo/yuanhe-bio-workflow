"""Validated registry mappings and conservative environment probes; no computation."""
import ast
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from .version import ROOT

STATUSES = {'active', 'experimental', 'deprecated', 'archived'}

def load(name='module_registry'):
    return json.loads((ROOT / 'configs' / (name + '.yaml')).read_text())

def validate():
    value = load(); errors = []; seen = set()
    for group, key in (('modules', 'module_id'), ('backends', 'backend_id')):
        for row in value[group]:
            ident = row[key]
            if ident in seen: errors.append('Duplicate identity: ' + ident)
            seen.add(ident)
            if row['status'] not in STATUSES: errors.append('Invalid status: ' + ident)
            name, _, symbol = row['entrypoint'].partition(':') if group == 'modules' else (row['implementation_path'], ':', row['entrypoint'])
            path = ROOT / (name or row['implementation_path'])
            if not path.exists(): errors.append('Missing implementation: ' + ident)
            elif symbol and path.is_file():
                definitions = {n.name for n in ast.parse(path.read_text()).body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))}
                definitions.update(a.asname or a.name for n in ast.parse(path.read_text()).body if isinstance(n, ast.ImportFrom) for a in n.names)
                if symbol not in definitions: errors.append('Missing symbol: ' + ident + ':' + symbol)
    bids = {b['backend_id'] for b in value['backends']}
    pids = {p['protocol_id'] + '@' + p['protocol_version'] for p in load('protocol_registry')['protocols']}
    for row in value['modules'] + value['backends']:
        if set(row.get('backend', [])) - bids: errors.append('Unknown backend: ' + row.get('module_id', ''))
        if set(row.get('supported_protocols', [])) - pids: errors.append('Unknown protocol')
        if row.get('parent') and row['parent'] not in bids: errors.append('Unknown parent')
    schema_ids={json.loads(p.read_text())['$id'] for p in (ROOT/'configs/schemas').glob('*.json')}
    for row in value['modules']+value['backends']:
        for field in ('input_schema','output_schema'):
            if row[field] not in schema_ids:errors.append('Unknown '+field+': '+row[field])
    if errors: raise ValueError('; '.join(errors))
    return {'valid': True, 'modules': len(value['modules']), 'backends': len(value['backends'])}

def probe_python(executable, packages):
    if not executable: return {'available': False, 'reason': 'Interpreter not configured', 'packages': {}}
    script = "import importlib, json; names=" + repr(packages) + "\ndef check(n):\n try: importlib.import_module(n); return True\n except Exception: return False\nprint(json.dumps({n:check(n) for n in names}))"
    try:
        p = subprocess.run([executable, '-c', script], capture_output=True, text=True, timeout=10)
        if p.returncode: return {'available': False, 'reason': 'Interpreter probe failed', 'packages': {}}
        values = json.loads(p.stdout)
        return {'available': all(values.values()), 'packages': values, 'executable': executable}
    except (OSError, subprocess.TimeoutExpired, ValueError):
        return {'available': False, 'reason': 'Interpreter unavailable', 'packages': {}}

def snapshot(config=None):
    config = config or {}; registry = copy.deepcopy(load()); probes = {}
    for row in registry['backends']:
        key = row['environment']; exe = sys.executable if key == 'current' else config.get(key, '')
        if key in ('chem_python', 'screen_python') and not exe: exe = config.get('chem_python') or sys.executable
        probe_key = (exe, tuple(row['dependencies']))
        if probe_key not in probes: probes[probe_key] = probe_python(exe, row['dependencies'])
        probe = probes[probe_key]
        binaries = {name: shutil.which(config.get(name) or name) for name in row['executables']}
        assets = {name: (ROOT / name).is_file() for name in row['required_files']}
        implemented = (ROOT / row['implementation_path']).exists()
        installed = implemented and probe['available'] and all(binaries.values()) and all(assets.values())
        # Import availability does not establish a working CUDA device, licensed model,
        # validated topology or reviewed run inputs. Only stdlib orchestration is ready here.
        runnable = bool(installed and row['backend_id'] == 'workbench' and row['status'] == 'active')
        if installed and row['backend_id']=='vina' and row['status']=='active':
            try:
                check=subprocess.run([binaries['vina'],'--version'],capture_output=True,text=True,timeout=5)
                import re
                match=re.search(r'Vina\s+v?(\d+\.\d+\.\d+)',check.stdout+' '+check.stderr,re.I)
                runnable=bool(check.returncode==0 and match and match[1]=='1.2.7')
                row['runtime_version_output']=(check.stdout+' '+check.stderr).strip()
            except (OSError,subprocess.TimeoutExpired):runnable=False
        row.update(implemented=implemented, installed=bool(installed), verified=False, runnable=runnable,
                   readiness='Available' if runnable else row['status'].title() if row['status'] != 'active' else 'Needs Validation' if installed else 'Not Installed',
                   dependency_probe=probe, executable_paths=binaries, assets_present=assets,
                   runtime_scope='Probe host only; execution-host validation and reviewed inputs are separate')
    for row in registry['modules']:
        children = [b for b in registry['backends'] if b['backend_id'] in row['backend']]
        row.update(name=row['display_name'], present=all(b['implemented'] for b in children),
                   installed=all(b['installed'] for b in children), verified=False,
                   runnable=all(b['runnable'] for b in children),
                   readiness='Available' if all(b['runnable'] for b in children) else 'Needs Validation' if all(b['installed'] for b in children) else 'Not Installed')
        row['capability_status'] = {c: {'implemented': row['present'], 'installed': row['installed'],
                                      'verified': False, 'runnable': row['runnable']} for c in row['capabilities']}
    return registry
