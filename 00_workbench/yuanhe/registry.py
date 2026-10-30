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
from .runtime import resolve_python, resolve_executable, runtime_env, vina_version, reinvent_environment_check

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
        for names in row.get('capability_backends', {}).values():
            if not names or set(names) - bids: errors.append('Unknown capability backend')
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
        p = subprocess.run([executable, '-c', script], capture_output=True, text=True, timeout=10,
                           env=runtime_env(executable))
        if p.returncode: return {'available': False, 'reason': 'Interpreter probe failed', 'packages': {}}
        values = json.loads(p.stdout)
        return {'available': all(values.values()), 'packages': values, 'executable': executable}
    except (OSError, subprocess.TimeoutExpired, ValueError):
        return {'available': False, 'reason': 'Interpreter unavailable', 'packages': {}}

def snapshot(config=None):
    config = config or {}; registry = copy.deepcopy(load()); probes = {}
    for row in registry['backends']:
        key = row['environment']; exe = resolve_python(config, key)
        probe_key = (exe, tuple(row['dependencies']))
        if probe_key not in probes: probes[probe_key] = probe_python(exe, row['dependencies'])
        probe = probes[probe_key]
        binaries = {name: resolve_executable(config, name, exe) for name in row['executables']}
        optional = {name: resolve_executable(config, name, exe) for name in row.get('optional_executables', [])}
        assets = {name: (ROOT / name).is_file() for name in row['required_files']}
        implemented = (ROOT / row['implementation_path']).exists()
        installed = implemented and probe['available'] and all(binaries.values()) and all(assets.values())
        configuration = reinvent_environment_check(config) if row['backend_id'] == 'reinvent4' else {'ok': True}
        # Import availability does not establish a working CUDA device, licensed model,
        # validated topology or reviewed run inputs. Only stdlib orchestration is ready here.
        runnable = bool(installed and row['backend_id'] == 'workbench' and row['status'] == 'active')
        if installed and row['backend_id']=='vina' and row['status']=='active':
            try:
                check=subprocess.run([binaries['vina'],'--version'],capture_output=True,text=True,timeout=5,
                                     env=runtime_env(exe))
                runnable=bool(check.returncode==0 and vina_version(check.stdout+' '+check.stderr)=='1.2.7')
                row['runtime_version_output']=(check.stdout+' '+check.stderr).strip()
            except (OSError,subprocess.TimeoutExpired):runnable=False
        row.update(implemented=implemented, installed=bool(installed), verified=False, runnable=runnable,
                   readiness='Available' if runnable else row['status'].title() if row['status'] != 'active' else
                   'Configuration Error' if installed and not configuration['ok'] else 'Needs Validation' if installed else 'Not Installed',
                   dependency_probe=probe, executable_paths=binaries, optional_executable_paths=optional,
                   assets_present=assets, resolved_environment={'key': key, 'python': exe},
                   configuration_check=configuration, configuration_valid=configuration['ok'],
                   runtime_scope='Probe host only; execution-host validation and reviewed inputs are separate')
    by_id = {b['backend_id']: b for b in registry['backends']}
    for row in registry['modules']:
        capabilities = {}
        for capability in row['capabilities']:
            alternatives = row.get('capability_backends', {}).get(capability)
            children = [by_id[name] for name in (alternatives or row['backend'])]
            combine = any if alternatives else all
            capabilities[capability] = {field: combine(b[field] for b in children)
                                       for field in ('implemented', 'installed', 'runnable', 'configuration_valid')}
            capabilities[capability].update(verified=False, backends=[b['backend_id'] for b in children])
        present = all(c['implemented'] for c in capabilities.values())
        installed = all(c['installed'] for c in capabilities.values())
        runnable = all(c['runnable'] for c in capabilities.values())
        partial = any(c['installed'] for c in capabilities.values())
        valid = all(c['configuration_valid'] for c in capabilities.values())
        row.update(name=row['display_name'], present=present, installed=installed, verified=False,
                   runnable=runnable, capability_status=capabilities,
                   configuration_valid=valid,
                   readiness='Available' if runnable else 'Configuration Error' if installed and not valid else 'Needs Validation' if installed else
                   'Partially Installed' if partial else 'Not Installed')
    return registry
