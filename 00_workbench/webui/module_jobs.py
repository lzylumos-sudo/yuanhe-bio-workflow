"""Independent UI jobs. Inputs and execution are isolated from research qualification."""
from __future__ import annotations
import base64
import binascii
import math
import os
from pathlib import Path
import re
import signal
import threading
import time
import uuid
from . import workflow as w

JOBS = Path(os.environ.get('YUANHE_MODULE_JOBS', str(w.UI/'module_jobs'))).resolve()
LIMIT = 64*1024*1024
SCOPE = 'Independent computation; no automatic admission to the qualified research pipeline'


def field(key, label, value, kind='text', **extra):
    return dict(key=key, label=label, default=value, type=kind, **extra)


COMMON = [field('seed', 'Random seed', 17, 'integer', min=1, max=2147483647)]
BOX = [field('center', 'Box center / angstrom (x, y, z)', [0, 0, 0], 'vector'),
       field('size', 'Box size / angstrom (x, y, z)', [20, 20, 20], 'vector', positive=True)]
MODULES = {
    'prepare': dict(label='Complex preparation', environment='chem_python', files={'complex': '.pdb'},
        fields=[field('receptor_chain', 'Receptor chain', 'A'), field('peptide_chain', 'Peptide chain', 'B'),
                field('first_residue', 'First receptor residue', 3, 'integer', min=-999, max=99999),
                field('last_residue', 'Last receptor residue', 500, 'integer', min=-999, max=99999),
                field('pocket_radius', 'Pocket radius / angstrom', 6, 'number', min=1, max=20)],
        description='Upload a protein/peptide complex and extract the pocket. No generation required.'),
    'peptide2mol': dict(label='Peptide2Mol generation', environment='chem_python', files={'pocket': '.pdb'},
        fields=COMMON+[field('budget', 'Sampling attempts', 1, 'integer', min=1, max=1000)],
        description='Upload an already prepared pocket PDB. Uses the pinned real model and counts actual attempts.'),
    'reinvent4': dict(label='REINVENT4 generation', environment='chem_python', files={},
        fields=COMMON+[field('scaffold', 'Reviewed scaffold SMILES with one [*:1] attachment', ''),
                       field('budget', 'Sampling attempts', 1, 'integer', min=1, max=1000)],
        description='Supply an explicit scaffold; no Peptide2Mol output is required.'),
    'risk': dict(label='Substructure risk screening', environment='chem_python', files={'library': '.smi,.sdf,.csv'},
        fields=[], description='Upload a molecule library. Existing rules and exclusions remain active.'),
    'availability': dict(label='Pre-docking availability screening', environment='chem_python', files={'library': '.smi,.sdf,.csv'},
        fields=[], description='Upload a molecule library. Runs both strict risk and availability gates, including synthesis-route evidence.'),
    'docking': dict(label='Candidate docking', environment='chem_python', files={'receptor': '.pdb', 'ligand': '.sdf'},
        fields=COMMON+BOX+[field('review_reason', 'Receptor, ligand state and box review rationale', '')],
        description='Upload one ligand and receptor. Recomputes both strict pre-docking gates before Vina; rejected candidates are not docked.'),
    'docking_control': dict(label='Public 1DRF docking control', environment='chem_python', files={},
        fields=[], description='Fixed public 1DRF/FOL reference with three seeds. Produces reference evidence, never candidate qualification.'),
    'screening': dict(label='Post-docking virtual screening', environment='screen_python',
        files={'receptor': '.pdb', 'poses': '.sdf'},
        fields=COMMON+BOX+[field('review_reason', 'Receptor, poses and original docking box review rationale', '')],
        description='Upload receptor and scored poses in the same coordinate frame, with the reviewed original docking box. Recomputes the existing risk, toxicity and synthesis-route gates before redocking; only accepted chemical states proceed.'),
    'interaction': dict(label='Interaction analysis', environment='chem_python', files={'complex': '.pdb', 'template': '.sdf?'},
        fields=[field('ligand', 'Ligand selector, e.g. BEN or A:BEN1', '')],
        description='Upload a complex and optional ligand chemical template. Calculates interactions without prior docking.'),
    'md': dict(label='MD: construct a membrane system', environment='md_python',
        files={'receptor': '.pdb', 'ligand': '.sdf'},
        fields=COMMON+[field('center', 'Membrane center / angstrom (x, y, z)', [0, 0, 0], 'vector'),
            field('normal', 'Membrane normal (x, y, z)', [0, 0, 1], 'vector'),
            field('orientation_source', 'Membrane orientation source and review rationale', ''),
            field('lipid', 'Lipid', 'POPE', 'choice', choices=['POPE', 'POPC', 'DPPC', 'DOPC']),
            field('production_ns', 'Production duration / ns (independent scope)', 500, 'number', min=.01, max=1000),
            field('equilibration_ns', 'Equilibration duration / ns', 2, 'number', min=.1, max=100),
            field('temperature_K', 'Temperature / K', 310, 'number', min=273, max=350)],
        description='Upload receptor and ligand directly. Build, inspect artifacts, review the system, then run and analyze with the existing MD engine.'),
    'md_tpr': dict(label='MD: run a prepared GROMACS input', environment='md_python',
        files={'production': '.tpr'},
        fields=[field('device', 'Execution device', 'gpu', 'choice', choices=['gpu', 'cpu']),
                field('threads', 'CPU threads', 4, 'integer', min=1, max=64),
                field('review_reason', 'System, force field and equilibration review rationale', '')],
        description='Upload a self-contained production TPR. Runs its recorded duration and verifies final step, checkpoint and completion; interrupted runs can resume.'),
}


def folder(ident):
    if not isinstance(ident, str) or not re.fullmatch(r'MOD-[a-f0-9]{32}', ident):
        raise ValueError('Invalid module job identifier')
    return JOBS/ident


def spec(kind):
    if kind not in MODULES:
        raise ValueError('Unknown independent module')
    return MODULES[kind]


def parameters(kind, supplied):
    if not isinstance(supplied, dict):
        raise ValueError('Parameters must be an object')
    fields = spec(kind)['fields']
    if set(supplied)-{f['key'] for f in fields}:
        raise ValueError('Unknown module parameter')
    result = {}
    for f in fields:
        v = supplied.get(f['key'], f['default']); t = f['type']; key = f['key']
        if t in {'integer', 'number'}:
            if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
                raise ValueError('Finite number required: '+key)
            if t == 'integer' and type(v) is not int:
                raise ValueError('Integer required: '+key)
            if not f['min'] <= v <= f['max']:
                raise ValueError('Parameter outside allowed range: '+key)
        elif t == 'vector':
            if not isinstance(v, list) or len(v) != 3 or any(type(n) not in {float, int} or not math.isfinite(n) for n in v):
                raise ValueError('Three finite numbers required: '+key)
            if f.get('positive') and any(n <= 0 or n > 100 for n in v):
                raise ValueError('Box dimensions must be positive and <= 100 angstrom')
            if key == 'normal' and sum(n*n for n in v) < 1e-12:
                raise ValueError('Membrane normal must be nonzero')
        elif not isinstance(v, str) or not v.strip() or len(v) > 4000:
            raise ValueError('Nonempty text required: '+key)
        if t == 'choice' and v not in f['choices']:
            raise ValueError('Unsupported choice: '+key)
        result[key] = v
    if kind == 'prepare' and result['first_residue'] > result['last_residue']:
        raise ValueError('Residue range is reversed')
    if kind == 'reinvent4':
        from rdkit import Chem
        mol = Chem.MolFromSmiles(result['scaffold'])
        points = [] if mol is None else [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
        if len(points) != 1 or mol.GetNumAtoms() <= 1 or not (points[0].GetIsotope() or points[0].GetAtomMapNum()):
            raise ValueError('Scaffold requires exactly one labelled attachment, for example [*:1]')
    if kind == 'interaction' and not re.fullmatch(r'[A-Za-z0-9_:.-]{1,40}', result['ligand']):
        raise ValueError('Invalid ligand selector')
    return result


def create(payload):
    kind = payload.get('module'); definition = spec(kind)
    if payload.get('acknowledged') is not True:
        raise ValueError('Acknowledge independent computation scope')
    p = parameters(kind, payload.get('parameters', {}))
    supplied = payload.get('files', {})
    if not isinstance(supplied, dict) or set(supplied)-set(definition['files']):
        raise ValueError('Unknown input file role')
    inputs = {}; total = 0
    for role, allowed in definition['files'].items():
        entry = supplied.get(role)
        if not entry and allowed.endswith('?'):
            continue
        if not isinstance(entry, dict):
            raise ValueError('Input file required: '+role)
        suffix = Path(str(entry.get('name', ''))).suffix.lower()
        if suffix not in allowed.rstrip('?').split(','):
            raise ValueError('Unsupported input format: '+role)
        try:
            content = base64.b64decode(entry.get('base64', ''), validate=True)
        except (ValueError, binascii.Error, TypeError) as e:
            raise ValueError('Invalid file encoding') from e
        if not content or len(content) > LIMIT:
            raise ValueError('Each input must be nonempty and <= 64 MiB')
        total += len(content)
        if total > 96*1024*1024:
            raise ValueError('Total input size exceeds 96 MiB')
        if suffix != '.tpr':
            text = content.decode('utf-8-sig')
            if '\0' in text:
                raise ValueError('Text input contains binary data')
            if suffix == '.pdb' and not any(x.startswith(('ATOM  ', 'HETATM')) for x in text.splitlines()):
                raise ValueError('PDB contains no atom records')
        inputs[role] = (suffix, content)
    ident = 'MOD-'+uuid.uuid4().hex; d = folder(ident); d.mkdir(parents=True)
    manifest = {}
    for role, (suffix, content) in inputs.items():
        path = d/'inputs'/(role+suffix); path.parent.mkdir(exist_ok=True); path.write_bytes(content)
        manifest[role] = {'path': path.relative_to(d).as_posix(), 'sha256': w.digest(path), 'size': len(content)}
    cfg = w.defaults()
    cfg.pop('research', None)
    from yuanhe_preflight.config import absolute_paths
    cfg['predocking'] = absolute_paths(cfg['predocking'], w.ROOT)
    request = {'schema': 'yuanhe.module-job/1', 'module': kind, 'parameters': p, 'inputs': manifest,
               'config': cfg, 'scope': SCOPE, 'candidate_qualified': False}
    w.atomic(d/'request.json', request)
    w.atomic(d/'source.json', w.source_snapshot())
    w.atomic(d/'job.json', {'id': ident, 'module': kind, 'created': time.time(),
        'request_sha256': w.digest(d/'request.json'), 'source_sha256': w.digest(d/'source.json')})
    return get(ident)


def verified(ident):
    d = folder(ident); meta = w.read(d/'job.json')
    if not meta:
        raise ValueError('Module job not found')
    if w.digest(d/'request.json') != meta['request_sha256'] or w.digest(d/'source.json') != meta['source_sha256']:
        raise ValueError('Module job provenance changed')
    req = w.read(d/'request.json')
    for row in req['inputs'].values():
        path = (d/row['path']).resolve()
        if not path.is_relative_to((d/'inputs').resolve()) or w.digest(path) != row['sha256']:
            raise ValueError('Module input changed')
    return d, req


def get(ident):
    d, req = verified(ident)
    ops = []
    for path in sorted((d/'operations').glob('*/operation.json'), reverse=True):
        rec = w.read(path)
        if rec['status'] == 'running' and (ident, 'module', path.parent.name) not in w._active and not w.live_record(rec, ident):
            rec = {**rec, 'status': 'interrupted', 'error': 'Worker no longer exists; inspect retained logs and checkpoints'}
        rec['id'] = path.parent.name; rec['result'] = w.read(path.parent/'result.json')
        log = path.parent/'process.log'
        if log.exists():
            with log.open('rb') as f:
                f.seek(max(0, log.stat().st_size-20000)); rec['log'] = f.read().decode('utf-8', 'replace')
        ops.append(rec)
    artifacts = [{'path': p.relative_to(d).as_posix(), 'size': p.stat().st_size}
                 for p in sorted(d.rglob('*')) if p.is_file() and not p.is_symlink() and '__pycache__' not in p.parts]
    return {**w.read(d/'job.json'), 'scope': SCOPE, 'candidate_qualified': False,
            'parameters': req['parameters'], 'inputs': req['inputs'], 'operations': ops,
            'artifacts': artifacts, 'systems': [w.read(p) for p in sorted((d/'md_systems').glob('*/system.json'))]}


def inventory():
    from .scheduling import exclusive_status
    return {'modules': {k: {**v, 'scope': SCOPE} for k, v in MODULES.items()},
            'resource_status': exclusive_status(w.ROOT),
            'jobs': sorted((w.read(p) for p in JOBS.glob('*/job.json')), key=lambda row: row['created'], reverse=True)[:100]}


def launch(ident, payload):
    from .scheduling import ResourceBusyError
    d, req = verified(ident)
    if w.read(d/'source.json') != w.source_snapshot():
        raise ValueError('Source changed; create a new independent job')
    kind = req['module']; action = payload.get('action', 'build' if kind == 'md' else 'run')
    if action not in ({'build', 'run', 'analyze'} if kind == 'md' else {'run'}):
        raise ValueError('Unsupported module action')
    system = payload.get('system_id')
    if kind == 'md' and action != 'build':
        from scientific.membrane import system_dir
        meta = w.read(system_dir(d, system)/'system.json', {})
        if meta.get('status') not in ({'completed'} if action == 'analyze' else {'built', 'completed'}):
            raise ValueError('A built system is required; analysis requires completed MD')
        if payload.get('system_reviewed') is not True or not str(payload.get('review_reason', '')).strip():
            raise ValueError('Review the system artifacts and record a rationale before running')
    with w._lock:
        if any(key[0] == ident for key in w._active):
            raise ResourceBusyError('Module job already has an active worker')
        history = get(ident)['operations']
        if any(x['status'] == 'running' for x in history):
            raise ResourceBusyError('A recorded worker is still running')
        if any(x['action'] == action and x['status'] == 'completed' for x in history):
            raise ValueError('This action completed; create a new job for a new experiment')
        if history and kind not in {'md', 'md_tpr'}:
            raise ValueError('Create a new job to preserve the previous experiment')
        if action == 'build' and history:
            raise ValueError('Create a new job after a failed build; partial artifacts are retained')
        lease = w.acquire_compute()
        out = None
        try:
            py = req['config'].get(spec(kind)['environment']) or req['config']['chem_python']
            if not Path(py).is_file():
                raise ValueError('Configure the module environment before execution')
            if kind == 'reinvent4':
                from yuanhe.runtime import reinvent_environment_check
                check = reinvent_environment_check(req['config'])
                if not check['ok']:
                    raise ValueError(check['detail'])
            op = str(time.time_ns()); out = d/'operations'/op; out.mkdir(parents=True)
            w.atomic(out/'operation.json', {'action': action, 'status': 'running', 'started': time.time(),
                      'system_id': system, 'system_reviewed': payload.get('system_reviewed') is True,
                      'review_reason': str(payload.get('review_reason', ''))[:4000]})
            from .readiness import runtime_env
            threads = req['parameters']['threads'] if kind == 'md_tpr' else 12 if kind == 'md' else 4
            environment = {**runtime_env(py), 'OMP_NUM_THREADS': str(threads),
                           'OPENMM_CPU_THREADS': str(threads), 'YUANHE_MD_THREADS': str(threads)}
            log = (out/'process.log').open('w', encoding='utf-8')
            w.spawn_job(ident, 'module', op, out/'operation.json',
                [py, str(w.UI/'module_worker.py'), '--run-id', ident, '--operation', op], lease, log,
                environment)
        except BaseException as exc:
            lease.close()
            if out:
                w.atomic(out/'operation.json', {**w.read(out/'operation.json'), 'status': 'failed', 'error': str(exc)})
            raise
    return get(ident)


def cancel(ident):
    d, _ = verified(ident)
    with w._lock:
        for path in (d/'operations').glob('*/operation.json'):
            rec = w.read(path); proc = w._active.get((ident, 'module', path.parent.name))
            pid = proc.pid if proc and proc.poll() is None else rec.get('pid') if w.live_record(rec, ident) else None
            if pid:
                identity = w.process_identity(pid)
                try:
                    os.killpg(pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                def finish(pid=pid, identity=identity):
                    time.sleep(30)
                    if w.live_record(identity, ident):
                        try:
                            os.killpg(pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                threading.Thread(target=finish, daemon=True).start()
            if rec['status'] == 'running':
                w.atomic(path, {**rec, 'status': 'cancelled', 'ended': time.time()})
    return get(ident)
