"""New-run provenance only. Never retrofit identities into historical results."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import uuid
from .version import ROOT, identity, environment
from .configuration import fingerprint, protocol_for

STATES = {'created', 'validated', 'queued', 'running', 'completed', 'failed', 'cancelled', 'blocked'}
REVIEW_STATES = {'not_reviewed', 'review_required', 'accepted', 'rejected'}

def now(): return datetime.now(timezone.utc).isoformat()

def new_run_id():
    return 'YH-' + datetime.now(timezone.utc).strftime('%Y%m%d') + '-' + uuid.uuid4().hex

def atomic(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    try:
        with tmp.open('w', encoding='utf-8') as f:
            json.dump(value, f, indent=2, allow_nan=False); f.write('\n'); f.flush(); os.fsync(f.fileno())
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()

def read(path, default=None):
    return json.loads(Path(path).read_text()) if Path(path).is_file() else default

def artifact(path, run_dir, stage, kind='unknown', content_hash=None):
    path = Path(path).resolve(); run_dir = Path(run_dir).resolve()
    if not path.is_relative_to(run_dir): raise ValueError('Artifact escapes run directory')
    name = path.relative_to(run_dir).as_posix(); sha = content_hash or digest(path)
    if kind=='unknown':
        kind={'.pdb':'structure','.sdf':'docking_pose' if stage=='docking' else 'ligand', '.top':'system_topology','.prmtop':'system_topology','.xtc':'trajectory','.dcd':'trajectory','.log':'execution_log'}.get(path.suffix,'file')
    return {'artifact_id': 'YHA-' + fingerprint([run_dir.name, name, sha]), 'artifact_type': kind,
            'run_id': run_dir.name, 'stage_id': stage, 'source': name, 'content_hash': sha,
            'schema_version': '1.0.0', 'validation_status': 'hash_recorded', 'storage_path': name}

def model_inventory():
    models = read(ROOT / '01_molecule_generation/configs/project.json', {}).get('models', {})
    result = {}
    sources=read(ROOT / 'deployment/model_dependencies.json', {})
    for key, value in models.items():
        name = value.get('checkpoint') or value.get('prior'); path = ROOT / '01_molecule_generation' / name
        # Hash installed model bytes; configured upstream versions alone are not verified versions.
        result[key] = {'source': name, 'source_url': sources.get(key,{}).get('source_url','unknown'), 'weights_source':sources.get(key,{}).get('weights_source','unknown'), 'declared_commit': value.get('repo_commit', 'unknown'),
                       'version': 'unknown', 'sha256': digest(path) if path.is_file() else 'unknown',
                       'expected_sha256': value.get('checkpoint_sha256', value.get('prior_sha256', 'unknown')),
                       'license': 'unknown', 'execution_verified': False}
    return result

def create(run_dir, config, requested, defaults):
    d = Path(run_dir)
    if (d / 'run_manifest.json').exists(): raise ValueError('Run manifest already exists')
    version = identity(); protocol = protocol_for(config); env = environment()
    inp = artifact(d / 'input.pdb', d, 'prepare', 'receptor')
    value = {'schema_version': '1.0.0', 'run_id': d.name, 'software_version': version['software_version'],
             'git_commit': version['git']['commit'], 'git_branch': version['git']['branch'], 'git_dirty': version['git']['dirty'],
             'protocol_id': protocol['protocol_id'], 'protocol_version': protocol['protocol_version'], 'protocol': protocol,
             'config_hash': fingerprint(config), 'target_id': config.get('target_profile', 'unknown'),
             'input_artifacts': [inp], 'model_versions': model_inventory(), 'tool_versions': {},
             'stage_status': {}, 'output_artifacts': [], 'start_time': None, 'end_time': None,
             'created_time': now(), 'status': 'created', 'scientific_review': 'not_reviewed',
             'environment': env, 'environment_hash': fingerprint(env),
             'resources': read(d / 'environment.json', {}).get('resources', {}),
             'data_versions': {'research_assets': read(d / 'research_assets.json', {}), 'target_input': inp}}
    from .contracts import validate_manifest
    validate_manifest(value)
    atomic(d / 'run_manifest.json', value)
    atomic(d / 'artifact_registry.json', {'schema_version': '1.0.0', 'artifacts': [inp]})
    atomic(d / 'resolved_config.yaml', {'schema_version': '1.0.0', 'config': config, 'sha256': fingerprint(config),
           'default': defaults, 'user_override': requested, 'protocol': protocol,
           'resolution_order': ['base_defaults', 'research_profile', 'task', 'stage_review'],
           'note': 'Existing engine-specific validators remain authoritative; stage resolved files record actual dispatch values.'})
    event(d, 'configuration_resolved', 'WORKBENCH', None, 'Run created; no computation started')
    return value

def event(d, event, module, stage, message, level='INFO'):
    path = Path(d) / 'events.jsonl'
    row = {'timestamp': now(), 'level': level, 'run_id': Path(d).name, 'stage_id': stage,
           'module_id': module, 'event': event, 'message': message}
    with path.open('a', encoding='utf-8') as f: f.write(json.dumps(row) + '\n')

def sync(d, stages):
    d = Path(d); value = read(d / 'run_manifest.json')
    if value is None: return None  # Historical runs retain their original identity.
    records = {s: read(d / s / 'stage.json', {}) for s in stages}
    value['stage_status'] = {s: r.get('status', 'created') for s, r in records.items()}
    states = set(value['stage_status'].values())
    value['status'] = ('running' if 'running' in states else 'failed' if 'failed' in states else
                       'cancelled' if 'cancelled' in states else 'blocked' if states & {'blocked', 'interrupted'} else
                       'completed' if states == {'completed'} else 'validated' if 'completed' in states else 'created')
    if any(r.get('started') for r in records.values()):
        value['start_time'] = datetime.fromtimestamp(min(r['started'] for r in records.values() if r.get('started')), timezone.utc).isoformat()
    ended = max((r.get('ended', 0) for r in records.values()), default=0)
    value['end_time'] = datetime.fromtimestamp(ended, timezone.utc).isoformat() if ended and value['status'] in {'completed', 'failed', 'cancelled', 'blocked'} else None
    artifacts = list(value['input_artifacts'])
    for stage, record in records.items():
        for name, sha in record.get('artifacts', {}).items():
            artifacts.append(artifact(d / name, d, stage, content_hash=sha))
    value['output_artifacts'] = artifacts[len(value['input_artifacts']):]
    value['tool_versions'] = {s: read(d / s / 'tool_versions.json') for s in stages if (d / s / 'tool_versions.json').is_file()}
    value['model_execution_records'] = {str(p.relative_to(d)): read(p) for s in ('peptide2mol','reinvent4') for p in sorted((d / s).glob('*/manifest.json'))}
    value['md_system_protocols'] = {p.parent.name: {'protocol': read(p).get('protocol', {'protocol_id':'unknown','protocol_version':'unknown'}), 'software_versions':read(p).get('software_versions',{}), 'source':str(p.relative_to(d)), 'sha256':digest(p)} for p in sorted((d / 'md_systems').glob('*/system.json'))}
    value['stage_execution_environments'] = {s: read(d / s / 'execution_environment.json') for s in stages if (d / s / 'execution_environment.json').is_file()}
    value['scientific_review'] = 'review_required' if value['status'] == 'completed' else 'not_reviewed'
    atomic(d / 'artifact_registry.json', {'schema_version': '1.0.0', 'artifacts': artifacts})
    from .contracts import validate_manifest
    validate_manifest(value)
    atomic(d / 'run_manifest.json', value)
    return value


def create_soluble(d, config, input_path):
    from .protocols import soluble
    from .contracts import validate_manifest
    d=Path(d);version=identity();protocol=soluble(config);env=environment()
    inp=artifact(input_path,d,'input','receptor')
    value={'schema_version':'1.0.0','run_id':d.name,'software_version':version['software_version'],
           'git_commit':version['git']['commit'],'git_branch':version['git']['branch'],'git_dirty':version['git']['dirty'],
           'protocol_id':protocol['protocol_id'],'protocol_version':protocol['protocol_version'],'protocol':protocol,
           'config_hash':fingerprint(config),'target_id':'unknown','input_artifacts':[inp],
           'model_versions':{},'tool_versions':{},'stage_status':{},'output_artifacts':[],
           'start_time':None,'end_time':None,'status':'created','scientific_review':'not_reviewed',
           'environment':env,'environment_hash':fingerprint(env),'created_time':now()}
    validate_manifest(value);atomic(d/'run_manifest.json',value)
    atomic(d/'resolved_config.yaml',{'config':config,'sha256':fingerprint(config),'protocol':protocol})
    atomic(d/'artifact_registry.json',{'schema_version':'1.0.0','artifacts':[inp]})


def sync_soluble(d,state):
    d=Path(d);value=read(d/'run_manifest.json')
    if value is None:return
    mapping={'NEW':'created','RUNNING':'running','COMPLETE':'completed','FAILED':'failed','CANCELLED':'cancelled',
             'WAITING_FOR_DECISION':'blocked','DEPENDENCY_BLOCKED':'blocked','PREPARED':'validated'}
    value['status']=mapping.get(state['status'],'blocked');value['native_status']=state['status']
    value['stage_status']={s:mapping.get(rec.get('status'),'blocked') for s,rec in state.get('stages',{}).items()}
    value['scientific_review']='review_required' if value['status']=='completed' else 'not_reviewed'
    if value['status']=='running' and not value['start_time']:value['start_time']=now()
    value['end_time']=now() if value['status'] in {'completed','failed','cancelled'} else None
    artifacts=[]
    for stage,rec in state.get('stages',{}).items():
        for name,sha in rec.get('files',{}).items():artifacts.append(artifact(d/name,d,stage,content_hash=sha))
    value['output_artifacts']=artifacts
    atomic(d/'artifact_registry.json',{'schema_version':'1.0.0','artifacts':value['input_artifacts']+artifacts})
    from .contracts import validate_manifest
    validate_manifest(value);atomic(d/'run_manifest.json',value)
