"""Audit a documented CPU10 -> CPU16 XML-state transfer and new RNG epoch.

This module never advances a simulation or modifies its provenance. The config
projection excludes only MD.cpu_threads; RNG continuity is audited separately
and must explicitly be RESEEDED, never represented as binary continuation.
"""
from __future__ import annotations
import copy
import hashlib
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET


class StateTransferAuditError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise StateTransferAuditError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def config_digest(config):
    # Exactly the existing modules.utils.core.config_hash serialization.
    return hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()


def science_projection(config):
    projected = copy.deepcopy(config)
    require(isinstance(projected.get('MD'), dict), 'Configuration has no MD mapping.')
    require('cpu_threads' in projected['MD'], 'Configuration has no cpu_threads.')
    del projected['MD']['cpu_threads']
    return config_digest(projected)


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def read_config(path):
    import yaml
    result = yaml.safe_load(Path(path).read_text(encoding='utf-8'))
    require(isinstance(result, dict), 'Configuration snapshot is not a mapping.')
    return result


def safe_path(module, value, within):
    require(isinstance(value, str) and value and '\\' not in value and ':' not in value,
            'Artifact paths must be module-relative POSIX paths.')
    relative = Path(value)
    require(not relative.is_absolute() and '..' not in relative.parts,
            'Artifact path is absolute or traverses parents.')
    resolved = (module / relative).resolve()
    require(resolved.is_relative_to(within.resolve()), 'Artifact escapes its transfer directory.')
    require(resolved.is_file(), 'Referenced migration artifact is missing: ' + value)
    return resolved


def local_artifact(directory, filename):
    require(isinstance(filename, str) and Path(filename).name == filename and
            filename not in ('', '.', '..') and '\\' not in filename,
            'Checkpoint/segment names must be plain filenames.')
    result = (directory / filename).resolve()
    require(result.is_relative_to(directory.resolve()) and result.is_file(),
            'Committed trajectory artifact is missing or escapes its replica.')
    return result


def timestep_ps(config):
    match = re.fullmatch(r'\s*([0-9]+(?:\.[0-9]+)?)\s*(fs|ps|ns)\s*', str(config['MD']['timestep']))
    require(match is not None, 'Invalid timestep in migration configuration.')
    value = float(match[1]) * {'fs': .001, 'ps': 1., 'ns': 1000.}[match[2]]
    require(math.isfinite(value) and value > 0, 'Invalid timestep in migration configuration.')
    return value


def physical_state(path):
    root = ET.parse(path).getroot()
    require(root.tag == 'State', 'XML snapshot is not an OpenMM State.')
    result = {'step': int(root.attrib['stepCount']), 'time_ps': float(root.attrib['time'])}
    require(math.isfinite(result['time_ps']), 'Nonfinite State time.')
    for name in ('PeriodicBoxVectors', 'Positions', 'Velocities'):
        element = root.find(name)
        require(element is not None and len(element) > 0, 'State lacks ' + name)
        values = tuple(tuple(float(row.attrib[axis]) for axis in ('x', 'y', 'z')) for row in element)
        require(all(math.isfinite(v) for row in values for v in row), 'Nonfinite State ' + name)
        result[name] = values
    require(len(result['Positions']) == len(result['Velocities']), 'State position/velocity count differs.')
    require(len(result['PeriodicBoxVectors']) == 3, 'State does not have three box vectors.')
    params = root.find('Parameters')
    require(params is not None, 'State lacks Context parameters.')
    result['parameters'] = {key: float(value) for key, value in params.attrib.items()}
    require(all(math.isfinite(v) for v in result['parameters'].values()), 'Nonfinite Context parameter.')
    require(result['parameters'].get('k_restraint') == 0., 'Transfer is not an unrestrained production state.')
    return result


def verify_state_transfers(run_folder, current_config, manifest, *, current_progress=None):
    """Strictly verify zero or one documented XML-state / reseeded RNG transfer.

    Binary compatibility and actual checkpoint reloading are tested by the
    separately checksummed probe and the existing final CPU16 checkpoint audit.
    Physical XML states and all pre-transfer committed files are checked here.
    """
    folder = Path(run_folder).resolve()
    require(folder.parent.name == 'results', 'Run folder must be below module/results.')
    module = folder.parent.parent
    run_id = manifest.get('run_id')
    require(run_id == folder.name, 'Manifest run identity differs from the requested run.')
    history = manifest.get('state_transfers', [])
    require(isinstance(history, list), 'state_transfers must be a list.')
    implementation_hash = digest(Path(__file__))
    if not history:
        require(current_progress is None, 'A candidate override requires explicit transfer history.')
        # Removing the history must not hide a changed runtime seed.
        for replica in range(1, current_config['replicates'] + 1):
            progress = read_json(folder / 'md' / f'replica_{replica:02d}' / 'progress.json')
            require(progress['seed'] == current_config['MD']['seed'] + replica,
                    'Runtime RNG seed changed without documented state-transfer history.')
        return {'status': 'NOT_APPLICABLE', 'verified_transfers': 0,
                'audit_implementation_sha256': implementation_hash, 'source_sha256': {}}
    require(len(history) == 1, 'This audit version supports exactly one explicit state transfer.')
    entry = history[0]
    require(isinstance(entry, dict), 'Invalid state-transfer history entry.')
    record_path = safe_path(module, entry.get('record'), folder / 'state_transfers')
    require(digest(record_path) == entry.get('sha256'), 'Transfer record checksum differs from manifest.')
    record = read_json(record_path)
    require(record.get('schema_version') == 1, 'Unsupported state-transfer schema.')
    require(record.get('mode') == 'XML_PHYSICAL_STATE_NEW_RNG_EPOCH', 'Undocumented binary/RNG transfer mode.')
    require(record.get('run_id') == run_id, 'Transfer references another run.')
    epoch = record.get('epoch_id')
    require(isinstance(epoch, str) and re.fullmatch(r'[A-Za-z0-9_-]+', epoch), 'Invalid epoch identifier.')
    transfer_dir = folder / 'state_transfers' / epoch
    require(record_path.parent == transfer_dir.resolve(), 'Record location differs from its epoch.')
    require(record.get('status') in ('CONFIG_APPLIED', 'RESUMED'), 'Transfer was not applied.')
    authorization = record.get('authorization', {})
    require(authorization.get('authorized') is True and authorization.get('source') == 'user' and
            type(authorization.get('cpu_threads')) is int and authorization['cpu_threads'] == 16,
            'Explicit user authorization for CPU16 is missing.')
    required_snapshots = ('config_before', 'config_after', 'manifest_before', 'progress_before',
                          'progress_after', 'source_checkpoint', 'source_xml', 'native16_checkpoint',
                          'native16_xml', 'probe_report')
    snapshots = record.get('snapshots', {})
    paths = {}
    source_hashes = {entry['record']: entry['sha256']}
    for name in required_snapshots:
        ref = snapshots.get(name, {})
        path = safe_path(module, ref.get('path'), transfer_dir)
        require(digest(path) == ref.get('sha256'), 'Snapshot checksum failed: ' + name)
        paths[name] = path
        source_hashes[ref['path']] = ref['sha256']
    require(len(set(paths.values())) == len(paths), 'Distinct snapshots alias the same file.')
    before_config, after_config = read_config(paths['config_before']), read_config(paths['config_after'])
    before_manifest = read_json(paths['manifest_before'])
    before_progress, after_progress = read_json(paths['progress_before']), read_json(paths['progress_after'])
    before, after = record.get('before', {}), record.get('after', {})
    require(before.get('threads') == 10 and after.get('threads') == 16, 'Only the documented CPU10 -> CPU16 transfer is supported.')
    require(type(before_config['MD']['cpu_threads']) is int and before_config['MD']['cpu_threads'] == 10 and
            type(after_config['MD']['cpu_threads']) is int and after_config['MD']['cpu_threads'] == 16,
            'Configuration thread counts disagree with transfer metadata.')
    require(before_config['MD']['platform'] == after_config['MD']['platform'] == 'CPU', 'Transfer must remain on CPU.')
    require(config_digest(before_config) == before.get('config_hash') == before_manifest.get('config_sha256'),
            'Original configuration hash was changed or misrepresented.')
    require(config_digest(after_config) == after.get('config_hash') == manifest.get('config_sha256') == config_digest(current_config),
            'Current configuration does not match the actual post-transfer hash.')
    projection = science_projection(before_config)
    require(projection == science_projection(after_config) == record.get('physical_protocol_hash'),
            'Configuration changed outside MD.cpu_threads.')
    require(before_manifest.get('run_id') == run_id and before_manifest.get('state_transfers', []) == [],
            'Original manifest lineage is not a first transfer of this run.')
    for field in ('input_sha256', 'md_implementation_sha256'):
        require(before_manifest.get(field) and before_manifest[field] == manifest.get(field),
                'Input or frozen implementation provenance changed during transfer: ' + field)
    for stage in ('prepare', 'build'):
        prior = before_manifest.get('stages', {}).get(stage, {})
        current = manifest.get('stages', {}).get(stage, {})
        require(prior.get('status') == current.get('status') == 'complete' and
                prior.get('files') == current.get('files'), 'Committed preparation/build history changed.')
    replica = record.get('replica', 1)
    require(type(replica) is int and replica == 1 and current_config['replicates'] == 1,
            'This transfer audit currently supports the single-replica public demo.')
    target = folder / 'md' / f'replica_{replica:02d}'
    candidate_only = current_progress is not None
    if candidate_only:
        require(current_progress == after_progress, 'Candidate progress differs from the frozen post-transfer snapshot.')
        for live, snapshot in ((folder / 'config.yaml', 'config_before'),
                               (folder / 'manifest.json', 'manifest_before'),
                               (target / 'progress.json', 'progress_before')):
            require(digest(live) == digest(paths[snapshot]),
                    'Candidate audit requires the original on-disk configuration/manifest/progress to remain unchanged.')
        progress = current_progress
    else:
        progress = read_json(target / 'progress.json')
    require(before_progress.get('phase') == after_progress.get('phase') == 'production',
            'State transfer must occur during production.')
    require(before_progress.get('platform') == after_progress.get('platform') == progress.get('platform') == 'CPU',
            'The trajectory must remain on the CPU platform.')
    changed = {key for key in set(before_progress) | set(after_progress)
               if before_progress.get(key) != after_progress.get(key)}
    require(changed <= {'checkpoint', 'checkpoint_sha256', 'seed', 'updated', 'state_transfer_epoch'},
            'Post-transfer progress changed outside checkpoint/RNG epoch metadata.')
    if 'state_transfer_epoch' in after_progress:
        require(after_progress['state_transfer_epoch'] == epoch, 'Progress epoch identifier disagrees.')
    boundary = before_progress.get('production_steps')
    require(type(boundary) is int and boundary > 0 and before_progress.get('current_step') == boundary and
            before.get('production_steps') == after.get('production_steps') == boundary,
            'Transfer production step boundary is inconsistent.')
    require(before_progress.get('seed') == before.get('seed') == before_config['MD']['seed'] + replica,
            'Original RNG seed does not match the original run.')
    new_seed = after.get('seed')
    require(type(new_seed) is int and 1 <= new_seed <= 2147383647 and new_seed != before['seed'],
            'New runtime seed must be distinct, positive and safe for seed+100000.')
    require(after_progress.get('seed') == progress.get('seed') == new_seed, 'New runtime RNG seed is undocumented or changed.')
    rng = record.get('rng', {})
    require(rng.get('continuity') == 'RESEEDED' and rng.get('barostat_adaptation') == 'RESET' and
            rng.get('integrator_seed') == new_seed and rng.get('barostat_seed') == new_seed + 100000,
            'RNG/adaptive-state reset is not explicitly and consistently disclosed.')
    committed = before_progress.get('segments')
    require(isinstance(committed, list) and committed, 'Transfer has no committed production segments.')
    require(all(type(s.get('steps')) is int and s['steps'] > 0 for s in committed) and
            sum(s['steps'] for s in committed) == boundary, 'Original committed steps do not sum to boundary.')
    require(progress.get('segments', [])[:len(committed)] == committed and progress.get('production_steps', -1) >= boundary,
            'Current trajectory does not preserve the original committed segment prefix.')
    for segment in committed:
        for kind in ('dcd', 'csv'):
            path = local_artifact(target, segment[kind])
            require(digest(path) == segment[kind + '_sha256'], 'A pre-transfer committed segment changed.')
            source_hashes[path.relative_to(module).as_posix()] = segment[kind + '_sha256']
    for label, state, binary_name, xml_name in (
            ('original', before_progress, 'source_checkpoint', 'source_xml'),
            ('new', after_progress, 'native16_checkpoint', 'native16_xml')):
        checkpoint = local_artifact(target, state['checkpoint'])
        xml = local_artifact(target, Path(state['checkpoint']).with_suffix('.xml').name)
        require(digest(checkpoint) == state['checkpoint_sha256'] == digest(paths[binary_name]),
                label + ' checkpoint was overwritten or differs from the immutable snapshot.')
        require(digest(xml) == digest(paths[xml_name]), label + ' XML state changed.')
        source_hashes[checkpoint.relative_to(module).as_posix()] = digest(checkpoint)
        source_hashes[xml.relative_to(module).as_posix()] = digest(xml)
    require(before_progress['checkpoint'] != after_progress['checkpoint'], 'Source checkpoint was overwritten by the new epoch.')
    require(before_progress.get('system_sha256') == after_progress.get('system_sha256') == progress.get('system_sha256') ==
            digest(folder / 'md/system/system.xml'), 'System definition changed across the transfer.')
    original_state, new_state = physical_state(paths['source_xml']), physical_state(paths['native16_xml'])
    require(original_state == new_state, 'XML physical positions/velocities/box/time/parameters changed at transfer.')
    require(original_state['step'] == boundary and math.isclose(original_state['time_ps'], boundary * timestep_ps(current_config),
            rel_tol=1e-10, abs_tol=1e-5), 'Transferred physical time/step disagrees with production boundary.')
    probe = read_json(paths['probe_report'])
    require(record.get('checkpoint_probe') == probe, 'Embedded checkpoint probe differs from its checksummed report.')
    for flag in ('passed', 'positions_exact', 'velocities_exact', 'box_exact', 'parameters_exact',
                 'new_seed_confirmed_after_XML_load', 'native16_roundtrip_verified',
                 'native16_roundtrip_exact', 'public_integrator_parameters_equal',
                 'old_binary_never_loaded_into_new_thread_count'):
        require(probe.get(flag) is True, 'Physical-state transfer probe is not verified: ' + flag)
    require(probe.get('actual_threads') == 16 and probe.get('native16_advanced_steps') == 0 and
            probe.get('native16_rng_diagnostic', {}).get('streams') == 16,
            'Native bridge thread count, zero-step boundary or per-thread RNG stream count is invalid.')
    require(probe.get('new_threads') == 16 and probe.get('step') == boundary and
            math.isclose(float(probe.get('time_ps', -1)), original_state['time_ps'], rel_tol=0, abs_tol=1e-7),
            'Checkpoint probe does not describe the transferred CPU16 boundary.')
    require(probe.get('source_checkpoint_sha256') == digest(paths['source_checkpoint']) and
            probe.get('source_xml_sha256') == digest(paths['source_xml']) and
            probe.get('native16_checkpoint_sha256') == digest(paths['native16_checkpoint']),
            'Checkpoint probe refers to different binary/XML artifacts.')
    require(probe.get('integrator_seed') == new_seed and probe.get('barostat_seed') == new_seed + 100000,
            'Actual new-context seeds differ from recorded RNG epoch.')
    return {'status': 'PASS', 'candidate_only': candidate_only, 'verified_transfers': 1, 'epoch_id': epoch, 'boundary_steps': boundary,
            'before_threads': 10, 'after_threads': 16, 'before_seed': before['seed'], 'after_seed': new_seed,
            'physical_protocol_hash': projection, 'source_sha256': source_hashes,
            'audit_implementation_sha256': implementation_hash,
            'scope': 'Continuous recorded physical state/time and unchanged force-field/integration settings; a new RNG epoch and reset barostat adaptive state are explicit. Not bitwise binary continuation.'}
