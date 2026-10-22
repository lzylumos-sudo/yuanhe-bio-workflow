"""Shared layer resolution and provenance; native scientific validators are retained."""
import copy
import hashlib
import json
from decimal import Decimal
from pathlib import Path
from .registry import load
from .version import ROOT


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()

def merge_layers(layers):
    result = {}; sources = {}; overrides = []
    def merge(dst, src, layer, prefix=''):
        if not isinstance(src, dict): raise ValueError('Configuration layer must be an object')
        for key, value in src.items():
            name = prefix + key
            if isinstance(value, dict) and isinstance(dst.get(key), dict):
                merge(dst[key], value, layer, name + '.')
            else:
                if key in dst and dst[key] != value:
                    overrides.append({'field': name, 'previous': copy.deepcopy(dst[key]), 'resolved': copy.deepcopy(value), 'source': layer})
                dst[key] = copy.deepcopy(value); sources[name] = layer
    for name, values in layers: merge(result, values, name)
    return {'config': result, 'sources': sources, 'overrides': overrides, 'sha256': fingerprint(result)}

def protocol_for(config):
    rc = config.get('research', {})
    version = '1.1.0' if rc.get('predocking',{}).get('engine') else '1.0.0'
    protocol = next(p for p in load('protocol_registry')['protocols'] if p['protocol_id'] == 'murj-discovery' and p['protocol_version']==version)
    def equivalent(key, expected, actual):
        # Freezing a project-relative location does not alter its scientific meaning.
        if key == 'predocking' and isinstance(actual, dict):
            expected, actual = copy.deepcopy(expected), copy.deepcopy(actual)
            for field in ('rules_path', 'cache_dir'):
                for value in (expected, actual):
                    if value.get(field): value[field] = str((ROOT / Path(value[field])).resolve())
        return expected == actual
    differences = {k: {'protocol': v, 'resolved': rc.get(k)} for k, v in protocol['locked'].items() if not equivalent(k, v, rc.get(k))}
    differences.update({k: {'protocol': v, 'resolved': config.get(k, v)} for k, v in protocol.get('workflow_locked', {}).items() if config.get(k, v) != v})
    return {'protocol_id': protocol['protocol_id'] if not differences else 'custom',
            'protocol_version': protocol['protocol_version'] if not differences else 'unknown',
            'base_protocol': protocol['protocol_id'] + '@' + protocol['protocol_version'],
            'conforms': not differences, 'differences': differences,
            'scope': 'orchestration recipe only; target, model and MD system protocols are recorded separately',
            'locked_fields': list(protocol['locked'])}

def integration_steps(duration_ns, timestep_fs):
    duration = Decimal(str(duration_ns)); step = Decimal(str(timestep_fs))
    if not duration.is_finite() or not step.is_finite() or duration <= 0 or step <= 0:
        raise ValueError('Duration and timestep must be finite and positive (ns, fs)')
    count = duration * Decimal(1000000) / step
    if count != count.to_integral_value(): raise ValueError('Duration is not an integer number of integration steps')
    return int(count)

def stage_config(config, reviews):
    # Existing stage review semantics are shallow top-level updates, not recursive merges.
    result = copy.deepcopy(config); sources = {k: 'run' for k in result}; overrides = []
    for stage, values in reviews:
        for key, value in values.items():
            if result.get(key) != value: overrides.append({'field': key, 'previous': result.get(key), 'resolved': value, 'source': 'review:' + stage})
            result[key] = copy.deepcopy(value); sources[key] = 'review:' + stage
    return {'config': result, 'sources': sources, 'overrides': overrides, 'sha256': fingerprint(result),
            'protocol': protocol_for(result)}

def write_resolved(path, record):
    # JSON syntax is valid YAML 1.2. Keeps the stdlib HTTP service dependency-free.
    from .records import atomic
    atomic(path, record)
