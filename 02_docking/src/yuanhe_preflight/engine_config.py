"""Versioned policy for new tasks; absence preserves historical v1 evidence."""
import copy
import math

DEFAULT = {
    'version': '3.0.0', 'mode': 'balanced', 'keep': 2000, 'missing_policy': 'review',
    'structure_enabled': True, 'catalog_actions': {'PAINS': 'review', 'BRENK': 'review'},
    'catalog_rule_actions': {},
    'toxicity': {'python': '', 'device': 'cpu', 'batch_size': 128, 'cpu_threads': 4,
                 'timeout_seconds': 1800, 'model_version': '2.0.1',
                 'endpoints': ['AMES', 'hERG', 'DILI', 'LD50_Zhu'],
                 'thresholds': {'AMES': 0.8, 'hERG': 0.8, 'DILI': 0.8, 'LD50_Zhu': 3.0}},
    'synthesis': {'sa_max': 6.0, 'sc_max': 4.0, 'scscore_path': '', 'scscore_sha256': ''},
    'retro_max_candidates': 10,

}


def defaults():
    return copy.deepcopy(DEFAULT)


def validate(value):
    if not isinstance(value, dict) or set(value) != set(DEFAULT):
        raise ValueError('Invalid screening engine fields')
    if value['version'] != '3.0.0' or value['mode'] not in ('discovery', 'balanced', 'delivery'):
        raise ValueError('Invalid screening engine version or mode')
    if value['missing_policy'] not in ('abort', 'review'):
        raise ValueError('Missing predictions require abort or review')
    for key, lo, hi in [('keep', 1, 1000000), ('retro_max_candidates', 0, 1000)]:
        if type(value[key]) is not int or not lo <= value[key] <= hi:
            raise ValueError('Invalid ' + key)
    if type(value['structure_enabled']) is not bool:
        raise ValueError('structure_enabled must be boolean')
    for name in ('toxicity', 'synthesis'):
        if not isinstance(value[name], dict) or set(value[name]) != set(DEFAULT[name]):
            raise ValueError('Invalid ' + name + ' fields')
    for name in ('catalog_actions', 'catalog_rule_actions'):
        if not isinstance(value[name], dict) or any(a not in ('review', 'exclude', 'warn') for a in value[name].values()):
            raise ValueError('Invalid catalog policy')
    if set(value['catalog_actions']) != {'PAINS', 'BRENK'}:
        raise ValueError('PAINS and BRENK policies required')
    tox = value['toxicity']
    if tox['device'] not in ('cpu', 'cuda', 'auto') or tox['model_version'] != '2.0.1':
        raise ValueError('Unsupported ADMET-AI version or device')
    if not isinstance(tox['endpoints'], list) or not tox['endpoints'] or len(set(tox['endpoints'])) != len(tox['endpoints']):
        raise ValueError('Distinct toxicity endpoints required')
    for ep in tox['endpoints']:
        threshold = tox['thresholds'].get(ep)
        if type(threshold) not in (int, float) or not math.isfinite(threshold):
            raise ValueError('Finite endpoint threshold required')
        if ep in ('AMES', 'hERG', 'DILI') and not 0 <= threshold <= 1:
            raise ValueError('Probability threshold must be 0..1')
    for key, lo, hi in [('batch_size', 1, 1000), ('cpu_threads', 1, 64), ('timeout_seconds', 1, 86400)]:
        if type(tox[key]) is not int or not lo <= tox[key] <= hi:
            raise ValueError('Invalid toxicity.' + key)
    for key, lo, hi in [('sa_max', 1, 10), ('sc_max', 1, 5)]:
        x = value['synthesis'][key]
        if type(x) not in (int, float) or not math.isfinite(x) or not lo <= x <= hi:
            raise ValueError('Invalid ' + key)
    return value
