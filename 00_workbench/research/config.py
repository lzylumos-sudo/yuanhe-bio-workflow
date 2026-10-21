"""Single definition used by UI/API/CLI and actual adapters."""
import copy, hashlib, json, math, platform, sys
from pathlib import Path
PATH = Path(__file__).resolve().parents[1] / 'configs/research_dual_blackwell.json'

def defaults(profile='research_dual_blackwell'):
    cfg = json.loads(PATH.read_text())
    if profile == 'smoke':
        cfg.update(profile='smoke', peptide_attempts=20, core_limit=2, reinvent_attempts=40,
                   generated_unique_limit=20, external_limit=20, initial_limit=40,
                   refine_limit=20, review_limit=8, md_limit=2, wetlab_min=1, wetlab_max=2,
                   generation_workers=1, vina_workers=1, vina_cpu=2, vina_thread_budget=2)
    elif profile != 'research_dual_blackwell':
        raise ValueError('Unknown research profile')
    return cfg

def resolve(value=None):
    value = value or {}
    cfg = defaults(value.get('profile', 'research_dual_blackwell'))
    if set(value) - set(cfg): raise ValueError('Unknown research configuration keys')
    for key, val in value.items():
        if key in ('initial', 'refine'):
            if not isinstance(val, dict) or set(val)-set(cfg[key]): raise ValueError('Unknown docking parameter')
            cfg[key].update(val)
        else: cfg[key] = copy.deepcopy(val)
    if cfg['version'] != '1.0.0': raise ValueError('Unsupported configuration version')
    for key in ('peptide_attempts','core_limit','reinvent_attempts','generated_unique_limit','external_limit',
                'initial_limit','refine_limit','review_limit','md_limit','md_replicas','md_concurrency',
                'wetlab_min','wetlab_max','generation_workers','vina_workers','vina_cpu','vina_thread_budget'):
        if type(cfg[key]) is not int or cfg[key] < 1: raise ValueError(key+' must be a positive integer')
    for key in ('md_production_ns','disk_reserve_gb'):
        if type(cfg[key]) not in (float,int) or not math.isfinite(cfg[key]) or cfg[key] <= 0: raise ValueError(key+' must be finite and positive')
    if cfg['md_replicas'] != 1: raise ValueError('One production trajectory per compound required')
    for low, high in [('core_limit','peptide_attempts'),('generated_unique_limit','reinvent_attempts'),
                      ('refine_limit','initial_limit'),('review_limit','refine_limit'),('md_limit','review_limit'),
                      ('wetlab_max','md_limit'),('wetlab_min','wetlab_max')]:
        if cfg[low] > cfg[high]: raise ValueError(low+' exceeds '+high)
    if cfg['initial_limit'] > cfg['generated_unique_limit']+cfg['external_limit']: raise ValueError('Initial pool exceeds source budgets')
    if cfg['vina_workers']*cfg['vina_cpu'] > cfg['vina_thread_budget']: raise ValueError('Vina CPU budget exceeded')
    if cfg['selection_mode'] not in ('diversity','score'): raise ValueError('Invalid selection mode')
    if cfg['reinvent_device'] not in ('auto','cpu','cuda'): raise ValueError('Invalid generator device')
    for phase in ('initial','refine'):
        for key, val in cfg[phase].items():
            if type(val) is not int or val < 1: raise ValueError(phase+'.'+key+' must be positive')
        if cfg[phase]['num_modes'] > 100: raise ValueError('Vina num_modes exceeds 100')
    if not isinstance(cfg['ligand_states'],dict):raise ValueError('ligand_states must map compound IDs to reviewed states')
    return cfg

def fingerprint(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

def plan(cfg):
    return {'parent_limits':{'initial':cfg['initial_limit'],'refine':cfg['refine_limit']},
            'docking_task_upper_bounds':{p:cfg['initial_limit' if p=='initial' else 'refine_limit'] *
                math.prod(cfg[p][k] for k in ('receptors','states','seeds')) for p in ('initial','refine')},
            'md_candidates':cfg['md_limit'],'trajectories_per_candidate':1,
            'production_ns_per_candidate':cfg['md_production_ns'],
            'total_production_ns':cfg['md_limit']*cfg['md_production_ns'],
            'gpu_concurrency_requested':cfg['md_concurrency'],'eta':None,
            'note':'Upper bounds, not promised qualified counts. ETA unknown until measured.'}

def environment():
    return {'python':sys.version,'platform':platform.platform(),'executable':sys.executable}
