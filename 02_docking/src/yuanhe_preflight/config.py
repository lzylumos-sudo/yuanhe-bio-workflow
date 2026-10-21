import sys
from pathlib import Path
from .toxicity_substructure_filter import DEFAULT_RULES

DEFAULT_PATH=Path(__file__).with_name('configs')/'predocking.example.json'

def defaults(root):
    from .common import read
    cfg=read(DEFAULT_PATH)
    cfg['rules_path']=str(DEFAULT_RULES.resolve());cfg['cache_dir']=str(Path(root)/'data/predocking_cache')
    cfg['retro']['python']=sys.executable
    return cfg

def validate(cfg):
    required={'rules_path','risk_policies','risk_exceptions','catalogs','stock_catalogs','evidence_max_age_days','retro','cache_dir','workers','batch_size','route_reviews','external_interfaces'}
    if not isinstance(cfg,dict) or set(cfg)!=required:raise ValueError('predocking_config_fields_mismatch')
    if cfg['external_interfaces']:raise ValueError('External interfaces are not implemented; use authorized local catalog imports. No structures are transmitted.')
    for key,lo,hi in [('workers',1,16),('batch_size',1,1000),('evidence_max_age_days',1,365)]:
        if type(cfg[key]) is not int or not lo<=cfg[key]<=hi:raise ValueError(key+' out of range')
    for key,lo,hi in [('time_limit_seconds',1,86400),('max_depth',1,30),('iteration_limit',1,1000000),('startup_grace_seconds',0,600)]:
        value=cfg['retro'].get(key)
        if type(value) is not int or not lo<=value<=hi:raise ValueError('retro.'+key+' out of range')
    for k in ('risk_exceptions','route_reviews','catalogs','stock_catalogs'):
        if not isinstance(cfg[k],list):raise ValueError(k+' must be an array')
    return cfg
