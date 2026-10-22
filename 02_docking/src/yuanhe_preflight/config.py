import json
import copy
import sys
from pathlib import Path

# Configuration must remain usable before the chemistry runtime is installed.
# Importing the rule evaluator here used to import RDKit in the HTTP process.
DEFAULT_RULES = Path(__file__).with_name('configs') / 'risk_rules.v1.json'

DEFAULT_PATH=Path(__file__).with_name('configs')/'predocking.example.json'

def absolute_paths(cfg, root):
    """Resolve deployment paths before workers change their working directory."""
    cfg=copy.deepcopy(cfg)
    def resolve(value):
        if isinstance(value,dict):
            for key,item in value.items():
                if key in {'rules_path','cache_dir','model_path','templates_path','stock_path','path','evidence_path','search_artifact','scscore_path'} and isinstance(item,str) and item:
                    value[key]=str((Path(root)/item).resolve())
                elif isinstance(item,(dict,list)):resolve(item)
        elif isinstance(value,list):
            for item in value:resolve(item)
    resolve(cfg)
    if '/' in cfg['retro'].get('python',''):
        cfg['retro']['python']=str((Path(root)/cfg['retro']['python']).absolute())
    if cfg.get('engine',{}).get('toxicity',{}).get('python') and '/' in cfg['engine']['toxicity']['python']:
        cfg['engine']['toxicity']['python']=str((Path(root)/cfg['engine']['toxicity']['python']).absolute())
    return cfg

def defaults(root):
    cfg=json.loads(DEFAULT_PATH.read_text(encoding='utf-8'))
    cfg['rules_path']=str(DEFAULT_RULES.resolve());cfg['cache_dir']=str(Path(root)/'data/predocking_cache')
    cfg['retro']['python']=sys.executable
    return cfg

def validate(cfg):
    required={'rules_path','risk_policies','risk_exceptions','retro','cache_dir','workers','batch_size','route_reviews'}
    if not isinstance(cfg,dict) or set(cfg)-{'engine'}!=required:raise ValueError('predocking_config_fields_mismatch')
    for key,lo,hi in [('workers',1,16),('batch_size',1,1000)]:
        if type(cfg[key]) is not int or not lo<=cfg[key]<=hi:raise ValueError(key+' out of range')
    for key,lo,hi in [('time_limit_seconds',1,86400),('max_depth',1,30),('iteration_limit',1,1000000),('startup_grace_seconds',0,600)]:
        value=cfg['retro'].get(key)
        if type(value) is not int or not lo<=value<=hi:raise ValueError('retro.'+key+' out of range')
    if 'engine' in cfg:
        from .engine_config import validate as validate_engine
        validate_engine(cfg['engine'])
    for k in ('risk_exceptions','route_reviews'):
        if not isinstance(cfg[k],list):raise ValueError(k+' must be an array')
    return cfg
