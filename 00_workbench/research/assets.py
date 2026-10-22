"""Freeze supplied scientific assets; absent optional assets stay explicitly absent."""
from pathlib import Path

def snapshot(root,cfg):
    from webui.workflow import digest
    rc=cfg['research'];paths=set()
    if rc.get('external_library',rc.get('supplier_library','')):paths.add(rc.get('external_library',rc.get('supplier_library','')))
    pc=rc['predocking'];paths.add(pc['rules_path'])
    sc=pc.get('engine',{}).get('synthesis',{}).get('scscore_path')
    if sc:paths.add(sc)
    for key in ('model_path','templates_path','stock_path'):
        if pc['retro'].get(key):paths.add(pc['retro'][key])
    def collect(value):
        if isinstance(value,dict):
            for key,val in value.items():
                if key in ('path','evidence_path') and isinstance(val,str) and val:paths.add(val)
                elif isinstance(val,(list,dict)):collect(val)
        elif isinstance(value,list):
            for item in value:collect(item)
    collect(rc['ligand_states']);collect(rc['receptor_ensemble'])
    result={}
    for value in sorted(paths):
        p=Path(value)
        if not p.is_absolute():p=root/p
        if not p.is_file():raise ValueError('Configured scientific asset is missing: '+str(p))
        result[str(p.resolve())]=digest(p)
    return result

def verify(d):
    from webui.workflow import read,digest
    for name,sha in read(d/'research_assets.json',{}).items():
        if not Path(name).is_file() or digest(name)!=sha:raise ValueError('Frozen research asset changed: '+name)
