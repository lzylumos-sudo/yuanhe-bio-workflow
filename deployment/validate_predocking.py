#!/usr/bin/env python3
"""Real CPU inference / 1000 unique NCI molecules; no simulated supply evidence.

Run in the ADMET-AI environment. Retro is optional and explicitly recorded.
This is execution acceptance, not a held-out clinical/experimental validation.
"""
import argparse,csv,json,resource,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'02_docking/src'))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--scscore',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--limit',type=int,default=1000)
    ap.add_argument('--retro-config',type=Path)
    ap.add_argument('--retro-python',default='')
    args=ap.parse_args()
    if not 1<=args.limit<=4500:raise ValueError('limit must be 1..4500')
    from rdkit import RDConfig
    from yuanhe_preflight.common import write,read,file_hash,canonical
    from yuanhe_preflight.config import defaults
    from yuanhe_preflight.pipeline import risk_stage,availability_stage
    from yuanhe_preflight.cli import import_candidates
    cfg=defaults(ROOT);out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    cfg['cache_dir']=str(out/'cache');engine=cfg['engine'];engine['toxicity']['python']=sys.executable
    engine['synthesis'].update(scscore_path=str(args.scscore.resolve()),scscore_sha256=file_hash(args.scscore))
    engine['retro_max_candidates']=0
    if args.retro_config:
        cfg['retro']=read(args.retro_config)
        cfg['retro']['python']=args.retro_python or sys.executable
        engine['retro_max_candidates']=2
    source=Path(RDConfig.RDDataDir)/'NCI/first_5K.smi';unique={}
    for line in source.read_text().splitlines():
        try:s=canonical(line.split()[0])
        except ValueError:continue
        unique.setdefault(s,line.split()[-1])
        if len(unique)==args.limit:break
    assert len(unique)==args.limit
    with (out/'input.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['id','smiles']);w.writeheader()
        w.writerows({'id':'NCI_'+n,'smiles':s} for s,n in unique.items())
    write(out/'config.json',cfg)
    started=time.monotonic();rows=import_candidates(out/'input.csv',out/'imported')
    a=risk_stage(rows,cfg,out/'risk',base=out,batch_id='NCI-real')
    b=availability_stage(read(out/'risk/records.json'),cfg,out/'availability')
    total=time.monotonic()-started
    report={'input':len(rows),'source':'RDKit Data/NCI/first_5K.smi','source_sha256':file_hash(source),
            'unique_valid':len(unique),'structure':a,'screening':b,'total_seconds':total,
            'molecules_per_second':len(rows)/total,'rss_unit':'bytes' if sys.platform=='darwin' else 'KiB',
            'max_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'child_max_rss':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
            'gpu':False,'docking_run':False,'md_run':False,'external_supplier_queries':False,
            'retro_enabled':args.retro_config is not None,'scientific_validation':'execution_only'}
    write(out/'acceptance.json',report)
    assert b['timings']['model_evaluated']>0 and b['timings']['model_failed']==0
    scored=[r for r in read(out/'availability/records.json') if r.get('screening')]
    assert all(r['synthesis']['prediction_status']=='SUCCESS' for r in scored)
    print(json.dumps({k:v for k,v in report.items() if k not in ('structure','screening')},indent=2))

if __name__=='__main__':main()
