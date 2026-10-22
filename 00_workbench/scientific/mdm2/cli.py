"""Repository-native benchmark lifecycle; labels are evaluation-only."""
import argparse
import json
import sys
from pathlib import Path
from .common import ROOT

def main(argv=None):
    for part in ('02_docking/src','01_molecule_generation'):
        sys.path.insert(0,str(ROOT/part))
    p=argparse.ArgumentParser(prog='yuanhe benchmark')
    p.add_argument('action',choices=['preflight','acquire','prepare-target','curate','prepare-library','throughput','freeze','run','evaluate','status'])
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--vina',type=Path)
    p.add_argument('--started-at',help='ISO8601 task start including timezone; existing run start is immutable')
    p.add_argument('--select-only',action='store_true',help='Review library matching before spending 3D preparation time')
    p.add_argument('--count',type=int,default=550)
    args=p.parse_args(argv)
    try:
        if args.action=='preflight':
            from .common import preflight
            if not args.vina:p.error('--vina required for preflight')
            result=preflight(args.root,args.vina,args.started_at)
        elif args.action=='acquire':
            from .sources import acquire
            result=acquire(args.root)
        elif args.action=='prepare-target':
            from .structure import prepare
            result=prepare(args.root)
        elif args.action=='curate':
            from .curation import curate
            result=curate(args.root)
        elif args.action=='prepare-library':
            from .library import prepare
            result=prepare(args.root,args.count,args.select_only)
        elif args.action in {'throughput','freeze','run'}:
            from .execution import dispatch
            result=dispatch(args.root,args.action)
        elif args.action=='evaluate':
            from .evaluation import evaluate
            result=evaluate(args.root)
        else:
            result={p.name:json.loads(p.read_text()) for p in args.root.glob('*.json') if p.name in {'run_manifest.json','progress.json','throughput.json','metrics.json'}}
        print(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False));return 0
    except Exception as exc:
        import traceback
        traceback.print_exc();return 2
