#!/usr/bin/env python3
"""Research CLI uses exactly the same task creation/start/review API as the UI."""
import argparse,json
from pathlib import Path
from research.config import resolve,plan
from webui import workflow

def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['preview','create','start','derive','status'])
    p.add_argument('--config',type=Path);p.add_argument('--run-id');p.add_argument('--stage');a=p.parse_args()
    cfg=json.loads(a.config.read_text()) if a.config else {}
    if a.action=='preview':result=plan(resolve(cfg.get('research',cfg)))
    elif a.action=='create':result=workflow.create(cfg)
    elif a.action=='start':result=workflow.start(a.run_id,a.stage)
    elif a.action=='derive':
        from research.lifecycle import derive
        result=derive(a.run_id,cfg)
    else:result=workflow.status(a.run_id)
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
