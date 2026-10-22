"""Read-only view of an explicitly configured, frozen benchmark run."""
import csv
import datetime
import json
import os
from pathlib import Path

FILES={'validation_report.md','metrics.json','reference_recovery.csv','all_rankings.csv',
       'molecule_stage_ledger.csv','assay_evidence.csv','scientific_top50.csv','scientific_top50.sdf',
       'deliverable_top50.csv','deliverable_top50.sdf','run_manifest.json','environment.json',
       'target_audit.json','recovery_metrics.png','property_matching.png','blinded_library.sdf','scaffold_recall.csv'}

def root():
    value=os.environ.get('YUANHE_BENCHMARK_ROOT')
    if not value:return None
    path=Path(value).resolve()
    if not path.is_dir() or not (path/'run_manifest.json').is_file():raise ValueError('Configured benchmark run is unavailable')
    return path

def read(path,default=None):
    return json.loads(path.read_text()) if path.is_file() else default

def ready(path):
    return (path/'metrics.json').is_file() and (path/'evaluation/ranking_seal.json').is_file() and (path/'validation_report.md').is_file()

def snapshot():
    path=root()
    if path is None:return {'configured':False}
    manifest=read(path/'run_manifest.json');progress=read(path/'progress.json',{})
    frozen=read(path/'inputs/preregistered_protocol.json',{})
    end=manifest.get('finished_at') or datetime.datetime.now(datetime.timezone.utc).isoformat()
    elapsed=(datetime.datetime.fromisoformat(end)-datetime.datetime.fromisoformat(manifest['started_at'])).total_seconds()/60
    result={'configured':True,'run_id':path.name,'status':manifest['status'],'progress':progress,
            'elapsed_minutes':elapsed,'hard_budget_minutes':manifest['hard_budget_minutes'],
            'hard_deadline':manifest['hard_deadline'],'population':frozen.get('input_parent_count'),
            'preparation_failures':len(frozen.get('preparation_failures',[])),
            'reused_jobs':read(path/'reuse_manifest.json',{}).get('imported_jobs',0),
            'scientific_git_commit':frozen.get('git_commit'),'results_ready':ready(path),
            'scope':'回顾性虚拟筛选验证；不进入MD；未知背景不是实验阴性'}
    # No labels, drug names, scores or recovery metrics are opened before evaluation.
    if result['results_ready']:
        result['metrics']=read(path/'metrics.json')
        with (path/'reference_recovery.csv').open() as fh:result['references']=list(csv.DictReader(fh))
        result['files']=sorted(name for name in FILES if (path/name).is_file())
    return result

def artifact(name):
    path=root()
    if path is None or not ready(path):raise ValueError('Results are not ready for unblinding')
    if name not in FILES:raise ValueError('Unknown benchmark artifact')
    target=(path/name).resolve()
    if target.parent!=path or not target.is_file():raise ValueError('Artifact unavailable')
    return target
