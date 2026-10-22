import json,sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from webui import benchmark as b

def test_ui_withholds_labels_until_sealed_report(tmp_path,monkeypatch):
    monkeypatch.setenv('YUANHE_BENCHMARK_ROOT',str(tmp_path))
    (tmp_path/'run_manifest.json').write_text(json.dumps({'status':'running','started_at':'2026-10-22T17:16:00+00:00','hard_budget_minutes':180,'hard_deadline':'2026-10-22T20:16:00+00:00'}))
    (tmp_path/'evaluation').mkdir();(tmp_path/'evaluation/labels.json').write_text('not readable as JSON')
    result=b.snapshot();assert not result['results_ready'] and 'references' not in result and 'metrics' not in result
    with pytest.raises(ValueError,match='not ready'):b.artifact('all_rankings.csv')
    (tmp_path/'metrics.json').write_text('{}');(tmp_path/'evaluation/ranking_seal.json').write_text('{}');(tmp_path/'validation_report.md').write_text('report')
    (tmp_path/'reference_recovery.csv').write_text('name,rank\nReference,1\n')
    assert b.snapshot()['references'][0]['rank']=='1'
    with pytest.raises(ValueError,match='Unknown'):b.artifact('../evaluation/labels.json')
    assert b.artifact('validation_report.md').read_text()=='report'

def test_ui_unconfigured_does_not_read_any_run(monkeypatch):
    monkeypatch.delenv('YUANHE_BENCHMARK_ROOT',raising=False)
    assert b.snapshot()=={'configured':False}
