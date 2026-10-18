import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from webui import workflow as w

@pytest.fixture
def stopped(tmp_path,monkeypatch):
 monkeypatch.setattr(w,'RUNS',tmp_path)
 d=tmp_path/'case';d.mkdir()
 w.atomic(d/'config.json',{'fixed':True})
 for name in ('prepare','peptide2mol','reinvent4','docking','screening'):
  w.atomic(d/name/'stage.json',{'status':'completed','config_sha256':w.digest(d/'config.json'),'artifacts':{}})
 initial={'01_molecule_generation/src/a.py':'same','03_virtual_screening/src/zvss/a.py':'old'}
 w.atomic(d/'source_snapshot.json',initial)
 current={**initial,'03_virtual_screening/src/zvss/a.py':'new'}
 monkeypatch.setattr(w,'source_snapshot',lambda:current)
 monkeypatch.setattr(w,'start',lambda rid,stage:{'id':rid,'stage':stage})
 return d,current

def test_revision_preserves_original_and_checks_upstream(stopped):
 d,current=stopped;original=w.digest(d/'source_snapshot.json');upstream=w.digest(d/'docking/stage.json')
 assert w.rescreen('case')['stage']=='screening'
 assert w.digest(d/'source_snapshot.json')==original
 assert w.digest(d/'docking/stage.json')==upstream
 assert len(list((d/'attempts').glob('screening_revision_*/stage.json')))==1
 w.verify_source(d)
 current['01_molecule_generation/src/a.py']='changed'
 with pytest.raises(ValueError):w.verify_source(d)

def test_revision_refuses_changed_inference(stopped):
 d,current=stopped;current['01_molecule_generation/src/a.py']='changed'
 with pytest.raises(ValueError,match='Upstream'):w.rescreen('case')
 assert (d/'screening/stage.json').is_file()

def test_revision_refuses_existing_md(stopped):
 d,_=stopped;w.atomic(d/'md_systems/MD_one/system.json',{'status':'built'})
 with pytest.raises(ValueError,match='Downstream'):w.rescreen('case')

def test_revision_refuses_tampered_upstream(stopped):
 d,_=stopped
 w.atomic(d/'docking/stage.json',{'status':'completed','config_sha256':'bad','artifacts':{}})
 with pytest.raises(RuntimeError):w.rescreen('case')

def test_rescreen_route_binds_run_id(monkeypatch):
 import json
 from webui import server
 monkeypatch.setattr(w,'rescreen',lambda rid:{'id':rid})
 status,_,payload=server.dispatch('POST','/api/project/run/case/rescreen',{}, {})
 assert status==202 and json.loads(payload)['id']=='case'

def test_md_continuation_preserves_all_upstream(stopped):
 d,current=stopped
 current['00_workbench/scientific/membrane.py']='repaired'
 before={k:w.digest(d/k/'stage.json') for k in ('prepare','peptide2mol','reinvent4','docking','screening')}
 result=w.continue_md('case')
 assert result['status']=='verified'
 assert before=={k:w.digest(d/k/'stage.json') for k in before}
 w.verify_source(d)
 current['02_docking/src/engine.py']='changed'
 with pytest.raises(ValueError):w.continue_md('case')

def test_md_continuation_rejects_tampered_screening(stopped):
 d,_=stopped
 w.atomic(d/'screening/stage.json',{'status':'completed','config_sha256':'changed'})
 with pytest.raises(RuntimeError):w.continue_md('case')
