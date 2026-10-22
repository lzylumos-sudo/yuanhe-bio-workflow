"""v2 engineering tests. Any injected predictions/vendor rows are test fixtures.

Real model and 1000-molecule acceptance is deployment/validate_predocking.py.
"""
import copy
import json
import sys
from pathlib import Path
import pytest
from rdkit import Chem

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from preflight.test_preflight import catalog, sdf
from yuanhe_preflight.config import defaults, validate
from yuanhe_preflight.engine_config import defaults as engine_defaults
from yuanhe_preflight.common import read, write, hash_obj, seal
from yuanhe_preflight.structure import screen, standardize
from yuanhe_preflight.toxicity_substructure_filter import RiskRules
from yuanhe_preflight.pipeline import risk_stage, availability_stage
from yuanhe_preflight.engine import rank
from yuanhe_preflight.gate import require_ligand


def test_new_default_and_legacy_schema(tmp_path):
    cfg=defaults(tmp_path)
    assert cfg['engine']['mode']=='balanced'
    validate(cfg)
    del cfg['engine']
    validate(cfg)


@pytest.mark.parametrize('change',[{'mode':'unsafe'},{'keep':0},{'missing_policy':'pass'},
                                   {'structure_enabled':'false'},{'retro_max_candidates':-1}])
def test_invalid_policies_fail_closed(tmp_path,change):
    cfg=defaults(tmp_path);cfg['engine'].update(change)
    with pytest.raises(ValueError):validate(cfg)


def test_full_catalog_review_not_toxicity_and_configurable():
    engine=engine_defaults();rules=RiskRules()
    r=screen('Oc1ccccc1O',rules,[],engine)
    assert r['filter_status']=='REVIEW'
    assert r['severity']=='medium'
    assert any(x['rule_id'].startswith('PAINS:') for x in r['alerts'])
    assert 'toxicity_risk' not in r['alert_types']
    engine['catalog_actions']={'PAINS':'warn','BRENK':'warn'}
    assert screen('Oc1ccccc1O',rules,[],engine)['decision']=='pass'
    assert screen('CC(=O)Cl',rules,[],engine)['filter_status']=='BLOCK'


def test_salt_parent_mapping_retains_exact_stereo():
    s=standardize('C[C@H]([NH3+])C(=O)O.[Cl-]')
    assert s['salt_changed'] and len(s['fragments'])==2
    assert '@' in s['canonical_smiles'] and '.' in s['exact_smiles']
    assert s['canonical_smiles']!=standardize('C[C@@H]([NH3+])C(=O)O.[Cl-]')['canonical_smiles']


def test_failed_model_retains_all_records_and_exports(tmp_path):
    cfg=defaults(tmp_path);cfg['engine']['toxicity']['python']='/nonexistent/chemistry-python'
    mol=sdf(tmp_path,'CCO')
    rows=[{'id':'one','smiles':'CCO','sdf':str(mol)}, {'id':'duplicate','smiles':'CCO','sdf':str(mol)},
          {'id':'bad','smiles':'invalid-smiles','sdf':str(mol)}]
    risk_stage(rows,cfg,tmp_path/'risk',base=tmp_path,batch_id='fixture')
    availability_stage(read(tmp_path/'risk/records.json'),cfg,tmp_path/'availability')
    records=read(tmp_path/'availability/records.json')
    assert len(records)==3 and read(tmp_path/'availability/qualified.json')==[]
    identity_index=read(tmp_path/'availability/identity_evidence.json')
    assert len(identity_index)==1
    for record in records:
        assert 'canonical_smiles' in record and 'ranking' in record
        for prediction in record['toxicity_predictions']:
            assert {'raw_prediction','prediction_unit','threshold','error_message'} <= prediction.keys()
        assert 'synthesis_difficulty' in record['synthesis']
    real=next(x for x in records if x.get('screening'))
    assert all(x['prediction_status']=='ERROR' and x['raw_prediction'] is None for x in real['toxicity_predictions'])
    assert real['synthesis']['sa_score']>0 and real['synthesis']['sc_score'] is None
    assert (tmp_path/'availability/screening.sqlite3').is_file()
    assert (tmp_path/'availability/passed.smi').read_text()==''
    cfg['engine']['missing_policy']='abort'
    with pytest.raises(RuntimeError,match='Required prediction'):
        availability_stage(read(tmp_path/'risk/records.json'),cfg,tmp_path/'aborted')
    assert len(read(tmp_path/'aborted/records.json'))==3


def test_unknown_endpoint_is_unsupported(tmp_path):
    from yuanhe_preflight.toxicity import predict_batch
    cfg=engine_defaults()['toxicity'];cfg.update(python='/missing',endpoints=['UNKNOWN'],thresholds={'UNKNOWN':1})
    r=predict_batch([{'compound_id':'unit','smiles':'CCO'}],cfg,str(tmp_path),tmp_path)
    assert r[0][0]['prediction_status']=='UNSUPPORTED' and r[0][0]['raw_prediction'] is None


def test_multimode_ranking_independent_scales():
    rows=[{'compound_id':str(i),'smiles':s,'synthesis':{'sa_score':sa,'sc_score':sc},'availability':{'status':a}}
          for i,(s,sa,sc,a) in enumerate([('c1ccccc1',2,3,'purchase_supported'),('Cc1ccccc1',3,2,'purchase_supported'),('C1CCCCC1',4,4,'route_supported')])]
    engine=engine_defaults();engine['mode']='discovery'
    assert [r['compound_id'] for r in rank(copy.deepcopy(rows),engine)]==['0','2','1']
    engine['mode']='delivery'
    assert [r['compound_id'] for r in rank(copy.deepcopy(rows),engine)]==['0','1','2']


def test_final_selection_gate_and_actual_docking_input(tmp_path,monkeypatch):
    """Injected model/supply FIXTURES isolate the orchestration from inference."""
    from yuanhe_preflight import engine
    import pipeline_bridge as bridge
    cfg=defaults(tmp_path);cfg['engine']['keep']=1
    cfg['catalogs']=[catalog(tmp_path,['CCO','CCCO'])]
    weight=tmp_path/'fixture-weights';weight.write_text('TEST ONLY NOT MODEL WEIGHTS')
    from yuanhe_preflight.common import file_hash
    cfg['engine']['synthesis'].update(scscore_path=str(weight),scscore_sha256=file_hash(weight))
    def tox(rows,*args):
        return [[{'endpoint':ep,'raw_prediction':0.1,'prediction_status':'SUCCESS','risk_level':'BELOW_THRESHOLD',
                  'threshold':cfg['engine']['toxicity']['thresholds'][ep],'provenance':{'weights':{'SIMULATED': 'TEST_ONLY'}}}
                 for ep in cfg['engine']['toxicity']['endpoints']] for row in rows]
    monkeypatch.setattr(engine,'predict_batch',tox)
    monkeypatch.setattr(engine,'score_batch',lambda rows,settings:[{'sa_score':2,'sc_score':2,'screening_status':'PASS','score_version':{'sc_weights_sha256':file_hash(weight)}} for row in rows])
    rows=[{'id':s,'smiles':s,'sdf':str(sdf(tmp_path,s,str(i)))} for i,s in enumerate(['CCO','CCCO','CN(C)N=O'])]
    risk_stage(rows,cfg,tmp_path/'risk',base=tmp_path,batch_id='UNIT_ONLY')
    availability_stage(read(tmp_path/'risk/records.json'),cfg,tmp_path/'availability')
    qualified=read(tmp_path/'availability/qualified.json');assert len(qualified)==1
    assert len(list((tmp_path/'availability/evidence').glob('*.json')))==1
    assert len(read(tmp_path/'availability/records.json'))==3
    selected=qualified[0];require_ligand(selected['source_sdf'],selected['qualification_evidence'])
    # Enter the real workbench adapter, stop at its executable boundary, and
    # assert the filter gate was called only for selected molecules.
    monkeypatch.setattr(bridge,'require',lambda *a:None)
    import yuanhe_preflight.gate as gate
    visited=[];original=gate.require_ligand
    monkeypatch.setattr(gate,'require_ligand',lambda path,proof: (visited.append(path),original(path,proof))[1])
    class BoundaryReached(Exception):pass
    monkeypatch.setattr(bridge,'executable',lambda *a: (_ for _ in ()).throw(BoundaryReached()))
    (tmp_path/'docking').mkdir()
    with pytest.raises(BoundaryReached):bridge.stage_docking(tmp_path,{'vina':'unit-fixture'},tmp_path/'docking')
    assert visited==[selected['source_sdf']]
    selection=tmp_path/'availability/selection.json';write(selection,{'qualified_ids':[]})
    with pytest.raises(ValueError,match='Selection evidence'):original(selected['source_sdf'],selected['qualification_evidence'])


def test_api_exposes_real_persisted_failures(tmp_path,monkeypatch):
    from webui import workflow,server
    run=tmp_path/'test';out=run/'availability';out.mkdir(parents=True)
    record={'id':'m1','molecule_id':'m1','toxicity_predictions':[{'prediction_status':'ERROR','raw_prediction':None}], 'decision':'review'}
    write(out/'records.json',[record]);(out/'molecules.csv').write_text('id\nm1\n')
    monkeypatch.setattr(workflow,'RUNS',tmp_path)
    body=json.loads(server.api_predocking_detail({'id':['m1']},'test','availability')[2])
    assert body['record']==record


def test_only_selected_2d_structures_receive_conformers(tmp_path):
    from preflight.test_preflight import eligible
    from yuanhe_preflight.ligand import prepare_selected
    cfg=defaults(tmp_path);cfg.pop('engine')
    row=eligible(tmp_path,cfg)
    source=Path(row['sdf'])
    mol=Chem.MolFromSmiles('CCO')
    from rdkit.Chem import AllChem
    AllChem.Compute2DCoords(mol)
    two_d=tmp_path/'selected-2d.sdf'
    with Chem.SDWriter(str(two_d)) as w:w.write(mol)
    row={**row,'sdf':str(two_d)}
    prepared=prepare_selected([row],tmp_path,tmp_path/'prepared',17)
    assert prepared[0]['compound_id']==row['compound_id']
    assert prepared[0]['state_id']==row['state_id']
    assert Chem.SDMolSupplier(prepared[0]['sdf'])[0].GetConformer().Is3D()
    assert not Chem.SDMolSupplier(str(two_d))[0].GetConformer().Is3D()
    assert source.is_file()


@pytest.mark.parametrize('bad', [{}, {'AMES':0.1}, {'AMES':float('nan'),'hERG':0.1,'DILI':0.1,'LD50_Zhu':2}, {'AMES':1.1,'hERG':0.1,'DILI':0.1,'LD50_Zhu':2}])
def test_incomplete_or_invalid_model_values_cannot_enter_success_cache(bad):
    from yuanhe_preflight.toxicity_worker import validate_raw
    with pytest.raises(ValueError):validate_raw(bad)


def test_regression_ld50_is_not_a_probability():
    from yuanhe_preflight.toxicity_worker import validate_raw
    validate_raw({'AMES':0.1,'hERG':0.2,'DILI':0.3,'LD50_Zhu':3.2})


def test_partial_retrosynthesis_is_not_route_found(tmp_path,monkeypatch):
    from yuanhe_preflight.availability_synthesis_filter import AvailabilityFilter
    cfg=defaults(tmp_path)
    f=AvailabilityFilter(cfg,'test')
    monkeypatch.setattr(f,'search',lambda smiles:{'status':'searched','routes':[{'type':'mol','smiles':smiles,'in_stock':False,'route_metadata':{'is_solved':False}}]})
    result=f.screen('CCO')
    assert result['route_status']=='NO_ROUTE_FOUND'
    assert result['decision']=='pending' and result['routes']
