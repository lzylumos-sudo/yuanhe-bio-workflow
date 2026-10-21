"""All vendor/route examples in this file are SIMULATED UNIT FIXTURES, not commercial evidence."""
import copy,csv,json,sys,subprocess
from datetime import date,timedelta
from pathlib import Path
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'02_docking/src'),str(ROOT/'00_workbench'),str(ROOT/'03_virtual_screening/src')]
from yuanhe_preflight.common import *
from yuanhe_preflight.config import defaults
from yuanhe_preflight.toxicity_substructure_filter import RiskRules
from yuanhe_preflight.catalog import Catalog
from yuanhe_preflight.availability_synthesis_filter import AvailabilityFilter,validate_route,assets,BUY_LABELS,ROUTE_LABEL
from yuanhe_preflight.pipeline import risk_stage,availability_stage
from yuanhe_preflight.gate import require_ligand,validate_evidence,QualificationError

@pytest.fixture
def cfg(tmp_path):
    return defaults(tmp_path)

def catalog(tmp_path,smiles=('CCO',),age=0,kind='in_stock',name='catalog'):
    d=tmp_path/name;d.mkdir(exist_ok=True);doc=d/'supply.txt';doc.write_text('SIMULATED UNIT TEST SUPPLY RECORD. Not a real vendor quote.')
    day=str(date.today()-timedelta(days=age))
    rows=[dict(smiles=s,vendor='SIMULATED_VENDOR',sku='TEST_'+str(n),supply_type=kind,checked_at=day,
               catalog_date=day,evidence_kind='stock_confirmation',evidence_file=doc.name,evidence_sha256=file_hash(doc)) for n,s in enumerate(smiles)]
    path=d/'catalog.csv'
    with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    return {'path':str(path),'version':'simulated-v1','source':'unit-test-only','license':'unit-test-fixture','trusted':True}

def sdf(tmp_path,smi,name='input'):
    mol=Chem.AddHs(Chem.MolFromSmiles(smi));p=AllChem.ETKDGv3();p.randomSeed=17;assert AllChem.EmbedMolecule(mol,p)==0
    file=tmp_path/(name+'.sdf')
    with Chem.SDWriter(str(file)) as w:w.write(mol)
    return file

def eligible(tmp_path,cfg,smi='CCO'):
    cfg['catalogs']=[catalog(tmp_path,[smi])]
    source=sdf(tmp_path,smi);row={'id':'original_1','smiles':smi,'sdf':str(source),'core':'C','seed_id':'seed1'}
    risk_stage([row],cfg,tmp_path/'risk',base=tmp_path,batch_id='SIMULATED_BATCH')
    availability_stage(read(tmp_path/'risk/records.json'),cfg,tmp_path/'availability')
    return read(tmp_path/'availability/qualified.json')[0]

@pytest.mark.parametrize('smiles,decision,category',[
 ('CN(C)N=O','exclude','toxicity_risk'),('CC(=O)Cl','exclude','nonspecific_reactivity'),
 ('CN=C=O','exclude','nonspecific_reactivity'),('CC1CO1','review','nonspecific_reactivity'),
 ('Oc1ccccc1O','pass','assay_interference'),('O=[N+]([O-])c1ccccc1','pass','medicinal_chemistry')])
def test_rule_examples_and_classification(smiles,decision,category):
    r=RiskRules().screen(smiles);assert r['decision']==decision;assert category in r['alert_types'];assert r['matched_atom_indices']
    if category=='assay_interference':assert 'toxicity_risk' not in r['alert_types']

def test_no_alert_wording():
    assert RiskRules().screen('CCO')['reason']=='未命中当前规则'

def test_rule_invalid_and_example_validation(tmp_path):
    rules=RiskRules().data;rules['rules'][0]['positive_examples']=['CCO'];p=tmp_path/'bad.json';write(p,rules)
    with pytest.raises(ValueError,match='rule_example'):RiskRules(p)
    rules['rules'][0]['SMARTS']='[';write(p,rules)
    with pytest.raises(ValueError,match='invalid_SMARTS'):RiskRules(p)

def test_whole_molecule_excluded_original_retained(tmp_path,cfg):
    s='CN(C)N=O';f=sdf(tmp_path,s);before=file_hash(f)
    risk_stage([{'id':'bad','smiles':s,'sdf':str(f)}],cfg,tmp_path/'risk',base=tmp_path,batch_id='unit')
    row=read(tmp_path/'risk/records.json')[0]
    assert row['decision']=='exclude' and row['smiles']==canonical(s) and row['original_smiles']==s
    assert file_hash(f)==before and read(tmp_path/'risk/raw_input.json')[0]['smiles']==s
    assert len(list(Chem.SDMolSupplier(str(tmp_path/'risk/excluded.sdf'))))==1
    availability_stage([row],cfg,tmp_path/'availability');assert read(tmp_path/'availability/qualified.json')==[]

def test_manual_exception_bound_to_structure_rules():
    rules=RiskRules();s=canonical('CN(C)N=O')
    ex={'rule_id':'TOX_N_NITROSAMINE','structure_smiles':s,'rules_sha256':rules.sha256,'reviewer':'Unit tester','reason':'SIMULATED exception','reviewed_at':'2026-10-21'}
    r=rules.screen(s,[ex]);assert r['decision']=='pass' and r['alerts'] and r['exceptions']==[ex]
    with pytest.raises(ValueError):rules.screen(s,[dict(ex,reason='')])
    with pytest.raises(ValueError):rules.screen(s,[dict(ex,structure_smiles='CCO')])

def test_exact_identity_stereo_isotope_salts(tmp_path):
    c=Catalog([catalog(tmp_path,['C[C@H](O)C(=O)O','[13CH3]CO','C[NH3+].[Cl-]'])])
    assert c.valid('C[C@H](O)C(=O)O')
    for s in ['C[C@@H](O)C(=O)O','CC(O)C(=O)O','CCO','CN','C[NH3+]','CCCO']:assert not c.valid(s)

def test_expired_and_catalog_only(tmp_path,cfg):
    c=Catalog([catalog(tmp_path,age=40)],30);assert c.match('CCO')[0]['status']=='supply_evidence_expired';assert not c.valid('CCO')
    c=Catalog([catalog(tmp_path,kind='catalog_only',name='unknown')]);assert not c.valid('CCO')
    assert c.match('CCO')[0]['status']=='catalog_hit_supply_unconfirmed'

def test_supply_document_tampering(tmp_path,cfg):
    row=eligible(tmp_path,cfg);assert require_ligand(row['source_sdf'],row['qualification_evidence'])
    (tmp_path/'catalog/supply.txt').write_text('modified')
    with pytest.raises(QualificationError,match='no_current_supply'):require_ligand(row['source_sdf'],row['qualification_evidence'])

def test_no_sa_fallback_and_missing_model(tmp_path,cfg):
    f=AvailabilityFilter(cfg,RiskRules().sha256);r=f.screen('CCO')
    assert r['status']=='model_missing' and r['decision']=='pending' and r['labels']==[]
    r['SA_score']=1
    assert r['decision']!='pass'

# A simulated one-step hydrolysis reverse tree, replayed by RDKit. Not a synthesis claim.
def tree():
    return {'type':'mol','smiles':'CCOC(C)=O','children':[{'type':'reaction','smiles':'CCOC(C)=O>>CC(=O)O.CCO',
      'metadata':{'template':'[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])O.[O:3][C:4]','template_source':'SIMULATED_UNIT_TEST'},
      'children':[{'type':'mol','smiles':'CC(=O)O','in_stock':True},{'type':'mol','smiles':'CCO','in_stock':True}]}]}

def test_complete_route_requires_all_terminal_supply(tmp_path):
    c=Catalog([catalog(tmp_path,['CC(=O)O','CCO'])]);r=validate_route(tree(),'CCOC(C)=O',c)
    assert r['status']=='complete' and r['step_count']==1 and len(r['terminal_materials'])==2
    assert validate_route(tree(),'CCOC(C)=O',Catalog())['status']=='needs_review'
    t=tree();t['children'][0]['children'][0]['smiles']='CCC';assert validate_route(t,'CCOC(C)=O',c)['status']=='needs_review'
    assert validate_route(tree(),'CCC',c)['status']=='identity_mismatch'

def test_stereo_route_unresolved(tmp_path):
    t={'type':'mol','smiles':'C[C@H](O)CCOC(C)=O','children':[{'type':'reaction','metadata':{'retro_template':'[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])O.[O:3][C:4]'},'children':[{'type':'mol','smiles':'CC(=O)O'},{'type':'mol','smiles':'C[C@H](O)CCO'}]}]}
    c=Catalog([catalog(tmp_path,['CC(=O)O','C[C@H](O)CCO'])]);r=validate_route(t,t['smiles'],c)
    assert 'stereochemistry_selectivity_requires_review' in r['problems'] and r['status']=='needs_review'
    review={'route_sha256':hash_obj(t),'target_smiles':canonical(t['smiles']),'reviewer':'SIMULATED','reason':'Unit fixture','reviewed_at':str(date.today()),'identity_stereo_selectivity_confirmed':True}
    assert validate_route(t,t['smiles'],c,route_review=review)['status']=='complete'

def fake_assets(tmp_path,cfg):
    for key in ('model','templates','stock'):
        p=tmp_path/(key+'.fixture');p.write_text('SIMULATED NOT REAL ASSET');cfg['retro'][key+'_path']=str(p)
        cfg['retro'][key+'_version']='SIMULATED';cfg['retro'][key+'_license']='Unit test fixture only'

def test_computed_route_label_and_gate(tmp_path,cfg,monkeypatch):
    fake_assets(tmp_path,cfg);cfg['stock_catalogs']=[catalog(tmp_path,['CC(=O)O','CCO'])]
    raw={'status':'searched','target_smiles':'CCOC(C)=O','routes':[tree()],'backend':'SIMULATED_UNIT_TEST'}
    manifest,_=assets(cfg['retro']);artifact=tmp_path/'raw.fixture.json';write(artifact,{'raw':raw,'raw_sha256':hash_obj(raw),'assets':manifest})
    monkeypatch.setattr(AvailabilityFilter,'search',lambda self,s:{**raw,'search_artifact':str(artifact)})
    source=sdf(tmp_path,'CCOC(C)=O');risk_stage([{'id':'route','smiles':'CCOC(C)=O','sdf':str(source)}],cfg,tmp_path/'risk',base=tmp_path,batch_id='UNIT')
    availability_stage(read(tmp_path/'risk/records.json'),cfg,tmp_path/'availability');row=read(tmp_path/'availability/qualified.json')[0]
    assert row['labels']==[ROUTE_LABEL];assert require_ligand(row['source_sdf'],row['qualification_evidence'])
    # Both acquisition types coexist without a new search.
    cfg['catalogs']=[catalog(tmp_path,['CCOC(C)=O'],name='product')]
    cfg['route_reviews']=[{'search_artifact':str(artifact)}]
    monkeypatch.setattr(AvailabilityFilter,'search',lambda *a:pytest.fail('Purchased targets must not trigger search'))
    both=AvailabilityFilter(cfg,RiskRules().sha256).screen('CCOC(C)=O')
    assert BUY_LABELS['in_stock'] in both['labels'] and ROUTE_LABEL in both['labels']

@pytest.mark.parametrize('status',['search_timeout','dependency_missing','no_route_found','stock_missing'])
def test_failures_never_pass(tmp_path,cfg,monkeypatch,status):
    monkeypatch.setattr(AvailabilityFilter,'search',lambda *a:{'status':status,'routes':[]})
    r=AvailabilityFilter(cfg,RiskRules().sha256).screen('CCO');assert r['status']==status and not r['labels'] and r['decision']=='pending'

def test_real_worker_missing_dependency(tmp_path,cfg):
    import importlib.util
    if importlib.util.find_spec('aizynthfinder'):pytest.skip('This environment has AiZynthFinder')
    fake_assets(tmp_path,cfg);cfg['retro']['python']=sys.executable
    r=AvailabilityFilter(cfg,RiskRules().sha256).screen('CCO');assert r['status']=='dependency_missing';assert not r['labels']

def test_timeout_process(tmp_path,cfg,monkeypatch):
    fake_assets(tmp_path,cfg)
    def timeout(*a,**k):raise subprocess.TimeoutExpired('SIMULATED',1)
    monkeypatch.setattr(subprocess,'run',timeout)
    r=AvailabilityFilter(cfg,RiskRules().sha256).screen('CCO');assert r['status']=='search_timeout' and not r['labels']

def test_gate_rejects_changed_structure_hash_and_missing_evidence(tmp_path,cfg):
    row=eligible(tmp_path,cfg);proof=read(row['qualification_evidence'])
    with pytest.raises(QualificationError,match='identity'):validate_evidence(proof,'CCCO')
    bad=copy.deepcopy(proof);bad['risk']['decision']='exclude'
    with pytest.raises(ValueError,match='hash'):validate_evidence(bad,'CCO')
    other=sdf(tmp_path,'CCCO','other')
    with pytest.raises(QualificationError):require_ligand(other,row['qualification_evidence'])
    with pytest.raises(QualificationError):require_ligand(other)

def test_public_vina_entry_refuses_before_process_launch(tmp_path,monkeypatch):
    from sbvs2.docking.vina_official import run_vina_job
    monkeypatch.setattr(subprocess,'run',lambda *a,**k:pytest.fail('Must reject before invoking tool'))
    with pytest.raises(QualificationError):run_vina_job('vina','receptor',sdf(tmp_path,'CCO'),None,None,17,tmp_path/'job')
    assert not (tmp_path/'job').exists()

def test_policy_and_source_changed_invalidates_gate(tmp_path,cfg):
    row=eligible(tmp_path,cfg);proof=read(row['qualification_evidence'])
    proof['config']['risk_policies']={'BRENK_NITRO':'review'};proof=seal(proof)
    with pytest.raises(QualificationError,match='risk'):validate_evidence(proof,'CCO')
    Path(row['source_sdf']).write_text('changed')
    with pytest.raises(ValueError):require_ligand(row['sdf'],row['qualification_evidence'])

def test_identity_recipe_matches_existing_zvss():
    from zvss.ids import compound_id,chemical_state_id
    s=canonical('C[C@H](O)C(=O)O');i=identity(s)
    assert i['compound_id']==compound_id(s);assert i['state_id']==chemical_state_id(s,0,'')

def test_stage_order_and_settings_roundtrip(tmp_path,cfg,monkeypatch):
    from webui import workflow as w
    keys=[x[0] for x in w.STAGES];assert keys[keys.index('reinvent4')+1:keys.index('docking')]==['risk','availability']
    monkeypatch.setattr(w,'RUNS',tmp_path);d=tmp_path/'case';write(d/'run.json',{'id':'case'})
    w.review('case','risk',{'predocking':cfg});assert read(d/'reviews/risk.json')['predocking']==cfg
    write(d/'risk/stage.json',{'status':'completed'});write(d/'reinvent4/candidates.json',[{'id':'unchanged'}]);before=file_hash(d/'reinvent4/candidates.json')
    w.revise_preflight('case','risk');assert not (d/'risk').exists();assert file_hash(d/'reinvent4/candidates.json')==before
    assert list((d/'attempts').glob('risk_evidence_revision_*/stage.json'))

def test_batch_dedup_and_provenance(tmp_path,cfg):
    cfg['catalogs']=[catalog(tmp_path)];f=sdf(tmp_path,'CCO')
    rows=[{'id':str(i),'smiles':s,'sdf':str(f),'core':str(i)} for i,s in enumerate(['CCO','OCC'])]
    c=risk_stage(rows,cfg,tmp_path/'risk',base=tmp_path,batch_id='UNIT');assert c['passed']==1 and c['duplicates']==1
    availability_stage(read(tmp_path/'risk/records.json'),cfg,tmp_path/'availability')
    row=read(tmp_path/'availability/qualified.json')[0];assert row['source_core']=='0' and row['generation_batch']=='UNIT'

def test_api_detail_has_highlight_and_raw_alerts(tmp_path,cfg,monkeypatch):
    from webui import workflow as w,server
    monkeypatch.setattr(w,'RUNS',tmp_path);d=tmp_path/'case';d.mkdir();s='CN(C)N=O';f=sdf(d,s)
    risk_stage([{'id':'test','smiles':s,'sdf':str(f)}],cfg,d/'risk',base=d,batch_id='UNIT')
    cid=read(d/'risk/records.json')[0]['id'];status,_,body=server.dispatch('GET','/api/project/run/case/predocking/risk',{'id':[cid]},None)
    payload=json.loads(body);assert status==200 and '<svg' in payload['svg'] and payload['record']['alert_ids']

def test_screening_import_keeps_evidence_and_identity(tmp_path,cfg):
    from types import SimpleNamespace
    from zvss.ingest.import_run import _import_single,ImportReport
    from zvss.ingest.poses import read_sdf_poses
    row=eligible(tmp_path,cfg)
    pose=tmp_path/'poses.sdf';pose.write_bytes(Path(row['source_sdf']).read_bytes());write(str(pose)+'.qualification.json',read(row['qualification_evidence']))
    raw=read_sdf_poses(SimpleNamespace(path=str(pose),relpath='poses.sdf',sha256=file_hash(pose)))[0][0]
    stored={}
    db=SimpleNamespace(upsert_compound=lambda r:stored.update(compound=r),upsert_state=lambda r:stored.update(state=r),upsert_job=lambda r:stored.update(job=r))
    report=ImportReport(run_id='UNIT');report.receptor={'version':'unit','sha256':'unit'}
    result=_import_single(SimpleNamespace(db=db,run_id='UNIT'),raw,SimpleNamespace(root=str(tmp_path)),{}, {},report)
    pose_row,_,compound,state,_=result
    assert compound['compound_id']==row['compound_id'] and state['state_id']==row['state_id']
    assert pose_row['meta']['predocking_evidence']['record_sha256']==row['qualification_sha256']
    assert compound['meta']['predocking']['availability']['labels']==row['labels']

def test_empty_qualified_pool_does_not_start_docking_or_screening(tmp_path,monkeypatch):
    import pipeline_bridge as bridge
    monkeypatch.setattr(bridge,'require',lambda *a:{})
    monkeypatch.setattr(bridge,'executable',lambda *a:pytest.fail('No tool should be launched'))
    write(tmp_path/'availability/qualified.json',[]);write(tmp_path/'risk/summary.json',{'input':10})
    out=tmp_path/'docking';out.mkdir()
    assert bridge.stage_docking(tmp_path,{},out)['input']==0
    screen=tmp_path/'screening';screen.mkdir()
    assert bridge.stage_screening(tmp_path,{},screen)['remaining']==0

def test_cross_task_cache_keys_and_resume(tmp_path,cfg,monkeypatch):
    fake_assets(tmp_path,cfg);calls=[]
    def simulated_process(argv,**kwargs):
        calls.append(argv)
        request=read(argv[argv.index('--request')+1]);target=argv[argv.index('--output')+1]
        write(target,{'status':'searched','target_smiles':request['smiles'],'routes':[], 'fixture':'SIMULATED CACHE UNIT TEST'})
    monkeypatch.setattr(subprocess,'run',simulated_process)
    first=AvailabilityFilter(cfg,'rules-v1').screen('CCO');assert len(calls)==1
    second=AvailabilityFilter(copy.deepcopy(cfg),'rules-v1').screen('OCC');assert len(calls)==1
    assert second['raw_search']['cache_hit'] and first['raw_search']['cache_key']==second['raw_search']['cache_key']
    AvailabilityFilter(cfg,'rules-v2').screen('CCO');assert len(calls)==2
    Path(cfg['retro']['model_path']).write_text('SIMULATED new model bytes')
    AvailabilityFilter(cfg,'rules-v2').screen('CCO');assert len(calls)==3
    cfg['catalogs']=[catalog(tmp_path,['CCC'])]
    AvailabilityFilter(cfg,'rules-v2').screen('CCO');assert len(calls)==4

def test_sdf_and_smiles_catalog_formats(tmp_path):
    spec=catalog(tmp_path);rows=list(csv.DictReader(Path(spec['path']).open()))
    path=tmp_path/'supplier.sdf';mol=Chem.MolFromSmiles('CCO')
    for k,v in rows[0].items():mol.SetProp(k,v)
    # SDF evidence reference is relative to this SDF directory.
    mol.SetProp('evidence_file','catalog/supply.txt')
    with Chem.SDWriter(str(path)) as w:w.write(mol)
    assert Catalog([{**spec,'path':str(path)}]).valid('CCO')
    tsv=tmp_path/'supplier.smi';rows[0]['evidence_file']='catalog/supply.txt'
    with tsv.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
    assert Catalog([{**spec,'path':str(tsv)}]).valid('CCO')

def test_export_includes_raw_vendor_supply_files(tmp_path,cfg):
    row=eligible(tmp_path,cfg)
    archive=row['source_evidence_archive']
    assert any(k.endswith('supply.txt') for k in archive)
    for value in archive.values():assert file_hash(tmp_path/'availability'/value['artifact'])==value['sha256']

def test_template_regioselectivity_requires_review(tmp_path):
    t={'type':'mol','smiles':'CCOC','children':[{'type':'reaction','metadata':{'template':'[C:1]-[O:2]>>[C:1].[O:2]'},'children':[{'type':'mol','smiles':'CC'},{'type':'mol','smiles':'CO'}]}]}
    stock=Catalog([catalog(tmp_path,['CC','CO','C','CCO'])]);r=validate_route(t,'CCOC',stock)
    assert 'template_selectivity_requires_review' in r['problems'] and r['status']=='needs_review'
