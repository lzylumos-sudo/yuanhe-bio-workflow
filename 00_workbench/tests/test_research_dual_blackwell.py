"""Target-host tests supplied but NOT executed in the static-only upgrade."""
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from research.config import defaults,resolve,plan,fingerprint

def test_default_production_budget():
    c=resolve();assert (c['md_limit'],c['md_replicas'],c['md_production_ns'],c['md_concurrency'])==(8,1,500,2)
    p=plan(c);assert p['total_production_ns']==4000
    assert p['docking_task_upper_bounds']=={'initial':100000,'refine':16000}

@pytest.mark.parametrize('change',[{'md_replicas':3},{'md_limit':51},{'refine_limit':100001},{'md_production_ns':float('nan')},{'vina_workers':25},{'version':'old'}])
def test_invalid_configuration(change):
    with pytest.raises(ValueError):resolve(change)

def test_nested_parameters_and_smoke_are_independent():
    c=resolve({'initial':{'exhaustiveness':16}})
    assert c['initial']['exhaustiveness']==16 and defaults()['initial']['exhaustiveness']==8
    assert resolve({'profile':'smoke'})['initial_limit']==40

def test_global_core_quotas():
    c=defaults();n=c['core_limit'];quotas=[c['reinvent_attempts']//n+(i<c['reinvent_attempts']%n) for i in range(n)]
    assert sum(quotas)==200000 and max(quotas)==1000

def test_stable_parent_identity_preserves_enantiomers():
    pytest.importorskip('rdkit')
    from research.qualification import identity
    assert identity('CC(=O)O')[0]==identity('CC(=O)[O-].[Na+]')[0]
    assert identity('N[C@H](C)C(=O)O')[0]!=identity('N[C@@H](C)C(=O)O')[0]

def test_no_evidence_never_qualifies(tmp_path):
    pytest.importorskip('rdkit')
    from research.qualification import acquisition
    assert acquisition(tmp_path,'CCO',[])['status']=='acquisition_pending'

def test_risk_requires_review_and_pains_cannot_be_toxicity_block():
    Chem=pytest.importorskip('rdkit.Chem')
    from research.qualification import risk
    rule={'id':'test','version':'unit','smarts':'[N+](=O)[O-]','rationale':'unit fixture','source':'unit fixture','action':'review'}
    review={'reviewer':'unit','date':'2026-10-21','rules_sha256':fingerprint([rule])}
    assert risk(Chem.MolFromSmiles('C[N+](=O)[O-]'),[rule],review)['status']=='risk_flagged_review'
    assert risk(Chem.MolFromSmiles('CCO'),[rule],review)['status']=='risk_screen_passed'
    with pytest.raises(ValueError):risk(Chem.MolFromSmiles('CCO'),[rule],None)
