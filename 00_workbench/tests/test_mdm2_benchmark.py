import json,sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'00_workbench'),str(ROOT/'01_molecule_generation'),str(ROOT/'02_docking/src')]
from scientific.mdm2.library import canonical_parent,ppi_pass
from scientific.mdm2.scoring import rank
from scientific.mdm2.evaluation import metrics_for
from scientific.peptide_contacts import contact_signature,contact_overlap

def test_enantiomers_remain_distinct_after_salt_stripping():
    assert canonical_parent('C[C@H](N)C(=O)O.Cl')==canonical_parent('C[C@H](N)C(=O)O')
    assert canonical_parent('C[C@H](N)C(=O)O')!=canonical_parent('C[C@@H](N)C(=O)O')

def test_full_library_failures_remain_misses():
    labels=[{'parent_id':str(i),'label':'active' if i<20 else 'unknown_background'} for i in range(550)]
    result=metrics_for([{'parent_id':'0'}],labels,550)
    assert result['top_k']['50']['Recall']==.05
    assert result['top_k']['50']['EF']==pytest.approx(.55)
    assert result['top_k']['50']['retrieved']==1

def test_ties_and_controls_do_not_use_labels():
    rows=[{'parent_id':'B','vina_affinity':-9.,'A1_recall':.1,'A2_recall':.7},
          {'parent_id':'A','vina_affinity':-8.,'A1_recall':.8,'A2_recall':.2}]
    assert [x['parent_id'] for x in rank(rows,'A0')]==['B','A']
    assert [x['parent_id'] for x in rank(rows,'A1')]==['A','B']
    assert [x['parent_id'] for x in rank(rows,'A2')]==['B','A']
    for r in rows:r.update(label='active',name='irrelevant')
    assert [x['parent_id'] for x in rank(rows,'A1')]==['A','B']

def test_contact_geometry_and_residue_identity():
    atoms=[(('A','19','PHE','CA'),[0,0,0]),(('A','23','TRP','CB'),[10,0,0])]
    sig=contact_signature(atoms,[[4,0,0]])
    assert sig=={'A:19:PHE'}
    assert contact_overlap(sig,sig)['recall']==1
    assert contact_signature(atoms,[[100,0,0]])==set()

def test_ppi_boundary_is_uniform():
    p={'MW':581.5,'cLogP':5,'TPSA':90,'RotB':5,'HBD':1,'HBA':6,'formal_charge':0}
    assert ppi_pass(p)
    assert not ppi_pass({**p,'MW':801})

def test_production_entry_still_requires_qualification(tmp_path,monkeypatch):
    from sbvs2.docking import vina_official as v
    called=[]
    monkeypatch.setattr(v,'_execute_vina_job',lambda *a,**k:called.append(True))
    with pytest.raises(Exception):v.run_vina_job('vina','receptor',tmp_path/'unknown.sdf',None,None,1,tmp_path/'job')
    assert called==[]

def test_benchmark_rejects_manifest_tamper_before_execution(tmp_path,monkeypatch):
    from sbvs2.docking import vina_official as v
    p=tmp_path/'manifest.json';p.write_text('{}');called=[]
    monkeypatch.setattr(v,'_execute_vina_job',lambda *a,**k:called.append(True))
    with pytest.raises(ValueError,match='manifest changed'):
        v.run_benchmark_job('vina','receptor','ligand',None,None,1,tmp_path/'job',manifest=p,manifest_sha256='incorrect')
    assert called==[]
