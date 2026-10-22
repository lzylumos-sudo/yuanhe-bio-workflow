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

def test_evaluation_exports_full_denominator_and_separate_product(tmp_path):
    # Synthetic engineering fixture, not scientific evidence or a docking run.
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from scientific import artifacts as a
    from scientific.mdm2.evaluation import evaluate
    from scientific.mdm2.common import now
    root=tmp_path/'fixture';inputs=root/'inputs';inputs.mkdir(parents=True)
    target=root/'target';target.mkdir()
    (target/'receptor_original.pdb').write_text('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n')
    a.atomic(target/'peptide_features.json',{'residues':['A:1:ALA']});a.atomic(target/'perturbed_features.json',{'residues':['A:2:ALA']})
    labels=[];library=[];prepared=[];filters=[]
    properties={'MW':550,'cLogP':4,'TPSA':90,'HBD':3,'HBA':4,'RotB':5,'formal_charge':0}
    for i in range(21):
        pid=f'M{i:04d}';library.append({'parent_id':pid,'smiles':'CCO','properties':properties})
        labels.append({'parent_id':pid,'label':'active' if i<20 else 'unknown_background',
                       'name':['Nutlin-3a','SAR405838','Brigimadlin'][i] if i<3 else pid,
                       'source_id':'fixture','scaffold':'fixture','properties':properties})
        prepared.append({'parent_id':pid,'status':'ok','state_id':pid+'_state'})
        filters.append({'parent_id':pid,'ppi_pass':True,'deliverable_pass':False,
                        'risk':{'decision':'pass','alert_ids':[]},'acquisition':{'decision':'pending','status':'model_missing'}})
    for name,data in [('library',library),('prepared',prepared),('filters',filters)]:a.atomic(inputs/(name+'.json'),data)
    a.atomic(root/'evaluation/labels.json',labels)
    mol=Chem.AddHs(Chem.MolFromSmiles('CCO'));AllChem.EmbedMolecule(mol,randomSeed=1)
    job=inputs/'job';job.mkdir()
    with Chem.SDWriter(str(job/'poses.sdf')) as w:w.write(mol)
    a.atomic(inputs/'job_records.json',[{'parent_id':'M0000','status':'ok','job_dir':str(job),
        'output_sha256':a.digest(job/'poses.sdf'),'poses':[{'vina_affinity':-7.}]}])
    (inputs/'blinded_library.sdf').write_text('fixture only\n')
    frozen={'files':{},'source_files':{},'label_seal_sha256':a.digest(root/'evaluation/labels.json'),'git_commit':'fixture','frozen_at':now()}
    a.atomic(inputs/'preregistered_protocol.json',frozen)
    (inputs/'preregistered_protocol.sha256').write_text(a.digest(inputs/'preregistered_protocol.json'))
    a.atomic(root/'environment.json',{'effective_cpu_quota':22})
    a.atomic(root/'run_manifest.json',{'started_at':now()})
    a.atomic(root/'throughput.json',{'elapsed_seconds':1,'mean_job_seconds':1,'failure_rate_including_preparation':0})
    result=evaluate(root)
    assert result['total_parents']==21 and result['docked']==1
    assert result['profiles']['ppi_benchmark']['A1']['top_k']['50']['Recall']==.05
    assert result['profiles']['current_product']['A1']['eligible_ranked_count']==0
    assert (root/'deliverable_top50.sdf').stat().st_size==0
    assert (root/'scientific_top50.sdf').stat().st_size>0
    assert len(a.read(root/'anonymous_rankings.json'))==126
    assert '回顾性' in (root/'validation_report.md').read_text()


def test_state_enumeration_cannot_invent_unassigned_carbon_stereo():
    from rdkit import Chem
    from scientific.mdm2.library import has_defined_stable_stereo
    assert not has_defined_stable_stereo(Chem.MolFromSmiles('CC(O)C(=O)O'))
    assert has_defined_stable_stereo(Chem.MolFromSmiles('C[C@H](O)C(=O)O'))
    assert has_defined_stable_stereo(Chem.MolFromSmiles('C[NH+](CC)CCC'))
