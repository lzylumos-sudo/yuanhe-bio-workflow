import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scientific import artifacts as a
from scientific.validation import replicas

def built(tmp_path):
 d=tmp_path/'run';folder=d/'md_systems/MD_original';folder.mkdir(parents=True)
 (folder/'system.gro').write_text('immutable input');(folder/'production.xtc').write_text('never copy trajectory')
 meta={'id':'MD_original','status':'completed','files':{'system.gro':a.digest(folder/'system.gro')},'parameters':{'production_ns':1},'compound_id':'CMP1'}
 a.atomic(folder/'system.json',meta)
 return d,folder

def test_replicas_preserve_only_inputs_and_distinct_seeds(tmp_path):
 d,src=built(tmp_path);before=a.digest(src/'system.json')
 result=replicas(d,{},dict(system_id='MD_original',reviewed=True,seeds=[17,41,97]),tmp_path)
 assert len(result['replicas'])==3
 for row in result['replicas']:
  dest=d/'md_systems'/row['system_id'];meta=a.read(dest/'system.json')
  assert meta['parameters']['replica_seed']==row['seed']
  assert meta['status']=='built'
  assert not (dest/'production.xtc').exists()
  assert a.digest(dest/'system.gro')==a.digest(src/'system.gro')
 assert a.digest(src/'system.json')==before

def test_replica_rejects_changed_parent(tmp_path):
 d,src=built(tmp_path);(src/'system.gro').write_text('tampered')
 with pytest.raises(ValueError,match='checksum'):replicas(d,{},dict(system_id='MD_original',reviewed=True),tmp_path)

@pytest.mark.parametrize('seeds',[[17,17],[True,41],[0,41],[17]])
def test_replica_rejects_invalid_seeds(tmp_path,seeds):
 d,src=built(tmp_path)
 with pytest.raises(ValueError):replicas(d,{},dict(system_id='MD_original',reviewed=True,seeds=seeds),tmp_path)

def test_effective_candidates_preserves_original_and_fails_on_changed_pose(tmp_path):
 from scientific.validation import effective_candidates
 d=tmp_path;src=d/'screening';src.mkdir()
 (src/'pose.sdf').write_text('fixed pose')
 rows=[dict(compound_id='CMP1',pose_id='POS1',sdf='screening/pose.sdf',qualification='preliminary')]
 a.atomic(src/'candidates.json',rows)
 base=src/'workspace/runs/screen/stages'
 a.atomic(base/'selection/summary.json',{});a.atomic(base/'redock/summary.json',{})
 op=d/'operations/a';a.atomic(op/'operation.json',{'status':'completed'})
 proof=dict(source_candidate_sha256=a.digest(src/'candidates.json'),source_selection_sha256=a.digest(base/'selection/summary.json'),source_redock_sha256=a.digest(base/'redock/summary.json'),qualification={'CMP1':dict(pose_id='POS1',reference_sha256=a.digest(src/'pose.sdf'),formal_qualified=True,requirements=[dict(status='pass')])})
 a.atomic(op/'result.json',proof);a.atomic(op/'screening_reassessment.json',proof)
 assert effective_candidates(d)[0]['qualification']=='formal'
 assert a.read(src/'candidates.json')[0]['qualification']=='preliminary'
 (src/'pose.sdf').write_text('modified')
 with pytest.raises(ValueError,match='pose mismatch'):effective_candidates(d)


def completed_replicas(tmp_path):
 from scientific.validation import MD_QUALITY_PROTOCOL
 ids=[]
 for seed in (17,41,97):
  ident=f'MD_{seed}';folder=tmp_path/'md_systems'/ident;folder.mkdir(parents=True)
  (folder/'production.xtc').write_bytes(b'fixture, not a real trajectory')
  a.atomic(folder/'system.json',dict(status='completed',compound_id='fixture',replica=dict(parent_system_sha256='same',velocity_seed=seed),acceptance_protocol=MD_QUALITY_PROTOCOL,parameters=dict(temperature_K=310)))
  a.atomic(folder/'analysis/summary.json',dict(source_hashes={'production.xtc':a.digest(folder/'production.xtc')},production_span_ns=5,quality_observations=dict(finite_trajectory=True,constraint_warnings=False),thermodynamics={'Temperature':{'mean':310}}))
  a.atomic(folder/'analysis/series.json',[dict(backbone_rmsd_A=1,ligand_rmsd_A=1,minimum_contact_A=3) for _ in range(10)])
  ids.append(ident)
 record_analyses(tmp_path,ids)
 return ids

def record_analyses(d,ids):
 a.atomic(d/'gromacs/systems.json',[{'system_id':ident,'files':{str(p.relative_to(d)):a.digest(p) for p in (d/'md_systems'/ident).rglob('*') if p.is_file()}} for ident in ids])


def test_quality_gate_does_not_assert_activity_or_convergence(tmp_path):
 from scientific.validation import summarize_replicas
 ids=completed_replicas(tmp_path)
 result=summarize_replicas(tmp_path,{},dict(system_ids=ids),tmp_path/'out')
 assert result['provisional_md_quality_accepted']
 assert result['md_validated_molecules']==0
 assert result['convergence']=='not_established'
 assert result['short_md_quality_passed_compounds']==['fixture']


def test_quality_gate_rejects_unstable_repeat_and_modified_raw_input(tmp_path):
 from scientific.validation import summarize_replicas
 ids=completed_replicas(tmp_path);folder=tmp_path/'md_systems'/ids[0]
 series=a.read(folder/'analysis/series.json');series[-1]['ligand_rmsd_A']=8
 a.atomic(folder/'analysis/series.json',series)
 with pytest.raises(ValueError,match='provenance'):summarize_replicas(tmp_path,{},dict(system_ids=ids),tmp_path/'untrusted')
 record_analyses(tmp_path,ids)
 result=summarize_replicas(tmp_path,{},dict(system_ids=ids),tmp_path/'out')
 assert not result['provisional_md_quality_accepted']
 assert any(not r['passed'] and 'ligand_pose_retention' in r['criterion'] for r in result['quality_criteria'])
 (folder/'production.xtc').write_bytes(b'modified')
 with pytest.raises(ValueError,match='provenance'):summarize_replicas(tmp_path,{},dict(system_ids=ids),tmp_path/'out2')


def test_ui_distinguishes_md_quality_from_experimental_validation(tmp_path,monkeypatch):
 from scientific.validation import summarize_replicas
 from webui import workflow
 ids=completed_replicas(tmp_path);op=tmp_path/'operations/comparison'
 result=summarize_replicas(tmp_path,{},dict(system_ids=ids),op)
 a.atomic(op/'result.json',result);a.atomic(op/'operation.json',{'status':'completed'})
 a.atomic(tmp_path/'config.json',{})
 a.atomic(tmp_path/'gromacs/stage.json',{'status':'completed','counts':{'scientifically_validated':0}})
 monkeypatch.setattr(workflow,'RUNS',tmp_path.parent)
 current=next(s for s in workflow.status(tmp_path.name)['stages'] if s['key']=='gromacs')
 assert current['counts']['md_quality_passed']==1
 assert current['counts']['scientifically_validated']==0
 assert 'md_quality_passed' not in a.read(tmp_path/'gromacs/stage.json')['counts']
 a.atomic(op/'result.json',{})
 with pytest.raises(ValueError,match='quality evidence'):workflow.status(tmp_path.name)
