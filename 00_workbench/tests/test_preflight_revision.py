"""Evidence editor regressions. All chemistry records are synthetic test fixtures."""
import copy
from contextlib import ExitStack
import json
from pathlib import Path
import sys
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from webui import workflow as w, preflight_revision as revision, server
from yuanhe_preflight.config import defaults
from yuanhe_preflight.structure import screen
from yuanhe_preflight.toxicity_substructure_filter import RiskRules
from yuanhe_preflight.common import canonical, hash_obj


@pytest.fixture
def case(tmp_path, monkeypatch):
    monkeypatch.setattr(w, 'RUNS', tmp_path)
    monkeypatch.setattr(w, '_active', {})
    monkeypatch.setattr(w, 'acquire_compute', ExitStack)
    d = tmp_path/'case'
    cfg = defaults(tmp_path)
    w.atomic(d/'config.json', {'predocking': cfg, 'research': {'predocking': cfg}})
    w.atomic(d/'run.json', {'id': 'case', 'config_sha256': w.digest(d/'config.json')})
    source = {'01_molecule_generation/src/model.py': 'original', '00_workbench/webui/server.py': 'fixed'}
    monkeypatch.setattr(w, 'source_snapshot', lambda: source.copy())
    w.atomic(d/'source_snapshot.json', {**source, '00_workbench/webui/server.py': 'old'})
    smi = canonical('CCCCCC')
    risk = screen(smi, RiskRules(), [], cfg['engine'])
    row = {'id': 'TEST_ONLY', 'compound_id': 'TEST_ONLY', 'smiles': smi, 'decision': 'review', 'risk': risk}
    for key in ('prepare', 'peptide2mol', 'reinvent4', 'risk', 'availability'):
        w.atomic(d/key/'records.json', [row] if key in ('risk', 'availability') else [])
        complete(d, key)
    return d, cfg, source


def complete(d, key):
    w.atomic(d/key/'stage.json', {'status': 'completed', 'config_sha256': w.digest(d/'config.json'),
                                 'artifacts': {f'{key}/records.json': w.digest(d/key/'records.json')}})


def payload(stage='risk'):
    return {'revision_token': revision.context('case', stage)['revision_token'],
            'reviewer': 'TEST_ONLY', 'reason': 'Synthetic regression fixture',
            'risk_reviews': [], 'route_reviews': []}


def risk_review():
    return {'compound_id': 'TEST_ONLY', 'rule_id': 'BRENK:Aliphatic_long_chain',
            'reviewer': 'TEST_ONLY', 'reason': 'Synthetic review fixture, not a scientific approval',
            'reviewed_at': '2026-10-28'}


def test_editor_get_and_legacy_empty_post_preserve_results(case):
    d, _, _ = case
    before = {str(p): p.read_bytes() for p in d.rglob('*') if p.is_file()}
    status, _, body = server.dispatch('GET', '/api/project/run/case/predocking/risk/revision', {}, {})
    assert status == 200 and json.loads(body)['candidates'][0]['compound_id'] == 'TEST_ONLY'
    with pytest.raises(ValueError):
        w.revise_preflight('case', 'risk', {})
    assert before == {str(p): p.read_bytes() for p in d.rglob('*') if p.is_file()}


def test_risk_save_preserves_policy_upstream_and_history(case):
    d, cfg, _ = case
    original = w.digest(d/'config.json')
    upstream = w.digest(d/'reinvent4/stage.json')
    req = payload(); req['risk_reviews'] = [risk_review()]
    result = w.revise_preflight('case', 'risk', req)
    assert result['screening_policy_unchanged']
    assert w.digest(d/'config.json') == original and w.digest(d/'reinvent4/stage.json') == upstream
    assert not (d/'risk').exists() and not (d/'availability').exists()
    folder = (d/result['revision']).parent
    assert (folder/'risk/records.json').is_file() and (folder/'availability/records.json').is_file()
    after = w.read(d/'reviews/risk.json')['predocking']
    assert {k: v for k, v in after.items() if k != 'risk_exceptions'} == {k: v for k, v in cfg.items() if k != 'risk_exceptions'}
    smi = w.read(folder/'risk/records.json')[0]['smiles']
    assert screen(smi, RiskRules(), after['risk_exceptions'], cfg['engine'])['decision'] == 'pass'
    assert w.read(folder/'risk/records.json')[0]['decision'] == 'review'
    w.verify_source(d)


@pytest.mark.parametrize('change', [
    lambda x: x.update(predocking={}), lambda x: x.update(reviewer=''),
    lambda x: x.update(revision_token='stale'), lambda x: x.update(risk_reviews={}),
    lambda x: x.update(risk_reviews=[{**risk_review(), 'compound_id': 'unknown'}]),
    lambda x: x.update(risk_reviews=[{**risk_review(), 'reason': ''}]),
    lambda x: x.update(route_reviews=[{}]),
])
def test_invalid_or_policy_edit_cannot_archive(case, change):
    d, _, _ = case
    req = payload(); change(req)
    before = w.digest(d/'risk/stage.json')
    with pytest.raises(ValueError): w.revise_preflight('case', 'risk', req)
    assert w.digest(d/'risk/stage.json') == before and not (d/'evidence_revisions').exists()


@pytest.mark.parametrize('target', ['docking', 'screening', 'gromacs', 'md_systems'])
def test_downstream_blocks_revision(case, target):
    d, _, _ = case
    req = payload(); (d/target).mkdir()
    with pytest.raises(ValueError, match='Downstream'): w.revise_preflight('case', 'risk', req)
    assert (d/'risk/records.json').exists()


def test_stale_editor_and_changed_upstream_rejected(case):
    d, _, source = case
    req = payload()
    source['01_molecule_generation/src/model.py'] = 'changed'
    with pytest.raises(ValueError, match='changed while'): w.revise_preflight('case', 'risk', req)
    with pytest.raises(ValueError, match='Upstream'): w.revise_preflight('case', 'risk', payload())


def test_missing_upstream_artifact_rejected(case):
    d, _, _ = case
    req = payload(); (d/'peptide2mol/records.json').unlink()
    with pytest.raises(RuntimeError, match='missing or modified'): w.revise_preflight('case', 'risk', req)
    assert (d/'risk/stage.json').is_file()


def test_mid_commit_failure_rolls_back(case, monkeypatch):
    d, _, _ = case
    req = payload(); original = w.atomic
    def fail_once(path, value):
        if Path(path).as_posix().endswith('reviews/risk.json'):
            raise OSError('SIMULATED disk failure')
        original(path, value)
    monkeypatch.setattr(w, 'atomic', fail_once)
    before = w.digest(d/'risk/records.json')
    with pytest.raises(OSError, match='disk failure'): w.revise_preflight('case', 'risk', req)
    assert w.digest(d/'risk/records.json') == before and (d/'availability/stage.json').exists()
    assert not (d/'active_source_revision.json').exists()
    assert not (d/'preflight_revision.pending.json').exists()


def test_crash_marker_blocks_execution(case):
    d, _, _ = case
    w.atomic(d/'preflight_revision.pending.json', {'record': 'incomplete'})
    with pytest.raises(ValueError, match='Incomplete evidence'): w.verify_source(d)
    with pytest.raises(ValueError, match='Incomplete evidence'): revision.context('case', 'risk')


def route_case(d, *, missing=False):
    # This tree is ONLY an API contract fixture, not a synthesis claim.
    tree = {'type': 'mol', 'smiles': 'CCCCCC', 'route_metadata': {'is_solved': not missing}}
    raw = {'target_smiles': 'CCCCCC', 'routes': [tree]}
    path = d/'availability/source_evidence/test.json'
    w.atomic(path, {'raw': raw, 'raw_sha256': hash_obj(raw), 'assets': {'fixture': True}})
    row = w.read(d/'availability/records.json')[0]
    row['availability'] = {'raw_search': {'search_artifact': str(path)}, 'asset_manifest': {'fixture': True},
        'routes': [{'tree': tree, 'route_sha256': hash_obj(tree), 'step_count': 1,
                    'problems': ['planner_did_not_solve_route'] if missing else ['template_selectivity_requires_review']}]}
    w.atomic(d/'availability/records.json', [row]); complete(d, 'availability')
    review = {k: v for k, v in risk_review().items() if k != 'rule_id'}
    review.update(route_sha256=hash_obj(tree), identity_stereo_selectivity_confirmed=True)
    req = payload('availability'); req['route_reviews'] = [review]
    return req


def test_route_source_preserved_after_archive(case):
    d, cfg, _ = case
    req = route_case(d)
    w.revise_preflight('case', 'availability', req)
    after = w.read(d/'reviews/availability.json')['predocking']
    assert Path(after['route_reviews'][0]['search_artifact']).is_file()
    assert {k: v for k, v in after.items() if k != 'route_reviews'} == {k: v for k, v in cfg.items() if k != 'route_reviews'}
    assert (d/'risk/stage.json').exists()


def test_incomplete_route_cannot_be_approved(case):
    d, _, _ = case
    with pytest.raises(ValueError, match='cannot replace'): w.revise_preflight('case', 'availability', route_case(d, missing=True))
    assert (d/'availability/stage.json').exists()


def test_risk_rerun_does_not_restore_stale_availability_settings(case):
    d, cfg, _ = case
    w.atomic(d/'reviews/availability.json', {'predocking': cfg})
    req = payload(); req['risk_reviews'] = [risk_review()]
    w.revise_preflight('case', 'risk', req)
    assert 'predocking' not in w.read(d/'reviews/availability.json')
    from yuanhe.configuration import stage_config
    resolved = stage_config(w.read(d/'config.json'), [('risk', w.read(d/'reviews/risk.json')), ('availability', w.read(d/'reviews/availability.json'))])
    assert len(resolved['config']['predocking']['risk_exceptions']) == 1


def test_hard_exclusion_cannot_be_waived(case):
    d, _, _ = case
    rows = w.read(d/'risk/records.json')
    for alert in rows[0]['risk']['alerts']: alert['action'] = 'exclude'
    w.atomic(d/'risk/records.json', rows); complete(d, 'risk')
    req = payload(); req['risk_reviews'] = [risk_review()]
    with pytest.raises(ValueError, match='hard exclusions'): w.revise_preflight('case', 'risk', req)


@pytest.mark.parametrize('field,value', [('structure_smiles', 'CCO'), ('rules_sha256', 'wrong'),
    ('catalog_policy_sha256', 'wrong'), ('reviewer', ''), ('rule_id', 'BRENK:unknown')])
def test_catalog_evidence_is_bound_to_identity_and_policy(case, field, value):
    d, cfg, _ = case
    row = w.read(d/'risk/records.json')[0]
    record = {**risk_review(), 'structure_smiles': row['smiles'], 'rules_sha256': row['risk']['rules_sha256'],
              'catalog_policy_sha256': row['risk']['catalog_policy_sha256']}
    record[field] = value
    with pytest.raises(ValueError): screen(row['smiles'], RiskRules(), [record], cfg['engine'])
