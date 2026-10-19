"""Decision invariants: invalid poses and missing evidence must not win."""
import math
from types import SimpleNamespace
import pytest
from zvss.scoring.aggregate import representative_pose
from zvss.select.pipeline import _pareto, _exclusion


def entry(cid, affinity=-8., mechanism=2., strain=5., tier='A', score=.5):
    return dict(compound_id=cid, pose_id='p'+cid, state_id='s'+cid, tier=tier,
                tier_reasons=[], evidence_score=score,
                raw=dict(affinity=affinity, mechanism=mechanism, strain=strain),
                components=dict(affinity=.5, mechanism=.5), components_missing=[],
                flags=dict(excessive_strain=False))


def test_valid_pose_not_hidden_by_failed_high_score():
    poses=[dict(pose_id='bad',assess_status='fail',scores={'minimizedAffinity':-20}),
           dict(pose_id='good',assess_status='pass',scores={'minimizedAffinity':-7})]
    assert representative_pose(poses)['pose_id']=='good'


@pytest.mark.parametrize('missing',[None,float('nan'),float('inf')])
def test_unknown_strain_is_not_pareto_advantage(missing):
    rows=[entry('known'),entry('unknown',strain=missing)]
    assert [r['compound_id'] for r in _pareto(rows)]==['known']


def test_pareto_uses_raw_objectives_not_pool_percentiles():
    rows=[entry('good',affinity=-9.,mechanism=3.),entry('dominated',affinity=-7.,mechanism=1.)]
    assert [r['compound_id'] for r in _pareto(rows)]==['good']


def test_invalid_tier_cannot_appear_on_pareto_front():
    assert not _pareto([entry('failed',tier='D',strain=0)])


def test_hard_strain_failure_cannot_be_rescued_by_tier_a():
    row=entry('failed');row['flags']['excessive_strain']=True
    assert _exclusion(row,{},None) is not None
