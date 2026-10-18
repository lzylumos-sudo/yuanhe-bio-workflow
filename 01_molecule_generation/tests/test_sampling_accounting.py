import json
import pytest
from src.adapters.peptide2mol import sampling_audit

def test_completed_zero_valid_is_a_negative_sampling_outcome(tmp_path):
 p=tmp_path/'audit.json';p.write_text(json.dumps(dict(requested_attempts=2,started_attempts=2,valid=0,failed=2,status='finished')))
 assert sampling_audit(p,2)['valid']==0

@pytest.mark.parametrize('change',[{'status':'running'},{'started_attempts':1},{'valid':1},{'failed':True},{'requested_attempts':3}])
def test_incomplete_or_inconsistent_attempts_never_pass(tmp_path,change):
 p=tmp_path/'audit.json';p.write_text(json.dumps(dict(dict(requested_attempts=2,started_attempts=2,valid=0,failed=2,status='finished'),**change)))
 with pytest.raises(ValueError):sampling_audit(p,2)
