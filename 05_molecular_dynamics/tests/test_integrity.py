"""Reject corrupted committed data before allocating an MD context."""
import json
import sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from modules.utils.core import atomic_json, digest, WorkflowError
from modules.md.openmm_engine import run_replica

@pytest.mark.parametrize('case,expected', [
    ('corrupt_segment', 'Committed trajectory'),
    ('wrong_step_sum', 'Committed segment steps'),
    ('corrupt_checkpoint', 'Checkpoint checksum'),
])
def test_reject_corrupt_history_before_context(tmp_path, case, expected):
    base = tmp_path
    (base/'system').mkdir()
    rep = base/'replica_01'
    rep.mkdir()
    # This deliberately invalid XML must never reach the OpenMM deserializer.
    (base/'system/system.xml').write_text('not an OpenMM System')
    for name in ('segment.dcd', 'segment.csv', 'segment.chk'):
        (rep/name).write_bytes(b'original committed bytes')
    segment = dict(dcd='segment.dcd', csv='segment.csv', steps=100,
                   dcd_sha256=digest(rep/'segment.dcd'), csv_sha256=digest(rep/'segment.csv'))
    state = dict(phase='production', production_steps=100, segments=[segment], seed=1,
                 system_sha256=digest(base/'system/system.xml'), current_step=100,
                 checkpoint='segment.chk', checkpoint_sha256=digest(rep/'segment.chk'))
    if case == 'corrupt_segment':
        (rep/'segment.csv').write_bytes(b'corrupt')
    elif case == 'wrong_step_sum':
        state['production_steps'] = 101
    else:
        (rep/'segment.chk').write_bytes(b'corrupt')
    atomic_json(rep/'progress.json', state)
    with pytest.raises(WorkflowError, match=expected):
        run_replica({'MD': {}}, {'md': base}, 1)
