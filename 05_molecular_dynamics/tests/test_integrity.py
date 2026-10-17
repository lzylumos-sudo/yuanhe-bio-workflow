"""Reject corrupted committed data before allocating an MD context."""
import json,sys,tempfile
from pathlib import Path
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root))
from modules.utils.core import atomic_json,digest,WorkflowError
from modules.md.openmm_engine import run_replica

with tempfile.TemporaryDirectory(prefix='nexus_integrity_') as temp:
    base=Path(temp);(base/'system').mkdir();rep=base/'replica_01';rep.mkdir()
    (base/'system/system.xml').write_text('not a System: must never be parsed in these corruption tests')
    for name in ['segment.dcd','segment.csv','segment.chk']:(rep/name).write_bytes(b'original committed bytes')
    segment={'dcd':'segment.dcd','csv':'segment.csv','steps':100,'dcd_sha256':digest(rep/'segment.dcd'),'csv_sha256':digest(rep/'segment.csv')}
    state={'phase':'production','production_steps':100,'segments':[segment],'seed':1,
           'system_sha256':digest(base/'system/system.xml'),'current_step':100,'checkpoint':'segment.chk','checkpoint_sha256':digest(rep/'segment.chk')}
    cases=[]
    for label,change,expected in [
        ('corrupt_segment',lambda s:(rep/'segment.csv').write_bytes(b'corrupt'),'Committed trajectory'),
        ('wrong_step_sum',lambda s:s.update(production_steps=101),'Committed segment steps'),
        ('corrupt_checkpoint',lambda s:(rep/'segment.chk').write_bytes(b'corrupt'),'Checkpoint checksum')]:
        for name in ['segment.csv','segment.chk']:(rep/name).write_bytes(b'original committed bytes')
        current=json.loads(json.dumps(state));change(current);atomic_json(rep/'progress.json',current)
        try:run_replica({'MD':{}},{'md':base},1)
        except WorkflowError as exc:
            if expected not in str(exc):raise
            cases.append({'case':label,'status':'PASS','rejected_before_context':True})
        else:raise AssertionError('Corrupt data accepted: '+label)
    result={'passed':True,'cases':cases}
    atomic_json(root/'installation_evidence/integrity_regression.json',result)
    print(json.dumps(result))
