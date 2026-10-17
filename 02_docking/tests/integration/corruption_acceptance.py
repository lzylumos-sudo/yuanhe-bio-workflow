from unittest.mock import patch
from sbvs.paths import ROOT,read_json,write_json
from sbvs.config import load
from sbvs.runner import create,run
source=ROOT/'data/demo/recovery_pair.smi';source.write_text('CCO first\nCCC second\n')
out=create(source,load('class_quick'));first=run(out);assert first['status']=='completed'
keys=list(first['completed_entries']);broken=keys[0];kept=keys[1];old_kept=first['completed_entries'][kept]
pose=out/'entries'/broken/'pose_001.sdf';pose.write_text('CORRUPTED_TEST_FIXTURE')
with patch('sbvs.engines.vina.dock',side_effect=RuntimeError('INJECTED_SINGLE_ENTRY_FAILURE')):
 result=run(out,resume=True)
assert result['status']=='partial_failed' and result['counts']['successful_entries']==1 and result['counts']['failed_entries']==1
assert result['completed_entries'][kept]==old_kept and broken not in result['completed_entries']
final=run(out,resume=True);assert final['status']=='completed' and final['counts']['successful_entries']==2 and final['completed_entries'][kept]==old_kept
write_json(ROOT/'reports/corruption_acceptance.json',{'status':'passed','run':str(out.relative_to(ROOT)),'verified_reused_entry':kept,'failed_entry_not_counted_success':True,'real_retry_completed':True});print(out)
