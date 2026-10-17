import json,time
from sbvs.paths import ROOT,write_json
from sbvs.config import load
from sbvs.runner import create,run
source=ROOT/'data/demo/cco.smi';source.write_text('CCO\n')
start=time.monotonic();out=create(source,load('class_quick'));result=run(out)
audit={'run':str(out.relative_to(ROOT)),'elapsed_seconds':time.monotonic()-start,'status':result['status']}
write_json(ROOT/'reports/cco_cli_smoke.json',audit);print(json.dumps(audit))
assert result['status']=='completed' and audit['elapsed_seconds']<180
