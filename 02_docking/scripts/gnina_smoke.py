import sys,json,time
from sbvs.paths import ROOT,write_json,read_json
from sbvs.config import load
from sbvs.runner import create,run
device=sys.argv[1];cfg=load('research_gnina');cfg.update(device=device,exhaustiveness=4,num_modes=3,per_ligand_timeout_seconds=180)
source=ROOT/'data/demo/cco.smi';source.write_text('CCO\n');start=time.monotonic();out=create(source,cfg);r=run(out)
audit={'device':device,'status':r['status'],'run':str(out.relative_to(ROOT)),'elapsed_seconds':time.monotonic()-start,'counts':r['counts'],'failures':r['failures']};write_json(ROOT/f'reports/gnina_{device}_smoke.json',audit)
caps=ROOT/'reports/gnina_capabilities.json';d=read_json(caps) if caps.exists() else {};d[device+'_passed']=r['status']=='completed';d[device+'_evidence']=audit['run'];write_json(caps,d);print(json.dumps(audit));sys.exit(0 if r['status']=='completed' else 2)
