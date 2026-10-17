import json
from sbvs.paths import ROOT,write_json
from sbvs.config import load
from sbvs.runner import create,run
from sbvs.export import read_rows
source=ROOT/'data/demo/state_validation.smi';source.write_text('CC(=O)O acetic_acid\nc1cc[nH]c(=O)c1 pyridone\n')
cfg=load('class_quick');cfg.update(ligand_states='enumerate',ph_low=6.4,ph_high=8.4,max_states=8,seeds=[11,23])
out=create(source,cfg);r=run(out);assert r['status']=='completed';states=read_rows(out/'states.csv');ranking=read_rows(out/'ranking.csv');assert len(ranking)==2 and all(int(s['seed_count'])==2 for s in states)
assert any('-' in __import__('json').loads((out/p['pose_file']).with_name('preparation.json').read_text())['identity'] for p in ranking)
write_json(ROOT/'reports/enumeration_acceptance.json',{'status':'passed','run':str(out.relative_to(ROOT)),'parents':len(ranking),'states':len(states),'seeds':[11,23],'tool':'molscrub 0.2.2','ph_range':[6.4,8.4]});print(out)
