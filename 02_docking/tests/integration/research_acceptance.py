import json,time,requests,numpy as np
from pathlib import Path
from rdkit import Chem
from Bio.PDB import PDBParser,MMCIFIO
from sbvs.paths import ROOT,read_json,write_json
from sbvs.export import read_rows
BASE='http://127.0.0.1:5902';HEAD={'X-SBVS-Request':'local-ui'};results={}
def post(path,**kw):
 r=requests.post(BASE+path,headers=HEAD,timeout=180,**kw);r.raise_for_status();return r.json()
with (ROOT/'data/raw/1DRF.pdb').open('rb') as f:
 plan=post('/api/receptors/plan',files={'receptor':('1DRF.pdb',f)},data={'chains':'A','model_index':'0'})
rid=plan['receptor_id'];policy={'confirmed':True,'box_confirmed':True,'keep_hetero':[],'repair_missing_sidechains':True,'altloc':'A','box':read_json(ROOT/'data/prepared/demo_1drf/box.json'),'reviewer':'AUTOMATED_PUBLIC_DEMO_ACCEPTANCE','test_only':True}
response=post('/api/receptors/'+rid+'/prepare',data={'policy':json.dumps(policy),'box_mode':'manual'})
assert response['status']=='prepared';results['custom_pdb_review']={'status':'passed','receptor_id':rid}
job=post('/api/jobs',data={'receptor_id':rid,'smiles':'CCO','preset':'class_quick','pocket_confirmed':'true','reviewed_box_hash':response['box_hash']})['id']
for _ in range(180):
 d=requests.get(BASE+'/api/jobs/'+job,timeout=10).json()
 if d['status'] in ('completed','failed','partial_failed'):break
 time.sleep(.5)
assert d['status']=='completed';results['custom_receptor_screen']={'status':'passed','run_id':job}
structure=PDBParser(QUIET=True).get_structure('demo',str(ROOT/'data/raw/1DRF.pdb'));io=MMCIFIO();io.set_structure(structure);cif=ROOT/'data/demo/1DRF_conversion_test.cif';io.save(str(cif))
with cif.open('rb') as f:plan=post('/api/receptors/plan',files={'receptor':('demo.cif',f)},data={'chains':'A'})
assert plan['chains']==['A'];results['mmcif_plan']={'status':'passed'}
with (ROOT/'data/raw/1DRF.pdb').open('rb') as f:
 r=requests.post(BASE+'/api/receptors/plan',headers=HEAD,files={'receptor':('demo.pdb',f)},data={'chains':'NONEXISTENT'},timeout=15)
assert r.status_code==400;results['missing_chain']={'status':'passed'}
for op,log in [('score_only','gnina_score_only.log'),('redock','gnina_redock.log')]:
 lines=(ROOT/'reports'/log).read_text().splitlines();out=Path(lines[-1]);m=read_json(out/'run_manifest.json');assert m['status']=='completed'
 rows=read_rows(out/'poses.csv');assert all(r['cnn_pose_score'] and r['cnn_affinity_pK'] and r['empirical_score_kcal_mol'] for r in rows)
 if op=='score_only':
  maximum=0
  for entry in m['completed_entries']:
   a=next(iter(Chem.SDMolSupplier(str(out/'entries'/entry/'input.sdf'),removeHs=False)));b=next(iter(Chem.SDMolSupplier(str(out/'entries'/entry/'poses.sdf'),removeHs=False)));a=Chem.RemoveHs(a);b=Chem.RemoveHs(b)
   match=b.GetSubstructMatch(a);assert len(match)==a.GetNumAtoms();delta=float(np.max(abs(a.GetConformer().GetPositions()-b.GetConformer().GetPositions()[list(match)])));maximum=max(maximum,delta)
  assert maximum<=.0011
  results[op]={'status':'passed','run':str(out.relative_to(ROOT)),'maximum_coordinate_change_A':maximum}
 else:results[op]={'status':'passed','run':str(out.relative_to(ROOT))}
gpu=[]
for p in (ROOT/'runs').glob('*/entries/*/gpu_samples.json'):
 if any(r.get('process_matches') for r in read_json(p)):gpu.append(str(p.relative_to(ROOT)))
results['actual_gpu_process_observed']={'status':'passed' if gpu else 'not_tested','evidence':gpu[:5]}
write_json(ROOT/'reports/research_acceptance.json',results);print(json.dumps(results,indent=2))
