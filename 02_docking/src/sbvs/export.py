import csv,html,json,time,uuid,shutil
from pathlib import Path
from rdkit import Chem
import numpy as np
from .paths import ROOT,rootpath,read_json,write_json,sha
def read_rows(path):
 with Path(path).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def report(run):
 run=rootpath(run);m=read_json(run/'run_manifest.json');ranking=read_rows(run/'ranking.csv');fail=read_rows(run/'failures.csv')
 esc=lambda x:html.escape(str(x));images=sorted((run/'images').glob('*.png')) if (run/'images').exists() else []
 body='<h1>Local SBVS Classroom report</h1><p>Run '+esc(m['run_id'])+'</p><h2>Actual run and environment</h2><pre>'+esc(json.dumps(m,ensure_ascii=False,indent=2))+'</pre>'
 body+='<h2>Receptor and box</h2><p>Bundled case source:RCSB 1DRF; For custom research inputs, consult configuration and preparation records.</p><pre>'+esc(json.dumps(read_json(rootpath(m['config']['box'])),indent=2))+'</pre>'
 body+='<h2>Parent compound Top10</h2><table><tr><th>Name</th><th>Vina kcal/mol</th><th>CNN pose</th><th>CNN affinity pK</th></tr>'
 for r in ranking[:10]:body+='<tr>'+''.join('<td>'+esc(r.get(k,''))+'</td>' for k in ('original_name','vina_affinity_kcal_mol','cnn_pose_score','cnn_affinity_pK'))+'</tr>'
 body+='</table><h2>Actual pose images</h2>'
 for p in images:body+='<img width="650" src="data:image/png;base64,'+__import__('base64').b64encode(p.read_bytes()).decode()+'">'
 if len(images)<2:body+='<p>Fewer than two actual pose screenshots; visual acceptance incomplete.</p>'
 body+='<h2>Failure list</h2><pre>'+esc(json.dumps(fail,ensure_ascii=False,indent=2))+'</pre><h2>Interpretation and limitations</h2><p>The demonstration library lacks activity labels; enrichment was not assessed.Vina is a docking score under the same protocol;CNN pose score is not a binding probability;CNN affinity is a predicted pK, and must not be added to kcal/mol . Pose QC is a geometric heuristic. Adding hydrogens does not establish the pH microstate. See preparation and raw commands in this run manifest. No MD. </p><h2>References</h2><p>AutoDock Vina: Eberhardt et al., 2021, doi:10.1021/acs.jcim.1c00203. GNINA: McNutt et al., 2021, doi:10.1186/s13321-021-00522-2. Meeko and RDKit Versions are recorded in environment locks and official sources.</p>'
 target=run/'report.html';target.write_text('<!doctype html><meta charset="utf-8"><style>body{font:16px sans-serif;max-width:1100px;margin:auto;padding:30px}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f3f5f7;padding:15px}td,th{border:1px solid #ccc;padding:8px}table{border-collapse:collapse}</style>'+body,encoding='utf-8');return target
def handoff(run,selection):
 run=rootpath(run);selection=read_json(rootpath(selection)) if isinstance(selection,(str,Path)) else selection
 if selection.get('human_confirmed') is not True or not selection.get('reviewer'):raise ValueError('Explicit human selection and reviewer required')
 rows=read_rows(run/'poses.csv');chosen=[r for r in rows if r['pose_id']==selection.get('pose_id')]
 if len(chosen)!=1:raise ValueError('Select exactly one real pose')
 row=chosen[0]
 if row['qc_flags'] and not selection.get('acknowledge_warnings'):raise ValueError('Review pose QC warnings before handoff')
 m=read_json(run/'run_manifest.json');cfg=m['config'];pose=(run/row['pose_file']).resolve()
 if not pose.is_relative_to(run):raise ValueError('Pose path escape')
 out=ROOT/'exports/md_handoff'/uuid.uuid4().hex;out.mkdir(parents=True)
 shutil.copyfile(rootpath(cfg['receptor_pdb']),out/'receptor_prepared.pdb');shutil.copyfile(pose,out/'ligand_selected.sdf');shutil.copyfile(rootpath(cfg['box']),out/'box.json')
 ligand=next(iter(Chem.SDMolSupplier(str(pose),removeHs=False)));exported=next(iter(Chem.SDMolSupplier(str(out/'ligand_selected.sdf'),removeHs=False)))
 delta=float(np.max(abs(ligand.GetConformer().GetPositions()-exported.GetConformer().GetPositions())))
 if delta>1e-6:raise ValueError('Export coordinate mismatch')
 rec=(out/'receptor_prepared.pdb').read_text();block=Chem.MolToPDBBlock(ligand);maxserial=max(int(l[6:11]) for l in rec.splitlines() if l.startswith(('ATOM  ','HETATM')));preview=[];mapping=[]
 for i,line in enumerate(l for l in block.splitlines() if l.startswith(('ATOM  ','HETATM'))):
  serial=maxserial+i+1
  if serial>99999:raise ValueError('PDB preview exceeds serial capacity')
  name=ligand.GetAtomWithIdx(i).GetSymbol()+str(i+1)
  if len(name)>4:raise ValueError('PDB preview atom name exceeds capacity')
  preview.append(f'HETATM{serial:5d} {name:>4} LIG L   1    '+line[30:]);mapping.append({'sdf_atom_index':i,'selected_pose_atom_index':i,'preview_atom_name':name,'preview_serial':serial})
 (out/'complex_preview.pdb').write_text('\n'.join(l for l in rec.splitlines() if not l.startswith(('END','CONECT')))+'\nTER\n'+'\n'.join(preview)+'\nEND\n')
 write_json(out/'atom_mapping.json',{'atoms':mapping,'max_sdf_coordinate_delta':delta,'preview_coordinate_precision_A':0.001})
 preparation=read_json(rootpath(cfg['preparation_record'])) if cfg.get('preparation_record') else read_json(ROOT/'data/prepared/demo_1drf/preparation.json')
 write_json(out/'selection_manifest.json',{'selection':selection,'run_id':m['run_id'],'pose':row,'protocol':cfg,'receptor_preparation':preparation,'selected_pose_sha256':sha(pose),'created':time.time(),'warning':'Preview only; no force field parameters, solvation, or MD performed. External receptor chemistry/protonation requires separate review.'})
 (out/'README.md').write_text('# Manual pose handoff\nligand_selected.sdf is the chemical master file.complex_preview.pdb is a same-frame preview, not a parameterized MD system. In a separate MD project, review protonation, metals,/cofactors, parameters and bonds.\n',encoding='utf-8');return out
