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
 body='<h1>本地 SBVS 课堂报告</h1><p>Run '+esc(m['run_id'])+'</p><h2>真实运行与环境</h2><pre>'+esc(json.dumps(m,ensure_ascii=False,indent=2))+'</pre>'
 body+='<h2>受体和盒子</h2><p>内置案例来源：RCSB 1DRF；科研自定义请以配置与准备记录为准。</p><pre>'+esc(json.dumps(read_json(rootpath(m['config']['box'])),indent=2))+'</pre>'
 body+='<h2>母体 Top10</h2><table><tr><th>名称</th><th>Vina kcal/mol</th><th>CNN pose</th><th>CNN affinity pK</th></tr>'
 for r in ranking[:10]:body+='<tr>'+''.join('<td>'+esc(r.get(k,''))+'</td>' for k in ('original_name','vina_affinity_kcal_mol','cnn_pose_score','cnn_affinity_pK'))+'</tr>'
 body+='</table><h2>实际姿态图</h2>'
 for p in images:body+='<img width="650" src="data:image/png;base64,'+__import__('base64').b64encode(p.read_bytes()).decode()+'">'
 if len(images)<2:body+='<p>不足两张实际姿态截图；图像验收尚未完成。</p>'
 body+='<h2>失败清单</h2><pre>'+esc(json.dumps(fail,ensure_ascii=False,indent=2))+'</pre><h2>解释与边界</h2><p>演示库没有实验活性标签，未评估富集。Vina 为同协议对接打分；CNN pose score 不是结合概率，CNN affinity 是预测 pK，不能与 kcal/mol 相加。姿态 QC 为几何启发式。加氢不等于确定 pH 微观状态。准备与原始命令详见本 run 清单。未进行 MD。</p><h2>引用</h2><p>AutoDock Vina: Eberhardt et al., 2021, doi:10.1021/acs.jcim.1c00203。GNINA: McNutt et al., 2021, doi:10.1186/s13321-021-00522-2。Meeko 与 RDKit 版本见环境锁文件及官方来源。</p>'
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
 (out/'README.md').write_text('# 人工姿态交付\nligand_selected.sdf 是化学主文件。complex_preview.pdb 仅为同坐标系预览，不是已参数化 MD 体系。请在独立 MD 项目审查质子化、金属/辅因子、参数及连接。\n',encoding='utf-8');return out
