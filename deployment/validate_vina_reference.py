#!/usr/bin/env python3
"""Real 1DRF/FOL redocking in the receptor frame; no ligand self-fitting."""
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'02_docking/src'))
def main():
 from rdkit import Chem
 from rdkit.Chem import rdMolAlign
 from sbvs2.docking.vina_official import run_vina_job
 from sbvs2.docking.box import Box
 from sbvs2.docking.protocol import Protocol
 p=argparse.ArgumentParser();p.add_argument('--vina',default='vina');p.add_argument('--output',required=True,type=Path);a=p.parse_args()
 base=ROOT/'02_docking/data';box=json.loads((base/'prepared/demo_1drf/box.json').read_text())
 reference=base/'reference/FOL_crystal.sdf';ref=Chem.RemoveHs(Chem.SDMolSupplier(str(reference),removeHs=False)[0])
 a.output.mkdir(parents=True,exist_ok=False);results=[]
 for seed in [17,41,97]:
  job=run_vina_job(a.vina,base/'prepared/demo_1drf/receptor_hydrogenated.pdb',reference,Box(tuple(box['center']),tuple(box['size']),'1DRF FOL reference'),Protocol.preset('standard',seeds=(seed,),exhaustiveness=32),seed,a.output/f'seed{seed}',cpu=4)
  if job['status']!='ok':raise RuntimeError(job)
  poses=[Chem.RemoveHs(m) for m in Chem.SDMolSupplier(str(a.output/f'seed{seed}'/'poses.sdf'))]
  rmsds=[float(rdMolAlign.CalcRMS(m,ref)) for m in poses]
  scores=[float(m.GetProp('vina_affinity_kcal_mol')) for m in poses]
  if scores!=sorted(scores):raise ValueError('Affinity ranking direction incorrect')
  results.append({'seed':seed,'best_affinity_kcal_mol':scores[0],'top_pose_rmsd_A':rmsds[0],'best_of_modes_rmsd_A':min(rmsds),'all_rmsd_A':rmsds,'all_affinity_kcal_mol':scores})
 result={'engine':'AutoDock Vina','version':'1.2.7','reference':'1DRF chain A crystallographic FOL','rmsd':'symmetry aware heavy atoms, unchanged receptor frame, no ligand self-fit','scope':'redocking engineering reference; no target efficacy inference','seeds':results}
 (a.output/'redocking.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
