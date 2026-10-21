#!/usr/bin/env python3
"""Read-only deployment checks; does not run a target or install software."""
import argparse,hashlib,json,os,platform,shutil,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def probe(argv,timeout=40,env=None):
 try:
  p=subprocess.run(argv,capture_output=True,text=True,timeout=timeout,env=env)
  return {'ok':p.returncode==0,'returncode':p.returncode,'output':(p.stdout+p.stderr)[-6000:]}
 except Exception as e:return {'ok':False,'error':str(e)}

def main():
 p=argparse.ArgumentParser();p.add_argument('--config');p.add_argument('--out');a=p.parse_args()
 cfg=json.loads(Path(a.config).read_text()) if a.config else {}
 checks={'linux_x86_64':{'ok':platform.system()=='Linux' and platform.machine() in ('x86_64','AMD64'),'platform':platform.platform()}}
 smi=shutil.which('nvidia-smi') or '/usr/lib/wsl/lib/nvidia-smi';checks['nvidia']=probe([smi,'--query-gpu=name,driver_version,memory.total','--format=csv,noheader'])
 candidates={'chem_python':'/opt/yuanhe/envs/chem/bin/python','peptide_python':'/opt/yuanhe/envs/peptide2mol/bin/python','reinvent_executable':'/opt/yuanhe/envs/reinvent4/bin/reinvent','vina':'vina','gmx':'gmx','prep_python':'/opt/yuanhe/envs/prep/bin/python','md_python':'/opt/yuanhe/envs/membrane/bin/python'}
 paths={k:cfg.get(k) or v for k,v in candidates.items()}
 checks['chemistry']=probe([paths['chem_python'],'-c','import rdkit,numpy,yaml,scipy,pandas,gemmi; print(rdkit.__version__)'])
 checks['screening']=probe([paths['chem_python'],'-c','import MDAnalysis,prolif,posebusters,fastapi; print("screening imports OK")'])
 checks['peptide_cuda']=probe([paths['peptide_python'],'-c','import torch,torch_geometric,torch_scatter,torch_cluster,lightning,hydra,rdkit,rootutils,lmdb; assert torch.cuda.is_available(), "CUDA unavailable"; x=torch.ones((8,8),device="cuda"); print(torch.__version__,torch.version.cuda,torch.cuda.get_device_name(0),float((x@x).sum()))'])
 checks['reinvent_cli']=probe([paths['reinvent_executable'],'--help'])
 sys.path.insert(0,str(ROOT/'00_workbench'))
 
 checks['preparation']=probe([paths['prep_python'],'-c','import pdbfixer,openmm; print(openmm.__version__)'])
 checks['membrane']=probe([paths['md_python'],'-c','import openmm,openmmforcefields,openff.toolkit,parmed,MDAnalysis,psi4,resp; print("membrane imports OK")'])
 checks['vina']=probe([paths['vina'],'--version'],env=os.environ.copy());checks['vina']['ok']=checks['vina']['ok'] and 'AutoDock Vina v1.2.7' in checks['vina'].get('output','');checks['gromacs']=probe([paths['gmx'],'--version'])
 gen=ROOT/'01_molecule_generation';model_cfg=json.loads((gen/'configs/project.json').read_text())['models']
 for name,pkey,hkey in [('peptide2mol','checkpoint','checkpoint_sha256'),('reinvent4','prior','prior_sha256')]:
  x=model_cfg[name];path=gen/x[pkey]
  actual=hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
  checks[name+'_weights']={'ok':actual==x[hkey],'sha256':actual,'path':str(path)}
 result={'checks':checks,'paths':paths,'all_runtime_checks_pass':all(x['ok'] for x in checks.values()),'real_model_inference':'not_tested','target_docking':'not_run','membrane_gromacs_preparation':'implemented_not_target_tested','second_machine_migration':'not_tested'}
 text=json.dumps(result,ensure_ascii=False,indent=2);print(text)
 if a.out:Path(a.out).write_text(text,encoding='utf-8')
 return 0 if result['all_runtime_checks_pass'] else 1
if __name__=='__main__':raise SystemExit(main())
