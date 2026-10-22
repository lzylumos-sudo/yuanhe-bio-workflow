"""Public engineering validation using real Yuanhe adapters; no private inputs."""
import os, sys, json, csv, re, time, shutil, subprocess, hashlib
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(sys.argv[1]).resolve()
sys.path[:0]=[str(ROOT/'02_docking/src'),str(ROOT/'00_workbench'),str(ROOT/'05_molecular_dynamics')]
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OPENMM_CPU_THREADS='1',MPLBACKEND='Agg')
def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(data,indent=2,allow_nan=False))
def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def checkpoint(name,values):
    write(OUT/(name+'.json'),values);print(name,json.dumps(values),flush=True)
def config():
    import yaml
    cfg=yaml.safe_load((ROOT/'05_molecular_dynamics/config/default.yaml').read_text())
    cfg.update(replicates=1,ligand={'resname':'LIG'},system_type='protein_ligand')
    cfg['MD'].update(production='0.2ps',equilibration='0.2ps',nvt='0.1ps',npt='0.1ps',heating_duration='0.1ps',
        timestep='1fs',segment_length='0.1ps',report_interval='0.05ps',checkpoint_interval='0.05ps',
        minimization_max_iterations=100,padding_nm=1.0,nonbonded_cutoff_nm=0.8,switch_distance_nm=0.7,
        platform='CUDA',require_gpu=True,cpu_threads=1,qm_threads=1,qm_memory='1 GB',seed=20261021)
    cfg['mmpbsa']['enabled']=False;cfg['mmgbsa']['enabled']=False
    cfg['analysis'].update(sasa_stride=1,sasa_points=20)
    return cfg
def paths():return {k:OUT/'dynamics'/k for k in ('md','analysis','reports')}
def prepare():
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from yuanhe_preflight.config import defaults
    from yuanhe_preflight.pipeline import risk_stage,availability_stage
    from modules.automation import backend
    inputs=OUT/'inputs';inputs.mkdir(exist_ok=True)
    mol=Chem.AddHs(Chem.MolFromSmiles('CCO'));p=AllChem.ETKDGv3();p.randomSeed=17;p.numThreads=1
    assert AllChem.EmbedMolecule(mol,p)==0
    AllChem.UFFOptimizeMolecule(mol,maxIters=50)
    sdf=inputs/'ethanol.sdf'
    with Chem.SDWriter(str(sdf)) as w:w.write(mol)
    # Same synthesis-only gate as production. Models live in independent environments.
    config_path=os.environ.get('YUANHE_PREDOCKING_CONFIG')
    if not config_path:raise ValueError('YUANHE_PREDOCKING_CONFIG must point to an installed real-model configuration')
    from yuanhe_preflight.config import validate,absolute_paths
    pre=absolute_paths(validate(read(Path(config_path))),ROOT)
    pre['cache_dir']=str(OUT/'predocking_cache');pre['workers']=1;pre['batch_size']=1
    risk_stage([{'id':'public_ethanol','smiles':'CCO','sdf':str(sdf),'source_core':'public-input'}],pre,OUT/'risk',base=OUT,batch_id='PUBLIC_ENGINEERING_ONLY')
    availability_stage(read(OUT/'risk/records.json'),pre,OUT/'availability')
    assert len(read(OUT/'availability/qualified.json'))==1
    cfg=config();pcfg=dict(cfg,ligand=None,system_type='protein_protein')
    plan={'protein_path':str(OUT/'sources/1CRN.pdb'),'variants':[],
          'decisions':{'construct_confirmed':True,'termini':'standard_charged'},
          'review_note':'Public 1CRN complete 46-residue chain; engineering-only solvent test, not a binding hypothesis.'}
    prep=backend.prepare(pcfg,paths(),plan)
    write(OUT/'md_config.json',cfg);write(OUT/'protein_plan.json',plan)
    checkpoint('preparation',{'status':'passed','protein':'1CRN','ligand':'ethanol','molecules':1,'synthesis_evidence':read(OUT/'availability/qualified.json')[0]['qualification_sha256'],'protein_preparation':prep})
def dock():
    import numpy as np
    from openmm import app,unit
    from sbvs2.docking.vina_official import run_vina_job
    from sbvs2.docking.protocol import Protocol
    from sbvs2.docking.box import Box
    receptor=paths()['md']/'prepared/protein.pdb'
    pdb=app.PDBFile(str(receptor));xyz=np.array(pdb.positions.value_in_unit(unit.angstrom))
    # Tiny explicit surface region; not an experimentally identified binding pocket.
    heavy=[a.index for a in pdb.topology.atoms() if a.element.symbol!='H']
    extreme=xyz[heavy][np.argmax(xyz[heavy,0])]
    center=extreme+np.array([2.,0.,0.]);box=Box(tuple(center.tolist()),(10.,10.,10.),'manual')
    q=read(OUT/'availability/qualified.json')[0]
    proto=Protocol.preset('quick',exhaustiveness=1,num_modes=1,seeds=(17,))
    result=run_vina_job(OUT/'sources/vina_1.2.7_linux_x86_64',receptor,Path(q['qualified_ligand']),box,proto,17,OUT/'docking',cpu=1,timeout=180,qualification_evidence=q['qualification_evidence'])
    checkpoint('docking_result',result)
    assert result['status']=='ok',result
    write(OUT/'box.json',{'center':list(box.center),'size':list(box.size),'source':'engineering-only public 1CRN surface region'})
def screen():
    import yaml
    folder=OUT/'screening_input';folder.mkdir(exist_ok=True)
    for src,dest in [(OUT/'docking/poses.sdf','poses.sdf'),(OUT/'docking/poses.sdf.qualification.json','poses.sdf.qualification.json'),
                     (paths()['md']/'prepared/protein.pdb','receptor.pdb'),(OUT/'box.json','box.json')]:shutil.copyfile(src,folder/dest)
    target=yaml.safe_load((ROOT/'03_virtual_screening/configs/targets/generic.yaml').read_text())
    target.update(target_id='generic',display_name='1CRN ethanol engineering only',native_length=46,
                  construct={'declared_length':46,'native_range':[1,46],'label_native_offset':0})
    target['structures']['primary']['pdb_id']='1CRN'
    (OUT/'target.yaml').write_text(yaml.safe_dump(target))
    budget=yaml.safe_load((ROOT/'03_virtual_screening/configs/budgets/quick.yaml').read_text())
    budget.update(deep_max_parents=1,md_max_parents=1,poses_per_parent_max=1)
    (OUT/'budget.yaml').write_text(yaml.safe_dump(budget))
    write(OUT/'host.json',{'paths':{'structures_dir':str(OUT/'sources')}})
    cmd=[sys.executable,'-m','zvss.cli','--root',str(OUT/'screening'),'run','--input',str(folder),'--target',str(OUT/'target.yaml'),'--budget',str(OUT/'budget.yaml'),'--host',str(OUT/'host.json'),'--run-id','public_minimal']
    with (OUT/'screening.log').open('w') as log:
        r=subprocess.run(cmd,env={**os.environ,'PYTHONPATH':str(ROOT/'03_virtual_screening/src')+':'+str(ROOT/'02_docking/src')},stdout=log,stderr=subprocess.STDOUT,timeout=240)
    assert r.returncode==0,'See screening.log'
    checkpoint('screening_result',{'status':'executed','formal_qualification':'Preserve actual screening result; no forced pass',
        'md_selection':'One public docked pose explicitly selected for engineering validation only, not a formal scientific candidate.'})
def build():
    from modules.automation import backend
    from rdkit import Chem
    sdf=OUT/'docking/poses.sdf';mol=Chem.SDMolSupplier(str(sdf),removeHs=False)[0]
    assert mol is not None and Chem.MolToSmiles(Chem.RemoveHs(mol))=='CCO'
    shutil.copyfile(sdf,paths()['md']/'prepared/ligand.sdf')
    write(OUT/'engineering_handoff.json',{'selected_pose':str(sdf),'sha256':sha(sdf),'selection_type':'explicit public engineering control',
        'formal_candidate':False,'scientific_claim':None,'receptor_sha256':sha(paths()['md']/'prepared/protein.pdb')})
    checkpoint('md_build',backend.build(config(),paths()))
def dynamics():
    from modules.automation import backend
    backend.run(config(),paths())
    checkpoint('md_execution',read(paths()['md']/'replica_01/progress.json'))
def analyze():
    from modules.automation import backend
    backend.analyze(config(),paths())
    audit=backend.audit_run(config(),paths(),reload_checkpoint=True)
    checkpoint('md_audit',audit)
    assert audit['passed'] and not audit['formal_duration_met']
if __name__=='__main__':
    phase=sys.argv[2];start=time.monotonic();globals()[phase]()
    write(OUT/(phase+'_timing.json'),{'seconds':time.monotonic()-start});print('COMPLETED',phase,flush=True)
