"""Gas-phase HF/6-31G* ESP and two-stage RESP, with a separate optimized conformer.

No AM1-BCC fallback. Input atom order and bound heavy-atom coordinates are immutable.
"""
from __future__ import annotations
import hashlib, importlib.metadata, json, os, subprocess, sys
from pathlib import Path

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path, value): Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')

def compute(mol, output, settings=None):
    import numpy as np
    from rdkit import Chem
    settings = dict(settings or {})
    if settings.get('charge_method','resp') != 'resp': raise ValueError('Formal parameterization requires RESP; no charge fallback')
    out = Path(output).resolve(); out.mkdir(parents=True,exist_ok=True)
    mol = Chem.AddHs(Chem.Mol(mol),addCoords=True)
    if len(Chem.GetMolFrags(mol))!=1 or not mol.GetNumConformers(): raise ValueError('One connected ligand with coordinates required')
    allowed={'H','C','N','O','F','P','S','Cl'}
    unsupported=sorted({a.GetSymbol() for a in mol.GetAtoms()}-allowed)
    if unsupported: raise ValueError('HF/6-31G* organic protocol does not qualify elements: '+','.join(unsupported))
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()) or settings.get('multiplicity',1)!=1:
        raise ValueError('Open-shell/special systems need a separately validated QM protocol')
    if any(flag=='?' for _,flag in Chem.FindMolChiralCenters(mol,includeUnassigned=True)):
        raise ValueError('Undefined ligand stereochemistry')
    net=Chem.GetFormalCharge(mol)
    if settings.get('net_charge',net)!=net: raise ValueError('Declared net charge differs from chemical structure')
    request={'schema':'hf-resp/2','method':'HF/6-31G*','optimization':'HF/6-31G* gas phase',
             'optimizer':{'opt_coordinates':'redundant','dynamic_level':1,'geom_maxiter':50},
             'net_charge':net,'multiplicity':1,'protonation':'explicit supplied chemical state; no state inference',
             'conformer':'single independent optimized copy of selected pose',
             'grid':{'vdw_scale_factors':[1.4,1.6,1.8,2.0],'points_per_A2':1.0},
             'stages':[{'a':0.0005,'b':0.1},{'a':0.001,'b':0.1}],
             'threads':int(settings.get('qm_threads',4)),'memory':settings.get('qm_memory','2 GB'),
             'atom_mapping':[{'index':a.GetIdx(),'element':a.GetSymbol(),'formal_charge':a.GetFormalCharge()} for a in mol.GetAtoms()]}
    source=out/'bound_input.sdf'
    if source.exists():
        # Reusing charges requires identical chemistry AND coordinates AND settings.
        old=Chem.SDMolSupplier(str(source),removeHs=False)[0]
        if old is None or Chem.MolToSmiles(old)!=Chem.MolToSmiles(mol) or not np.allclose(old.GetConformer().GetPositions(),mol.GetConformer().GetPositions(),atol=1e-4,rtol=0):
            raise ValueError('Existing RESP input differs; use a new parameter directory')
    else:
        with Chem.SDWriter(str(source)) as w:w.write(mol)
    request['input_sha256']=sha(source)
    if (out/'charges.json').exists():
        result=json.loads((out/'charges.json').read_text())
        if result['request']!=request: raise ValueError('RESP cache parameters differ')
        for name,expected in result['files'].items():
            if sha(out/name)!=expected:raise ValueError('RESP cached evidence changed: '+name)
    else:
        write(out/'request.json',request)
        with (out/'driver.log').open('w') as log:
            subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker',str(out)],stdout=log,stderr=subprocess.STDOUT,check=True)
        result=json.loads((out/'charges.json').read_text())
    q=np.asarray(result['charges_e'],dtype=float)
    if len(q)!=mol.GetNumAtoms() or not np.isfinite(q).all() or abs(q.sum()-net)>1e-6: raise ValueError('RESP charge vector invalid')
    return q,result

def worker(out):
    import numpy as np
    import psi4, resp
    from rdkit import Chem
    out=Path(out).resolve();os.chdir(out)
    req=json.loads(Path('request.json').read_text());mol=Chem.SDMolSupplier('bound_input.sdf',removeHs=False)[0]
    psi4.set_num_threads(req['threads']);psi4.set_memory(req['memory']);psi4.core.set_output_file('psi4.log',False)
    xyz=mol.GetConformer().GetPositions()
    lines=[f"{req['net_charge']} {req['multiplicity']}"]+[f'{a.GetSymbol()} {p[0]:.12f} {p[1]:.12f} {p[2]:.12f}' for a,p in zip(mol.GetAtoms(),xyz)]
    qm=psi4.geometry('\n'.join(lines+['units angstrom','no_reorient','no_com','symmetry c1']))
    qm.set_name('ligand');psi4.set_options({'basis':'6-31G*','scf_type':'pk','e_convergence':1e-10,'d_convergence':1e-8,'g_convergence':'gau_tight'})
    # OptKing's documented recovery switches coordinates/trust radius when
    # internal-coordinate back transformations fail. Tight convergence and the
    # electronic method remain mandatory; exhaustion still raises an error.
    psi4.set_options(req.get('optimizer', {}))
    energy=psi4.optimize('hf/6-31G*',molecule=qm)
    optimized=Chem.Mol(mol);points=np.asarray(qm.geometry())*psi4.constants.bohr2angstroms
    for i,p in enumerate(points):optimized.GetConformer().SetAtomPosition(i,p.tolist())
    with Chem.SDWriter('qm_optimized.sdf') as w:w.write(optimized)
    options={'METHOD_ESP':'hf','BASIS_ESP':'6-31G*','VDW_SCALE_FACTORS':[1.4,1.6,1.8,2.0],
             'VDW_POINT_DENSITY':1.0,'RESP_A':0.0005,'RESP_B':0.1,'IHFREE':True}
    first=resp.resp([qm],options)[1]
    # Graph-defined alkyl CH2/CH3 groups; never infer bonds from QM distances.
    groups=[];free=set()
    for atom in mol.GetAtoms():
        hs=[n.GetIdx()+1 for n in atom.GetNeighbors() if n.GetAtomicNum()==1]
        if atom.GetAtomicNum()==6 and atom.GetHybridization()==Chem.HybridizationType.SP3 and len(hs)>=2:
            groups.append(hs);free.update([atom.GetIdx()+1,*hs])
    second=dict(options, RESP_A=0.001,
                CONSTRAINT_CHARGE=[[float(first[i]),[i+1]] for i in range(len(first)) if i+1 not in free],
                CONSTRAINT_GROUP=groups, GRID=['1_ligand_grid.dat'],ESP=['1_ligand_grid_esp.dat'])
    q=resp.resp([qm],second)[1] if free else first.copy()
    for name in ['1_ligand_grid.dat','1_ligand_grid_esp.dat']:
        if not Path(name).is_file() or not np.isfinite(np.loadtxt(name)).all(): raise ValueError('Invalid ESP artifact: '+name)
    if abs(sum(q)-req['net_charge'])>1e-6: raise ValueError('RESP net charge failed')
    np.savetxt('resp_charges.dat',q,fmt='%.10f')
    write('charges.json',{'request':req,'charges_e':q.tolist(),'stage1_charges_e':first.tolist(),
        'stage2_constraints':{'equivalent_hydrogens_1based':groups,'fixed_atom_indices_1based':[i+1 for i in range(len(first)) if i+1 not in free]},
        'stage2_status':'fitted' if free else 'all_atoms_fixed_to_stage1',
        'optimized_energy_hartree':float(energy),'charge_sum_e':float(sum(q)),
        'versions':{'psi4':psi4.__version__,**{n:importlib.metadata.version(n) for n in ['resp','rdkit','numpy']}},
        'files':{n:sha(n) for n in ['bound_input.sdf','qm_optimized.sdf','1_ligand_grid.dat','1_ligand_grid_esp.dat','resp_charges.dat']}})

def gaff_topology(mol, output, settings=None):
    """Real antechamber/parmchk2/tleap chain with imported RESP charges."""
    import numpy as np, parmed
    settings=dict(settings or {});out=Path(output).resolve()
    q,report=compute(mol,out,settings)
    ff=settings.get('gaff','gaff2')
    if ff not in {'gaff','gaff2'}: raise ValueError('Explicit GAFF or GAFF2 required')
    commands=[['antechamber','-i','bound_input.sdf','-fi','sdf','-o','ligand.mol2','-fo','mol2','-at',ff,'-c','rc','-cf','resp_charges.dat','-nc',str(report['request']['net_charge']),'-rn','LIG','-pf','y'],
              ['parmchk2','-i','ligand.mol2','-f','mol2','-o','ligand.frcmod','-s','2' if ff=='gaff2' else '1']]
    for command in commands:
        with (out/(command[0]+'.log')).open('w') as log: subprocess.run(command,cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True)
    m=parmed.load_file(str(out/'ligand.mol2'))
    from rdkit import Chem
    source=Chem.SDMolSupplier(str(out/'bound_input.sdf'),removeHs=False)[0]
    if [a.atomic_number for a in m.atoms]!=[a.GetAtomicNum() for a in source.GetAtoms()] or not np.allclose(m.coordinates,source.GetConformer().GetPositions(),atol=1e-3,rtol=0): raise ValueError('Antechamber changed atom order or bound pose')
    if not np.allclose([a.charge for a in m.atoms],q,atol=1e-6,rtol=0): raise ValueError('Antechamber overwrote RESP charges')
    if 'ATTN' in (out/'ligand.frcmod').read_text(): raise ValueError('GAFF missing parameter; inspect frcmod')
    (out/'tleap.in').write_text(f'source leaprc.{ff}\nLIG = loadmol2 ligand.mol2\nloadamberparams ligand.frcmod\ncheck LIG\nsaveamberparm LIG ligand.prmtop ligand.inpcrd\nquit\n')
    with (out/'tleap.log').open('w') as log: subprocess.run(['tleap','-f','tleap.in'],cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True)
    topology=parmed.load_file(str(out/'ligand.prmtop'),str(out/'ligand.inpcrd'))
    if not np.allclose([a.charge for a in topology.atoms],q,atol=1e-6,rtol=0): raise ValueError('Final topology charge mismatch')
    topology.save(str(out/'ligand.top'),overwrite=True)
    converted=parmed.load_file(str(out/'ligand.top'))
    if not np.allclose([a.charge for a in converted.atoms],q,atol=1e-6,rtol=0): raise ValueError('GROMACS conversion changed RESP charges')
    write(out/'topology_validation.json',{'forcefield':ff,'charge_method':'HF/6-31G* two-stage RESP','atoms':len(q),'charge_sum_e':float(sum(q)),
           'bound_pose_preserved':True,'charge_max_error_e':float(max(abs(a.charge-b) for a,b in zip(topology.atoms,q))),
           'files':{n:sha(out/n) for n in ['ligand.prmtop','ligand.inpcrd','ligand.top','ligand.mol2','ligand.frcmod']}})
    return q,report

if __name__=='__main__':
    if len(sys.argv)==3 and sys.argv[1]=='--worker':worker(sys.argv[2])
    else:
        import argparse
        from rdkit import Chem
        p=argparse.ArgumentParser();p.add_argument('--ligand',required=True);p.add_argument('--output',required=True);p.add_argument('--gaff',default='gaff2',choices=['gaff','gaff2']);a=p.parse_args()
        gaff_topology(Chem.SDMolSupplier(a.ligand,removeHs=False)[0],a.output,{'gaff':a.gaff})
