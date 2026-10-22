"""Source curation/preparation. Only anonymous chemistry is exported to compute."""
import json
from pathlib import Path
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
from scientific import artifacts as a
from .common import ROOT, sha, write_csv


def canonical_parent(smiles):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None:raise ValueError('invalid_smiles')
    mol=rdMolStandardize.LargestFragmentChooser(preferOrganic=True).choose(mol)
    Chem.SanitizeMol(mol);Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
    return Chem.MolToSmiles(mol,isomericSmiles=True)


def scaffold(smiles):
    mol=Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol),isomericSmiles=False)


def properties(smiles):
    from src.chemistry import calculate_properties
    return calculate_properties(smiles)


def all_activities(root):
    return [row for p in sorted((Path(root)/'evaluation/raw').glob('activity_*.json'))
            if '.source.' not in p.name for row in json.loads(p.read_text())['activities']]


def reference_smiles(root):
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from .structure import category
    result={}
    for ccd,name in [('NUT','Nutlin-3a'),('7HC','SAR405838'),('G7I','Brigimadlin')]:
        path=Path(root)/f'evaluation/ccd/{ccd}.cif';d=MMCIF2Dict(str(path))
        desc=category(d,'_pdbx_chem_comp_descriptor')
        candidates=[r for r in desc if r['type']=='SMILES_CANONICAL' and 'OpenEye' in r['program']]
        if not candidates:raise ValueError('Missing stereo-aware CCD descriptor '+ccd)
        result[name]={'smiles':canonical_parent(candidates[0]['descriptor']),'ccd':ccd,'ccd_sha256':a.digest(path)}
    # CCD 7HC's isomeric descriptor omits trans-1,4-cyclohexanol relative
    # stereochemistry. The primary Figure 1a and ChEMBL2381408 specify it.
    sar={r['canonical_smiles'] for r in all_activities(root) if r['molecule_chembl_id']=='CHEMBL2381408'}
    if len(sar)!=1:raise ValueError('Ambiguous SAR405838 source structure')
    s=canonical_parent(next(iter(sar)))
    flat=lambda x:Chem.MolToSmiles(Chem.MolFromSmiles(x),isomericSmiles=False)
    if flat(s)!=flat(result['SAR405838']['smiles']):raise ValueError('SAR405838 CCD connectivity mismatch')
    result['SAR405838'].update(smiles=s,stereo_source='CHEMBL2381408 + primary Figure 1a, trans-cyclohexanol',
                              ccd_descriptor_limitation='relative cyclohexanol stereo absent; not accepted as complete identity')
    return result


def background_pool(root):
    from src.chemistry import calculate_properties
    activities=all_activities(root)
    known_ids={a['molecule_chembl_id'] for a in activities}
    known_smiles=set();known_scaffolds=set()
    for r in activities:
        try:
            s=canonical_parent(r['canonical_smiles']);known_smiles.add(s);known_scaffolds.add(scaffold(s))
        except Exception:pass
    result={}
    for p in sorted((Path(root)/'evaluation/raw').glob('background_*.json')):
        if '.source.' in p.name:continue
        for r in json.loads(p.read_text())['molecules']:
            try:
                if r['molecule_chembl_id'] in known_ids:continue
                raw=(r.get('molecule_structures') or {}).get('canonical_smiles')
                s=canonical_parent(raw)
                if s in known_smiles or s in result:continue
                m=Chem.MolFromSmiles(s)
                if any(x[1]=='?' for x in Chem.FindMolChiralCenters(m,includeUnassigned=True)):continue
                if any(atom.GetAtomicNum() not in {1,6,7,8,9,15,16,17,35,53} for atom in m.GetAtoms()):continue
                props=calculate_properties(s)
                if not ppi_pass(props):continue
                result[s]={'smiles':s,'source_id':r['molecule_chembl_id'],'original_smiles':raw,
                           'properties':props,'scaffold':scaffold(s),'source_file':str(p.relative_to(root)),
                           'known_mdm2_scaffold':scaffold(s) in known_scaffolds}
            except Exception:continue
    return list(result.values())


def ppi_pass(p):
    return (150<=p['MW']<=800 and -1<=p['cLogP']<=7.5 and p['TPSA']<=180 and p['RotB']<=18
            and p['HBD']<=6 and p['HBA']<=15 and -2<=p['formal_charge']<=2)


STATE_RULE={'algorithm':'molscrub 0.2.2 acid/base at pH 7.4 and tautomer enumeration',
            'selection':'highest RDKit tautomer score; lowest absolute charge; canonical isomeric SMILES',
            'max_states':1,'conformers':1,'etkdg_seed':20261023,'forcefield':'MMFF94s',
            'max_forcefield_iterations':400,'crystal_coordinates_as_input':False}


def prepare_one(row,folder):
    from molscrub import Scrub
    from sbvs2.docking.pdbqt import ligand_to_pdbqt
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    identifier=row['parent_id'];record_path=folder/(identifier+'.json')
    if record_path.exists():
        saved=a.read(record_path)
        if saved.get('input_smiles')!=row['smiles'] or saved.get('state_rule')!=STATE_RULE:raise ValueError('Preparation cache configuration changed')
        if saved.get('status')=='ok' and a.digest(folder/saved['sdf'])!=saved['sha256']:raise ValueError('Preparation cache SDF changed')
        return saved
    result={'parent_id':identifier,'input_smiles':row['smiles'],'state_rule':STATE_RULE,'status':'failed'}
    try:
        molecule=Chem.MolFromSmiles(row['smiles'])
        enumerated=Scrub(ph_low=7.4,ph_high=7.4,skip_gen3d=True)(molecule)
        if not enumerated:raise ValueError('no_protonation_state')
        enumerator=rdMolStandardize.TautomerEnumerator()
        selected=sorted(enumerated,key=lambda m:(-enumerator.ScoreTautomer(m),abs(Chem.GetFormalCharge(m)),Chem.MolToSmiles(m,isomericSmiles=True)))[0]
        state_smiles=Chem.MolToSmiles(selected,isomericSmiles=True)
        gen=Scrub(skip_acidbase=True,skip_tautomers=True,numconfs=1,etkdg_rng_seed=20261023,max_ff_iter=400)
        prepared=gen(Chem.MolFromSmiles(state_smiles))[0]
        # Strip all data-source names and molecule properties before compute.
        for name in list(prepared.GetPropNames(includePrivate=True,includeComputed=False)):prepared.ClearProp(name)
        prepared.SetProp('_Name',identifier);prepared.SetProp('parent_id',identifier)
        state_id='s_'+sha([identifier,state_smiles])[:20];prepared.SetProp('state_id',state_id)
        sdf=folder/(identifier+'.sdf');temporary=folder/(identifier+'.sdf.tmp')
        with Chem.SDWriter(str(temporary)) as writer:writer.write(prepared)
        temporary.replace(sdf)
        text,mapping=ligand_to_pdbqt(prepared,0);(folder/(identifier+'.pdbqt')).write_text(text)
        result.update(status='ok',sdf=sdf.name,sha256=a.digest(sdf),state_id=state_id,state_smiles=state_smiles,
                      enumerated_state_count=len(enumerated),selected_state_count=1,mapping=mapping)
    except Exception as exc:result['error']=type(exc).__name__+': '+str(exc)
    a.atomic(record_path,result);return result


def prepare_development(root):
    root=Path(root);out=root/'development';out.mkdir(exist_ok=True)
    if (out/'prepared.json').exists():return a.read(out/'prepared.json')
    pool=background_pool(root);refs=reference_smiles(root)
    generator=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
    fps=[generator.GetFingerprint(Chem.MolFromSmiles(r['smiles'])) for r in refs.values()]
    eligible=[]
    for row in pool:
        if row['known_mdm2_scaffold']:continue
        fp=generator.GetFingerprint(Chem.MolFromSmiles(row['smiles']))
        if max(DataStructs.BulkTanimotoSimilarity(fp,fps))>=.4:continue
        if 350<=row['properties']['MW']<=700:eligible.append(row)
    eligible.sort(key=lambda r:sha([20261023,r['smiles']]))
    picked=[];seen=set()
    for row in eligible:
        if row['scaffold'] in seen:continue
        seen.add(row['scaffold']);picked.append(row)
        if len(picked)==50:break
    if len(picked)!=50:raise ValueError('Insufficient independent development background molecules')
    for i,r in enumerate(picked):r['parent_id']=f'D{i+1:04d}'
    a.atomic(out/'source_selection.json',picked)
    rows=[]
    for row in picked:
        result=prepare_one(row,out/'ligands');rows.append(result)
        a.atomic(out/'preparation_progress.json',{'completed':len(rows),'input':50})
    a.atomic(out/'prepared.json',rows);return rows


def prepare(root,count=550):
    root=Path(root);out=root/'inputs';out.mkdir(exist_ok=True)
    if (out/'preregistered_protocol.json').exists():raise ValueError('Formal library already frozen; create a new run to change inputs')
    curated=a.read(root/'evaluation/curated_actives.json')
    if not curated or not all(r.get('primary_source_verified') for r in curated):
        raise ValueError('Primary-source verified curation required')
    development=a.read(root/'development/source_selection.json')
    if not development:raise ValueError('Independent development selection required')
    # Matching and labels are performed by the curator, not by the screener.
    seen=set();actives=[]
    for row in curated:
        s=canonical_parent(row['smiles'])
        if s in seen:raise ValueError('Curated parent duplicate')
        seen.add(s);actives.append({**row,'smiles':s,'properties':properties(s),'scaffold':scaffold(s),'label':'active'})
    if not 20<=len(actives)<=30:raise ValueError('Quick benchmark requires 20-30 verified actives')
    banned={r['smiles'] for r in development};banned_scaffolds={r['scaffold'] for r in development}
    if any(r['smiles'] in banned or r['scaffold'] in banned_scaffolds for r in actives):raise ValueError('Development/test chemical series overlap')
    pool=[r for r in background_pool(root) if r['smiles'] not in seen|banned and r['scaffold'] not in banned_scaffolds]
    names=['MW','cLogP','TPSA','HBD','HBA','RotB','formal_charge'];scale=np.array([100,1.5,30,1,2,3,1])
    values=np.array([[r['properties'][k] for k in names] for r in pool])/scale
    targets=np.array([[r['properties'][k] for k in names] for r in actives])/scale
    rng=np.random.default_rng(20261023);order=rng.permutation(len(actives));selected=[];used=set();scaffold_counts={}
    for i in range(max(0,count-len(actives))):
        target=targets[order[i%len(order)]];dist=((values-target)**2).sum(1)
        candidate=None
        for idx in np.argsort(dist,kind='stable'):
            idx=int(idx);s=pool[idx]['scaffold']
            if idx not in used and scaffold_counts.get(s,0)<3:candidate=idx;break
        if candidate is None:break
        used.add(candidate);r=pool[candidate];scaffold_counts[r['scaffold']]=scaffold_counts.get(r['scaffold'],0)+1
        selected.append({**r,'name':r['source_id'],'label':'unknown_background','matching_distance':float(dist[candidate])})
    if len(selected)<50:raise ValueError('Insufficient independent property-matched background pool')
    rows=actives+selected;rows.sort(key=lambda r:sha([20261023,r['smiles']]))
    labels=[];anonymous=[]
    for i,row in enumerate(rows,1):
        identifier=f'M{i:04d}';labels.append({**row,'parent_id':identifier,'run_id':root.name,'profile':'all','arm':'all','seed':20261023})
        anonymous.append({'parent_id':identifier,'smiles':row['smiles'],'properties':row['properties']})
    a.atomic(root/'evaluation/labels.json',labels);write_csv(root/'evaluation/labels.csv',labels)
    write_csv(root/'assay_evidence.csv',[r for r in labels if r['label']=='active'])
    a.atomic(out/'library.json',anonymous)
    with Chem.SDWriter(str(out/'blinded_library.sdf')) as writer:
        for row in anonymous:
            mol=Chem.MolFromSmiles(row['smiles']);mol.SetProp('_Name',row['parent_id']);mol.SetProp('parent_id',row['parent_id']);writer.write(mol)
    # Membership is fixed before preparation. Failures remain in all denominators.
    def summary(items):
        matrix=np.array([[r['properties'][k] for k in names] for r in items])
        return {'n':len(items),'mean':dict(zip(names,matrix.mean(0).tolist())),
                'std':dict(zip(names,matrix.std(0).tolist()))}
    ap,bp=summary(actives),summary(selected)
    diagnostic={'active':ap,'unknown_background':bp,'pool_size':len(pool),'selection_seed':20261023,
                'algorithm':'round-robin nearest properties, fixed scales, maximum 3 per exact Murcko scaffold',
                'scales':dict(zip(names,scale.tolist())),
                'standardized_mean_difference':{k:(ap['mean'][k]-bp['mean'][k])/max(1e-8,((ap['std'][k]**2+bp['std'][k]**2)/2)**.5) for k in names},
                'unknown_background_is_experimental_negative':False,'development_parent_count':50,
                'development_exact_scaffolds_excluded':True,'formal_parent_count':len(rows)}
    a.atomic(root/'evaluation/matching_diagnostics.json',diagnostic)
    prepared=[]
    for row in anonymous:
        prepared.append(prepare_one(row,out/'ligands'))
        a.atomic(out/'preparation_progress.json',{'completed':len(prepared),'input':len(anonymous)})
    a.atomic(out/'prepared.json',prepared)
    return {'input_parents':len(rows),'verified_actives':len(actives),'unknown_background':len(selected),
            'prepared':sum(r['status']=='ok' for r in prepared),'preparation_failures':sum(r['status']!='ok' for r in prepared)}
