"""Label-free computation. This module never reads evaluation data."""
from pathlib import Path
from rdkit import Chem
from scientific import artifacts as a
from scientific.peptide_contacts import read_heavy_atoms, contact_signature, contact_overlap
from .library import ppi_pass, coordinates_retain_state
from .common import write_csv

PROFILES={'current_product':{'risk':'existing defaults, exclude/review block','acquisition':'existing default evidence gate',
    'property_preferences':'advisory project preferences; not an implemented hard MW filter'},
    'ppi_benchmark':{'MW':[150,800],'cLogP':[-1,7.5],'TPSA_max':180,'RotB_max':18,'HBD_max':6,'HBA_max':15,
    'formal_charge':[-2,2],'risk':'annotate existing alerts','acquisition':'annotate unknown; required for deliverable list'}}

def filter_library(root):
    from yuanhe_preflight.config import defaults
    from yuanhe_preflight.toxicity_substructure_filter import RiskRules
    from yuanhe_preflight.availability_synthesis_filter import AvailabilityFilter
    root=Path(root);cfg=defaults(root);rules=RiskRules(cfg['rules_path'],cfg['risk_policies'])
    availability=AvailabilityFilter(cfg,rules.sha256)
    a.atomic(root/'inputs/product_filter_config.json',cfg);a.atomic(root/'inputs/filter_profiles.json',PROFILES)
    result=[]
    for row in a.read(root/'inputs/library.json'):
        risk=rules.screen(row['smiles']);supply=availability.screen(row['smiles'])
        result.append({'parent_id':row['parent_id'],'ppi_pass':ppi_pass(row['properties']),
                       'risk':risk,'acquisition':supply,'product_pass':risk['decision']=='pass' and supply['decision']=='pass',
                       'deliverable_pass':risk['decision']=='pass' and supply['decision']=='pass'})
    a.atomic(root/'inputs/filters.json',result)
    return result

def rank(rows,arm):
    if arm=='A0':key=lambda x:(x['vina_affinity'],x['parent_id'])
    elif arm in ('A1','A2'):key=lambda x:(-x[arm+'_recall'],x['vina_affinity'],x['parent_id'])
    else:raise ValueError('Unknown arm')
    return [{**r,'rank':i,'arm':arm} for i,r in enumerate(sorted(rows,key=key),1)]

def score(root):
    root=Path(root);receptor=read_heavy_atoms(root/'target/receptor_original.pdb')
    true=set(a.read(root/'target/peptide_features.json')['residues'])
    perturbed=set(a.read(root/'target/perturbed_features.json')['residues'])
    prepared={r['parent_id']:r for r in a.read(root/'inputs/prepared.json')}
    jobs={r['parent_id']:r for r in a.read(root/'inputs/job_records.json',[])}
    rows=[]
    for row in a.read(root/'inputs/library.json'):
        pid=row['parent_id'];prep=prepared[pid];job=jobs.get(pid,{})
        result={'parent_id':pid,'state_id':prep.get('state_id'),'preparation_status':prep['status'],
                'docking_status':job.get('status','not_started'),'error':prep.get('error') or job.get('error'),
                'vina_affinity':None,'pose_id':None,'pose_sdf':None,'pose_index':None,'A1_recall':None,'A2_recall':None}
        if job.get('status')=='ok':
            path=Path(job['job_dir'])/'poses.sdf'
            if a.digest(path)!=job['output_sha256']:raise ValueError('Docking output modified: '+pid)
            index=min(range(len(job['poses'])),key=lambda i:job['poses'][i]['vina_affinity'])
            mol=Chem.SDMolSupplier(str(path),removeHs=True)[index]
            if mol is None:raise ValueError('Invalid output molecule '+pid)
            result.update(raw_vina_status='ok',raw_vina_affinity=job['poses'][index]['vina_affinity'],pose_id=pid+'_pose_'+str(index+1),pose_sdf=str(path),pose_index=index,pose_sha256=a.digest(path))
            result['selected_pose_stereo_pass']=coordinates_retain_state(mol,prep['state_smiles'])
            if not result['selected_pose_stereo_pass']:
                result.update(docking_status='failed_pose_stereo_qc',error='Selected minimum-affinity pose changes specified stereochemistry; no alternative pose substitution')
                rows.append(result);continue
            contacts=contact_signature(receptor,mol.GetConformer().GetPositions())
            result.update(vina_affinity=job['poses'][index]['vina_affinity'],pose_id=pid+'_pose_'+str(index+1),
                          pose_sdf=str(path),pose_index=index,pose_sha256=a.digest(path),
                          A1_recall=contact_overlap(true,contacts)['recall'],A2_recall=contact_overlap(perturbed,contacts)['recall'],
                          contact_residues=sorted(contacts))
        rows.append(result)
    a.atomic(root/'screening_scores.json',rows);return rows
