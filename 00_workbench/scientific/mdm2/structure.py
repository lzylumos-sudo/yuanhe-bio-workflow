"""Entity-aware mmCIF extraction and peptide-only pocket definition."""
import csv
import json
from pathlib import Path
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from scientific import artifacts as a
from scientific.peptide_contacts import read_heavy_atoms, contact_signature


def category(data, prefix):
    fields = [k for k in data if k.startswith(prefix + '.')]
    if not fields:
        return []
    return [{k.split('.', 1)[1]: data[k][i] for k in fields} for i in range(len(data[fields[0]]))]


def atom_rows(data):
    return [r for r in category(data, '_atom_site')
            if r.get('pdbx_PDB_model_num', '1') == '1' and r['type_symbol'] not in {'H','D'}
            and r['label_alt_id'] in {'.', 'A', '?'}]


def write_pdb(rows, path):
    lines = []
    for i, r in enumerate(rows, 1):
        name = r['auth_atom_id']; name = (' ' + name).ljust(4) if len(name)<4 else name[:4]
        ins = r.get('pdbx_PDB_ins_code', '?'); ins = ' ' if ins in {'?', '.'} else ins
        lines.append(f"ATOM  {i:5d} {name} {r['auth_comp_id']:>3s} {r['auth_asym_id']:1s}{int(r['auth_seq_id']):4d}{ins:1s}   "
                     f"{float(r['Cartn_x']):8.3f}{float(r['Cartn_y']):8.3f}{float(r['Cartn_z']):8.3f}"
                     f"{float(r['occupancy']):6.2f}{float(r['B_iso_or_equiv']):6.2f}          {r['type_symbol']:>2s}\n")
    Path(path).write_text(''.join(lines) + 'END\n')


def audit_cif(path):
    d = MMCIF2Dict(str(path))
    return {'file': str(path), 'sha256': a.digest(path), 'entities': category(d,'_entity'),
            'polymers': category(d,'_entity_poly'), 'source_gen': category(d,'_entity_src_gen'),
            'source_natural': category(d,'_entity_src_nat'), 'references': category(d,'_struct_ref'),
            'sequence_alignment': category(d,'_struct_ref_seq'),
            'sequence_differences': category(d,'_struct_ref_seq_dif'),
            'unobserved_residues': category(d,'_pdbx_unobs_or_zero_occ_residues'),
            'nonpolymers': category(d,'_pdbx_entity_nonpoly'), 'citations': category(d,'_citation'),
            'resolution_A': d.get('_refine.ls_d_res_high'),
            'audit_status': 'entity_and_sequence_tables_extracted; review sequence differences explicitly'}


def peptide_chemical_features(peptide_path,receptor_atoms):
    """Descriptive RDKit chemical features in the experimental peptide frame.

    These are exported for interface interpretation. The existing benchmark
    ranking remains the frozen residue-contact method, not a fitted pharma score.
    """
    from rdkit import Chem, RDConfig
    from rdkit.Chem import ChemicalFeatures
    mol=Chem.MolFromPDBFile(str(peptide_path),removeHs=False,sanitize=True)
    if mol is None:raise ValueError('Cannot assign peptide chemistry')
    factory=ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir)/'BaseFeatures.fdef'))
    xyz=np.array([p for _,p in receptor_atoms]);coords=mol.GetConformer().GetPositions();features=[]
    for feature in factory.GetFeaturesForMol(mol):
        if feature.GetFamily() not in {'Donor','Acceptor','NegIonizable','PosIonizable','Aromatic','Hydrophobe','LumpedHydrophobe'}:continue
        ids=feature.GetAtomIds()
        if np.linalg.norm(coords[list(ids),None,:]-xyz[None,:,:],axis=2).min()>4.5:continue
        atoms=[]
        for i in ids:
            info=mol.GetAtomWithIdx(i).GetPDBResidueInfo()
            atoms.append({'chain':info.GetChainId().strip(),'residue':info.GetResidueNumber(),'resname':info.GetResidueName().strip(),'atom':info.GetName().strip()})
        pos=feature.GetPos();features.append({'family':feature.GetFamily(),'type':feature.GetType(),'pos':[pos.x,pos.y,pos.z],'atoms':atoms})
    if not features:raise ValueError('No peptide chemical interface features')
    return features


def prepare(root):
    root = Path(root); out = root/'target'; out.mkdir(exist_ok=True)
    source = root/'sources/1YCR.cif'; d = MMCIF2Dict(str(source))
    refs = category(d,'_struct_ref')
    receptor_ids = {r['entity_id'] for r in refs if r.get('pdbx_db_accession')=='Q00987'}
    peptide_ids = {r['entity_id'] for r in refs if r.get('pdbx_db_accession')=='P04637'}
    if len(receptor_ids)!=1 or len(peptide_ids)!=1:
        raise ValueError('Cannot uniquely identify MDM2 and p53 from UniProt entity mappings')
    rows = atom_rows(d)
    def entity_atoms(ids):
        selected = [r for r in rows if r['label_entity_id'] in ids and r['group_PDB']=='ATOM']
        chains = sorted({r['label_asym_id'] for r in selected})
        if len(chains)!=1:
            raise ValueError('Main structure requires an unambiguous entity instance')
        return selected
    rec, pep = entity_atoms(receptor_ids), entity_atoms(peptide_ids)
    write_pdb(rec,out/'receptor_original.pdb'); write_pdb(pep,out/'peptide.pdb')
    mapping={}
    for r in rec+pep:
        key=(r['label_asym_id'],r['label_seq_id'])
        mapping[key]={k:r[k] for k in ['label_entity_id','label_asym_id','auth_asym_id','label_seq_id','auth_seq_id','label_comp_id','auth_comp_id']}
    with (out/'residue_mapping.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(next(iter(mapping.values()))));w.writeheader();w.writerows(mapping.values())
    rec_atoms=read_heavy_atoms(out/'receptor_original.pdb'); pep_atoms=read_heavy_atoms(out/'peptide.pdb')
    xyz=np.asarray([p for _,p in rec_atoms]); pxyz=np.asarray([p for _,p in pep_atoms])
    distances=np.linalg.norm(pxyz[:,None,:]-xyz[None,:,:],axis=2)
    interface=pxyz[distances.min(1)<=4.5]
    if not len(interface):raise ValueError('No peptide interface contacts')
    center=(interface.min(0)+interface.max(0))/2; size=np.maximum(interface.max(0)-interface.min(0)+10,20)
    if size.max()>40:raise ValueError('Unexpectedly large peptide box; inspect entity and geometry')
    box={'center':center.tolist(),'size':size.tolist(),'source':'1YCR peptide contacting heavy atoms only',
         'rule':'peptide heavy atoms <=4.5 A from receptor; bounding box +5 A per face, minimum edge 20 A'}
    a.atomic(out/'box.json',box)
    density={}
    for (key,_),ds in zip(pep_atoms,distances):
        item=density.setdefault(':'.join(key[:3]),{'contact_atom_pairs':0,'contacting_peptide_atoms':0})
        item['contact_atom_pairs']+=int((ds<=4.5).sum());item['contacting_peptide_atoms']+=int(ds.min()<=4.5)
    # Existing peptide contact overlap is used, without MurJ anchor or membrane terms.
    features={'implementation':'scientific.peptide_contacts.contact_signature/contact_overlap',
              'kind':'residue_contact_fingerprint; not a validated chemical pharmacophore equivalence model',
              'automatic':True,'cutoff_A':4.5,'residues':sorted(contact_signature(rec_atoms,pxyz)),
              'peptide_residue_contact_density':density,
              'features':[{'element':r['type_symbol'],'residue':r['auth_comp_id'],'author_residue':r['auth_seq_id'],
                           'atom':r['auth_atom_id'],'pos':[float(r[x]) for x in ('Cartn_x','Cartn_y','Cartn_z')]}
                          for r,ds in zip(pep,distances) if ds.min()<=4.5]}
    a.atomic(out/'peptide_features.json',features)
    # Same contact atoms/elements; one predeclared rigid perturbation selected only
    # by pocket geometry. Reject severe receptor overlap and empty/all-identical contacts.
    from scipy.spatial.transform import Rotation
    rng=np.random.default_rng(20261023); base=interface; pivot=base.mean(0)
    original=contact_signature(rec_atoms,base); perturb=None
    for attempt in range(1000):
        axis=rng.normal(size=3);axis/=np.linalg.norm(axis)
        rotation=Rotation.from_rotvec(axis*rng.uniform(.15,.45)).as_matrix()
        shift=rng.uniform(-1.5,1.5,3); moved=(base-pivot)@rotation.T+pivot+shift
        nearest=np.linalg.norm(moved[:,None,:]-xyz[None,:,:],axis=2).min(1)
        sig=contact_signature(rec_atoms,moved)
        if (nearest.min()>=1.5 and (nearest<2).mean()<=.10 and np.all(np.abs(moved-center)<=size/2)
            and .75*len(original)<=len(sig)<=1.25*len(original) and sig!=original):
            perturb={'seed':20261023,'attempt':attempt,'rotation':rotation.tolist(),'translation':shift.tolist(),
                     'pivot':pivot.tolist(),'residues':sorted(sig),'minimum_receptor_distance_A':float(nearest.min()),
                     'fraction_under_2A':float((nearest<2).mean()),'coordinates':moved.tolist(),
                     'same_feature_types_and_count':True,'algorithm':'bounded rigid perturbation, geometric rejection only'};break
    if perturb is None:raise ValueError('Could not construct a valid A2 geometric control')
    chemical=peptide_chemical_features(out/'peptide.pdb',rec_atoms)
    a.atomic(out/'chemical_features.json',{'used_for_ranking':False,'true':chemical,'perturbed':[{**f,'pos':((np.array(f['pos'])-pivot)@rotation.T+pivot+shift).tolist()} for f in chemical]})
    a.atomic(out/'perturbed_features.json',perturb)
    from sbvs2.receptor.prepare import prepare_receptor, ReceptorPolicy
    # PDBFixer's unconstrained missing-terminal placement can put OXT close
    # enough to CA to trigger a spurious Meeko proximity bond. Complete only
    # missing terminal oxygen with ideal planar carboxylate geometry, keeping
    # every observed heavy atom fixed. This is recorded, not a hidden deletion.
    last_id=rec[-1]['auth_seq_id'];terminal={r['auth_atom_id']:r for r in rec if r['auth_seq_id']==last_id}
    completed=list(rec);completion=[]
    if 'OXT' not in terminal and all(k in terminal for k in ('N','CA','C','O')):
        point=lambda row:np.array([float(row[k]) for k in ('Cartn_x','Cartn_y','Cartn_z')])
        carbon=point(terminal['C']);v1=point(terminal['CA'])-carbon;v2=point(terminal['O'])-carbon
        direction=-(v1/np.linalg.norm(v1)+v2/np.linalg.norm(v2));direction/=np.linalg.norm(direction)
        pos=carbon+1.25*direction
        added={**terminal['O'],'auth_atom_id':'OXT','label_atom_id':'OXT','occupancy':'1.0','B_iso_or_equiv':'0.0'}
        for key,value in zip(('Cartn_x','Cartn_y','Cartn_z'),pos):added[key]=str(value)
        completed.append(added)
        completion=[{'residue':last_id,'atom':'OXT','rule':'planar trigonal carboxylate, C-OXT 1.25 A; experimental coordinates unchanged','pos':pos.tolist()}]
    write_pdb(completed,out/'receptor_completed.pdb')
    manifest=prepare_receptor(out/'receptor_completed.pdb',ReceptorPolicy(ph=7.4,ph_confirmed=True,minimize=False),out/'prepared')
    manifest['terminal_completion']=completion
    from sbvs2.docking.pdbqt import receptor_to_pdbqt
    receptor_to_pdbqt(out/'prepared/receptor_prepared.pdb',out/'receptor')
    result={'main':audit_cif(source),'references':[audit_cif(p) for p in sorted((root/'evaluation/structures').glob('*.cif'))],
            'target_type':'soluble_PPI','box':box,'preparation':manifest,
            'membrane_depth':'not_applicable','transmembrane_direction':'not_applicable','legacy_target_residues':'not_applicable',
            'reference_crystal_coordinates_used_for_search':False,
            'p53_checks':{r:any(key.split(':')[1:]==[num,name] and value['contact_atom_pairs']>0 for key,value in density.items()) for r,num,name in [('Phe19','19','PHE'),('Trp23','23','TRP'),('Leu26','26','LEU')]}}
    # Verify the receptor coordinates actually used by Meeko stayed in the
    # peptide frame; added hydrogens/OXT are excluded from this comparison.
    original=dict(rec_atoms);prepared_atoms=dict(read_heavy_atoms(out/'prepared/receptor_prepared.pdb'))
    if not set(original)<=set(prepared_atoms):raise ValueError('Observed receptor atom mapping lost')
    frame_delta=np.array([np.array(prepared_atoms[k])-v for k,v in original.items()])
    frame_rms=float(np.sqrt(np.mean(np.sum(frame_delta**2,axis=1))))
    if frame_rms>.01:raise ValueError('Prepared receptor shifted from peptide frame')
    result['observed_heavy_atom_frame_rmsd_A']=frame_rms
    for entry in [result['main']]+result['references']:
        differences=entry['sequence_differences']
        entry['review']={'substitutions_recorded':[r for r in differences if 'tag' not in r.get('details','').lower()],
                         'expression_tag_records':sum('tag' in r.get('details','').lower() for r in differences),
                         'main_receptor':Path(entry['file']).stem=='1YCR',
                         'crystal_ligand_pose_used_for_screening':False,
                         'primary_structure_citation_available':any(r.get('pdbx_database_id_DOI') not in (None,'?','.') for r in entry['citations'])}
        entry['audit_status']='entity, sequence differences, source identity and citation reviewed; expression tags distinguished from substitutions'
    a.atomic(root/'target_audit.json',result)
    return result
