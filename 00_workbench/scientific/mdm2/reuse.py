"""Rebind anonymous IDs only when the actual Vina input bytes are identical."""
import hashlib
import shutil
from pathlib import Path
from rdkit import Chem
from scientific import artifacts as a
from .common import sha,now


def reuse(root,source_root):
    from sbvs2.docking.pdbqt import ligand_to_pdbqt
    root,source_root=Path(root),Path(source_root)
    fp=root/'inputs/preregistered_protocol.json';sp=source_root/'inputs/preregistered_protocol.json'
    for path in (fp,sp):
        if a.digest(path)!=path.with_suffix('.sha256').read_text().strip():raise ValueError('Frozen manifest changed')
    target,source=a.read(fp),a.read(sp)
    keys=['receptor_sha256','engine_sha256','protocol_hash','box_hash','cpu_per_worker']
    if any(target[k]!=source[k] for k in keys) or target['seeds']!=source['seeds']:raise ValueError('Docking identity differs; no reuse allowed')
    for rel,digest in target['files'].items():
        if a.digest(root/rel)!=digest:raise ValueError('Target input changed: '+rel)
    candidates={}
    for path in sorted((source_root/'inputs/jobs').glob('*/*/attempt_*/job_record.json')):
        record=a.read(path)
        if record['status']!='ok' or record.get('output_sha256')!=a.digest(path.parent/'poses.sdf'):continue
        if (record['receptor_sha256']!=target['receptor_sha256'] or record['protocol_hash']!=target['protocol_hash'] or
            record['box_hash']!=target['box_hash'] or record['seed'] not in target['seeds'] or record['cpu']!=target['cpu_per_worker']):continue
        candidates.setdefault(a.digest(path.parent/'ligand.pdbqt'),path)
    imported=[]
    for ligand in target['ligands']:
        if a.digest(ligand['path'])!=ligand['sha256']:raise ValueError('Frozen ligand changed')
        mol=Chem.SDMolSupplier(ligand['path'],removeHs=False)[0]
        mol=Chem.AddHs(mol,addCoords=True);text,_=ligand_to_pdbqt(mol,0)
        pdbqt_sha=hashlib.sha256(text.encode()).hexdigest();path=candidates.get(pdbqt_sha)
        if path is None:continue
        identity={k:target[k] for k in keys};identity.update(ligand_sha256=ligand['sha256'],seed=20261023,preparation=target['state_rule'])
        cache_key=sha(identity);dest=root/'inputs/jobs'/ligand['parent_id']/cache_key/'attempt_1'
        if dest.exists():continue
        dest.parent.mkdir(parents=True,exist_ok=True);shutil.copytree(path.parent,dest)
        # Preserve original raw command, timestamps, IDs and hashes untouched.
        # The binding records why an anonymous-name-only difference is reusable.
        binding={'reused_at':now(),'source_job_record':str(path),'source_job_record_sha256':a.digest(path),
            'source_manifest_sha256':a.digest(sp),'target_manifest_sha256':a.digest(fp),
            'target_parent_id':ligand['parent_id'],'target_state_id':ligand['state_id'],
            'target_sdf_sha256':ligand['sha256'],'identical_vina_pdbqt_sha256':pdbqt_sha,
            'criteria':'exact Vina ligand bytes + receptor + engine + box + protocol + seed + CPU; no score-based selection',
            'original_raw_record_preserved':True}
        a.atomic(dest/'reuse_evidence.json',binding);imported.append(binding)
    result={'source_root':str(source_root),'target_root':str(root),'imported_jobs':len(imported),'jobs':imported}
    a.atomic(root/'reuse_manifest.json',result);return {'imported_jobs':len(imported)}
