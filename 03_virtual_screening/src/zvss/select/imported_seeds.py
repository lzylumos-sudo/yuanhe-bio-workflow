"""Verify real upstream Vina seeds against the exact selected pose, without redocking."""
from pathlib import Path
from ..ids import sha256_file
from ..util import read_json

def evaluate(ctx,entry,proof,planned,supporting_required):
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    from ..chem.rdio import iter_sdf
    path=Path(proof['path'])
    if sha256_file(str(path))!=proof['sha256']:raise ValueError('Imported seed evidence changed')
    data=read_json(str(path));pose=ctx.db.one('SELECT coords_path, compound_id FROM pose WHERE pose_id=?',(entry['pose_id'],))
    if not pose or not pose.get('coords_path'):return None
    reference=next(iter_sdf(pose['coords_path']),None)
    if reference is None:return None
    mol=Chem.RemoveHs(reference[0]);canon=Chem.MolToSmiles(mol,isomericSmiles=True)
    attempted=set();succeeded=set();support=set();matching_reference=False
    for record in data['jobs']:
        if record['smiles']!=canon:continue
        attempted.add(record['seed'])
        if record['status']!='ok':continue
        source=Path(record['poses'])
        if sha256_file(str(source))!=record['poses_sha256']:raise ValueError('Imported Vina poses changed')
        succeeded.add(record['seed'])
        for candidate,_ in iter_sdf(str(source)):
            candidate=Chem.RemoveHs(candidate)
            if Chem.MolToSmiles(candidate,isomericSmiles=True)!=canon:raise ValueError('Vina seed chemical identity differs')
            rms=float(rdMolAlign.CalcRMS(candidate,mol,maxMatches=10000))
            if rms<1e-4:matching_reference=True
            if rms<=data['cluster_rmsd_A']:support.add(record['seed'])
    observed={'seeds_attempted':len(attempted),'seeds_succeeded':len(succeeded),'supporting_seeds_same_cluster':len(support),
              'selected_pose_found':matching_reference,'rmsd_threshold_A':data['cluster_rmsd_A'],'receptor_sha256':data['receptor_sha256']}
    return {'passed':matching_reference and len(succeeded)>=planned and len(support)>=supporting_required,'observed':observed}
