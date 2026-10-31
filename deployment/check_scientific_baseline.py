#!/usr/bin/env python3
"""Strict protected-source hashes with explicitly scoped, versioned revisions."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASELINE=json.loads((ROOT/'docs/architecture/tracked_file_audit.json').read_text())
PREFIXES=('01_molecule_generation/src/adapters/peptide2mol.py','01_molecule_generation/src/adapters/reinvent4.py',
          '01_molecule_generation/src/chemistry.py','01_molecule_generation/configs/',
          '02_docking/src/sbvs2/docking/','02_docking/src/sbvs2/ligands/','02_docking/src/sbvs2/receptor/',
          '03_virtual_screening/src/zvss/membrane/','03_virtual_screening/src/zvss/stages.py','03_virtual_screening/configs/',
          '04_interaction_analysis/src/', '05_molecular_dynamics/modules/md/',
          '05_molecular_dynamics/modules/parameterization/','05_molecular_dynamics/modules/analysis/',
          '05_molecular_dynamics/config/','00_workbench/configs/','data/','archive/')

def check():
    files=[f for f in BASELINE['files'] if f['path'].startswith(PREFIXES)]
    revision_file=ROOT/'docs/architecture/presentation_revisions.json'
    revisions=json.loads(revision_file.read_text()) if revision_file.is_file() else {'schema':'yuanhe.presentation-revisions/1','revisions':[]}
    if revisions['schema']!='yuanhe.presentation-revisions/1':raise ValueError('Unknown presentation revision schema')
    allowed={'04_interaction_analysis/src/zlis/diagram2d.py'}
    original={f['path']:f['sha256'] for f in files};revised={}
    for revision in revisions['revisions']:
        path=revision['path']
        if path not in allowed or path in revised or revision['scope']!='presentation_only':raise ValueError('Invalid presentation revision: '+path)
        if revision['baseline_sha256']!=original.get(path):raise ValueError('Presentation revision does not match preserved baseline: '+path)
        revised[path]=revision['revised_sha256']
    presentation_count=len(revised)
    optimizer_file=ROOT/'docs/architecture/optimizer_revisions.json'
    optimizers=json.loads(optimizer_file.read_text()) if optimizer_file.is_file() else {'schema':'yuanhe.optimizer-revisions/1','revisions':[]}
    if optimizers['schema']!='yuanhe.optimizer-revisions/1':raise ValueError('Unknown optimizer revision schema')
    for revision in optimizers['revisions']:
        path=revision['path']
        if path!='05_molecular_dynamics/modules/parameterization/hf_resp.py' or path in revised or revision['scope']!='optimizer_numerics' or revision['request_schema']!='hf-resp/2':
            raise ValueError('Invalid optimizer revision: '+path)
        if revision['baseline_sha256']!=original.get(path):raise ValueError('Optimizer revision does not match preserved baseline: '+path)
        revised[path]=revision['revised_sha256']
    errors=[f['path'] for f in files if not (ROOT/f['path']).is_file() or hashlib.sha256((ROOT/f['path']).read_bytes()).hexdigest()!=revised.get(f['path'],f['sha256'])]
    if errors:raise ValueError('Scientific baseline drift: '+', '.join(errors))
    return {'strictly_equal_files':len(files)-len(revised),'versioned_presentation_files':presentation_count,
            'versioned_optimizer_files':len(revised)-presentation_count,'baseline_commit':BASELINE['baseline_commit']}
if __name__=='__main__':print(json.dumps(check()))
