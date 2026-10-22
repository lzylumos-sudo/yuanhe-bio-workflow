#!/usr/bin/env python3
"""Strict equality for protected scientific implementations and default configurations."""
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
    errors=[f['path'] for f in files if not (ROOT/f['path']).is_file() or hashlib.sha256((ROOT/f['path']).read_bytes()).hexdigest()!=f['sha256']]
    if errors:raise ValueError('Scientific baseline drift: '+', '.join(errors))
    return {'strictly_equal_files':len(files),'baseline_commit':BASELINE['baseline_commit']}
if __name__=='__main__':print(json.dumps(check()))
