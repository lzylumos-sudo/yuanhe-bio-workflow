import os,json,csv,hashlib,shutil,subprocess,sys
from pathlib import Path
import numpy as np
from rdkit import Chem
from sbvs.paths import ROOT,read_json,write_json,sha
from sbvs.export import read_rows,handoff
from sbvs.preparation.ligands import identity
results={};run=ROOT/read_json(ROOT/'reports/class20_result.json')['run'];rows=read_rows(run/'poses.csv')
assert len(rows)>=20
for row in rows:
 m=next(iter(Chem.SDMolSupplier(str(run/row['pose_file']),removeHs=False)));prep=read_json((run/row['pose_file']).parent/'preparation.json');assert identity(m)==prep['identity'];assert np.isfinite(m.GetConformer().GetPositions()).all()
results['vina_sdf_chemistry']={'status':'passed','pose_count':len(rows)}
ref=next(iter(Chem.SDMolSupplier(str(ROOT/'data/reference/FOL_crystal.sdf'),removeHs=False)));mapping=read_json(ROOT/'data/reference/FOL_mapping.json');xyz={}
for l in (ROOT/'data/raw/1DRF.pdb').read_text().splitlines():
 if l.startswith('HETATM') and l[17:20]=='FOL':xyz[l[12:16].strip()]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
for name,i in mapping['ccd_atom_to_index'].items():assert np.allclose(ref.GetConformer().GetPositions()[i],xyz[name],atol=1e-6)
box=read_json(ROOT/'data/prepared/demo_1drf/box.json');assert np.all(abs(ref.GetConformer().GetPositions()-np.array(box['center']))<=np.array(box['size'])/2)
results['reference_coordinate_and_box']={'status':'passed','heavy_atoms':len(xyz)}
try:handoff(run,{'pose_id':rows[0]['pose_id'],'human_confirmed':False})
except ValueError:pass
else:raise AssertionError('Unconfirmed handoff was accepted')
# Explicitly labelled test fixture; this is not a scientific selection by the user.
out=handoff(run,{'pose_id':rows[0]['pose_id'],'human_confirmed':True,'reviewer':'AUTOMATED_ACCEPTANCE_FIXTURE_NOT_USER_SELECTION','acknowledge_warnings':True,'test_only':True})
assert sha(out/'ligand_selected.sdf')==sha(run/rows[0]['pose_file']);assert all((out/f).exists() for f in ['receptor_prepared.pdb','ligand_selected.sdf','complex_preview.pdb','box.json','atom_mapping.json','selection_manifest.json','README.md'])
results['md_handoff']={'status':'passed','test_only':True,'path':str(out.relative_to(ROOT)),'no_md_performed':True}
write_json(ROOT/'reports/artifact_acceptance.json',results);print(json.dumps(results,indent=2))
