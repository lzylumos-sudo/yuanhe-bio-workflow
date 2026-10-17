from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from zvss.stages import stage_redock
from zvss.strain.strain import _strain

def test_disabled_redock_never_probes_or_launches_engine():
 for disabled,limit in [(['redock'],20),([],0)]:
  ctx=SimpleNamespace(config=SimpleNamespace(budget=SimpleNamespace(not_run=disabled,redock_max_parents=limit)))
  with patch('zvss.redock.gnina.run_redock',side_effect=AssertionError('must not run')):
   assert stage_redock(ctx)['status']=='not_run'

def test_strain_preserves_input_and_records_convergence():
 mol=Chem.AddHs(Chem.MolFromSmiles('CCO'));AllChem.EmbedMolecule(mol,randomSeed=17)
 before=np.array(mol.GetConformer().GetPositions())
 result=_strain(mol)
 assert np.array_equal(before,mol.GetConformer().GetPositions())
 assert result['hydrogen_minimization_converged'] and result['reference_minimization_converged']
 assert result['strain_kcal_mol']>=-0.01
 assert result['protocol']=='heavy_atoms_fixed_H_relaxed_minus_local_minimum'


def test_strain_metadata_does_not_falsely_claim_nonconvergence():
 from zvss.strain.strain import run
 from unittest.mock import Mock
 mol=Chem.AddHs(Chem.MolFromSmiles('CCO'));AllChem.EmbedMolecule(mol,randomSeed=17)
 ctx=SimpleNamespace(config=SimpleNamespace(budget=SimpleNamespace()),db=Mock())
 pose={'compound_id':'c','pose_id':'p'}
 with patch('zvss.data.pose_rows',return_value=[pose]), patch('zvss.data.pose_mol',return_value=mol), patch('zvss.scoring.aggregate.representative_pose',return_value=pose):
  run(ctx)
 for call in ctx.db.record_assessment.call_args_list:
  if call.args[2]=='strain' and call.args[3] in {'force_field','pose_id','protocol'}:
   assert call.args[5].value=='pass' and call.args[6] is None
