"""Offline pipeline regression using archived 1DRF output, without new docking."""
from pathlib import Path
import json,os,shutil,subprocess,sys
import pytest

PROJECT=Path(__file__).resolve().parents[1]
DATA=Path(__file__).parent/'data/1drf'

@pytest.fixture(scope='module')
def workspace(tmp_path_factory):
 root=tmp_path_factory.mktemp('portable-screening');inputs=root/'inputs';inputs.mkdir()
 shutil.copy2(DATA/'receptor_prepared.pdb',inputs/'receptor.pdb')
 shutil.copy2(DATA/'box.json',inputs/'box.json')
 shutil.copy2(DATA/'poses.sdf',inputs/'poses.sdf')
 structures=root/'structures';structures.mkdir();shutil.copy2(DATA/'receptor_prepared.pdb',structures/'1DRF.pdb')
 host=root/'host.json';host.write_text(json.dumps({'paths':{'structures_dir':str(structures)}}))
 for run_id in ('e2e1','e2e2'):
  args=[sys.executable,'-m','zvss.cli','--root',str(root),'run','--input',str(inputs),'--target',str(PROJECT/'configs/targets/generic.yaml'),'--budget',str(PROJECT/'configs/budgets/quick.yaml'),'--host',str(host),'--run-id',run_id]
  result=subprocess.run(args,capture_output=True,text=True,env={**os.environ,'PYTHONPATH':str(PROJECT/'src')},timeout=180)
  assert result.returncode==0,result.stdout+result.stderr
 return root

def read(root,run,path):return json.loads((root/'runs'/run/path).read_text())

def test_pipeline_imports_every_pose_and_preserves_the_real_exclusion(workspace):
 summary=read(workspace,'e2e1','stages/import/summary.json')
 assert summary['poses_imported']==3 and summary['poses_failed']==0
 selection=read(workspace,'e2e1','stages/selection/summary.json')
 assert selection['selection_version'].startswith('SEL')
 assert not selection['primary']
 assert selection['candidates_considered']==1
 assert selection['excluded']==1 and selection['eligible']==0
 assert selection['exclusions'][0]['reason_code']=='excessive_strain'
 assert not selection['preliminary'] and selection['formal_qualified']==0
 assert (workspace/'runs/e2e1/reports/report.md').is_file()

def test_second_run_reproduces_selection_and_candidates(workspace):
 selections=[read(workspace,r,'stages/selection/summary.json')['selection_version'] for r in ('e2e1','e2e2')]
 assert selections[0]==selections[1]
 assert (workspace/'runs/e2e1/reports/candidates.csv').read_bytes()==(workspace/'runs/e2e2/reports/candidates.csv').read_bytes()

def test_zero_formal_does_not_silently_create_md_handoff(workspace):
 summary=read(workspace,'e2e1','stages/md/md_summary.json')
 assert summary['status']=='blocked_no_formally_qualified_candidates'
 assert summary['count']==0 and summary['packages']==[]
