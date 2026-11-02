"""Docking result/pose server contracts: job_id is the authoritative key."""
import json
from pathlib import Path
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'00_workbench'))
from webui import server, workflow


def _sdf(path, smiles='CCO'):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=17) == 0
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)


def _run(tmp_path):
    run = tmp_path/'run'
    job = run/'docking/jobs/DOCKJOB'
    job.mkdir(parents=True)
    _sdf(job/'poses.sdf')
    (run/'docking/results.json').write_text(json.dumps([
        {'id': 'M1', 'job_id': 'DOCKJOB', 'seed': 17, 'status': 'ok',
         'vina_affinity_kcal_mol': -5.2, 'qualification_evidence': {'nested': 'x'*200000}},
    ]))
    (run/'docking/receptor.pdb').write_text(
        'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n')
    return run


def test_docking_pose_resolves_by_job_id(tmp_path, monkeypatch):
    run = _run(tmp_path)
    monkeypatch.setattr(server.workflow, 'run_dir', lambda rid: run)
    status, ctype, body = server.api_project_pose({'id': ['DOCKJOB']}, 'rid', 'docking')
    data = json.loads(body)
    assert status == 200
    assert data['atoms']
    assert 'AutoDock Vina' in data['caption']


def test_docking_rows_use_authoritative_json_for_large_evidence(tmp_path, monkeypatch):
    run = _run(tmp_path)
    monkeypatch.setattr(server.workflow, 'run_dir', lambda rid: run)
    status, ctype, body = server.api_project_rows({'offset': ['0']}, 'rid', 'docking')
    data = json.loads(body)
    assert status == 200
    assert data['total'] == 1
    assert data['rows'][0]['job_id'] == 'DOCKJOB'
