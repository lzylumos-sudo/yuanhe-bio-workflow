"""Real geometry checks across legacy and receptor-specific screening layouts."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scientific import artifacts as a, peptide_contacts as contacts
from scientific.validation import peptide_evidence
from webui import workflow


def pdb(path, offset=0, name='ALA'):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(''.join(
        f'ATOM  {i:5d}  CA  {name} A{i:4d}    {x+offset:8.3f}{0:8.3f}{0:8.3f}  1.00  0.00           C\n'
        for i, x in enumerate((0, 10, 20), 1)))


def fixture(d, research=True, name='ALA', affinity=-4):
    Chem = pytest.importorskip('rdkit.Chem')
    pdb(d/'prepare/receptor_full.pdb')
    pdb(d/'prepare/peptide_reference.pdb', offset=1)
    receptor = d/(name+'.pdb'); pdb(receptor, name=name)
    sha = a.digest(receptor)
    folder = d/'screening'/sha[:16] if research else d/'screening'
    if research:
        target = folder/'docking_input/receptor.pdb'
    else:
        target = d/'docking/receptor.pdb'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(receptor.read_bytes())
    sdf = folder/'selected_poses/pose.sdf'; sdf.parent.mkdir(parents=True, exist_ok=True)
    mol = Chem.MolFromSmiles('C'); conf = Chem.Conformer(1)
    conf.SetAtomPosition(0, (1., 0., 0.)); mol.AddConformer(conf)
    with Chem.SDWriter(str(sdf)) as writer:
        writer.write(mol)
    row = dict(compound_id='same_compound', pose_id='same_pose', sdf=str(sdf.relative_to(d)))
    if research:
        row['research_receptor'] = {'sha256':sha, 'path':str(receptor)}
    a.atomic(folder/'candidates.json', [row])
    a.atomic(folder/'workspace/runs/screen/stages/selection/summary.json',
             {'primary':[{'compound_id':row['compound_id'], 'raw':{'affinity':affinity}}]})
    return row, folder


def test_legacy_layout(tmp_path):
    row, folder = fixture(tmp_path, research=False)
    result = contacts.calculate(tmp_path, {}, {}, tmp_path/'result')
    assert result['candidates'][0]['recall'] > 0
    assert result['controls']['self']['jaccard'] == 1
    assert result['controls']['translated_100A_beyond_receptor']['jaccard'] == 0


def test_automatic_screening_before_aggregate_exists(tmp_path):
    row, folder = fixture(tmp_path)
    assert not (tmp_path/'screening/candidates.json').exists()
    cfg = {'research':{'profile':'smoke'}, '_screen_receptor':row['research_receptor']}
    result = contacts.calculate(tmp_path, cfg, {}, folder/'peptide_evidence')
    assert result['candidates'][0]['affinity_kcal_mol'] == -4
    assert (folder/'peptide_evidence/peptide_evidence.json').is_file()


def test_standalone_operation_uses_receptor_workspace(tmp_path, monkeypatch):
    row, folder = fixture(tmp_path)
    a.atomic(tmp_path/'screening/candidates.json', [row])
    checked = []
    monkeypatch.setattr(workflow, 'require', lambda d, stage: checked.append(stage))
    result = peptide_evidence(tmp_path, {'research':{'profile':'smoke'}}, {}, tmp_path/'operation')
    assert checked == ['screening', 'prepare']
    assert result['candidates'][0]['receptor_sha256'] == row['research_receptor']['sha256']
    assert result['receptor_results'][0]['result']['candidates'][0]['affinity_kcal_mol'] == -4


def test_receptors_with_same_candidate_ids_keep_separate_scores(tmp_path):
    row1, _ = fixture(tmp_path, affinity=-4)
    # A harmless record after the atoms changes identity without changing geometry.
    source = tmp_path/'second.pdb'; source.write_bytes((tmp_path/'ALA.pdb').read_bytes()+b'END\n')
    sha2 = a.digest(source); folder2 = tmp_path/'screening'/sha2[:16]
    row2 = {**row1, 'research_receptor':{'sha256':sha2, 'path':str(source)}}
    copied = folder2/'docking_input/receptor.pdb'; copied.parent.mkdir(parents=True)
    copied.write_bytes(source.read_bytes())
    a.atomic(folder2/'workspace/runs/screen/stages/selection/summary.json',
             {'primary':[{'compound_id':row2['compound_id'], 'raw':{'affinity':-8}}]})
    a.atomic(tmp_path/'screening/candidates.json', [row1, row2])
    result = contacts.calculate(tmp_path, {'research':True}, {}, tmp_path/'result')
    assert [r['affinity_kcal_mol'] for r in result['candidates']] == [-4, -8]
    assert len({r['receptor_sha256'] for r in result['candidates']}) == 2
    for rec in result['receptor_results']:
        assert a.digest(tmp_path/'result'/rec['artifact']) == rec['sha256']


def test_changed_receptor_is_rejected(tmp_path):
    row, folder = fixture(tmp_path)
    a.atomic(tmp_path/'screening/candidates.json', [row])
    (folder/'docking_input/receptor.pdb').write_text('changed')
    with pytest.raises(ValueError, match='checksum'):
        contacts.calculate(tmp_path, {'research':True}, {}, tmp_path/'result')


def test_research_does_not_fall_back_to_legacy_on_missing_identity(tmp_path):
    fixture(tmp_path, research=False)
    with pytest.raises(ValueError, match='receptor SHA-256'):
        contacts.calculate(tmp_path, {'research':True}, {}, tmp_path/'result')


def test_empty_research_pool(tmp_path):
    a.atomic(tmp_path/'screening/candidates.json', [])
    result = contacts.calculate(tmp_path, {'research':True}, {}, tmp_path/'result')
    assert result['status'] == 'no_candidates' and result['candidates'] == []


def test_candidate_path_escape_is_rejected(tmp_path):
    row, folder = fixture(tmp_path)
    row['sdf'] = '../outside.sdf'
    a.atomic(tmp_path/'screening/candidates.json', [row])
    with pytest.raises(ValueError, match='artifact|escapes'):
        contacts.calculate(tmp_path, {'research':True}, {}, tmp_path/'result')
