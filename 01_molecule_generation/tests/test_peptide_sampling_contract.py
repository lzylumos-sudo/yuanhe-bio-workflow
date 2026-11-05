"""Execute the actual hash-pinned patched loop with synthetic model outputs, no GPU."""
import ast
import hashlib
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace
import traceback

import pytest
from src.adapters import peptide2mol

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def source():
    path = Path(os.environ.get('PEPTIDE2MOL_TEST_SOURCE', str(
        ROOT / '01_molecule_generation/vendor/Peptide2Mol/src/models/components/moldiff_test/model.py')))
    if not path.is_file():
        pytest.skip('Fetch the pinned CPU loop fixture with deployment/fetch_peptide_test_fixture.py')
    row = json.loads((ROOT / 'deployment/vendor_patches/manifest.json').read_text())[0]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == row['after_sha256']
    return path


class AttrDict(dict):
    __getattr__ = dict.__getitem__


class Tensor:
    def cpu(self): return self
    def numpy(self): return self


def exercise(source, tmp_path, mode, *, budget=3, batch_size=1):
    tree = ast.parse(source.read_text())
    method = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == 'sample')
    code = compile(ast.fix_missing_locations(ast.Module(body=[method], type_ignores=[])), str(source), 'exec')
    class Fake:
        sample_batch_size = batch_size
        num_mols = budget
        sample_max = 20
        guidance = None
        log_dir = str(tmp_path / 'model')
        pdb_dir = str(tmp_path)
        calls = 0
        counts = []
        def make_mydata_placeholder(self, **kw):
            return dict(batch_node=Tensor(), halfedge_index=Tensor(), batch_halfedge=Tensor(),
                        ref_data=SimpleNamespace(name=['test']), n_nodes_list=[])
        def sample_model(self, **kw):
            self.calls += 1; self.counts.append(kw['n_graphs'])
            if mode == 'device_error': raise RuntimeError('CUDA fixture failure')
            return {}
        def seperate_outputs(self, outputs, n_graphs, *args):
            if mode == 'separation_error': raise ValueError('separation fixture')
            if mode == 'empty': return []
            return [dict(pred=[None, None, None], halfedge_index=None)] * n_graphs
        def seperate_outputs_origin(self, *args): raise ValueError('origin fixture')
        def decode_output(self, **kw):
            if mode == 'decode_error': raise ValueError('decode fixture')
            return {}
        def reconstruct_from_generated_with_edges(self, mol):
            if mode == 'reconstruction_error': raise ValueError('reconstruction fixture')
            return 'C.C' if mode == 'disconnected' else 'CCO'
    fake = Fake(); fake.model = SimpleNamespace(sample=fake.sample_model)
    target = Path(fake.log_dir + '_SDF'); target.mkdir()
    if mode != 'origin_error': (target / 'test_origin.pdb').write_text('fixture')
    chem = SimpleNamespace(MolToSmiles=lambda mol: mol,
                           MolToMolFile=lambda mol, p: Path(p).write_text(mol))
    namespace = dict(EasyDict=AttrDict, os=os, json=json, traceback=traceback, Chem=chem)
    exec(code, namespace)
    if mode == 'device_error':
        with pytest.raises(RuntimeError, match='CUDA fixture'):
            namespace['sample'](fake, SimpleNamespace(pos=SimpleNamespace(device='cpu')))
    else:
        namespace['sample'](fake, SimpleNamespace(pos=SimpleNamespace(device='cpu')))
    return fake, target


@pytest.mark.parametrize('mode', ['empty', 'origin_error', 'separation_error', 'decode_error',
                                  'reconstruction_error', 'disconnected', 'valid'])
def test_exact_attempt_budget_and_complete_failure_accounting(source, tmp_path, mode):
    fake, target = exercise(source, tmp_path, mode)
    audit = peptide2mol.sampling_audit(target / 'attempt_budget.json', 3)
    assert fake.calls == 3 and sum(fake.counts) == 3
    assert audit['valid'] == (3 if mode == 'valid' else 0)
    assert audit['valid'] + audit['failed'] == 3
    assert len(list(target.glob('*.sdf'))) == audit['valid']


def test_partial_last_batch_does_not_overshoot_budget(source, tmp_path):
    fake, target = exercise(source, tmp_path, 'valid', budget=3, batch_size=2)
    assert fake.counts == [2, 1]
    assert peptide2mol.sampling_audit(target / 'attempt_budget.json', 3)['valid'] == 3


def test_device_failure_is_not_a_completed_sampling_run(source, tmp_path):
    fake, target = exercise(source, tmp_path, 'device_error')
    audit = json.loads((target / 'attempt_budget.json').read_text())
    assert fake.calls == 1 and audit['status'] == 'failed'
    with pytest.raises(ValueError): peptide2mol.sampling_audit(target / 'attempt_budget.json', 3)


def test_unpatched_model_rejected_before_any_prepare_or_inference(tmp_path, monkeypatch):
    repo = tmp_path / 'vendor'; (repo / 'notebooks').mkdir(parents=True); (repo / 'src').mkdir()
    (repo / 'notebooks/deal_with_mol_test_from_pdb_pal_5A.py').write_text('')
    (repo / 'src/eval.py').write_text('')
    checkpoint = tmp_path / 'model'; checkpoint.write_bytes(b'fixture')
    pocket = tmp_path / 'poc.pdb'; pocket.write_text('')
    monkeypatch.setattr(peptide2mol.platform, 'system', lambda: 'Linux')
    monkeypatch.setattr(peptide2mol, '_run', lambda *a: pytest.fail('must fail before launching model'))
    with pytest.raises(RuntimeError, match='bounded-sampling patch'):
        peptide2mol.run(pocket_pdb=pocket, checkpoint=checkpoint,
                       checkpoint_sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
                       output_dir=tmp_path / 'run', budget=2, seed=1, repo=repo, executable=sys.executable)
    assert json.loads((tmp_path / 'run/manifest.json').read_text())['commands'] == []


def test_verified_patch_is_required_by_readiness(source, tmp_path, monkeypatch):
    repo = tmp_path / 'vendor'; (repo / 'notebooks').mkdir(parents=True); (repo / 'src').mkdir()
    (repo / 'notebooks/deal_with_mol_test_from_pdb_pal_5A.py').write_text('')
    (repo / 'src/eval.py').write_text('')
    path = repo / 'src/models/components/moldiff_test/model.py'
    path.parent.mkdir(parents=True); path.write_bytes(source.read_bytes())
    monkeypatch.setattr(peptide2mol.platform, 'system', lambda: 'Linux')
    assert peptide2mol.platform_readiness(repo, sys.executable)['ready']
    path.write_text(path.read_text() + '# changed\n')
    assert not peptide2mol.platform_readiness(repo, sys.executable)['ready']
