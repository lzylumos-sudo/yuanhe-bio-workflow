"""The same official SA scorer must work in both supported RDKit distributions."""
import importlib
from pathlib import Path
from types import SimpleNamespace

import pytest
from yuanhe_preflight.sa_score import load_sascorer


@pytest.fixture(autouse=True)
def clear_cache():
    load_sascorer.cache_clear()
    yield
    load_sascorer.cache_clear()


def test_installed_scorer_loads_real_fragment_data():
    from rdkit import Chem
    scorer = load_sascorer()
    assert 1 <= scorer.calculateScore(Chem.MolFromSmiles('CCO')) <= 10


def test_conda_contrib_layout_without_python_package(tmp_path, monkeypatch):
    from rdkit import RDConfig
    source = tmp_path / 'Contrib/SA_Score/sascorer.py'; source.parent.mkdir(parents=True)
    source.write_text('def calculateScore(mol): return 2.5\n')
    source.with_name('fpscores.pkl.gz').write_bytes(b'synthetic fragment fixture')
    monkeypatch.setattr(RDConfig, 'RDContribDir', str(source.parent.parent))
    original = importlib.import_module
    def missing(name, *a, **kw):
        if name == 'rdkit.Contrib.SA_Score.sascorer':
            raise ModuleNotFoundError("No module named 'rdkit.Contrib'", name='rdkit.Contrib')
        return original(name, *a, **kw)
    monkeypatch.setattr(importlib, 'import_module', missing)
    scorer = load_sascorer()
    assert Path(scorer.__file__) == source
    assert scorer.calculateScore(None) == 2.5


def test_missing_fragment_data_fails_explicitly(tmp_path, monkeypatch):
    module = SimpleNamespace(__file__=str(tmp_path / 'sascorer.py'))
    monkeypatch.setattr(importlib, 'import_module', lambda name: module)
    with pytest.raises(RuntimeError, match='fragment data is missing'):
        load_sascorer()


def test_internal_missing_dependency_is_not_disguised(monkeypatch):
    def broken(name): raise ModuleNotFoundError('missing dependency', name='unexpected_dependency')
    monkeypatch.setattr(importlib, 'import_module', broken)
    with pytest.raises(ModuleNotFoundError) as error:
        load_sascorer()
    assert error.value.name == 'unexpected_dependency'
