"""Load the official SA scorer from either the PyPI or Conda RDKit layout."""
from functools import lru_cache
import importlib
import importlib.util
from pathlib import Path


@lru_cache(maxsize=1)
def load_sascorer():
    try:
        scorer = importlib.import_module('rdkit.Contrib.SA_Score.sascorer')
    except ModuleNotFoundError as exc:
        # Do not disguise a broken dependency inside an otherwise present scorer.
        if exc.name not in {'rdkit.Contrib', 'rdkit.Contrib.SA_Score',
                            'rdkit.Contrib.SA_Score.sascorer'}:
            raise
        from rdkit import RDConfig
        source = Path(RDConfig.RDContribDir) / 'SA_Score' / 'sascorer.py'
        if not source.is_file():
            raise RuntimeError('RDKit SA_Score is missing from ' + str(source) +
                               '; install the complete RDKit Contrib data in this environment') from exc
        spec = importlib.util.spec_from_file_location('_yuanhe_rdkit_sascorer', source)
        scorer = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(scorer)
    fragments = Path(scorer.__file__).with_name('fpscores.pkl.gz')
    if not fragments.is_file():
        raise RuntimeError('RDKit SA_Score fragment data is missing: ' + str(fragments))
    return scorer
