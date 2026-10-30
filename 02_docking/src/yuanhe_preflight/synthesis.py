"""Ertl SA (RDKit implementation) and Coley SCScore, with real external weights.

SCScore equations/fingerprints follow connorcoley/scscore (MIT). Only the
official 1024-bit boolean model is accepted; never unpickle downloaded files.
Both scores increase with difficulty, but are thresholded independently.
"""
import gzip
import importlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator

from .common import file_hash


def _load_sa_scorer():
    """Use the installed RDKit scorer in either wheel or Conda layout."""
    try:
        return importlib.import_module('rdkit.Contrib.SA_Score.sascorer')
    except ModuleNotFoundError as exc:
        if exc.name not in {'rdkit.Contrib', 'rdkit.Contrib.SA_Score',
                            'rdkit.Contrib.SA_Score.sascorer'}:
            raise
    from rdkit import RDConfig
    source=Path(RDConfig.RDContribDir)/'SA_Score'/'sascorer.py'
    if not source.is_file() or not source.with_name('fpscores.pkl.gz').is_file():
        raise ModuleNotFoundError('Installed RDKit SA_Score implementation or fragment data is missing')
    spec=importlib.util.spec_from_file_location('yuanhe_preflight._rdkit_sascorer',source)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


sascorer = _load_sa_scorer()


class SCScore:
    def __init__(self, path, expected_sha256):
        path = Path(path)
        if not expected_sha256 or file_hash(path) != expected_sha256:
            raise ValueError('SCScore weight SHA256 mismatch or absent')
        with gzip.open(path, 'rt') as handle:
            self.weights = [np.asarray(x, dtype=np.float64) for x in json.load(handle)]
        dims = [1024, 300, 300, 300, 300, 300, 1]
        if len(self.weights) != 12:
            raise ValueError('Expected official 1024bool SCScore network')
        for i, (n, m) in enumerate(zip(dims, dims[1:])):
            if self.weights[2*i].shape != (n, m) or self.weights[2*i+1].shape != (m,):
                raise ValueError('SCScore network shape mismatch')
        if not all(np.isfinite(x).all() for x in self.weights):
            raise ValueError('Nonfinite model weights')
        self.sha256 = expected_sha256
        self.fp = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024, includeChirality=True)

    def predict(self, molecules):
        if not molecules:
            return []
        x = np.asarray([self.fp.GetFingerprintAsNumPy(m) for m in molecules], dtype=np.float64)
        for i in range(0, len(self.weights), 2):
            x = x @ self.weights[i] + self.weights[i+1]
            if i < len(self.weights)-2:
                x = np.maximum(x, 0)
        return (1 + 4 / (1 + np.exp(-np.clip(x[:, 0], -700, 700)))).tolist()


def score_batch(rows, config):
    error = None
    try:
        model = SCScore(config['scscore_path'], config['scscore_sha256'])
    except Exception as exc:
        model = None
        error = type(exc).__name__ + ': ' + str(exc)
    mols = [Chem.MolFromSmiles(r['smiles']) for r in rows]
    scores = model.predict(mols) if model else [None] * len(rows)
    records = []
    for row, mol, sc in zip(rows, mols, scores):
        sa = float(sascorer.calculateScore(mol))
        failed = not math.isfinite(sa) or sc is None
        hard = sa > config['sa_max'] or (sc is not None and sc > config['sc_max'])
        records.append({
            'molecule_id': row['compound_id'], 'sa_score': sa, 'sc_score': sc,
            'score_version': {'sa': 'Ertl-2009/RDKit-' + rdBase.rdkitVersion,
                              'sc': 'Coley-2018/1024bool', 'sc_weights_sha256': config['scscore_sha256'],
                              'sa_fragments_sha256': file_hash(Path(sascorer.__file__).with_name('fpscores.pkl.gz'))},
            'synthesis_difficulty': 'UNKNOWN' if failed else 'HIGH' if hard else 'WITHIN_THRESHOLDS',
            'screening_status': 'REVIEW' if failed or hard else 'PASS',
            'prediction_status': 'ERROR' if failed else 'SUCCESS',
            'rejection_reason': error or ('independent_synthesis_threshold_exceeded' if hard else ''),
            'thresholds': {'sa_max': config['sa_max'], 'sc_max': config['sc_max']},
            'interpretation': 'SA 1..10 and SC 1..5: larger is harder; neither proves a complete synthesis route'
        })
    return records
