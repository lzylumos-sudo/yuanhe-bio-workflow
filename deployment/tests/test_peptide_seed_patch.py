"""The deployed inference entry point must apply the requested seed first."""
import ast
import hashlib
import json
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

DEPLOYMENT = Path(__file__).resolve().parents[1]
FIXTURE = Path(__file__).with_name('fixtures') / 'peptide2mol' / 'eval.py.txt'


@pytest.fixture
def evaluate(tmp_path):
    patches = DEPLOYMENT / 'vendor_patches'
    row = next(row for row in json.loads((patches / 'manifest.json').read_text())
               if row['vendor'] == 'Peptide2Mol' and row['file'] == 'src/eval.py')
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() == row['before_sha256']
    destination = tmp_path / row['file']
    destination.parent.mkdir(parents=True)
    shutil.copyfile(FIXTURE, destination)
    subprocess.run(['patch', '--batch', '-p1', '-i', str(patches / row['patch'])],
                   cwd=tmp_path, check=True, capture_output=True, text=True)
    assert hashlib.sha256(destination.read_bytes()).hexdigest() == row['after_sha256']
    tree = ast.parse(destination.read_text())
    assert any(isinstance(node, ast.ImportFrom) and node.module == 'lightning'
               and any(alias.name == 'seed_everything' for alias in node.names)
               for node in tree.body)
    function = next(node for node in tree.body
                    if isinstance(node, ast.FunctionDef) and node.name == 'evaluate')
    function.decorator_list = []
    function.returns = None
    for argument in function.args.args:
        argument.annotation = None
    code = compile(ast.Module(body=[function], type_ignores=[]), str(destination), 'exec')
    return code


@pytest.mark.parametrize('seed', [0, 17, 41, None])
def test_seed_precedes_data_model_and_trainer_initialization(evaluate, seed):
    calls = []
    trainer = SimpleNamespace(callback_metrics={},
                              test=lambda **kwargs: calls.append(('test',)))

    def instantiate(config, **kwargs):
        calls.append(('instantiate', config._target_))
        return trainer if config._target_ == 'trainer' else object()

    namespace = {
        'seed_everything': lambda value, *, workers: calls.append(('seed', value, workers)),
        'hydra': SimpleNamespace(utils=SimpleNamespace(instantiate=instantiate)),
        'log': SimpleNamespace(info=lambda *args: None),
        'instantiate_loggers': lambda config: [],
        'load_trusted_checkpoint': lambda *args: calls.append(('checkpoint',)),
    }
    exec(evaluate, namespace)
    config = SimpleNamespace(
        ckpt_path='public-test-checkpoint', seed=seed,
        data=SimpleNamespace(_target_='data'),
        model=SimpleNamespace(_target_='model'),
        trainer=SimpleNamespace(_target_='trainer'),
    )
    config.get = lambda name, default=None: getattr(config, name, default)
    namespace['evaluate'](config)
    expected = [('instantiate', 'data'), ('instantiate', 'model'),
                ('instantiate', 'trainer'), ('checkpoint',), ('test',)]
    if seed is not None:
        expected.insert(0, ('seed', seed, True))
    assert calls == expected
