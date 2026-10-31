"""Optimizer boundary and cache contracts; no quantum calculation is mocked as success."""
import json
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from modules.parameterization import hf_resp


class WorkerBoundary(RuntimeError):
    pass


@pytest.fixture
def recorded_request(tmp_path, monkeypatch):
    mol = Chem.AddHs(Chem.MolFromSmiles('CO'))
    conformer = Chem.Conformer(mol.GetNumAtoms())
    points = [(0, 0, 0), (1.42, 0, 0), (-.36, 1.03, 0),
              (-.36, -.51, .89), (-.36, -.51, -.89), (1.75, .9, 0)]
    for index, point in enumerate(points):
        conformer.SetAtomPosition(index, point)
    mol.AddConformer(conformer)
    monkeypatch.setattr(hf_resp.subprocess, 'run', Mock(side_effect=WorkerBoundary))
    with pytest.raises(WorkerBoundary):
        hf_resp.compute(mol, tmp_path)
    return mol, json.loads((tmp_path/'request.json').read_text())


def test_request_records_initial_hessian_and_unchanged_protocol(recorded_request):
    _, request = recorded_request
    assert request['schema'] == 'hf-resp/4'
    assert request['optimizer'] == {
        'opt_coordinates': 'cartesian', 'dynamic_level': 0,
        'full_hess_every': 0, 'geom_maxiter': 50, 'g_convergence': 'gau_tight',
    }
    assert request['method'] == 'HF/6-31G*'
    assert request['stages'] == [{'a': .0005, 'b': .1}, {'a': .001, 'b': .1}]
    assert request['grid']['vdw_scale_factors'] == [1.4, 1.6, 1.8, 2.0]


def test_worker_forwards_optimizer_and_stops_before_resp_on_failure(
        tmp_path, monkeypatch, recorded_request):
    _, request = recorded_request
    psi4 = SimpleNamespace(
        set_num_threads=Mock(), set_memory=Mock(), set_options=Mock(),
        core=SimpleNamespace(set_output_file=Mock()),
        geometry=Mock(return_value=SimpleNamespace(set_name=Mock())),
        optimize=Mock(side_effect=WorkerBoundary),
    )
    resp = SimpleNamespace(resp=Mock())
    monkeypatch.setitem(sys.modules, 'psi4', psi4)
    monkeypatch.setitem(sys.modules, 'resp', resp)
    monkeypatch.chdir(tmp_path)
    with pytest.raises(WorkerBoundary):
        hf_resp.worker(tmp_path)
    assert psi4.optimize.call_args.args == ('hf/6-31G*',)
    assert psi4.optimize.call_args.kwargs['optimizer_keywords'] == request['optimizer']
    assert psi4.set_options.call_args.args[0]['scf_type'] == 'pk'
    resp.resp.assert_not_called()
    assert not (tmp_path/'charges.json').exists()


def test_previous_optimizer_cache_is_rejected(tmp_path, recorded_request):
    mol, request = recorded_request
    request['schema'] = 'hf-resp/3'
    request['optimizer'].pop('full_hess_every')
    (tmp_path/'charges.json').write_text(json.dumps({'request': request, 'files': {}}))
    with pytest.raises(ValueError, match='RESP cache parameters differ'):
        hf_resp.compute(mol, tmp_path)
