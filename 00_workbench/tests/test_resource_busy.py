"""Real lock contention and API recovery; no scientific engine is executed."""
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from webui import module_jobs, scheduling, server, workflow as w


@pytest.fixture
def isolated_job(tmp_path, monkeypatch):
    monkeypatch.setattr(w, 'ROOT', tmp_path)
    monkeypatch.setattr(w, '_active', {})
    monkeypatch.setattr(w, 'source_snapshot', lambda: {'synthetic_source': 'unchanged'})
    monkeypatch.setattr(module_jobs, 'JOBS', tmp_path/'module_jobs')
    ident = 'MOD-' + 'a'*32
    d = module_jobs.folder(ident)
    d.mkdir(parents=True)
    w.atomic(d/'source.json', w.source_snapshot())
    w.atomic(d/'request.json', {'schema': 'yuanhe.module-job/1', 'module': 'docking_control',
                              'inputs': {}, 'parameters': {}, 'config': {'chem_python': sys.executable}})
    w.atomic(d/'job.json', {'id': ident, 'module': 'docking_control', 'created': 0,
                          'source_sha256': w.digest(d/'source.json'),
                          'request_sha256': w.digest(d/'request.json')})
    return ident, d


def dispatch_start(ident):
    status, _, body = server.dispatch('POST', '/api/modules/jobs/'+ident+'/start', {}, {})
    return status, json.loads(body)


@pytest.mark.skipif(sys.platform != 'linux', reason='Linux flock contract')
def test_busy_api_leaves_job_unstarted_and_can_be_started_after_release(isolated_job, tmp_path, monkeypatch):
    ident, d = isolated_job
    before = {name: (d/name).read_bytes() for name in ['job.json', 'source.json', 'request.json']}
    spawned = []

    def synthetic_worker(rid, kind, op, record, argv, lease, log, env):
        spawned.append((rid, kind, op))
        w.atomic(record, {**w.read(record), 'status': 'completed'})
        log.close()
        lease.close()

    monkeypatch.setattr(w, 'spawn_job', synthetic_worker)
    with scheduling.acquire(tmp_path):
        assert module_jobs.inventory()['resource_status']['available'] is False
        status, body = dispatch_start(ident)
        assert status == 409
        assert body['error_code'] == 'RESOURCE_BUSY'
        assert body['recoverable'] is True
        assert body['details']['queued'] is False
        assert body['run_id'] == ident
        assert 'not been queued' in body['suggested_action']
        assert body['error'] == body['message']
        assert not spawned
        assert not (d/'operations').exists()
        assert all((d/name).read_bytes() == content for name, content in before.items())
    assert module_jobs.inventory()['resource_status']['available'] is True
    status, body = dispatch_start(ident)
    assert status == 200
    assert len(spawned) == 1
    assert body['operations'][0]['status'] == 'completed'


@pytest.mark.skipif(sys.platform != 'linux', reason='Linux flock contract')
def test_probe_preserves_exclusive_and_lane_policies_and_releases_its_lock(tmp_path):
    assert scheduling.exclusive_status(tmp_path)['available'] is True
    with scheduling.acquire(tmp_path, tmp_path/'run', 'compute'):
        assert scheduling.exclusive_status(tmp_path)['available'] is False
        with scheduling.acquire(tmp_path, tmp_path/'run', 'receptor'):
            with pytest.raises(scheduling.ResourceBusyError):
                scheduling.acquire(tmp_path)
    assert scheduling.exclusive_status(tmp_path)['available'] is True
    with scheduling.acquire(tmp_path):
        with pytest.raises(scheduling.ResourceBusyError):
            scheduling.acquire(tmp_path, tmp_path/'run', 'receptor')


def test_active_module_reports_recoverable_conflict(isolated_job, monkeypatch):
    ident, d = isolated_job
    monkeypatch.setattr(w, '_active', {(ident, 'module', 'synthetic'): object()})
    status, body = dispatch_start(ident)
    assert status == 409
    assert body['error_code'] == 'RESOURCE_BUSY'
    assert not (d/'operations').exists()


def test_provenance_guard_is_not_reclassified_as_busy(isolated_job):
    ident, d = isolated_job
    (d/'request.json').write_text('{}')
    status, body = dispatch_start(ident)
    assert status == 400
    assert body['error_code'] == 'INPUT_ERROR'
    assert body['recoverable'] is False
    assert not (d/'operations').exists()


def test_genuine_runtime_failure_still_returns_500(isolated_job, monkeypatch):
    ident, _ = isolated_job

    def fail(*args):
        raise RuntimeError('synthetic execution failure')

    monkeypatch.setattr(module_jobs, 'launch', fail)
    status, body = dispatch_start(ident)
    assert status == 500
    assert body['error_code'] == 'RUNTIME_ERROR'
    assert body['recoverable'] is False


def test_unavailable_probe_is_unknown_not_an_assertion_of_free_capacity(monkeypatch, tmp_path):
    def fail(*args):
        raise OSError('synthetic lock filesystem failure')

    monkeypatch.setattr(scheduling, 'acquire', fail)
    status = scheduling.exclusive_status(tmp_path)
    assert status['available'] is None
    assert status['queued'] is False
