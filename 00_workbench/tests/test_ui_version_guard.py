"""HTTP-level release mismatch protection; dispatch is synthetic and never computes."""
import http.client
import json
from pathlib import Path
import sys
import threading

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from webui import server


@pytest.fixture
def http_service(monkeypatch):
    calls = []

    def synthetic_dispatch(method, path, query, body):
        calls.append((method, path, body))
        return server._json({'dispatched': True})

    monkeypatch.setattr(server, 'dispatch', synthetic_dispatch)
    monkeypatch.setattr(server, 'VERSION', 'test-release-1')
    service = server.Server(('127.0.0.1', 0), server.Handler)
    thread = threading.Thread(target=service.serve_forever, kwargs={'poll_interval': .01}, daemon=True)
    thread.start()

    def request(method, path, version=None, *, browser=None):
        headers = {'Content-Type': 'application/json'}
        if version is not None:
            headers['X-Yuanhe-UI-Version'] = version
        if browser in ('origin', 'both'):
            headers['Origin'] = 'http://127.0.0.1:' + str(service.server_address[1])
        if browser in ('fetch', 'both'):
            headers['Sec-Fetch-Site'] = 'same-origin'
        connection = http.client.HTTPConnection('127.0.0.1', service.server_address[1], timeout=5)
        try:
            connection.request(method, path, body=b'{}' if method == 'POST' else None, headers=headers)
            response = connection.getresponse()
            return response.status, json.loads(response.read())
        finally:
            connection.close()

    yield request, calls
    service.shutdown()
    service.server_close()
    thread.join(timeout=2)


@pytest.mark.parametrize('path', [
    '/api/modules/jobs',
    '/api/modules/jobs/MOD-' + 'a'*32 + '/start',
    '/api/project/runs',
    '/api/project/run/example/review/gromacs',
    '/api/studio/environment',
    '/api/studio/upload',
])
def test_stale_browser_writes_are_refused_before_dispatch(http_service, path):
    request, calls = http_service
    status, payload = request('POST', path, 'test-release-old', browser='both')
    assert status == 409
    assert payload['error_code'] == 'UI_VERSION_MISMATCH'
    assert payload['recoverable'] is True
    assert payload['details'] == {'ui_version': 'test-release-old', 'service_version': 'test-release-1'}
    assert 'No change was applied' in payload['suggested_action']
    assert not calls


@pytest.mark.parametrize('browser', ['origin', 'fetch', 'both'])
@pytest.mark.parametrize('path', [
    '/api/modules/jobs', '/api/modules/jobs/MOD-' + 'a'*32 + '/start',
    '/api/studio/environment', '/api/studio/upload',
])
def test_cached_browser_without_version_is_refused_before_dispatch(http_service, browser, path):
    request, calls = http_service
    status, payload = request('POST', path, browser=browser)
    assert status == 409
    assert payload['error_code'] == 'UI_VERSION_MISMATCH'
    assert payload['recoverable'] is True
    assert payload['details'] == {'ui_version': None, 'service_version': 'test-release-1'}
    assert not calls


@pytest.mark.parametrize('browser', ['origin', 'fetch', 'both'])
def test_current_browser_with_version_is_accepted(http_service, browser):
    request, calls = http_service
    assert request('POST', '/api/modules/jobs', 'test-release-1', browser=browser)[0] == 200
    assert len(calls) == 1


def test_page_that_was_compatible_is_blocked_after_service_upgrade(http_service, monkeypatch):
    request, calls = http_service
    assert request('POST', '/api/modules/jobs', 'test-release-1')[0] == 200
    monkeypatch.setattr(server, 'VERSION', 'test-release-2')
    status, payload = request('POST', '/api/modules/jobs', 'test-release-1')
    assert status == 409
    assert payload['details']['service_version'] == 'test-release-2'
    assert len(calls) == 1
    assert request('POST', '/api/modules/jobs', 'test-release-2')[0] == 200
    assert len(calls) == 2


@pytest.mark.parametrize('path', [
    '/api/modules/jobs/MOD-' + 'a'*32 + '/cancel',
    '/api/project/run/example/cancel',
    '/api/jobs/abc123/cancel',
    '/api/studio/run/example/cancel-operation',
])
@pytest.mark.parametrize('version', [None, 'test-release-old'])
def test_stale_browser_can_still_stop_jobs(http_service, path, version):
    request, calls = http_service
    assert request('POST', path, version, browser='both')[0] == 200
    assert len(calls) == 1


def test_headerless_cli_and_reads_remain_compatible(http_service):
    request, calls = http_service
    assert request('POST', '/api/modules/jobs')[0] == 200
    assert request('GET', '/api/v1/versions', 'test-release-old')[0] == 200
    assert request('GET', '/api/v1/versions', browser='both')[0] == 200
    assert len(calls) == 3


def test_cancel_suffix_does_not_exempt_unrelated_writes(http_service):
    request, calls = http_service
    assert request('POST', '/api/studio/environment/cancel', 'test-release-old')[0] == 409
    assert not calls
