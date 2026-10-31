"""Transient transport failures cannot publish incomplete model assets."""
import hashlib
from http.client import IncompleteRead
import io
from pathlib import Path
import ssl
import sys
from urllib.error import HTTPError, URLError

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import install_predocking_assets as assets

CONTENT = b'verified public asset fixture'
ENTRY = dict(name='fixture', url='https://example.invalid/model',
             sha256=hashlib.sha256(CONTENT).hexdigest())


@pytest.mark.parametrize('failure', [
    URLError(ssl.SSLEOFError(8, 'TLS EOF')),
    ConnectionResetError('connection reset'),
    TimeoutError('read timed out'),
    HTTPError(ENTRY['url'], 503, 'temporary service failure', {}, None),
])
def test_transient_failure_retries_before_atomic_verified_publish(tmp_path, monkeypatch, failure):
    target = tmp_path / 'model'
    calls, waits = [], []
    def open_asset(*args, **kwargs):
        assert not target.exists()
        calls.append(args)
        if len(calls) == 1:
            raise failure
        return io.BytesIO(CONTENT)
    monkeypatch.setattr(assets, 'urlopen', open_asset)
    monkeypatch.setattr(assets.time, 'sleep', waits.append)
    assets.download_verified(ENTRY, target)
    assert target.read_bytes() == CONTENT
    assert len(calls) == 2 and waits == [1]
    assert not target.with_suffix('.part').exists()


def test_partial_transfer_is_discarded_before_retry(tmp_path, monkeypatch):
    class Interrupted(io.BytesIO):
        def read(self, size=-1):
            if self.tell():
                raise IncompleteRead(b'', len(CONTENT))
            return super().read(3)
    responses = iter([Interrupted(b'bad partial bytes'), io.BytesIO(CONTENT)])
    monkeypatch.setattr(assets, 'urlopen', lambda *a, **kw: next(responses))
    monkeypatch.setattr(assets.time, 'sleep', lambda _: None)
    target = tmp_path / 'model'
    assets.download_verified(ENTRY, target)
    assert target.read_bytes() == CONTENT
    assert not target.with_suffix('.part').exists()


def test_persistent_transport_failure_is_bounded(tmp_path, monkeypatch):
    calls, waits = [], []
    def fail(*args, **kwargs):
        calls.append(args)
        raise URLError(ssl.SSLEOFError(8, 'persistent TLS EOF'))
    monkeypatch.setattr(assets, 'urlopen', fail)
    monkeypatch.setattr(assets.time, 'sleep', waits.append)
    target = tmp_path / 'model'
    with pytest.raises(URLError):
        assets.download_verified(ENTRY, target)
    assert len(calls) == 3 and waits == [1, 2]
    assert not target.exists() and not target.with_suffix('.part').exists()


@pytest.mark.parametrize('failure', [
    URLError(ssl.SSLCertVerificationError(1, 'certificate rejected')),
    HTTPError(ENTRY['url'], 403, 'forbidden', {}, None),
    PermissionError('destination denied'),
])
def test_permanent_failure_is_not_retried(tmp_path, monkeypatch, failure):
    calls = []
    def fail(*args, **kwargs):
        calls.append(args)
        raise failure
    monkeypatch.setattr(assets, 'urlopen', fail)
    monkeypatch.setattr(assets.time, 'sleep', lambda _: pytest.fail('must not retry'))
    target = tmp_path / 'model'
    with pytest.raises(type(failure)):
        assets.download_verified(ENTRY, target)
    assert len(calls) == 1
    assert not target.exists() and not target.with_suffix('.part').exists()


def test_hash_mismatch_is_not_retried_or_published(tmp_path, monkeypatch):
    monkeypatch.setattr(assets, 'urlopen', lambda *a, **kw: io.BytesIO(b'wrong bytes'))
    monkeypatch.setattr(assets.time, 'sleep', lambda _: pytest.fail('must not retry'))
    target = tmp_path / 'model'
    with pytest.raises(ValueError, match='hash mismatch'):
        assets.download_verified(ENTRY, target)
    assert not target.exists() and not target.with_suffix('.part').exists()
