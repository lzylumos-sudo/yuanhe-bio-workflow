import hashlib
import io
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import fetch_generation_models as fetch
from mirror_lock import mirrored_lock


@pytest.mark.parametrize('range_supported', [True, False])
def test_partial_download_resumes_or_safely_restarts(tmp_path, monkeypatch, range_supported):
    content = b'0123456789'
    entry = dict(name='fixture', url='https://example.invalid/model', bytes=len(content),
                 sha256=hashlib.sha256(content).hexdigest())
    target = tmp_path / 'model'; target.with_name('model.part').write_bytes(content[:3])
    def request(req, **kw):
        assert req.get_header('Range') == 'bytes=3-'
        result = io.BytesIO(content[3:] if range_supported else content)
        result.status = 206 if range_supported else 200
        result.headers = {'Content-Range': 'bytes 3-9/10'}
        return result
    monkeypatch.setattr(fetch.urllib.request, 'urlopen', request)
    fetch.download(entry, target)
    assert target.read_bytes() == content


def test_wrong_range_cannot_be_promoted(tmp_path, monkeypatch):
    target = tmp_path / 'model'; target.with_name('model.part').write_bytes(b'abc')
    response = io.BytesIO(b'def'); response.status = 206
    response.headers = {'Content-Range': 'bytes 0-2/6'}
    monkeypatch.setattr(fetch.urllib.request, 'urlopen', lambda *a, **kw: response)
    with pytest.raises(ValueError, match='range'):
        fetch.download(dict(name='fixture', url='https://example.invalid', bytes=6, sha256='0'*64), target)
    assert not target.exists() and target.with_name('model.part').read_bytes() == b'abc'


def test_resume_requires_matching_pins_but_allows_explicit_https_mirror(tmp_path):
    dest = tmp_path / 'bundle'; record = fetch.plan(dest)
    dest.mkdir(); receipt = dest / 'model-bundle.json'; receipt.write_text(json.dumps(record))
    resumed = fetch.plan(dest, resume=True, urls={'reinvent4': 'https://example.invalid/prior'})
    assert resumed['models'][1]['sha256'] == record['models'][1]['sha256']
    record['models'][1]['sha256'] = '0'*64; receipt.write_text(json.dumps(record))
    with pytest.raises(ValueError, match='identities'): fetch.plan(dest, resume=True)


def test_mirror_changes_location_only_and_retains_artifact_identity():
    lock = '# pinned\n@EXPLICIT\nhttps://conda.anaconda.org/conda-forge/linux-64/pkg.conda#abcd\n'
    result = mirrored_lock(lock, 'https://example.invalid/conda-forge')
    assert result == lock.replace('https://conda.anaconda.org/conda-forge', 'https://example.invalid/conda-forge')
    with pytest.raises(ValueError): mirrored_lock(lock, 'http://example.invalid')
