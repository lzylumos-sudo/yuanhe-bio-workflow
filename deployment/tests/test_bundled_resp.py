import importlib.util
import json
from pathlib import Path

import pytest


spec = importlib.util.spec_from_file_location('verify_resp_source', Path(__file__).resolve().parents[1]/'public_minimal/verify_resp_source.py')
source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)


def test_bundled_resp_retains_original_source_pin():
    path, proof = source.verify()
    assert path.is_file()
    assert proof['commit'] == 'e0da367e9bef96eb469b5696ff7d6ea58e45d638'
    assert proof['files'] == 26
    assert proof['upstream_content_compared'] is False


def test_corrupt_resp_bundle_is_rejected(tmp_path, monkeypatch):
    (tmp_path/'resp-source.json').write_text(json.dumps({'archive': 'bad.zip', 'sha256': '0'*64}))
    (tmp_path/'bad.zip').write_bytes(b'corrupt')
    monkeypatch.setattr(source, 'ROOT', tmp_path)
    with pytest.raises(ValueError, match='hash mismatch'):
        source.verify()
