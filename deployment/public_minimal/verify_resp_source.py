"""Verify the bundled, unchanged RESP source before installation."""
import argparse
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import ssl
import time
from urllib.error import URLError
from urllib.request import urlopen
import zipfile

ROOT = Path(__file__).resolve().parent


def content(archive):
    result = {}
    with zipfile.ZipFile(io.BytesIO(archive)) as bundle:
        for info in bundle.infolist():
            if info.is_dir():
                continue
            parts = PurePosixPath(info.filename).parts
            if len(parts) < 2 or '..' in parts or parts[0] == '/':
                raise ValueError('Unexpected RESP archive path')
            name = '/'.join(parts[1:])
            if name in result or (info.external_attr >> 16) & 0o170000 == 0o120000:
                raise ValueError('Duplicate path or symlink in RESP source')
            result[name] = hashlib.sha256(bundle.read(info)).hexdigest()
    if 'LICENSE' not in result or 'resp/espfit.py' not in result:
        raise ValueError('Incomplete RESP source archive')
    return result


def verify(upstream=False):
    record = json.loads((ROOT/'resp-source.json').read_text())
    archive = ROOT/record['archive']
    data = archive.read_bytes()
    if hashlib.sha256(data).hexdigest() != record['sha256']:
        raise ValueError('Bundled RESP archive hash mismatch')
    files = content(data)
    if hashlib.sha256((ROOT/'vendor/RESP-LICENSE').read_bytes()).hexdigest() != files['LICENSE']:
        raise ValueError('RESP license differs from upstream source')
    if upstream:
        for attempt in range(3):
            try:
                with urlopen(record['upstream_archive_url'], timeout=120) as response:
                    official = response.read()
                break
            except (URLError, TimeoutError, ConnectionError) as exc:
                if isinstance(getattr(exc, 'reason', exc), ssl.SSLCertVerificationError) or attempt == 2:
                    raise
                time.sleep(attempt+1)
        if hashlib.sha256(official).hexdigest() != record['upstream_archive_sha256']:
            raise ValueError('Original upstream RESP archive hash mismatch')
        if files != content(official):
            raise ValueError('Bundled RESP files differ from the original pinned upstream source')
    return archive.resolve(), {'commit': record['commit'], 'sha256': record['sha256'],
                               'files': len(files), 'upstream_content_compared': upstream}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upstream', action='store_true', help='CI-only comparison to original pinned archive')
    parser.add_argument('--path', action='store_true', help='Print verified archive path for the installer')
    args = parser.parse_args()
    path, proof = verify(args.upstream)
    print(path if args.path else json.dumps(proof))
