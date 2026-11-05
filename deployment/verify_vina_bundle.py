#!/usr/bin/env python3
"""Verify the unchanged official Vina executable without executing it."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def verified_binary(root=ROOT):
    spec = json.loads((root / 'deployment/vina-bundle.json').read_text())
    archive = (root / spec['archive']).read_bytes()
    if hashlib.sha256(archive).hexdigest() != spec['archive_sha256']:
        raise ValueError('Vina archive hash mismatch')
    license_text = (root / spec['license_path']).read_bytes()
    if hashlib.sha256(license_text).hexdigest() != spec['license_sha256']:
        raise ValueError('Vina license hash mismatch')
    binary = gzip.decompress(archive)
    if (len(binary) != spec['binary_bytes'] or
            hashlib.sha256(binary).hexdigest() != spec['binary_sha256']):
        raise ValueError('Official Vina binary size/hash mismatch')
    return binary, spec


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upstream', action='store_true')
    args = parser.parse_args()
    binary, spec = verified_binary()
    if args.upstream:
        with tempfile.TemporaryDirectory() as directory:
            for key, expected in [('upstream_url', binary),
                                  ('upstream_license_url', (ROOT / spec['license_path']).read_bytes())]:
                path = Path(directory) / key
                subprocess.run(['curl', '-fL', '--retry', '3', '--retry-all-errors',
                                '--connect-timeout', '20', '--max-time', '120',
                                '--retry-max-time', '180', spec[key], '-o', str(path)], check=True)
                if path.read_bytes() != expected:
                    raise ValueError('Bundled Vina differs from official upstream: ' + key)
    print(json.dumps({'version': spec['version'], 'binary_sha256': spec['binary_sha256'],
                      'bytes': len(binary), 'upstream_compared': args.upstream,
                      'executable_run': False}))


if __name__ == '__main__':
    main()
