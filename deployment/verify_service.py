#!/usr/bin/env python3
"""Read-only verification that the expected deployment actually serves a URL."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]


def verify(base_url, expected_version, expected_commit, fetch=None, fetch_bytes=None):
    if fetch is None:
        def fetch(path):
            with urlopen(base_url.rstrip('/') + path, timeout=30) as response:
                return json.load(response)
        if fetch_bytes is None:
            def fetch_bytes(path):
                with urlopen(base_url.rstrip('/') + path, timeout=30) as response:
                    return response.read()
    health = fetch('/api/health')
    versions = fetch('/api/v1/versions')
    if health.get('ok') is not True or health.get('version') != expected_version or versions.get('software_version') != expected_version:
        raise ValueError('Wrong service version; inspect service startup and port ownership')
    source = versions.get('source', {})
    git = versions.get('git', {})
    commit = git.get('commit')
    if commit == 'unknown':
        clean_source = source.get('verified') is True and git.get('dirty') is not True
        commit = source.get('commit') if clean_source else None
    else:
        clean_source = git.get('dirty') is False
    if commit != expected_commit or not clean_source:
        raise ValueError('Wrong or modified source commit; version text alone is insufficient')
    assets = {}
    if fetch_bytes is not None:
        static = ROOT / '00_workbench/webui/static'
        for name in ('index.html', 'modules.html', 'js/platform-version.js', 'js/api.js',
                     'js/studio.js', 'js/studio-tools.js', 'js/module-jobs.js'):
            expected = hashlib.sha256((static / name).read_bytes()).hexdigest()
            actual = hashlib.sha256(fetch_bytes('/' + name)).hexdigest()
            if actual != expected:
                raise ValueError('Frontend differs from the expected source: ' + name)
            assets[name] = actual
    diagnostics = fetch('/api/studio/diagnostics')
    return {'service_verified': True, 'software_version': expected_version,
            'source_commit': commit, 'diagnostics_endpoint': diagnostics.get('status'),
            'readiness_passed': diagnostics.get('all_runtime_checks_pass', False),
            'frontend_verified': bool(assets), 'frontend_sha256': assets,
            'scientific_acceptance': 'not_run'}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--url', required=True)
    ap.add_argument('--expected-commit', required=True)
    ap.add_argument('--expected-version', default=(ROOT / 'VERSION').read_text().strip())
    ap.add_argument('--out', type=Path)
    args = ap.parse_args()
    try:
        report = verify(args.url, args.expected_version, args.expected_commit)
    except Exception as exc:
        report = {'service_verified': False, 'error': str(exc), 'scientific_acceptance': 'not_run'}
    text = json.dumps(report, indent=2)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + '\n')
    print(text)
    return 0 if report['service_verified'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
