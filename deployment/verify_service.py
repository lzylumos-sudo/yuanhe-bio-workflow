#!/usr/bin/env python3
"""Read-only verification that the expected deployment actually serves a URL."""
import argparse
import json
from pathlib import Path
import sys
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]


def verify(base_url, expected_version, expected_commit, fetch=None):
    if fetch is None:
        def fetch(path):
            with urlopen(base_url.rstrip('/') + path, timeout=30) as response:
                return json.load(response)
    health = fetch('/api/health')
    versions = fetch('/api/v1/versions')
    if health.get('ok') is not True or health.get('version') != expected_version or versions.get('software_version') != expected_version:
        raise ValueError('Wrong service version; inspect service startup and port ownership')
    source = versions.get('source', {})
    commit = versions.get('git', {}).get('commit')
    if commit == 'unknown':
        commit = source.get('commit') if source.get('verified') else None
    if commit != expected_commit or versions.get('git', {}).get('dirty') is True:
        raise ValueError('Wrong or modified source commit; version text alone is insufficient')
    diagnostics = fetch('/api/studio/diagnostics')
    return {'service_verified': True, 'software_version': expected_version,
            'source_commit': commit, 'diagnostics_endpoint': diagnostics.get('status'),
            'readiness_passed': diagnostics.get('all_runtime_checks_pass', False),
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
