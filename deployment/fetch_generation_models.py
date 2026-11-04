#!/usr/bin/env python3
"""Acquire pinned public generation models into a new, reusable offline bundle."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from http.client import IncompleteRead
from urllib.error import HTTPError, URLError
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
DOWNLOADS = {
    'peptide2mol': ('Peptide2Mol', 'checkpoint',
        'https://github.com/BLUE-Flowing/Peptide2Mol/releases/download/v1.0/PMT_major.ckpt', 67462479),
    'reinvent4': ('REINVENT4', 'prior',
        'https://zenodo.org/api/records/20701824/files/libinvent.prior/content', 94756131),
}


def plan(destination, root=ROOT):
    destination = Path(destination)
    if not destination.is_absolute() or destination.exists() or destination.is_symlink():
        raise ValueError('Choose an absolute NEW bundle directory; existing paths are preserved')
    pins = json.loads((root / 'deployment/model_dependencies.json').read_text())
    entries = []
    for key, (name, field, url, size) in DOWNLOADS.items():
        model = pins[key]
        entries.append(dict(name=name, repository=model['source_url'], commit=model['repo_commit'],
                            path=model[field], sha256=model[field + '_sha256'], url=url, bytes=size))
    return dict(schema='yuanhe.model-bundle/1', destination=str(destination), models=entries,
                status='planned', inference='not_run')


def download(entry, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_name(target.name + '.part')
    if target.exists() or target.is_symlink():
        raise ValueError('Refusing to overwrite model weights: ' + str(target))
    if partial.exists() or partial.is_symlink():
        raise FileExistsError(partial)
    for attempt in range(4):
        url = entry['url']
        if attempt % 2 and url == DOWNLOADS['reinvent4'][2]:
            url = 'https://zenodo.org/records/20701824/files/libinvent.prior?download=1'
        try:
            download_once(entry, partial, url)
            break
        except (HTTPError, URLError, TimeoutError, ConnectionError, IncompleteRead) as exc:
            reason = getattr(exc, 'reason', exc)
            retry = (exc.code in (408, 429, 500, 502, 503, 504) if isinstance(exc, HTTPError)
                     else isinstance(reason, (TimeoutError, ConnectionError, IncompleteRead)))
            if not retry or attempt == 3:
                raise
            partial.unlink(missing_ok=True)
            print(f"Transient download failure for {entry['name']}: {exc}; retrying", flush=True)
            time.sleep(2 ** attempt)
    partial.replace(target)


def download_once(entry, partial, url):
    digest, size = hashlib.sha256(), 0
    request = urllib.request.Request(url, headers={'User-Agent': 'Yuanhe-model-bootstrap/1'})
    with urllib.request.urlopen(request, timeout=120) as response, partial.open('xb') as output:
        while chunk := response.read(1 << 20):
            size += len(chunk)
            if size > entry['bytes']:
                raise ValueError('Unexpected download size: ' + entry['name'])
            digest.update(chunk)
            output.write(chunk)
    if size != entry['bytes'] or digest.hexdigest() != entry['sha256']:
        raise ValueError('Model size or SHA-256 mismatch: ' + entry['name'])


def execute(record):
    if not shutil.which('git'):
        raise ValueError('Install Git first; see DEPLOY.md')
    destination = Path(record['destination'])
    destination.mkdir(parents=True, exist_ok=False)
    state = {**record, 'status': 'fetching', 'verified': []}
    receipt = destination / 'model-bundle.json'
    try:
        for entry in record['models']:
            print('Fetching pinned source and weights: ' + entry['name'], flush=True)
            source = destination / 'vendor' / entry['name']
            source.mkdir(parents=True)
            commands = [['git', 'init', str(source)],
                        ['git', '-C', str(source), 'config', 'core.autocrlf', 'false'],
                        ['git', '-C', str(source), 'remote', 'add', 'origin', entry['repository']],
                        ['git', '-C', str(source), 'fetch', '--depth=1', 'origin', entry['commit']],
                        ['git', '-C', str(source), 'checkout', '--detach', 'FETCH_HEAD']]
            for argv in commands:
                subprocess.run(argv, check=True, timeout=900)
            actual = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
            if actual != entry['commit'] or not (source / 'LICENSE').is_file():
                raise ValueError('Source revision or upstream LICENSE mismatch: ' + entry['name'])
            download(entry, destination / entry['path'])
            state['verified'].append(entry['name'])
        state['status'] = 'verified'
    except BaseException as exc:
        state.update(status='failed', error=str(exc))
        raise
    finally:
        receipt.write_text(json.dumps(state, indent=2) + '\n', encoding='utf-8')
    print('Bundle verified: ' + str(destination))
    print('Use this directory as model_source. Models were not installed or executed.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', required=True, type=Path)
    parser.add_argument('--execute', action='store_true', help='Download after reviewing upstream terms in DEPLOY.md')
    args = parser.parse_args()
    try:
        record = plan(args.destination)
        if args.execute:
            execute(record)
        else:
            print(json.dumps(record, indent=2))
            print('Plan only; add --execute to download. See upstream terms in DEPLOY.md.')
    except (ValueError, OSError, subprocess.SubprocessError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
