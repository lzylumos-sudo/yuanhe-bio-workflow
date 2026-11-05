#!/usr/bin/env python3
"""Acquire pinned public generation models into a new, reusable offline bundle."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import ssl
import time
from http.client import IncompleteRead
from urllib.error import HTTPError, URLError
import urllib.request
import urllib.error
from source_checkout import checked_checkout

ROOT = Path(__file__).resolve().parents[1]
DOWNLOADS = {
    'peptide2mol': ('Peptide2Mol', 'checkpoint',
        'https://github.com/BLUE-Flowing/Peptide2Mol/releases/download/v1.0/PMT_major.ckpt', 67462479),
    'reinvent4': ('REINVENT4', 'prior',
        'https://zenodo.org/api/records/20701824/files/libinvent.prior/content', 94756131),
}


def plan(destination, root=ROOT, *, resume=False, urls=None):
    destination = Path(destination)
    if not destination.is_absolute() or destination.is_symlink() or (destination.exists() and not resume):
        raise ValueError('Choose an absolute NEW bundle directory; existing paths are preserved')
    pins = json.loads((root / 'deployment/model_dependencies.json').read_text())
    entries = []
    for key, (name, field, url, size) in DOWNLOADS.items():
        url = (urls or {}).get(key) or url
        if not url.startswith('https://'):
            raise ValueError('Model download overrides require HTTPS; pinned hashes remain mandatory')
        model = pins[key]
        entries.append(dict(name=name, repository=model['source_url'], commit=model['repo_commit'],
                            path=model[field], sha256=model[field + '_sha256'], url=url, bytes=size))
    if resume:
        saved = json.loads((destination / 'model-bundle.json').read_text())
        def identity(rows):
            return [{k: v for k, v in row.items() if k != 'url'} for row in rows]
        if saved.get('schema') != 'yuanhe.model-bundle/1' or identity(saved['models']) != identity(entries):
            raise ValueError('Existing bundle does not match the pinned source/weight identities')
    return dict(schema='yuanhe.model-bundle/1', destination=str(destination), models=entries,
                status='planned', inference='not_run')


def download(entry, target, retries=4):
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_name(target.name + '.part')
    if target.exists() or target.is_symlink():
        raise ValueError('Refusing to overwrite model weights: ' + str(target))
    if partial.is_symlink():
        raise ValueError('Refusing a symlink partial download')
    for attempt in range(retries):
        try:
            offset = partial.stat().st_size if partial.exists() else 0
            if offset > entry['bytes']:
                raise ValueError('Partial download exceeds pinned size')
            if offset < entry['bytes']:
                headers = {'User-Agent': 'Yuanhe-model-bootstrap/1'}
                if offset:
                    headers['Range'] = f'bytes={offset}-'
                url = entry['url']
                if attempt % 2 and url == DOWNLOADS['reinvent4'][2]:
                    url = 'https://zenodo.org/records/20701824/files/libinvent.prior?download=1'
                request = urllib.request.Request(url, headers=headers)
                with urllib.request.urlopen(request, timeout=120) as response:
                    status = getattr(response, 'status', 200)
                    if status == 206:
                        content_range = response.headers.get('Content-Range', '')
                        if not content_range.startswith(f'bytes {offset}-') or not content_range.endswith('/' + str(entry['bytes'])):
                            raise ValueError('Invalid partial response range')
                    elif status == 200:
                        offset = 0  # Servers may ignore Range; restart without appending duplicate bytes.
                    else:
                        raise ValueError('Unexpected download HTTP status: ' + str(status))
                    with partial.open('ab' if offset else 'wb') as output:
                        size = offset
                        while chunk := response.read(1 << 20):
                            size += len(chunk)
                            if size > entry['bytes']:
                                raise ValueError('Unexpected download size: ' + entry['name'])
                            output.write(chunk)
            if partial.stat().st_size != entry['bytes']:
                if attempt + 1 == retries:
                    raise ValueError('Model size mismatch after retries; partial download retained: ' + entry['name'])
                raise OSError('Incomplete model download; partial bytes retained')
            digest = hashlib.sha256()
            with partial.open('rb') as handle:
                for chunk in iter(lambda: handle.read(1 << 20), b''):
                    digest.update(chunk)
            if digest.hexdigest() != entry['sha256']:
                raise ValueError('Model size or SHA-256 mismatch: ' + entry['name'])
            partial.replace(target)
            return
        except (OSError, urllib.error.URLError, IncompleteRead) as exc:
            reason = getattr(exc, 'reason', exc)
            transient = (exc.code in (408, 429, 500, 502, 503, 504) if isinstance(exc, HTTPError)
                         else isinstance(reason, (OSError, IncompleteRead)))
            if isinstance(reason, ssl.SSLCertVerificationError) or not transient or attempt + 1 == retries:
                raise
            print(f"Transient download failure for {entry['name']}: {exc}; partial bytes retained", flush=True)
            time.sleep(2 ** attempt)


def execute(record, *, resume=False):
    if not shutil.which('git'):
        raise ValueError('Install Git first; see DEPLOY.md')
    destination = Path(record['destination'])
    destination.mkdir(parents=True, exist_ok=resume)
    state = {**record, 'status': 'fetching', 'verified': []}
    receipt = destination / 'model-bundle.json'
    try:
        for entry in record['models']:
            print('Fetching pinned source and weights: ' + entry['name'], flush=True)
            source = destination / 'vendor' / entry['name']
            source.mkdir(parents=True, exist_ok=resume)
            commands = [['git', 'init', str(source)],
                        ['git', '-C', str(source), 'config', 'core.autocrlf', 'false'],
                        ['git', '-C', str(source), 'remote', 'add', 'origin', entry['repository']],
                        ['git', '-C', str(source), 'fetch', '--depth=1', 'origin', entry['commit']]]
            if resume and (source / '.git/config').is_file():
                remotes = subprocess.check_output(['git', '-C', str(source), 'remote'], text=True).split()
                if 'origin' in remotes:
                    commands[2] = ['git', '-C', str(source), 'remote', 'set-url', 'origin', entry['repository']]
            for argv in commands:
                subprocess.run(argv, check=True, timeout=900)
            checked_checkout(source, 'FETCH_HEAD')
            actual = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
            if actual != entry['commit'] or not (source / 'LICENSE').is_file():
                raise ValueError('Source revision or upstream LICENSE mismatch: ' + entry['name'])
            weight = destination / entry['path']
            if resume and weight.is_file() and not weight.is_symlink():
                from restore_models import sha256
                if weight.stat().st_size != entry['bytes'] or sha256(weight) != entry['sha256']:
                    raise ValueError('Existing model weight differs; refusing overwrite: ' + str(weight))
            else:
                download(entry, weight)
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
    parser.add_argument('--resume', action='store_true', help='Resume a matching recorded bundle; reverify all source and weight hashes')
    parser.add_argument('--peptide2mol-url', help='Explicit HTTPS mirror for the same hash-pinned checkpoint')
    parser.add_argument('--libinvent-url', help='Explicit HTTPS mirror for the same hash-pinned prior')
    args = parser.parse_args()
    try:
        record = plan(args.destination, resume=args.resume,
                      urls={'peptide2mol': args.peptide2mol_url, 'reinvent4': args.libinvent_url})
        if args.execute:
            execute(record, resume=args.resume)
        else:
            print(json.dumps(record, indent=2))
            print('Plan only; add --execute to download. See upstream terms in DEPLOY.md.')
    except (ValueError, OSError, subprocess.SubprocessError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
