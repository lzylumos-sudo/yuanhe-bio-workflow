# -*- coding: utf-8 -*-
"""webui/capture_fixtures.py -- snapshot the JSON API into webui/testdata/.

The snapshots are used by ``node webui/test_pages.mjs`` so the front-end can be
exercised offline (no browser needed).  The server is started on an ephemeral
port inside this process, so the fixtures are exactly what the API returns.

    python webui/capture_fixtures.py [--wait-pocket]

Both the fixtures and the pocket-voxel cache are written under webui/.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import threading
import time
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from webui import server                                            # noqa: E402

ENDPOINTS = {
    'meta': '/api/meta',
    'overview': '/api/overview',
    'hotspots': '/api/hotspots',
    'pharmacophore': '/api/pharmacophore',
    'pocket': '/api/pocket',
    'scaffolds': '/api/scaffolds',
    'filters': '/api/filters',
    'poses': '/api/poses',
    'evidence_report': '/api/evidence/report',
    'evidence_priors': '/api/evidence/priors',
    'evidence_census': '/api/evidence/census',
    'evidence_interface': '/api/evidence/interface',
    'docs': '/api/docs',
    'jobs': '/api/jobs',
    'scene_cj3': '/api/scene?ref=SglCJ3&ligand=0',
    'scene_m': '/api/scene?ref=SglM&ligand=1',
    'pocket_points': '/api/pocket/points',
    'docs_01': '/api/docs/' + urllib.parse.quote('01_Molecular scaffold design algorithm.md'),
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--timeout', type=float, default=180.0,
                    help='how long to wait for the pocket voxel cache')
    a = ap.parse_args()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'testdata')
    os.makedirs(out, exist_ok=True)

    httpd = server.Server(('127.0.0.1', 0), server.Handler)
    port = httpd.server_address[1]
    threading.Thread(target=httpd.serve_forever, kwargs=dict(poll_interval=0.2),
                     daemon=True).start()
    base = f'http://127.0.0.1:{port}'
    print(f'[fixtures] local server {base}')

    deadline = time.time() + a.timeout
    for name, path in ENDPOINTS.items():
        while True:
            try:
                with urllib.request.urlopen(base + path, timeout=120) as resp:
                    body = resp.read()
                break
            except Exception as exc:                                # 202 for pocket points
                if name == 'pocket_points' and time.time() < deadline:
                    time.sleep(2.0)
                    continue
                print(f'[fixtures] {name}: {exc}')
                body = b'{}'
                break
        with open(os.path.join(out, f'{name}.json'), 'wb') as fh:
            fh.write(body)
        print(f'[fixtures] {name:<20} {len(body):>9,} B')

    payload = json.load(open(os.path.join(out, 'pocket_points.json'), encoding='utf-8'))
    if payload.get('status') != 'ready':
        print('[fixtures] warning: pocket voxels not ready (front-end test stubs it)')
    httpd.shutdown()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
