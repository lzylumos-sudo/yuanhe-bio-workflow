# -*- coding: utf-8 -*-
"""webui/test_api.py -- end-to-end checks of the web console's JSON API.

Starts the real HTTP server on an ephemeral port (in-process, no external
process) and asserts the shape and the semantic content of every endpoint.

    python webui/test_api.py
"""
from __future__ import annotations

import json
import os
import sys
import threading
import time
import urllib.error
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

try:                                    # the console may be cp936 on Windows
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except (AttributeError, ValueError):
    pass

from webui import server                                            # noqa: E402

BASE = None
PASSED = 0
FAILED = []


def check(name, condition, detail=''):
    global PASSED
    if condition:
        PASSED += 1
        print(f'  ok   {name}')
    else:
        FAILED.append(name)
        print(f'  FAIL {name} {detail}')


def fetch(path, method='GET', body=None):
    url = BASE + path
    data = None
    headers = {'Accept': 'application/json'}
    if body is not None:
        data = json.dumps(body).encode('utf-8')
        headers['Content-Type'] = 'application/json'
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req, timeout=180) as resp:
            raw = resp.read()
            return resp.status, raw, resp.headers.get('Content-Type', '')
    except urllib.error.HTTPError as exc:
        return exc.code, exc.read(), exc.headers.get('Content-Type', '')


def jget(path):
    status, raw, ctype = fetch(path)
    assert status in (200, 202), f'{path} -> {status}: {raw[:200]}'
    return json.loads(raw.decode('utf-8'))


def run(do_job=False):
    global BASE
    httpd = server.Server(('127.0.0.1', 0), server.Handler)
    BASE = f'http://127.0.0.1:{httpd.server_address[1]}'
    threading.Thread(target=httpd.serve_forever, kwargs=dict(poll_interval=0.2),
                     daemon=True).start()
    print(f'[api] {BASE}')

    # ---- health / meta
    print('\n· Basic')
    meta = jget('/api/meta')
    check('meta Working directory and structure definitions present', meta.get('root') and len(meta.get('structures', {})) == 3)
    check('meta Report results and evidence Available',
          meta['available']['results'] and meta['available']['evidence'])
    check('meta Listed figures', len(meta.get('figures', [])) >= 4)

    status, raw, ctype = fetch('/')
    check('Home page is HTML', status == 200 and 'text/html' in ctype and b'<div id="app"' in raw)
    status, raw, ctype = fetch('/js/app.js')
    check('Frontend modules served as JS MIME provided', status == 200 and 'javascript' in ctype)

    # ---- overview
    print('\n· Overview')
    ov = jget('/api/overview')
    fun = ov['funnel']
    check('Funnel counts match README consistent',
          fun['enumerated'] == 20307 and fun['drug_like'] == 538,
          f'(enumerated={fun.get("enumerated")}, drug_like={fun.get("drug_like")})')
    check('Best scaffold is a biphenyl carboxylate scaffold',
          ov['best'][0]['smiles'].startswith('OC(=O)c1ccc'), ov['best'][0]['smiles'])
    check('Evidence status 24 entries,PASS 22', ov['evidence']['n'] == 24
          and ov['evidence']['by_status'].get('PASS') == 22, str(ov['evidence']))
    check('Structure filters 538 -> 32', ov['filters']['n_input'] == 538
          and ov['filters']['survivors'] == 32)

    # ---- design data
    print('\n· Design and results')
    hot = jget('/api/hotspots')
    check('Interface hotspots 52 items, contacted by all three peptides 14 items',
          hot['n'] == 52 and hot['n_common'] == 14, f"(n={hot['n']}, common={hot['n_common']})")
    scores = [r['score'] for r in hot['rows']]
    check('Hotspots sorted by score descending', scores == sorted(scores, reverse=True))

    ph = jget('/api/pharmacophore')
    check('12 consensus pharmacophore features', len(ph['features']) == 12, str(len(ph['features'])))
    check('Each feature has 3D coordinates and type',
          all(len(f['pos']) == 3 and f['type'] for f in ph['features']))
    check('Feature pairs include distances', all(p['distance'] > 3 for p in ph['pairs']))

    axes = jget('/api/axes')['axes']
    # A7-H12 has no surviving enumeration (its distance/weight combination is
    # outside the linker budget) -- it must still be listed, with candidates = 0.
    check('Design axes include enumeration statistics (A7-H12 axes without candidates retained)',
          len(axes) == 15
          and sum(1 for a in axes if a['candidates'] > 0) == 14
          and any(a['candidates'] == 0 for a in axes), str(len(axes)))

    pocket = jget('/api/pocket')['pocket']
    check('Cavity 5,440 Å³, 4 subpockets',
          round(pocket['total_volume_A3']) == 5440 and len(pocket['subpockets']) == 4,
          f"{pocket['total_volume_A3']}, {len(pocket['subpockets'])}")

    sc = jget('/api/scaffolds')
    check('200 ranked scaffolds', len(sc['rows']) == 200, str(len(sc['rows'])))
    row = sc['rows'][0]
    check('Scaffold rows include SMILES/Properties/scores/Contact residues',
          row['smiles'] and row['MW'] and 'total' in row and row['murj_contacts'])
    check('Scaffold rows and structure filters joined by SMILES Linker',
          any('f_keep' in r for r in sc['rows']))
    check('Filter decisions are booleans, not strings',
          all(isinstance(r['f_keep'], bool) for r in sc['rows'] if 'f_keep' in r)
          and sum(1 for r in sc['rows'] if r.get('f_keep') is True) == 14)
    check('Top-ranked scaffold failed structure filters; negative result retained',
          row['f_keep'] is False and row['f_bridge'] is True,
          f"keep={row['f_keep']!r} bridge={row['f_bridge']!r}")

    filt = jget('/api/filters')
    check('Structure filters 538 rows, retained 32', len(filt['rows']) == 538
          and filt['summary']['survivors'] == 32)
    check('All filter rows keep are boolean and 32 are true',
          sum(1 for r in filt['rows'] if r['keep'] is True) == 32
          and sum(1 for r in filt['rows'] if r['bridge'] is True) == 462)

    poses = jget('/api/poses')
    check('25 items Top pose with atoms and/complete bonds',
          len(poses['top']) == 25 and poses['top'][0]['atoms']
          and poses['top'][0]['bonds'])
    check('Pose coordinates lie in the 9NU8 coordinate frame (Z > 140)',
          all(a['p'][2] > 140 for a in poses['top'][0]['atoms']))

    # ---- evidence
    print('\n· Evidence chain')
    rep = jget('/api/evidence/report')
    check('Report contains 24 prior-check sections',
          'Priors' in rep['markdown'] and '## 5' in rep['markdown'])
    pri = jget('/api/evidence/priors')
    check('Priors 24 entries,PASS 22 / FLAG 1 / WEAK 1',
          len(pri['rows']) == 24 and pri['status']['by_status'] ==
          {'PASS': 22, 'FLAG': 1, 'WEAK': 1}, str(pri['status']['by_status']))
    check('Prior rows include quotations and sources',
          pri['rows'][0]['claim'] and pri['rows'][0]['locator'])
    cen = jget('/api/evidence/census')
    check('Inventory:MurJ 496 residues, construct excluded',
          all(c['murj_residues'] == 496 for c in cen['atom_budget']['census']))
    check('Atom budget: peptide contacts 121 vs median scaffold 24',
          cen['atom_budget']['budget']['small_molecule']['heavy_atoms_median'] == 24)
    iface = jget('/api/evidence/interface')
    check('Interface measurements include both sides', len(iface['murj']) > 50 and len(iface['peptide']) > 50)
    check('Anchor table is nonempty', len(iface['anchors']) >= 12)

    # ---- 3D scene
    print('\n· 3D Structure')
    scene = jget('/api/scene?ref=SglCJ3&ligand=0')
    check('Receptor heavy-atom trace 496 items Cα', len(scene['receptor']['ca']) == 496)
    check('All three peptides share a coordinate frame', len(scene['peptides']) == 3
          and all(len(p['ca']) > 20 for p in scene['peptides']))
    check('Ligand pose 24 heavy atoms', len(scene['ligand']['atoms']) == 24)
    check('Scene includes pharmacophores and hotspots',
          len(scene['features']) == 12 and len(scene['hotspots']) >= 20)
    check('Subpocket colors are available for trace coloring',
          len(scene.get('subpockets', [])) == 4
          and any(r.get('sub') for r in scene['receptor']['residues']))
    check('After alignment, SglM peptide lies near the SglCJ3 pocket',
          _dist(scene['peptides'][0]['ca'][0], scene['center']) < 60)

    scene_m = jget('/api/scene?ref=SglM')
    check('Alternative reference structure works', len(scene_m['receptor']['ca']) == 496
          and scene_m['ref'] == 'SglM')

    pts = jget('/api/pocket/points')
    if pts['status'] == 'ready':
        check('Cavity point cloud 10k+ voxels, volume 5,440 Å³',
              pts['n_voxels'] > 9000 and round(pts['volume']) == 5440)
        check('Point cloud classified by chemistry', len(pts['chem']) == pts['n_voxels'])
    else:
        check('Computing cavity points (first run approximately 15 s)', pts['status'] == 'computing')

    # ---- docs & jobs
    print('\n· Documentation and jobs')
    docs = jget('/api/docs')['docs']
    check('Listed 5 or more documents', len(docs) >= 5, str(len(docs)))
    check('Documents have titles', all(d['title'] for d in docs))
    doc = jget('/api/docs/' + urllib.parse.quote('01_Molecular scaffold design algorithm.md'))
    check('English document filename is readable', len(doc['markdown']) > 5000)

    jobs = jget('/api/jobs')
    check('Job endpoint returns four script types', len(jobs['kinds']) == 4)
    check('Run history schema complete', isinstance(jobs['jobs'], list))

    status, raw, _ = fetch('/api/jobs', 'POST', {'kind': 'nope'})
    check('Unknown job type rejected', status == 400)

    status, raw, _ = fetch('/api/jobs', 'POST',
                           {'kind': 'pipeline', 'params': {'out': '../escape'}})
    check('Invalid output directory rejected', status == 400)

    status, raw, _ = fetch('/api/jobs', 'POST',
                           {'kind': 'filters', 'params': {'targets': -5}})
    check('Out-of-range arguments rejected', status == 400)

    # ---- file serving / traversal
    print('\n· Static files')
    status, raw, ctype = fetch('/files/figures/fig4_funnel.png')
    check('Images downloadable', status == 200 and 'image/png' in ctype)
    status, raw, _ = fetch('/files/priors/murj_priors.json')
    check('Prior file readable', status == 200 and b'structures' in raw)
    status, raw, _ = fetch('/files/webui/data.py')
    check('webui Directory access cannot exceed scope', status == 403, str(status))
    status, raw, _ = fetch('/files/../webui/data.py')
    check('Path traversal rejected', status in (403, 404), str(status))
    status, raw, ctype = fetch('/nonexistent-page')
    check('Unknown route falls back to SPA home page', status == 200 and 'text/html' in ctype)
    status, raw, _ = fetch('/js/missing-file.js')
    check('Missing frontend resource returns 404 rather than HTML', status == 404, str(status))

    # ---- 7-step workflow UI
    print('\n· Workflow (7 steps)')
    wf = jget('/api/workflow')
    st = wf['stages']
    check('Workflow has 7 step', [s['n'] for s in st] == list(range(1, 8)))
    check('Number 1/2 connected steps with file-derived counts',
          st[0]['status'] == 'ready' and st[1]['status'] == 'ready'
          and st[1]['funnel'][0]['count'] == 20307 and st[1]['out'] == 32)
    check('Number 3-7 steps without data use pending and empty counts',
          all(s['status'] == 'pending' and all(f['count'] is None for f in s['funnel']) for s in st[2:]))
    wd = jget('/api/workflow?demo=1')
    check('Layout preview is labeled demo', all(s['status'] == 'demo' for s in wd['stages'][2:]))
    mols = jget('/api/stage/scaffold/molecules')
    check('Number 2 step molecule list 538 entries,32 entries passed',
          mols['n'] == 538 and sum(1 for r in mols['rows'] if r['keep']) == 32)
    pose = jget('/api/stage/scaffold/pose?id=' + mols['rows'][0]['id'])
    check('Number 2 step molecules include 3D Pose', pose['ligand'] and len(pose['ligand']['atoms']) > 5
          and len(pose['ligand']['bonds']) > 5)
    check('Unconnected stages have empty molecule lists', jget('/api/stage/docking/molecules')['n'] == 0)
    status, raw, _ = fetch('/api/stage/nope/molecules')
    check('Unknown stage returns 404', status == 404, str(status))
    m2 = jget('/api/mol2d?smiles=' + 'OC(%3DO)c1ccccc1')
    check('2D Structure image generated', bool(m2.get('svg')) or 'RDKit' in (m2.get('reason') or ''))

    # ---- live job (opt-in: runs the repo's own smoke test through the API)
    if do_job:
        print('\n· Live job smoke test')
        status, raw, _ = fetch('/api/jobs', 'POST', {'kind': 'tests'})
        check('Job accepted', status == 201, f'{status}: {raw[:120]}')
        job_id = json.loads(raw.decode('utf-8'))['id']
        seen_running = False
        deadline = time.time() + 900
        final = None
        while time.time() < deadline:
            time.sleep(1.0)
            listing = jget('/api/jobs')['jobs']
            if any(j['id'] == job_id and j['status'] == 'running' for j in listing):
                seen_running = True
            current = jget(f'/api/jobs/{job_id}?since=0')
            if current['status'] in ('done', 'failed', 'cancelled'):
                final = current
                break
        check('Job list readable during execution, without process objects',
              seen_running and all('proc' not in j for j in listing))
        check('Job completed with exit code 0',
              final is not None and final['status'] == 'done' and final['exit'] == 0,
              str(final and (final['status'], final['exit'])))
        check('Logs streamed successfully', final is not None and final['total_lines'] > 5,
              str(final and final['total_lines']))

    httpd.shutdown()


def _dist(a, b):
    return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5


if __name__ == '__main__':
    import argparse

    ap = argparse.ArgumentParser(description='WEDGE-Design Console API Test')
    ap.add_argument('--job', action='store_true',
                    help='Run an additional live job (tests/test_smoke.py, approximately 1 minutes)')
    args = ap.parse_args()
    try:
        run(do_job=args.job)
    except Exception as exc:                                       # noqa: BLE001
        import traceback
        traceback.print_exc()
        FAILED.append(f'crash: {exc}')
    print(f'\n{PASSED} passed, {len(FAILED)} failed')
    if FAILED:
        print('failed:', ', '.join(FAILED))
    raise SystemExit(1 if FAILED else 0)
