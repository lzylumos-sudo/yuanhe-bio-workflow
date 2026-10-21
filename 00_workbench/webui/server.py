# -*- coding: utf-8 -*-
"""webui.server -- dependency-free HTTP server for the WEDGE-Design web console.

    python run_webui.py                # http://127.0.0.1:8777
    python run_webui.py --port 9000 --no-browser

Serves the single-page front-end from ``webui/static`` and a small JSON API over
the files that the pipeline writes (``results/``, ``evidence/``).
Only the Python standard library is used, so the console runs on the same
laptop CPU as the pipeline itself.
"""
from __future__ import annotations

import json
from pathlib import Path
import mimetypes
import os
import posixpath
import re
import socket
import sys
import threading
import time
import urllib.parse
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from webui import data, jobs, stages, workflow, studio                                   # noqa: E402

STATIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static')
ROOT = data.ROOT
FILE_DIRS = ('results', 'evidence', 'figures', 'docs', 'priors')
VERSION = '0.4.0-rc1'

MIME = {
    '.html': 'text/html; charset=utf-8', '.css': 'text/css; charset=utf-8',
    '.js': 'text/javascript; charset=utf-8', '.mjs': 'text/javascript; charset=utf-8',
    '.json': 'application/json; charset=utf-8', '.md': 'text/markdown; charset=utf-8',
    '.png': 'image/png', '.jpg': 'image/jpeg', '.svg': 'image/svg+xml',
    '.ico': 'image/x-icon', '.csv': 'text/csv; charset=utf-8',
    '.sdf': 'chemical/x-mdl-sdfile', '.pdb': 'chemical/x-pdb',
    '.woff2': 'font/woff2', '.txt': 'text/plain; charset=utf-8',
}


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _json(payload, status=200):
    body = json.dumps(payload, ensure_ascii=False, default=str).encode('utf-8')
    return status, 'application/json; charset=utf-8', body


def _error(message, status=400, **extra):
    return _json(dict(error=message, status=status, **extra), status)


class HttpError(Exception):
    def __init__(self, message, status=400):
        super().__init__(message)
        self.status = status


# ---------------------------------------------------------------------------
# API
# ---------------------------------------------------------------------------
def api_health(q):
    return _json(dict(ok=True, version=VERSION, time=time.time(), jobs=len(jobs.list_jobs(999))))


def api_meta(q):
    return _json(data.meta())


def api_overview(q):
    summary = data.summary() or {}
    return _json(dict(
        summary=summary,
        funnel=summary.get('funnel', {}),
        filters=data.filter_summary(),
        evidence=data.evidence_status(),
        n_scaffolds=summary.get('n_scaffolds'),
        best=summary.get('best', []),
        features=summary.get('features', []),
        axes=data.axes(),
        hotspots=(data.hotpsots_raw() or [])[:12],
        poses=len(data.poses_top()),
        pocket=data.pocket(),
        meta=data.meta(),
    ))


def api_hotspots(q):
    rows = data.hotpsots_raw()
    return _json(dict(rows=rows, n=len(rows),
                      n_common=sum(1 for r in rows if float(r.get('convergence') or 0) >= 1.0)))


def api_pharmacophore(q):
    ph = data.pharmacophore() or {}
    return _json(dict(features=ph.get('features', []), pairs=ph.get('pairs', []),
                      axes=data.axes(), pocket=data.pocket()))


def api_pocket(q):
    return _json(dict(pocket=data.pocket(), status=data.pocket_status()))


def api_pocket_points(q):
    step = float(q.get('step', ['0.8'])[0] or 0.8)
    step = max(0.5, min(2.0, step))
    payload = data.pocket_points(step)
    return _json(payload, 200 if payload['status'] == 'ready' else 202)


def api_axes(q):
    return _json(dict(axes=data.axes()))


def api_scaffolds(q):
    rows = data.scaffold_rows()
    return _json(dict(rows=rows, n=len(rows), filters=data.filter_summary(),
                      poses=len(data.poses_top())))


def api_filters(q):
    return _json(dict(rows=data.filter_rows(), summary=data.filter_summary(),
                      n=len(data.filter_rows())))


def api_poses(q):
    return _json(dict(top=data.poses_top(), diverse=data.poses_diverse()))


def api_scene(q):
    ref = (q.get('ref', [data.REF_TAG])[0] or data.REF_TAG)
    ligand = q.get('ligand', [None])[0]
    radius = float(q.get('radius', ['13'])[0] or 13)
    radius = max(6.0, min(30.0, radius))
    return _json(data.scene(ref=ref, ligand=ligand, radius=radius))


def api_evidence_report(q):
    return _json(dict(markdown=data.evidence_report(),
                      status=data.evidence_status(),
                      generated=data._mtime(data.rel('evidence/report.md')),  # noqa: SLF001
                      atom_budget=data.atom_budget()))


def api_evidence_priors(q):
    return _json(dict(rows=data.priors_check(), status=data.evidence_status(),
                      resolved=data.priors_resolved()))


def api_evidence_census(q):
    return _json(dict(census=data.census(), atom_budget=data.atom_budget(),
                      status=data.evidence_status()))


def api_evidence_interface(q):
    return _json(dict(murj=data.interface_murj(), peptide=data.interface_peptide(),
                      anchors=data.anchors(),
                      hotspots=data.hotpsots_raw()))


def api_evidence_anchors(q):
    return _json(dict(anchors=data.anchors()))


def api_docs(q):
    return _json(dict(docs=data.docs_list()))


def api_doc(q, name):
    text = data.doc(name)
    if not text:
        raise HttpError(f'Document not found: {name}', 404)
    return _json(dict(name=name, markdown=text))


def api_job_kinds(q):
    return _json(dict(kinds=jobs.list_kinds()))


def api_jobs_get(q):
    return _json(dict(jobs=jobs.list_jobs(40), kinds=jobs.list_kinds()))


def api_job_one(q, job_id):
    since = int(q.get('since', ['0'])[0] or 0)
    job = jobs.get_job(job_id, since=since)
    if job is None:
        raise HttpError(f'Job not found: {job_id}', 404)
    return _json(job)


def api_jobs_start(q, body):
    kind = (body or {}).get('kind')
    payload = (body or {}).get('params') or {}
    try:
        job_id = jobs.start(kind, payload)
    except ValueError as exc:
        raise HttpError(str(exc), 400) from exc
    except RuntimeError as exc:
        raise HttpError(str(exc), 409) from exc
    return _json(dict(id=job_id), 201)


def api_job_cancel(q, body, job_id):
    try:
        ok = jobs.cancel(job_id)
    except KeyError as exc:
        raise HttpError(f'Job not found: {job_id}', 404) from exc
    return _json(dict(cancelled=ok, id=job_id))


def _demo(q):
    return (q.get('demo', ['0'])[0] or '0') in ('1', 'true')


def api_workflow(q):
    return _json(stages.workflow(_demo(q)))


def _stage_key(key):
    if key not in stages.BY_KEY:
        raise HttpError(f'Unknown stage: {key}', 404)
    return key


def api_stage_molecules(q, key):
    return _json(stages.molecules(_stage_key(key), _demo(q)))


def api_stage_pose(q, key):
    return _json(stages.pose(_stage_key(key), q.get('id', [''])[0], _demo(q)))


def api_mol2d(q):
    return _json(stages.mol2d(q.get('smiles', [''])[0][:400]))


# ---------------------------------------------------------------------------
# routing table: (method, regex, handler)
# ---------------------------------------------------------------------------

def api_project(q):
    return _json(workflow.inventory())

def api_research_preview(q,body):
    from research.config import resolve,plan
    from research.resources import detect
    rc=resolve(body.get('research',body))
    return _json({'resolved':rc,'plan':plan(rc),'resources':detect(workflow.ROOT)})

def api_research_derive(q,body,rid):
    from research.lifecycle import derive
    return _json(derive(rid,body),201)

def api_research_handoff(q,body,rid):
    from research.lifecycle import handoff
    return _json(handoff(rid,body))

def api_project_create(q, body):
    return _json(workflow.create(body), 201)

def api_project_run(q, rid):
    return _json(workflow.status(rid))

def api_project_start(q, body, rid, key):
    return _json(workflow.start(rid,key), 202)

def api_project_review(q, body, rid, key):
    return _json(workflow.review(rid,key,body))

def api_project_rescreen(q, body, rid):
    return _json(workflow.rescreen(rid),202)

def api_project_cancel(q, body, rid):
    return _json(workflow.cancel(rid))

def api_project_log(q, rid, key):
    return _json(workflow.log(rid,key))

def api_project_rows(q, rid, key):
    import csv
    if key not in [x[0] for x in workflow.STAGES]:
        raise HttpError('Unknown stage',404)
    path=workflow.run_dir(rid)/key/'molecules.csv'
    if not path.exists():return _json({'rows':[],'total':0})
    with path.open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    if key=='screening':
        from scientific.validation import effective_candidates
        rows=effective_candidates(workflow.run_dir(rid))
    if key in ('risk','availability'):
        rows=workflow.read(workflow.run_dir(rid)/key/'records.json',[])
    query=q.get('search',[''])[0].lower()
    if query:rows=[r for r in rows if query in str(r).lower()]
    offset=max(0,int(q.get('offset',['0'])[0]))
    return _json({'rows':rows[offset:offset+100],'total':len(rows),'offset':offset})

def api_project_pose(q, rid, key):
    from rdkit import Chem
    d=workflow.run_dir(rid);ident=q.get('id',[''])[0]
    path=None;caption='Molecular coordinates; not aligned to the protein coordinate frame'
    if key=='peptide2mol':
        rows=workflow.read(d/key/'seeds.json',[])
        row=next((x for x in rows if x['id']==ident),None)
        if row:path=d/row['raw_sdf']
    elif key in ('risk','availability'):
        rows=workflow.read(d/key/'records.json',[])
        row=next((x for x in rows if x.get('id')==ident),None)
        if row and row.get('source_sdf'):path=Path(row['source_sdf'])
    elif key in ('reinvent4','docking'):
        rows=workflow.read(d/'reinvent4/candidates.json',[])
        row=next((x for x in rows if x['id']==ident),None)
        if row:path=d/row['sdf'] if key=='reinvent4' else d/'docking/jobs'/ident/'poses.sdf'
    elif key=='screening':
        rows=workflow.read(d/'screening/candidates.json',[])
        row=next((x for x in rows if x.get('pose_id')==ident or x.get('compound_id')==ident),None)
        if row and row.get('sdf'):path=d/row['sdf']
    if path is None or not path.resolve().is_relative_to(d.resolve()) or not path.is_file():
        return _json({'atoms':[],'caption':'No readable 3D artifact for this record'})
    mol=next((m for m in Chem.SDMolSupplier(str(path),removeHs=True) if m is not None),None)
    if mol is None or not mol.GetNumConformers():return _json({'atoms':[],'caption':'SDF Unreadable coordinates'})
    atoms=[{'p':list(mol.GetConformer().GetAtomPosition(a.GetIdx())),'chain':'LIG','name':'CA','element':a.GetSymbol()} for a in mol.GetAtoms()]
    if key in ('docking','screening'):
        caption='GNINA Original first pose and actual docking receptor; not automatically the best pose'
        if key=='screening':caption='Screening-selected pose and actual docking receptor; see candidate table for qualification'
        receptor=d/'docking/receptor.pdb'
        for line in receptor.read_text().splitlines():
            if line.startswith('ATOM  '):
                xyz=[float(line[30:38]),float(line[38:46]),float(line[46:54])]
                if any(sum((x-y)**2 for x,y in zip(xyz,a['p']))<100 for a in atoms if a['chain']=='LIG'):
                    atoms.append({'p':xyz,'chain':'REC','name':line[12:16].strip()})
    return _json({'atoms':atoms,'caption':caption})

def api_predocking_revise(q,body,rid,key):
    return _json(workflow.revise_preflight(rid,key))

def api_predocking_detail(q,rid,key):
    if key not in ('risk','availability'):raise HttpError('Unknown evidence stage',404)
    d=workflow.run_dir(rid);ident=q.get('id',[''])[0]
    rows=workflow.read(d/key/'records.json',[])
    row=next((r for r in rows if r.get('id')==ident),None)
    if row is None:raise HttpError('Molecule not found',404)
    svg=d/'risk'/(row['id']+'.svg')
    return _json({'record':row,'svg':svg.read_text() if svg.is_file() else '',
                  'note':'No match to current rules does not establish non-toxicity; a computed route is not completed synthesis. Final purchasing or commissioned synthesis requires separate review.'})

def api_project_structure(q, rid):
    d=workflow.run_dir(rid);path=d/'input.pdb'
    atoms=[]
    if path.exists():
        for line in path.read_text().splitlines():
            if line.startswith('ATOM  '):
                try:atoms.append({'p':[float(line[30:38]),float(line[38:46]),float(line[46:54])],'chain':line[21:22], 'name':line[12:16].strip(),'residue':line[17:27].strip()})
                except ValueError:pass
    return _json({'atoms':atoms})

def api_studio(q,body,action):
    if action=='upload':return _json(studio.upload(body),201)
    if action=='environment':return _json(studio.environment_save(body))
    raise ValueError('Unknown operation')

def api_studio_run(q,body,rid,action):
    if action=='continue-md':return _json(workflow.continue_md(rid))
    if action=='export':return _json(studio.export_run(rid,body))
    if action=='core':return _json(studio.save_core(rid,body))
    if action=='box':return _json(studio.box(rid,body))
    if action=='membrane-preview':return _json(studio.membrane_preview(rid,body))
    if action=='receptor-apply':return _json(studio.apply_receptor(rid,body))
    return _json(studio.launch(rid,action,body),202)

def api_studio_get(q,rid,action):
    if action=='operations':return _json({'operations':studio.operations(rid)})
    if action=='seeds':return _json({'rows':workflow.read(workflow.run_dir(rid)/'peptide2mol/seeds.json',[])})
    if action=='core':return _json(studio.core(rid,q.get('id',[''])[0]))
    if action=='md':return _json({'systems':studio.md_list(rid)})
    if action=='md-structure':return _json(studio.md_structure(rid,q.get('id',[''])[0]))
    if action=='trajectory':return _json(studio.md_data(rid,q.get('id',[''])[0]))
    if action=='candidates':
        import csv
        path=workflow.run_dir(rid)/'screening/molecules.csv'
        if not path.exists():return _json({'rows':[]})
        workflow.require(workflow.run_dir(rid),'screening')
        orientation=workflow.read(workflow.run_dir(rid)/'screening/workspace/runs/screen/stages/membrane/summary.json',{})
        from scientific.validation import effective_candidates
        return _json({'rows':effective_candidates(workflow.run_dir(rid)), 'orientation':orientation.get('plane') if orientation.get('orientation_verified') else None})
    raise ValueError('Unknown read operation')

ROUTES = [
    ('POST', r'^/api/studio/(?P<action>upload|environment)$', api_studio),
    ('POST', r'^/api/studio/run/(?P<rid>[A-Za-z0-9_-]+)/(?P<action>[a-z_-]+)$', api_studio_run),
    ('GET', r'^/api/studio/run/(?P<rid>[A-Za-z0-9_-]+)/(?P<action>[a-z_-]+)$', api_studio_get),
    ('GET', r'^/api/project$', api_project),
    ('POST', r'^/api/project/runs$', api_project_create),
    ('POST', r'^/api/project/run/(?P<rid>[A-Za-z0-9_-]+)/rescreen$', api_project_rescreen),
    ('GET', r'^/api/project/run/(?P<rid>[A-Za-z0-9_-]+)$', api_project_run),
    ('POST', r'^/api/project/run/(?P<rid>[A-Za-z0-9_-]+)/predocking/(?P<key>risk|availability)/revise$', api_predocking_revise),
    ('GET', r'^/api/project/run/(?P<rid>[A-Za-z0-9_-]+)/predocking/(?P<key>risk|availability)$', api_predocking_detail),
    ('GET', r'^/api/project/run/(?P<rid>[A-Za-z0-9_-]+)/structure$', api_project_structure),
    ('GET', r'^/api/project/run/(?P<rid>[A-Za-z0-9_-]+)/pose/(?P<key>[a-z0-9]+)$', api_project_pose),
    ('GET', r'^/api/project/run/(?P<rid>[A-Za-z0-9_-]+)/rows/(?P<key>[a-z0-9]+)$', api_project_rows),
    ('GET', r'^/api/project/run/(?P<rid>[A-Za-z0-9_-]+)/log/(?P<key>[a-z0-9]+)$', api_project_log),
    ('POST', r'^/api/project/run/(?P<rid>[A-Za-z0-9_-]+)/start/(?P<key>[a-z0-9]+)$', api_project_start),
    ('POST', r'^/api/project/run/(?P<rid>[A-Za-z0-9_-]+)/cancel$', api_project_cancel),
    ('POST', r'^/api/project/run/(?P<rid>[A-Za-z0-9_-]+)/review/(?P<key>[a-z0-9]+)$', api_project_review),
    ('GET', r'^/api/health$', api_health),
    ('GET', r'^/api/meta$', api_meta),
    ('GET', r'^/api/overview$', api_overview),
    ('GET', r'^/api/workflow$', api_workflow),
    ('POST', r'^/api/research/preview$', api_research_preview),
    ('POST', r'^/api/research/run/(?P<rid>[^/]+)/derive$', api_research_derive),
    ('POST', r'^/api/research/run/(?P<rid>[^/]+)/handoff$', api_research_handoff),
    ('GET', r'^/api/stage/(?P<key>[a-z_]+)/molecules$', api_stage_molecules),
    ('GET', r'^/api/stage/(?P<key>[a-z_]+)/pose$', api_stage_pose),
    ('GET', r'^/api/mol2d$', api_mol2d),
    ('GET', r'^/api/hotspots$', api_hotspots),
    ('GET', r'^/api/pharmacophore$', api_pharmacophore),
    ('GET', r'^/api/pocket$', api_pocket),
    ('GET', r'^/api/pocket/points$', api_pocket_points),
    ('GET', r'^/api/axes$', api_axes),
    ('GET', r'^/api/scaffolds$', api_scaffolds),
    ('GET', r'^/api/filters$', api_filters),
    ('GET', r'^/api/poses$', api_poses),
    ('GET', r'^/api/scene$', api_scene),
    ('GET', r'^/api/evidence/report$', api_evidence_report),
    ('GET', r'^/api/evidence/priors$', api_evidence_priors),
    ('GET', r'^/api/evidence/census$', api_evidence_census),
    ('GET', r'^/api/evidence/interface$', api_evidence_interface),
    ('GET', r'^/api/evidence/anchors$', api_evidence_anchors),
    ('GET', r'^/api/docs$', api_docs),
    ('GET', r'^/api/docs/(?P<name>[^/]+)$', api_doc),
    ('GET', r'^/api/jobs$', api_jobs_get),
    ('GET', r'^/api/jobs/kinds$', api_job_kinds),
    ('GET', r'^/api/jobs/(?P<job_id>[0-9a-f]+)$', api_job_one),
    ('POST', r'^/api/jobs$', api_jobs_start),
    ('POST', r'^/api/jobs/(?P<job_id>[0-9a-f]+)/cancel$', api_job_cancel),
]
COMPILED = [(m, re.compile(p), h) for m, p, h in ROUTES]


def dispatch(method, path, query, body):
    for m, rx, handler in COMPILED:
        if m != method:
            continue
        match = rx.match(path)
        if not match:
            continue
        kwargs = match.groupdict()
        try:
            if method == 'POST':
                return handler(query, body, **kwargs) if kwargs else handler(query, body)
            return handler(query, **kwargs) if kwargs else handler(query)
        except HttpError:
            raise
        except (ValueError, KeyError) as exc:
            return _error(f'{type(exc).__name__}: {exc}', 400)
        except OSError as exc:
            return _error(f'File read failed: {exc}', 500)
        except Exception as exc:                            # noqa: BLE001
            import traceback
            traceback.print_exc()
            return _error(f'{type(exc).__name__}: {exc}', 500)
    raise HttpError(f'Unknown endpoint: {method} {path}', 404)


# ---------------------------------------------------------------------------
# static files
# ---------------------------------------------------------------------------
def _safe_join(base, rel_url):
    rel_url = urllib.parse.unquote(rel_url).split('?')[0].split('#')[0]
    rel_url = rel_url.replace('\\', '/')
    norm = posixpath.normpath('/' + rel_url).lstrip('/')
    full = os.path.realpath(os.path.join(base, norm))
    if os.path.commonpath([full, os.path.realpath(base)]) != os.path.realpath(base):
        raise HttpError('Invalid path', 403)
    return full


class Handler(BaseHTTPRequestHandler):
    server_version = f'WEDGE-Design-UI/{VERSION}'
    protocol_version = 'HTTP/1.1'

    def _valid_host(self):
        host = self.headers.get('Host', '').split(':')[0].lower()
        return host in {'localhost', '127.0.0.1', self.server.server_address[0]}

    def log_message(self, fmt, *args):                       # quieter log
        if os.environ.get('WEDGE_UI_VERBOSE'):
            sys.stderr.write('%s - %s\n' % (self.address_string(), fmt % args))

    def _send(self, status, ctype, body, extra=None):
        self.send_response(status)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        for key, value in (extra or {}).items():
            self.send_header(key, value)
        self.end_headers()
        if self.command != 'HEAD':
            try:
                self.wfile.write(body)
            except (BrokenPipeError, ConnectionResetError):
                pass

    def do_OPTIONS(self):                                    # noqa: N802
        self.send_response(204)
        self.send_header('Allow', 'GET, POST, HEAD, OPTIONS')
        self.send_header('Content-Length', '0')
        self.end_headers()

    def do_HEAD(self):                                       # noqa: N802
        self.do_GET()

    def do_GET(self):                                        # noqa: N802
        if not self._valid_host():
            self._send(*_error('Host mismatch', 403)); return
        parsed = urllib.parse.urlsplit(self.path)
        path = urllib.parse.unquote(parsed.path)
        query = urllib.parse.parse_qs(parsed.query)
        try:
            if path.startswith('/api/'):
                status, ctype, body = dispatch('GET', path, query, None)
                self._send(status, ctype, body)
                return
            if path.startswith('/artifacts/'):
                full = _safe_join(str(workflow.RUNS), path[len('/artifacts/'):])
                if not os.path.isfile(full): raise HttpError('Artifact not found',404)
                self._send_file(full, {'Content-Disposition': 'attachment', 'X-Content-Type-Options': 'nosniff'})
                return
            if path.startswith('/files/'):
                self._serve_data_file(path[len('/files/'):], query)
                return
            self._serve_static(path)
        except HttpError as exc:
            status, ctype, body = _error(str(exc), exc.status)
            self._send(status, ctype, body)
        except Exception as exc:                             # noqa: BLE001
            import traceback
            traceback.print_exc()
            status, ctype, body = _error(f'{type(exc).__name__}: {exc}', 500)
            self._send(status, ctype, body)

    def do_POST(self):                                       # noqa: N802
        if not self._valid_host():
            self.close_connection = True
            self._send(*_error('Host mismatch', 403)); return
        parsed = urllib.parse.urlsplit(self.path)
        query = urllib.parse.parse_qs(parsed.query)
        origin = self.headers.get('Origin')
        host = self.headers.get('Host', '')
        if origin and origin != 'http://' + host:
            self.close_connection = True
            self._send(*_error('Cross-site write request denied', 403)); return
        if self.headers.get('Sec-Fetch-Site') == 'cross-site':
            self.close_connection = True
            self._send(*_error('Cross-site write request denied', 403)); return
        try:
            length = int(self.headers.get('Content-Length') or 0)
            if not 0 <= length <= ((12 if parsed.path == '/api/studio/upload' else 1) * 1024 * 1024):
                raise ValueError()
        except ValueError:
            self.close_connection = True
            self._send(*_error('Invalid request body length or limit exceeded: 12 MB', 413)); return
        raw = self.rfile.read(length) if length else b'{}'
        try:
            body = json.loads(raw.decode('utf-8') or '{}')
            if not isinstance(body, dict):
                raise ValueError('JSON body must be an object')
        except ValueError:
            self._send(*_error('Request body is not valid JSON', 400))
            return
        try:
            status, ctype, payload = dispatch('POST', parsed.path, query, body)
            self._send(status, ctype, payload)
        except HttpError as exc:
            self._send(*_error(str(exc), exc.status))
        except Exception as exc:                             # noqa: BLE001
            import traceback
            traceback.print_exc()
            self._send(*_error(f'{type(exc).__name__}: {exc}', 500))

    # -- static ------------------------------------------------------------
    def _serve_static(self, path):
        if path in ('', '/'):
            path = '/index.html'
        full = _safe_join(STATIC, path)
        if not os.path.isfile(full):
            # a missing asset is an error (a 200 with index.html would surface as
            # "Unexpected token '<'" in the browser); only extension-less paths
            # fall back to the SPA shell.
            if os.path.splitext(path)[1]:
                raise HttpError(f'Frontend resource missing: {path}', 404)
            full = os.path.join(STATIC, 'index.html')
            if not os.path.isfile(full):
                raise HttpError('Frontend runtime files missing webui/static/index.html', 500)
        self._send_file(full)

    def _serve_data_file(self, rel_path, query):
        rel_path = urllib.parse.unquote(rel_path)
        top = rel_path.split('/')[0]
        if top not in FILE_DIRS:
            raise HttpError(f'Downloads are restricted to files under {"/".join(FILE_DIRS)} this directory', 403)
        full = _safe_join(os.path.join(ROOT, top), rel_path[len(top) + 1:])
        if not os.path.isfile(full):
            raise HttpError(f'File not found: {rel_path}', 404)
        name = os.path.basename(full)
        extra = {}
        if query.get('download'):
            extra['Content-Disposition'] = (
                "attachment; filename*=UTF-8''" + urllib.parse.quote(name))
        self._send_file(full, extra)

    def _send_file(self, full, extra=None):
        stat = os.stat(full)
        ext = os.path.splitext(full)[1].lower()
        ctype = MIME.get(ext) or mimetypes.guess_type(full)[0] or 'application/octet-stream'
        self.send_response(200)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(stat.st_size))
        self.send_header('Last-Modified', self.date_time_string(stat.st_mtime))
        self.send_header('Cache-Control', 'no-cache')
        for key, value in (extra or {}).items():
            self.send_header(key, value)
        self.end_headers()
        if self.command == 'HEAD':
            return
        with open(full, 'rb') as fh:
            try:
                while True:
                    chunk = fh.read(262144)
                    if not chunk:
                        break
                    self.wfile.write(chunk)
            except (BrokenPipeError, ConnectionResetError):
                pass


class Server(ThreadingHTTPServer):
    daemon_threads = True
    allow_reuse_address = True

    def handle_error(self, request, client_address):
        # Browser refresh,/tab closure and preconnection probes can terminate the RST connection; these interruptions are not errors;
        # do not print traceback(to avoid flooding the console with WinError 10054 messages).
        exc = sys.exc_info()[1]
        if isinstance(exc, (ConnectionResetError, ConnectionAbortedError, BrokenPipeError)):
            return
        super().handle_error(request, client_address)


def _pick_port(host, port, tries=20):
    for offset in range(tries):
        candidate = port + offset
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            try:
                sock.bind((host, candidate))
            except OSError:
                continue
        return candidate
    raise OSError(f'{host}:{port}-{port + tries} all occupied')


def serve(host='127.0.0.1', port=8777, open_browser=True, quiet=False):
    if not os.path.isdir(STATIC):
        print(f'[webui] Frontend directory missing: {STATIC}', file=sys.stderr)
        return 2
    port = _pick_port(host, port)
    httpd = Server((host, port), Handler)
    url = f'http://{host}:{port}/'
    if not quiet:
        print('=' * 62)
        print(f'  Yuanhe Bio Research Studio · {VERSION}')
        print(f'  {url}')
        print(f'  Working directory {ROOT}')
        print(f'  Python {sys.version.split()[0]}   Scientific tasks use separate environments')
        print('  Ctrl+C End')
        print('=' * 62, flush=True)
    if open_browser:
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever(poll_interval=0.3)
    except KeyboardInterrupt:
        print('\n[webui] Stopped')
    finally:
        httpd.server_close()
    return 0
