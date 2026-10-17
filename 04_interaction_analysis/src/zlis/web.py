"""Local web UI (standard library only).

Security: binds to 127.0.0.1 only and requires a per-launch random token (in the
URL printed at start-up). There is no user account system; do not expose this
port to a network.
"""
from __future__ import annotations

import cgi
import json
import secrets
import tempfile
import threading
import traceback
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from .pipeline import PROJECT_ROOT

RESULTS = PROJECT_ROOT / "results"
TOKEN = secrets.token_urlsafe(16)

PAGE = """<!doctype html><html lang="zh"><head><meta charset="utf-8"><title>ZLIS</title>
<style>body{font-family:Helvetica,Arial,sans-serif;max-width:900px;margin:24px auto;padding:0 12px}
fieldset{margin:12px 0}label{display:block;margin:6px 0}pre{background:#f4f4f4;padding:8px;overflow:auto}
button{padding:6px 14px}.ok{color:#176f2c}.bad{color:#b00020}</style></head><body>
<h2>Zhiyu Ligand Interaction Studio</h2>
<p>Local service, listening only on 127.0.0.1, Access requires a random token; there is no account system. Do not expose publicly.</p>
<form id="f" method="post" action="/run?token=__TOKEN__" enctype="multipart/form-data">
<fieldset><legend>Input</legend>
<label>Complex PDB file <input type="file" name="pdb" accept=".pdb,.ent" required></label>
<label>Ligand (e.g. A:BEN1 or BEN, Leave empty to analyze all candidates) <input name="ligand"></label>
<label>Ligand template SDF/SMILES file (Optional; chemistry only, coordinates unchanged) <input type="file" name="template"></label>
<label>Net charge (Optional without a template) <input name="charge" size="4"></label>
<label><input type="checkbox" name="offline" value="1"> Offline (Do not access RCSB)</label>
<label><input type="checkbox" name="nopng" value="1"> Skip PyMOL Ray tracing PNG (Faster)</label>
</fieldset><button>Run analysis</button></form>
<h3>Existing results</h3><ul>__LIST__</ul>
</body></html>"""


def _list_results() -> str:
    if not RESULTS.exists():
        return "<li>(none)</li>"
    items = []
    for d in sorted(p for p in RESULTS.iterdir() if p.is_dir() and not p.name.startswith(".")):
        links = " ".join(f'<a href="/files/{d.name}/{f}?token={TOKEN}" target="_blank">{f}</a>'
                         for f in ("interaction_2d_editor.html", "interaction_2d_audit.html", "interactions.xlsx",
                                   "binding_site_main.png") if (d / f).exists())
        items.append(f"<li>{d.name}: {links}</li>")
    return "".join(items) or "<li>(none)</li>"


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *a):  # quiet
        pass

    def _ok(self, ctype, data: bytes, code=200):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _auth(self) -> bool:
        q = parse_qs(urlparse(self.path).query)
        if q.get("token", [""])[0] != TOKEN:
            self._ok("text/plain; charset=utf-8", b"403 token required", 403)
            return False
        return True

    def do_GET(self):
        if not self._auth():
            return
        u = urlparse(self.path)
        if u.path == "/":
            html = PAGE.replace("__TOKEN__", TOKEN).replace("__LIST__", _list_results())
            return self._ok("text/html; charset=utf-8", html.encode())
        if u.path.startswith("/files/"):
            rel = Path(u.path[len("/files/"):])
            target = (RESULTS / rel).resolve()
            if RESULTS.resolve() not in target.parents or not target.is_file():
                return self._ok("text/plain", b"not found", 404)
            ctype = {".html": "text/html; charset=utf-8", ".png": "image/png", ".svg": "image/svg+xml",
                     ".json": "application/json", ".pdf": "application/pdf"}.get(target.suffix, "application/octet-stream")
            return self._ok(ctype, target.read_bytes())
        self._ok("text/plain", b"not found", 404)

    def do_POST(self):
        if not self._auth():
            return
        if urlparse(self.path).path != "/run":
            return self._ok("text/plain", b"not found", 404)
        form = cgi.FieldStorage(fp=self.rfile, headers=self.headers,
                                environ={"REQUEST_METHOD": "POST", "CONTENT_TYPE": self.headers["Content-Type"]})
        try:
            from .cli import run_all
            with tempfile.TemporaryDirectory() as td:
                pdb_item = form["pdb"]
                name = Path(pdb_item.filename).name or "upload.pdb"
                pdb = Path(td) / name
                pdb.write_bytes(pdb_item.file.read())
                tpl = None
                if "template" in form and getattr(form["template"], "filename", ""):
                    tpl = Path(td) / Path(form["template"].filename).name
                    tpl.write_bytes(form["template"].file.read())
                charge = form.getvalue("charge")
                res = run_all(str(pdb), str(RESULTS), ligand=(form.getvalue("ligand") or "").strip() or None,
                              template=str(tpl) if tpl else None, charge=int(charge) if charge else None,
                              offline=bool(form.getvalue("offline")), no_png=bool(form.getvalue("nopng")))
            body = "<p class='ok'>Completed</p><pre>" + json.dumps(res, indent=1, ensure_ascii=False) + "</pre>"
            body += f'<p><a href="/?token={TOKEN}">Back</a></p>'
            self._ok("text/html; charset=utf-8", ("<meta charset=utf-8>" + body).encode())
        except Exception as ex:
            msg = "<p class='bad'>Failed</p><pre>" + "".join(traceback.format_exception_only(type(ex), ex)) + "</pre>"
            self._ok("text/html; charset=utf-8", ("<meta charset=utf-8>" + msg).encode(), 500)


def serve(port: int = 8765, open_browser: bool = True) -> None:
    RESULTS.mkdir(exist_ok=True)
    srv = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    url = f"http://127.0.0.1:{port}/?token={TOKEN}"
    print("ZLIS web UI:", url, flush=True)
    if open_browser:
        threading.Timer(0.8, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
