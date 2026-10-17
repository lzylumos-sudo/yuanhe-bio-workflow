"""Local web UI (127.0.0.1 only) over the run database and report files."""

from __future__ import annotations

import os
from typing import Any

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse
from fastapi.templating import Jinja2Templates

from . import __version__
from .util import read_json

TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), "web", "templates")


def create_app(root: str | None = None) -> FastAPI:
    from .runctx import load_context

    app = FastAPI(title="ZVSS", version=__version__)
    templates = Jinja2Templates(directory=TEMPLATE_DIR)
    base = os.path.abspath(root) if root else None

    def _ctx(run_id: str = "inspect"):
        return load_context(run_id=run_id, root=base, create=False)

    @app.get("/", response_class=HTMLResponse)
    def index(request: Request) -> HTMLResponse:
        ctx = _ctx()
        runs = ctx.db.list_runs(limit=50)
        return templates.TemplateResponse(request, "index.html", {"runs": runs, "version": __version__})

    @app.get("/run/{run_id}", response_class=HTMLResponse)
    def run_view(request: Request, run_id: str) -> HTMLResponse:
        ctx = _ctx(run_id)
        run = ctx.db.get_run(run_id)
        if not run:
            raise HTTPException(status_code=404, detail="run not found")
        stages = {}
        directory = os.path.join(ctx.run_dir, "stages")
        if os.path.isdir(directory):
            for name in sorted(os.listdir(directory)):
                summary = os.path.join(directory, name, "summary.json")
                error = os.path.join(directory, name, "error.json")
                if os.path.isfile(error):
                    stages[name] = {"status": "error", "detail": read_json(error, {})}
                elif os.path.isfile(summary):
                    stages[name] = {"status": "done", "detail": read_json(summary, {})}
        scoring = (stages.get("scoring") or {}).get("detail") or {}
        selection = (stages.get("selection") or {}).get("detail") or {}
        return templates.TemplateResponse(
            request,
            "run.html",
            {
                "run": run,
                "run_id": run_id,
                "stages": stages,
                "ranking": (scoring.get("ranking") or [])[:40],
                "primary": selection.get("primary") or [],
                "exploratory": selection.get("exploratory") or [],
                "exclusions": selection.get("exclusions") or [],
                "selection_version": selection.get("selection_version"),
            },
        )

    @app.get("/api/runs")
    def api_runs() -> JSONResponse:
        return JSONResponse(_ctx().db.list_runs(limit=50))

    @app.get("/api/run/{run_id}")
    def api_run(run_id: str) -> JSONResponse:
        ctx = _ctx(run_id)
        run = ctx.db.get_run(run_id)
        if not run:
            raise HTTPException(status_code=404, detail="run not found")
        payload: dict[str, Any] = {"run": run, "stages": {}}
        directory = os.path.join(ctx.run_dir, "stages")
        if os.path.isdir(directory):
            for name in sorted(os.listdir(directory)):
                path = os.path.join(directory, name, "summary.json")
                if os.path.isfile(path):
                    payload["stages"][name] = read_json(path, {})
        return JSONResponse(payload)

    @app.get("/report/{run_id}", response_class=PlainTextResponse)
    def report_text(run_id: str) -> PlainTextResponse:
        ctx = _ctx(run_id)
        path = os.path.join(ctx.reports_dir(), "report.md")
        if not os.path.exists(path):
            raise HTTPException(status_code=404, detail="report not found")
        with open(path, "r", encoding="utf-8") as handle:
            return PlainTextResponse(handle.read())

    @app.get("/api/health")
    def health() -> JSONResponse:
        return JSONResponse({"status": "ok", "version": __version__})

    return app
