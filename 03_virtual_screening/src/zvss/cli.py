"""Command line interface."""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Any

from . import __version__


def _json_print(payload: Any) -> None:
    print(json.dumps(payload, ensure_ascii=False, indent=2, default=str))


def cmd_doctor(args: argparse.Namespace) -> int:
    from . import compat
    from .redock import vina

    report: dict[str, Any] = {"version": __version__, "python": sys.version.split()[0]}
    report["scratch"] = {"root": compat.scratch_root(), "ascii": compat.ascii_scratch_available()}
    try:
        compat.bootstrap()
        report["rdkit_data"] = compat.rdkit_data_dir()
    except Exception as exc:  # noqa: BLE001
        report["rdkit_data"] = f"unavailable: {exc}"
    checks: dict[str, Any] = {}
    for name in ("rdkit", "numpy", "scipy", "pandas", "gemmi", "MDAnalysis", "prolif", "posebusters", "fastapi", "yaml"):
        try:
            module = __import__(name)
            checks[name] = getattr(module, "__version__", "imported")
        except Exception as exc:  # noqa: BLE001
            checks[name] = f"missing: {type(exc).__name__}"
    report["packages"] = checks
    if args.with_tools:
        from .runctx import load_context

        ctx = load_context(root=args.root, target_profile=args.target, create=False)
        report["vina"] = vina.detect_vina(ctx).to_dict()
    _json_print(report)
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    from .runctx import load_context
    from .stages import run_all

    ctx = load_context(
        run_id=args.run_id,
        root=args.root,
        host_config=args.host,
        target_profile=args.target,
        budget_config=args.budget,
        stage_config=args.stages,
    )
    inputs = {
        "input_dirs": [os.path.abspath(p) for p in args.input] if args.input else [],
        "manifest": os.path.abspath(args.manifest) if args.manifest else None,
        "recursive": not args.no_recursive,
        "hints": json.loads(args.hints) if args.hints else {},
    }
    print(f"run_id: {ctx.run_id}")
    print(f"run_dir: {ctx.run_dir}")
    print(f"config_digest: {ctx.config.digest()}")

    def progress(stage: str, payload: dict[str, Any]) -> None:
        print(f"  [{payload.get('status')}] {stage}")

    state = run_all(ctx, inputs, stages=args.only.split(",") if args.only else None, resume=args.resume, progress=progress)
    if state.get("errors"):
        print("errors:", json.dumps(state["errors"], ensure_ascii=False, indent=2)[:4000])
        return 1
    report = (state.get("stages", {}).get("report") or {})
    print("report:", report.get("files", {}).get("report"))
    return 0


def cmd_serve(args: argparse.Namespace) -> int:
    import uvicorn

    from .api import create_app

    app = create_app(root=args.root)
    uvicorn.run(app, host="127.0.0.1", port=args.port, log_level=args.log_level)
    return 0


def cmd_runs(args: argparse.Namespace) -> int:
    from .runctx import load_context

    ctx = load_context(run_id="inspect", root=args.root, target_profile=args.target, create=False)
    _json_print(ctx.db.list_runs(limit=args.limit))
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    from .runctx import load_context
    from .report import render
    from .util import read_json

    ctx = load_context(run_id=args.run_id, root=args.root, target_profile=args.target, create=False)
    state = read_json(os.path.join(ctx.run_dir, "run_state.json"), {}) or {}
    _json_print(render.write_all(ctx, state))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="zvss", description="Zhiyu Virtual Screening Studio")
    parser.add_argument("--version", action="version", version=f"zvss {__version__}")
    parser.add_argument("--root", default=None, help="project root (default: repository root)")
    sub = parser.add_subparsers(dest="command", required=True)

    doctor = sub.add_parser("doctor", help="environment and dependency report")
    doctor.add_argument("--target", default=None)
    doctor.add_argument("--with-tools", action="store_true", help="also probe external tools (slow)")
    doctor.set_defaults(func=cmd_doctor)

    run = sub.add_parser("run", help="run the full screening pipeline")
    run.add_argument("--input", nargs="+", required=True, help="directories containing docking output")
    run.add_argument("--manifest", default=None, help="YAML/JSON input manifest")
    run.add_argument("--target", default=None, help="target profile yaml")
    run.add_argument("--budget", default=None, help="budget yaml")
    run.add_argument("--host", default=None, help="host config yaml")
    run.add_argument("--stages", default=None, help="stage settings yaml")
    run.add_argument("--run-id", default=None)
    run.add_argument("--hints", default=None, help="JSON hints (receptor, seed, box)")
    run.add_argument("--only", default=None, help="comma separated stage subset")
    run.add_argument("--resume", action="store_true", help="reuse cached stage results")
    run.add_argument("--no-recursive", action="store_true")
    run.set_defaults(func=cmd_run)

    serve = sub.add_parser("serve", help="start the local web UI")
    serve.add_argument("--port", type=int, default=8765)
    serve.add_argument("--log-level", default="info")
    serve.set_defaults(func=cmd_serve)

    runs = sub.add_parser("runs", help="list recorded runs")
    runs.add_argument("--limit", type=int, default=20)
    runs.add_argument("--target", default=None)
    runs.set_defaults(func=cmd_runs)

    report = sub.add_parser("report", help="re-render the report of an existing run")
    report.add_argument("--run-id", required=True)
    report.add_argument("--target", default=None)
    report.set_defaults(func=cmd_report)
    return parser


def main(argv: list[str] | None = None) -> int:
    from . import compat

    compat.configure_temp()
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args) or 0)
    except Exception as exc:  # noqa: BLE001
        print(f"error: {type(exc).__name__}: {exc}", file=sys.stderr)
        if os.environ.get("ZVSS_DEBUG"):
            raise
        return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
