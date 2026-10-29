"""Local public-data browser acceptance server with genuine isolated module workers."""
import argparse
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'00_workbench'))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--root', type=Path, required=True)
    ap.add_argument('--port', type=int, default=0); ap.add_argument('--gmx', default='')
    ap.add_argument('--vina', default=''); ap.add_argument('--reinvent-executable', default='')
    args = ap.parse_args(); root = args.root.resolve(); root.mkdir(parents=True, exist_ok=True)
    os.environ['YUANHE_MODULE_JOBS'] = str(root/'jobs')
    from webui import module_jobs as jobs, workflow as w, server
    jobs.JOBS = root/'jobs'; w.RUNS = root/'formal_runs'
    cfg = w.defaults()
    cfg.update(chem_python=sys.executable, screen_python=sys.executable, md_python=sys.executable,
               gmx=args.gmx, vina=args.vina or 'vina')
    if args.reinvent_executable:
        cfg['reinvent_executable'] = args.reinvent_executable
    # No synthetic engine result: unavailable environments must fail visibly.
    w.defaults = lambda: cfg
    if args.gmx:
        from build_module_md_fixture import build
        build(args.gmx, root/'md_fixture')
    httpd = server.Server(('127.0.0.1', args.port), server.Handler)
    print('FIXTURE_URL=http://127.0.0.1:'+str(httpd.server_port)+'/', flush=True)
    try:
        httpd.serve_forever(poll_interval=.1)
    finally:
        httpd.server_close()


if __name__ == '__main__':
    main()
