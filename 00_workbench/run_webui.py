# -*- coding: utf-8 -*-
"""Start the Yuanhe workbench; serving the UI does not launch calculations.

Use --no-browser for a headless local server. Default: 127.0.0.1:8777.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def main():
    ap = argparse.ArgumentParser(description='Yuanhe Bio Research Studio')
    ap.add_argument('--host', default='127.0.0.1')
    ap.add_argument('--port', type=int, default=8777)
    ap.add_argument('--no-browser', action='store_true')
    ap.add_argument('--quiet', action='store_true')
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    runtime = root / 'deployment/runtime.json'
    try:
        cfg = json.loads(runtime.read_text()) if runtime.is_file() else {}
        selected = os.environ.get('YUANHE_UI_PYTHON') or cfg.get('chem_python') or sys.executable
        selected = shutil.which(selected) or selected
        if not os.path.isabs(selected):
            selected = str(root / selected)
        if not os.path.isfile(selected) or not os.access(selected, os.X_OK):
            raise ValueError('UI Python is missing or not executable: ' + selected)
        if os.path.abspath(selected) != os.path.abspath(sys.executable):
            env = dict(os.environ, YUANHE_UI_PYTHON=selected, PYTHONNOUSERSITE='1')
            env.pop('PYTHONHOME', None); env.pop('PYTHONPATH', None)
            os.execve(selected, [selected, '-B', str(Path(__file__).resolve()), *sys.argv[1:]], env)
        from webui.readiness import ui_check
        check = ui_check(sys.executable, root)
        if not check['ok']:
            raise ValueError('UI Python cannot import RDKit/NumPy/drawing dependencies: ' +
                             str(check.get('error') or check.get('output')))
    except (OSError, ValueError, TypeError) as exc:
        print(str(exc) + '\nAdministrator: follow deployment/WORKBENCH_DEPLOYMENT.md. '
              'Configure chem_python and restart the service. For dependency reports run '
              'python3 deployment/doctor.py --config deployment/runtime.json.', file=sys.stderr)
        return 2
    from webui.server import serve
    return serve(host=args.host, port=args.port,
                 open_browser=not args.no_browser, quiet=args.quiet)


if __name__ == '__main__':
    raise SystemExit(main())
