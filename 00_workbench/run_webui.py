# -*- coding: utf-8 -*-
"""Start the Yuanhe workbench; serving the UI does not launch calculations.

Use --no-browser for a headless local server. Default: 127.0.0.1:8777.
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from webui.server import serve  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description='Yuanhe Bio Research Studio')
    ap.add_argument('--host', default='127.0.0.1')
    ap.add_argument('--port', type=int, default=8777)
    ap.add_argument('--no-browser', action='store_true')
    ap.add_argument('--quiet', action='store_true')
    args = ap.parse_args()
    return serve(host=args.host, port=args.port,
                 open_browser=not args.no_browser, quiet=args.quiet)


if __name__ == '__main__':
    raise SystemExit(main())
