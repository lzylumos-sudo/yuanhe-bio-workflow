# -*- coding: utf-8 -*-
"""Run the pepwedge scaffold-design pipeline.

    python run_pipeline.py [--root .] [--out results]
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pepwedge.pipeline import run  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', default='.')
    ap.add_argument('--out', default='results')
    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args()
    run(a.root, a.out, verbose=not a.quiet)


if __name__ == '__main__':
    main()
