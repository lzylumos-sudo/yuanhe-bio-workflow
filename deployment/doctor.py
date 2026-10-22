#!/usr/bin/env python3
"""Complete workbench readiness report without requiring an existing task."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / '00_workbench'))
from webui import workflow
from webui.readiness import collect


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--ui-python', help='Actual service Python, if different from chem_python')
    args = parser.parse_args()
    try:
        cfg = workflow.defaults()
        if args.config:
            cfg.update(json.loads(args.config.read_text()))
        result = collect(cfg, ROOT, ui_python=args.ui_python or cfg.get('chem_python') or sys.executable)
    except Exception as exc:
        result = {'status': 'blocked', 'all_runtime_checks_pass': False, 'error': str(exc)}
    text = json.dumps(result, indent=2, ensure_ascii=False)
    print(text)
    if args.out:
        workflow.atomic(args.out, result)
    return 0 if result['all_runtime_checks_pass'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
