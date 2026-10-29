"""Isolated independent-module worker, launched by the UI job manager."""
import argparse
import os
import time
import traceback
from webui import workflow as w
from webui import module_jobs as jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-id', required=True)
    ap.add_argument('--operation', required=True)
    args = ap.parse_args()
    if not args.operation.isdigit():
        raise ValueError('Invalid operation')
    d, req = jobs.verified(args.run_id)
    out = d/'operations'/args.operation
    rec = w.read(out/'operation.json')
    if not rec or rec['status'] != 'running':
        return 1
    from scientific import worker_processes
    worker_processes.install()
    try:
        if w.read(d/'source.json') != w.source_snapshot():
            raise ValueError('Source changed before module execution')
        from scientific.module_adapters import execute
        result = execute(d, req, rec, out)
        w.atomic(out/'result.json', {'scope': jobs.SCOPE, 'candidate_qualified': False, 'result': result})
        if w.read(out/'operation.json')['status'] == 'running':
            w.atomic(out/'operation.json', {**w.read(out/'operation.json'), 'status': 'blocked' if result.get('status') == 'blocked_by_strict_preflight' else 'completed', 'ended': time.time()})
        return 0
    except Exception as exc:
        traceback.print_exc()
        if w.read(out/'operation.json')['status'] == 'running':
            w.atomic(out/'operation.json', {**w.read(out/'operation.json'), 'status': 'failed',
                                           'error': str(exc), 'ended': time.time()})
        return 1
    finally:
        worker_processes.stop_children(worker_processes.descendants(os.getpid()))


if __name__ == '__main__':
    raise SystemExit(main())
