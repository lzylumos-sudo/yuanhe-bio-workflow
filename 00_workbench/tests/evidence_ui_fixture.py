"""Local-only UI regression server. Synthetic records; never run scientific jobs."""
import argparse
from contextlib import ExitStack
import csv
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from webui import workflow as w, server
from yuanhe_preflight.common import hash_obj
from yuanhe_preflight.structure import screen
from yuanhe_preflight.toxicity_substructure_filter import RiskRules


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--port', type=int, default=8879)
    args = parser.parse_args()
    args.root.mkdir(parents=True, exist_ok=False)
    w.RUNS = args.root.resolve()
    w.acquire_compute = ExitStack  # No jobs are allowed in this local fixture server.
    def no_jobs(*args, **kwargs):
        raise ValueError('Synthetic UI fixture: scientific execution forbidden')
    w.start = no_jobs
    cfg = w.defaults()
    for name in ('UI_RISK_FIXTURE', 'UI_ROUTE_FIXTURE'):
        d = w.RUNS/name
        w.atomic(d/'config.json', cfg)
        w.atomic(d/'run.json', {'id': name, 'config_sha256': w.digest(d/'config.json')})
        w.atomic(d/'source_snapshot.json', w.source_snapshot())
        risk = screen('CCCCCC', RiskRules(), [], cfg['predocking']['engine'])
        tree = {'type': 'mol', 'smiles': 'CCCCCC', 'route_metadata': {'is_solved': False}}
        row = {'id': 'SYNTHETIC_ONLY', 'compound_id': 'SYNTHETIC_ONLY', 'smiles': 'CCCCCC',
               'decision': 'review', 'reason': 'SYNTHETIC UI FIXTURE - NOT A SCIENTIFIC RESULT', 'risk': risk,
               'availability': {'routes': [{'route_sha256': hash_obj(tree), 'tree': tree,
                    'step_count': 1, 'problems': ['planner_did_not_solve_route']}]},
               'large_fixture_evidence': 'x' * 150000}
        for key in ('prepare', 'peptide2mol', 'reinvent4', 'risk', 'availability'):
            w.atomic(d/key/'records.json', [row] if key in ('risk', 'availability') else [])
            w.atomic(d/key/'stage.json', {'status': 'completed', 'counts': {},
                'config_sha256': w.digest(d/'config.json'),
                'artifacts': {f'{key}/records.json': w.digest(d/key/'records.json')}})
            if key in ('risk', 'availability'):
                with (d/key/'molecules.csv').open('w', newline='', encoding='utf-8') as f:
                    writer = csv.writer(f); writer.writerow(['id', 'evidence']); writer.writerow(['SYNTHETIC_ONLY', 'x'*150000])
    w.inventory = lambda: {'version': {'software_version': 'UI-REGRESSION-FIXTURE', 'git': {'commit': 'synthetic'}},
        'notice': 'Synthetic local regression fixtures. No inference, docking or MD.', 'registered_modules': [],
        'runs': ['UI_RISK_FIXTURE', 'UI_ROUTE_FIXTURE'], 'defaults': cfg}
    httpd = server.Server(('127.0.0.1', args.port), server.Handler)
    print('FIXTURE_URL=http://127.0.0.1:' + str(httpd.server_port) + '/', flush=True)
    try:
        httpd.serve_forever(poll_interval=.1)
    finally:
        httpd.server_close()


if __name__ == '__main__': main()
