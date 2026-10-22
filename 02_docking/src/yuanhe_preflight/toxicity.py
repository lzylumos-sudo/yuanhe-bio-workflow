"""Batch subprocess adapter with explicit endpoint semantics and failure records."""
import subprocess
import sys
import time
from pathlib import Path
from .common import write, read
from .toxicity_worker import SUPPORTED


def predict_batch(rows, config, cache_dir, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    request = out/'toxicity.request.json'
    result = out/'toxicity.raw.json'
    write(request, {'smiles': [r['smiles'] for r in rows], 'config': config, 'cache_dir': cache_dir})
    # Never consume an output left by an earlier failed attempt.
    if result.exists():
        result.unlink()
    error = None
    if rows:
        try:
            with (out/'toxicity.log').open('w') as log:
                proc = subprocess.run([config['python'] or sys.executable, str(Path(__file__).with_name('toxicity_worker.py')),
                                       '--request', str(request), '--output', str(result)],
                                      stdout=log, stderr=log, timeout=config['timeout_seconds'], check=False)
            if proc.returncode:
                error = 'Model worker failed; see toxicity.log (exit '+str(proc.returncode)+')'
        except (OSError, subprocess.TimeoutExpired) as exc:
            error = type(exc).__name__ + ': ' + str(exc)
    try:
        raw = read(result) if result.is_file() else {'records': {}}
    except (ValueError, OSError):
        raw = {'records': {}}
    outputs = []
    for row in rows:
        saved = raw['records'].get(row['smiles'], {})
        details = []
        for endpoint in config['endpoints']:
            value = saved.get('raw', {}).get(endpoint)
            status = 'UNSUPPORTED' if endpoint not in SUPPORTED else 'SUCCESS' if saved.get('status') == 'SUCCESS' and value is not None else 'ERROR'
            high = value is not None and value >= config['thresholds'][endpoint]
            details.append({'molecule_id': row['compound_id'], 'model_name': 'ADMET-AI',
                            'model_version': config['model_version'], 'endpoint': endpoint,
                            'raw_prediction': value, 'prediction_unit': SUPPORTED.get(endpoint),
                            'risk_level': ('HIGH' if high else 'BELOW_THRESHOLD') if status == 'SUCCESS' else 'UNKNOWN',
                            'threshold': config['thresholds'][endpoint], 'prediction_status': status,
                            'error_message': '' if status == 'SUCCESS' else saved.get('error') or error or 'Endpoint unavailable',
                            'applicability_domain': 'UNSUPPORTED', 'cache_hit': saved.get('cache_hit', False),
                            'provenance': saved.get('provenance', {})})
        outputs.append(details)
    return outputs
