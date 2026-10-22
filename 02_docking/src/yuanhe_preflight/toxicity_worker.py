"""Isolated ADMET-AI 2.0.1 inference. Successful raw predictions are resumable.

Called by absolute filename from a separate environment. No network/model
download occurs during execution. CUDA respects the scheduler's visible lease.
"""
import argparse
import importlib.metadata
import json
import math
import os
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from yuanhe_preflight.common import hash_obj, file_hash, write, read, seal, verify_seal

SUPPORTED = {'AMES': 'probability', 'hERG': 'probability', 'DILI': 'probability',
             'LD50_Zhu': 'log10(1/(mol/kg))'}


def validate_raw(values):
    if set(values) != set(SUPPORTED):
        raise ValueError('Required model endpoints incomplete')
    for ep, value in values.items():
        if type(value) not in (int, float) or not math.isfinite(value):
            raise ValueError('Nonfinite or nonnumeric model prediction')
        if SUPPORTED[ep] == 'probability' and not 0 <= value <= 1:
            raise ValueError('Model probability outside 0..1')


def run(request, output):
    config = request['config']
    if config['device'] == 'cpu':
        os.environ['CUDA_VISIBLE_DEVICES'] = ''
    os.environ.setdefault('MPLCONFIGDIR',str(Path(output).parent/'mpl-cache'))
    import torch
    import admet_ai
    from admet_ai import ADMETModel
    from admet_ai.constants import DEFAULT_MODELS_DIR
    from filelock import FileLock

    version = importlib.metadata.version('admet-ai')
    if version != config['model_version']:
        raise ValueError('ADMET-AI version mismatch: ' + version)
    if config['device'] == 'cuda' and not torch.cuda.is_available():
        raise RuntimeError('CUDA requested but unavailable')
    torch.set_num_threads(config['cpu_threads'])
    torch.manual_seed(17)
    torch.use_deterministic_algorithms(True)
    manifest = {str(p.relative_to(DEFAULT_MODELS_DIR)): file_hash(p)
                for p in sorted(Path(DEFAULT_MODELS_DIR).rglob('*.pt'))}
    if not manifest:
        raise ValueError('ADMET-AI weights missing')
    provenance = {'model_name': 'ADMET-AI', 'model_version': version, 'weights': manifest,
                  'implementation_sha256': file_hash(Path(admet_ai.__file__).with_name('admet_model.py')),
                  'worker_sha256': file_hash(__file__),
                  'packages': {k: importlib.metadata.version(k) for k in ('torch', 'chemprop', 'rdkit', 'numpy')},
                  'device': 'cuda' if torch.cuda.is_available() else 'cpu',
                  'applicability_domain': 'UNSUPPORTED', 'seed': 17,
                  'inference': {k:config[k] for k in ('device','cpu_threads','batch_size')}}
    cache = Path(request['cache_dir']) / 'toxicity-v2'
    cache.mkdir(parents=True, exist_ok=True)
    records = {}
    pending = []
    for smi in dict.fromkeys(request['smiles']):
        key = hash_obj({'smiles': smi, 'provenance': provenance})
        path = cache / (key + '.json')
        try:
            rec = read(path)
            verify_seal(rec)
            validate_raw(rec['raw'])
            if rec.get('status') != 'SUCCESS' or rec['smiles'] != smi or rec['provenance'] != provenance:
                raise ValueError('Prediction cache identity mismatch')
            records[smi] = {**rec, 'cache_hit': True}
        except (OSError, ValueError, KeyError):
            pending.append((smi, path))
    model = ADMETModel(include_physchem=False, drugbank_path=None, num_workers=0) if pending else None
    started = time.monotonic()

    def predict(batch):
        try:
            data = model.predict([x[0] for x in batch])
            for smi, path in batch:
                values = {ep: float(data.loc[smi, ep]) for ep in SUPPORTED if ep in data.columns}
                validate_raw(values)
                rec = seal({'smiles': smi, 'raw': values, 'provenance': provenance, 'status': 'SUCCESS'})
                with FileLock(str(path)+'.lock'):
                    write(path, rec)
                records[smi] = {**rec, 'cache_hit': False}
        except Exception as exc:
            # A bad batch is retried individually; failed predictions are never cached.
            if len(batch) > 1:
                for item in batch:
                    predict([item])
            else:
                records[batch[0][0]] = {'smiles': batch[0][0], 'raw': {}, 'provenance': provenance,
                                      'status': 'ERROR', 'error': type(exc).__name__+': '+str(exc)}
    for start in range(0, len(pending), config['batch_size']):
        predict(pending[start:start+config['batch_size']])
        write(output, {'records': records, 'provenance': provenance, 'completed': len(records),
                       'elapsed_seconds': time.monotonic()-started})
        write(Path(output).parent/'progress.json',{'stage':'toxicity','completed':len(records),'input':len(request['smiles'])})
        print(json.dumps({'stage': 'toxicity', 'completed': len(records), 'total': len(request['smiles'])}), flush=True)
    write(output, {'records': records, 'provenance': provenance, 'completed': len(records),
                   'elapsed_seconds': time.monotonic()-started})


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--request', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    run(read(args.request), Path(args.output))
