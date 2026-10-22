"""Pre-Docking Screening Engine: layered decisions, independent objectives.

Called inside the existing availability stage before qualification sidecars are
issued. Structural exclusions, duplicates and failed inference remain exported.
"""
import itertools
import time
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from .common import seal, hash_obj, verify_seal
from .toxicity import predict_batch
from .synthesis import score_batch


def evaluate(rows, cfg, out):
    settings = cfg['engine']
    eligible = [r for r in rows if r['decision'] == 'pass']
    started = time.monotonic()
    predictions = predict_batch(eligible, settings['toxicity'], cfg['cache_dir'], out)
    tox_seconds = time.monotonic()-started
    started = time.monotonic()
    synthesis = score_batch(eligible, settings['synthesis']) if eligible else []
    sa_seconds = time.monotonic()-started
    scored = {}
    for row, tox, syn in zip(eligible, predictions, synthesis):
        reasons = []
        if any(x['prediction_status'] != 'SUCCESS' for x in tox):
            reasons.append('toxicity_incomplete')
        if any(x['risk_level'] == 'HIGH' for x in tox):
            reasons.append('predicted_toxicity_threshold_exceeded')
        if syn['screening_status'] != 'PASS':
            reasons.append(syn['rejection_reason'])
        record = seal({'smiles': row['smiles'], 'molecule_id': row['compound_id'],
                       'toxicity_predictions': tox, 'synthesis': syn,
                       'policy_sha256': hash_obj(settings), 'version': '2.0.0',
                       'decision': 'review' if reasons else 'pass', 'reasons': reasons})
        scored[row['compound_id']] = {**row, 'screening': record, 'toxicity_predictions': tox,
                                      'synthesis': syn, 'decision': record['decision'],
                                      'reason': '; '.join(reasons) if reasons else row['reason']}
    results = [scored[r['compound_id']] if r['decision'] == 'pass' else dict(r) for r in rows]
    return results, {'toxicity_seconds': tox_seconds, 'synthesis_seconds': sa_seconds,
                     'model_evaluated': len(eligible),
                     'model_failed': sum(any(p['prediction_status'] != 'SUCCESS' for p in t) for t in predictions)}


def validate_screening(record, settings, smiles):
    verify_seal(record)
    if record['smiles'] != smiles or record['policy_sha256'] != hash_obj(settings):
        raise ValueError('Screening identity or policy changed')
    tox = record['toxicity_predictions']
    if [x['endpoint'] for x in tox] != settings['toxicity']['endpoints']:
        raise ValueError('Toxicity endpoint evidence incomplete')
    import math
    for p in tox:
        value = p['raw_prediction']
        if p['prediction_status'] != 'SUCCESS' or not p.get('provenance', {}).get('weights'):
            raise ValueError('Real model provenance and successful inference required')
        if not isinstance(value, (int, float)) or not math.isfinite(value) or value >= settings['toxicity']['thresholds'][p['endpoint']]:
            raise ValueError('Toxicity threshold not satisfied')
    syn = record['synthesis']
    for key, threshold in [('sa_score', 'sa_max'), ('sc_score', 'sc_max')]:
        value = syn.get(key)
        if not isinstance(value, (int, float)) or not math.isfinite(value) or value > settings['synthesis'][threshold]:
            raise ValueError('Synthesis threshold not satisfied')
    if syn['score_version']['sc_weights_sha256'] != settings['synthesis']['scscore_sha256'] or record['decision'] != 'pass':
        raise ValueError('Synthesis model evidence changed or screening blocked')
    from .common import file_hash
    if file_hash(settings['synthesis']['scscore_path']) != settings['synthesis']['scscore_sha256']:
        raise ValueError('SCScore weights changed')


def rank(rows, settings):
    """Independent threshold tiers, lexicographic objectives, scaffold diversity.

    No addition of SA, SC or endpoint raw scores. Deterministic ties use IDs.
    Discovery alternates scaffolds; balanced alternates within acquisition tier;
    delivery prefers verified acquisition and independent difficulty objectives.
    """
    def key(r):
        syn = r.get('synthesis', {})
        a = r.get('availability', {})
        acquisition = 0 if a.get('status') == 'purchase_supported' else 1 if a.get('status') == 'route_supported' else 2
        tox = r.get('toxicity_predictions', [])
        tox_tier = sum(x.get('raw_prediction', 0) >= x.get('threshold', 1) * .75 for x in tox)
        return (acquisition, tox_tier, syn.get('sa_score') or 10, syn.get('sc_score') or 5, r['compound_id'])
    ordered = sorted(rows, key=key)
    for row in ordered:
        row['ranking'] = {'mode': settings['mode'], 'objectives': list(key(row)),
                          'scaffold': MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(row['smiles']), includeChirality=True) or row['smiles'],
                          'method': 'lexicographic-independent-thresholds-scaffold-round-robin'}
    if settings['mode'] != 'delivery':
        groups = {}
        for row in ordered:
            tier = key(row)[0] if settings['mode'] == 'balanced' else 0
            groups.setdefault(tier, {}).setdefault(row['ranking']['scaffold'], []).append(row)
        ordered = [r for tier in sorted(groups) for layer in itertools.zip_longest(*groups[tier].values()) for r in layer if r is not None]
    for i, row in enumerate(ordered):
        row['ranking']['rank'] = i + 1
    return ordered
