"""Evidence-only preflight revisions. Reading/cancelling never invalidates results."""
import copy
import hashlib
import json
import time
from pathlib import Path
from . import workflow as w


def source_allowed(name):
    return name in {
        '00_workbench/webui/server.py', '00_workbench/webui/workflow.py',
        '00_workbench/webui/preflight_revision.py',
        '02_docking/src/yuanhe_preflight/structure.py',
    }


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def context(rid, stage):
    if stage not in ('risk', 'availability'):
        raise ValueError('Invalid preflight stage')
    d = w.run_dir(rid)
    cfg = w.read(d/'config.json')
    if not cfg or w.read(d/'run.json', {}).get('config_sha256') != w.digest(d/'config.json'):
        raise ValueError('Run missing or frozen configuration changed')
    if (d/'preflight_revision.pending.json').exists():
        raise ValueError('Incomplete evidence revision: administrator recovery required before continuing')
    if any((d/k).exists() for k in ('docking', 'screening', 'gromacs', 'md_systems')):
        raise ValueError('Downstream artifacts exist; create a separate evidence-review run')
    for key, _, _ in w.STAGES:
        if w.read(d/key/'stage.json', {}).get('status') in ('running', 'interrupted'):
            raise ValueError('Stop running or interrupted work before revising evidence')
    if any(w.read(p, {}).get('status') in ('running', 'interrupted') for p in (d/'operations').glob('*/operation.json')):
        raise ValueError('Stop background preparation before revising evidence')
    upstream = ['prepare', 'peptide2mol', 'reinvent4'] + (['risk'] if stage == 'availability' else [])
    for key in upstream:
        w.require(d, key)
    affected = ['risk', 'availability'] if stage == 'risk' else ['availability']
    for key in affected:
        if w.read(d/key/'stage.json', {}).get('status') == 'completed':
            w.require(d, key)
    predocking = copy.deepcopy(cfg['predocking'])
    for key in ('risk', 'availability'):
        predocking = copy.deepcopy(w.read(d/'reviews'/f'{key}.json', {}).get('predocking', predocking))
    paths = ['config.json', 'source_snapshot.json', 'active_source_revision.json',
             'reviews/risk.json', 'reviews/availability.json']
    paths += [f'{key}/stage.json' for key in upstream + affected]
    paths += [f'{key}/records.json' for key in affected]
    state = {p: w.digest(d/p) if (d/p).is_file() else None for p in paths}
    state['current_source'] = w.source_snapshot()
    previous = w.read(d/'active_source_revision.json')
    if previous:
        record = d/'source_revisions'/previous['file']
        if (record.parent != d/'source_revisions' or not record.is_file() or
            w.digest(record) != previous['sha256'] or
            w.read(record).get('initial_source_sha256') != w.digest(d/'source_snapshot.json')):
            raise ValueError('Previous source revision evidence changed')
    rows = w.read(d/stage/'records.json', [])
    candidates = []
    for row in rows:
        candidates.append({
            'compound_id': row.get('compound_id', row.get('id')),
            'decision': row.get('decision'), 'reason': row.get('reason'),
            'alerts': [{'rule_id': a['rule_id'], 'action': a['action']}
                       for a in row.get('risk', {}).get('alerts', [])],
            'routes': [{'route_sha256': r.get('route_sha256'), 'status': r.get('status'),
                        'problems': r.get('problems', []),
                        'reviewable': bool(r.get('step_count') and
                            r.get('tree', {}).get('route_metadata', {}).get('is_solved') is True and
                            not set(r.get('problems', [])) - {'stereochemistry_selectivity_requires_review', 'template_selectivity_requires_review'})}
                       for r in row.get('availability', {}).get('routes', [])],
        })
    return {'revision_token': fingerprint(state), 'stage': stage, 'candidates': candidates,
            'existing_reviews': {k: predocking.get(k, []) for k in ('risk_exceptions', 'route_reviews')},
            'policy': 'Thresholds, rules, models, budgets and complete-route requirements remain frozen.'}


def _reviewer(item):
    for key in ('reviewer', 'reason', 'reviewed_at'):
        if not isinstance(item.get(key), str) or not item[key].strip():
            raise ValueError('Each evidence record requires reviewer, reason and reviewed_at')


def submit(rid, stage, payload):
    if not isinstance(payload, dict) or set(payload) != {'revision_token', 'reviewer', 'reason', 'risk_reviews', 'route_reviews'}:
        raise ValueError('Only evidence reviews are editable; execution and screening settings are frozen')
    for key in ('reviewer', 'reason'):
        if not isinstance(payload[key], str) or not payload[key].strip():
            raise ValueError('Revision requires a reviewer and reason')
    for key in ('risk_reviews', 'route_reviews'):
        if not isinstance(payload[key], list):
            raise ValueError(key + ' must be an array')
    if (stage == 'risk' and payload['route_reviews']) or (stage == 'availability' and payload['risk_reviews']):
        raise ValueError('Review evidence in its own stage')
    d = w.run_dir(rid)
    with w._lock:
        if w._active:
            raise ValueError('Wait for active work to finish')
        lease = w.acquire_compute()
        try:
            current = context(rid, stage)
            if payload['revision_token'] != current['revision_token']:
                raise ValueError('Evidence changed while the editor was open; reopen and review the latest records')
            cfg = w.read(d/'config.json')
            predocking = copy.deepcopy(cfg['predocking'])
            for key in ('risk', 'availability'):
                predocking = copy.deepcopy(w.read(d/'reviews'/f'{key}.json', {}).get('predocking', predocking))
            policy = lambda p: {k: v for k, v in p.items() if k not in ('risk_exceptions', 'route_reviews')}
            if policy(predocking) != policy(cfg['predocking']):
                raise ValueError('Existing review changed frozen screening settings; create a new run')
            initial = w.read(d/'source_snapshot.json')
            source = w.source_snapshot()
            if initial is None:
                raise ValueError('Source snapshot missing; create a new run')
            changed = {k for k in set(initial) | set(source) if initial.get(k) != source.get(k)}
            if any(not source_allowed(k) for k in changed):
                raise ValueError('Upstream or scientific policy code changed; create a new run')
            rows = {r.get('compound_id', r.get('id')): r for r in w.read(d/stage/'records.json', [])}
            sources = []
            metadata = {'reviewer', 'reason', 'reviewed_at'}
            for item in payload['risk_reviews']:
                if not isinstance(item, dict) or set(item) != metadata | {'compound_id', 'rule_id'}:
                    raise ValueError('Risk reviews require compound_id, rule_id and review metadata')
                _reviewer(item)
                row = rows.get(item['compound_id'], {})
                risk = row.get('risk', {})
                alert = next((a for a in risk.get('alerts', []) if a['rule_id'] == item['rule_id']), None)
                if not alert or alert['action'] != 'review':
                    raise ValueError('Only a recorded review alert can be reviewed; hard exclusions cannot be waived')
                if risk.get('rules_sha256') != w.digest(predocking['rules_path']):
                    raise ValueError('Risk rules changed')
                record = {k: item[k] for k in metadata | {'rule_id'}}
                record.update(structure_smiles=row['smiles'], rules_sha256=risk['rules_sha256'])
                if item['rule_id'].startswith(('PAINS:', 'BRENK:')):
                    record['catalog_policy_sha256'] = risk['catalog_policy_sha256']
                if any(x['rule_id'] == record['rule_id'] and x['structure_smiles'] == record['structure_smiles'] for x in predocking['risk_exceptions']):
                    raise ValueError('This alert already has a review; preserve it and create a separate run to replace it')
                predocking['risk_exceptions'].append(record)
            for item in payload['route_reviews']:
                if not isinstance(item, dict) or set(item) != metadata | {'compound_id', 'route_sha256', 'identity_stereo_selectivity_confirmed'}:
                    raise ValueError('Route reviews require compound_id, route_sha256, chemistry confirmation and review metadata')
                _reviewer(item)
                row = rows.get(item['compound_id'], {})
                available = row.get('availability', {})
                route = next((r for r in available.get('routes', []) if r.get('route_sha256') == item['route_sha256']), None)
                permitted = {'stereochemistry_selectivity_requires_review', 'template_selectivity_requires_review'}
                if (not route or item['identity_stereo_selectivity_confirmed'] is not True or
                    not route.get('step_count') or route.get('tree', {}).get('route_metadata', {}).get('is_solved') is not True or
                    set(route.get('problems', [])) - permitted):
                    raise ValueError('A review cannot replace missing steps, templates, identity checks or planning-pool materials')
                if any(x.get('route_sha256') == item['route_sha256'] for x in predocking['route_reviews']):
                    raise ValueError('This route already has a review')
                path = Path(available.get('raw_search', {}).get('search_artifact', '')).resolve()
                roots = [d.resolve(), Path(predocking['cache_dir']).resolve()]
                if not any(path.is_relative_to(root) for root in roots) or not path.is_file():
                    raise ValueError('Recorded local route search artifact is missing')
                evidence = w.read(path)
                if (evidence.get('raw_sha256') != fingerprint(evidence.get('raw')) or
                    evidence.get('assets') != available.get('asset_manifest') or
                    evidence.get('raw', {}).get('target_smiles') != row['smiles'] or
                    not any(fingerprint(t) == item['route_sha256'] for t in evidence.get('raw', {}).get('routes', []))):
                    raise ValueError('Route search evidence changed')
                record = {k: item[k] for k in metadata | {'route_sha256', 'identity_stereo_selectivity_confirmed'}}
                record['target_smiles'] = row['smiles']
                predocking['route_reviews'].append(record)
                sources.append((path, w.digest(path), record))
            from yuanhe_preflight.config import validate
            validate(predocking)
            return _commit(d, stage, payload, predocking, sources, source, changed, initial)
        finally:
            lease.close()


def _commit(d, stage, payload, predocking, sources, source, changed, initial):
    """Rollback handled I/O failures; a process crash leaves a fail-closed marker."""
    stamp = str(time.time_ns())
    folder = d/'evidence_revisions'/stamp
    affected = ['risk', 'availability'] if stage == 'risk' else ['availability']
    paths = ['reviews/risk.json', 'reviews/availability.json', 'active_source_revision.json']
    before = {p: (d/p).read_bytes() if (d/p).is_file() else None for p in paths}
    moved = []
    pending = d/'preflight_revision.pending.json'
    try:
        folder.mkdir(parents=True)
        for src, sha, record in sources:
            dest = folder/'sources'/(sha + '.json')
            dest.parent.mkdir(exist_ok=True)
            dest.write_bytes(src.read_bytes())
            if w.digest(dest) != sha:
                raise ValueError('Route artifact changed during revision')
            record['search_artifact'] = str(dest.resolve())
        proof = {'scope': 'preflight', 'stage': stage, 'reviewer': payload['reviewer'], 'reason': payload['reason'],
                 'recorded_at': time.time(), 'request': payload, 'original_config_sha256': w.digest(d/'config.json'),
                 'before_reviews': {p: data.decode('utf-8') if data is not None else None for p, data in before.items()},
                 'predocking': predocking, 'source': source, 'initial_source_sha256': w.digest(d/'source_snapshot.json'),
                 'changes': {k: {'before': initial.get(k), 'after': source.get(k)} for k in sorted(changed)},
                 'upstream_stage_records': {k: w.digest(d/k/'stage.json') for k in ('prepare', 'peptide2mol', 'reinvent4')},
                 'archives': {k: str((folder/k).relative_to(d)) for k in affected if (d/k).exists()}}
        w.atomic(folder/'revision.json', proof)
        w.atomic(pending, {'record': str((folder/'revision.json').relative_to(d))})
        for key in affected:
            if (d/key).exists():
                (d/key).rename(folder/key)
                moved.append(key)
        review = w.read(d/'reviews'/f'{stage}.json', {})
        review['predocking'] = predocking
        w.atomic(d/'reviews'/f'{stage}.json', review)
        if stage == 'risk':
            later = w.read(d/'reviews/availability.json', {})
            later.pop('predocking', None)
            w.atomic(d/'reviews/availability.json', later)
        revision = d/'source_revisions'/('preflight_' + stamp + '.json')
        w.atomic(revision, proof)
        w.atomic(d/'active_source_revision.json', {'file': revision.name, 'sha256': w.digest(revision)})
        pending.unlink()
        return {'status': 'pending', 'stage': stage, 'revision': str((folder/'revision.json').relative_to(d)),
                'generation_unchanged': True, 'screening_policy_unchanged': True}
    except Exception:
        for p, data in before.items():
            if data is None:
                (d/p).unlink(missing_ok=True)
            else:
                (d/p).parent.mkdir(parents=True, exist_ok=True)
                (d/p).write_bytes(data)
        for key in reversed(moved):
            (folder/key).rename(d/key)
        if folder.exists():
            w.atomic(folder/'aborted.json', {'status': 'rolled_back'})
        pending.unlink(missing_ok=True)
        raise
