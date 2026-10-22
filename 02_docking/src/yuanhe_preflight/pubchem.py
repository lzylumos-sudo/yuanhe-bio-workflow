"""Public, opt-in PubChem supplier-deposited catalog matching (no API key).

Depositor identity, SKU and product URL are evidence of a catalog match, not
current stock. Remote hits remain UNCONFIRMED until dated supplier evidence is
imported through Catalog. No price, lead time or stock is inferred.
"""
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen
from urllib.error import HTTPError
from filelock import FileLock
from .common import canonical, hash_obj, file_hash, read, write, seal, verify_seal

BASE = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/'


class PubChem:
    def __init__(self, config, cache_dir):
        self.cfg = config
        self.cache = Path(cache_dir)/'pubchem-v1'
        self.cache.mkdir(parents=True, exist_ok=True)

    def query(self, smiles):
        smiles = canonical(smiles)
        key = hash_obj({'smiles': smiles, 'config': self.cfg})
        path = self.cache/(key+'.json')
        with FileLock(str(self.cache/'service.lock'), timeout=600):
            if path.exists():
                try:
                    saved = read(path); verify_seal(saved)
                    if saved['status'] == 'SUCCESS' and time.time()-saved['queried_timestamp'] < 86400:
                        return {**saved, 'cache_hit': True}
                except (ValueError, KeyError, OSError):
                    pass
            evidence = []
            def get(suffix):
                # Below PubChem's five requests/second policy, including all workers.
                time.sleep(.25)
                url = BASE + suffix
                req = Request(url, headers={'User-Agent': 'YuanheBioPredocking/2.0'})
                for attempt in range(3):
                    try:
                        with urlopen(req, timeout=self.cfg['timeout_seconds']) as response:
                            raw = json.load(response)
                        break
                    except HTTPError as exc:
                        if exc.code not in (429, 500, 502, 503, 504) or attempt == 2:
                            raise
                        time.sleep(attempt+1)
                artifact = self.cache/(hash_obj({'url': url, 'raw': raw})+'.response.json')
                write(artifact, {'url': url, 'body': raw})
                evidence.append({'url': url, 'path': str(artifact), 'sha256': file_hash(artifact)})
                return raw
            records = []; status = 'SUCCESS'; error = ''; cid = None
            try:
                props = get('compound/smiles/'+quote(smiles, safe='')+'/property/IsomericSMILES/JSON')['PropertyTable']['Properties']
                exact = [p for p in props if canonical(p.get('SMILES', p.get('IsomericSMILES', ''))) == smiles]
                if not exact:
                    status = 'STRUCTURE_MISMATCH'
                    error = 'PubChem standardization differs; salt/stereo/tautomer equivalence not accepted'
                else:
                    cid = exact[0]['CID']
                    sids = get(f'compound/cid/{cid}/sids/JSON?sid_type=source')['InformationList']['Information'][0].get('SID', [])
                    sids = sids[:self.cfg['max_sids']]
                    if sids:
                        substances = get('substance/sid/'+','.join(map(str, sids))+'/JSON')['PC_Substances']
                        for sub in substances:
                            db = sub.get('source', {}).get('db', {})
                            vendor = db.get('name', '')
                            sku = db.get('source_id', {}).get('str')
                            product = next((x['sburl'] for x in sub.get('xref', []) if x.get('sburl', '').startswith(('http://', 'https://'))), None)
                            linked = any(c.get('id', {}).get('id', {}).get('cid') == cid for c in sub.get('compound', []))
                            if vendor in self.cfg['vendors'] and sku and product and linked:
                                records.append({'vendor': vendor, 'sku': sku, 'purchase_url': product,
                                                'sid': sub['sid']['id'], 'cid': cid, 'exact_smiles': smiles,
                                                'status': 'catalog_hit_supply_unconfirmed', 'source': 'PubChem supplier deposition',
                                                'stereochemistry_match': True, 'salt_difference': False,
                                                'availability_label': 'PURCHASABLE_UNCONFIRMED'})
            except HTTPError as exc:
                status = 'NOT_MATCHED' if exc.code == 404 and cid is None else 'ERROR'
                error = 'HTTP '+str(exc.code)
            except Exception as exc:
                status = 'ERROR'; error = type(exc).__name__+': '+str(exc)
            result = seal({'smiles': smiles, 'status': status, 'records': records, 'evidence': evidence,
                           'query_time': datetime.now(timezone.utc).isoformat(), 'queried_timestamp': time.time(),
                           'error_message': error, 'cid': cid, 'cache_hit': False,
                           'coverage': 'First '+str(self.cfg['max_sids'])+' depositor records; non-exhaustive',
                           'label': 'PURCHASABLE_UNCONFIRMED' if records else 'NOT_MATCHED' if status in ('SUCCESS', 'NOT_MATCHED') else 'AVAILABILITY_UNKNOWN'})
            write(path, result)
            return result
