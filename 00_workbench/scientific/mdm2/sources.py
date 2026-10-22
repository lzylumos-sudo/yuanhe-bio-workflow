"""Official source acquisition. Raw activities remain unverified until reviewed."""
import concurrent.futures
import datetime
import hashlib
import json
import time
import urllib.request
from pathlib import Path

BASE = 'https://www.ebi.ac.uk/chembl/api/data/'


def fetch(url, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.with_suffix(path.suffix + '.source.json').exists():
        return path
    last = None
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, headers={'User-Agent': 'YuanheBio-research/MDM2-benchmark'})
            with urllib.request.urlopen(request, timeout=90) as response:
                content = response.read()
            if not content:
                raise ValueError('Empty official source')
            if path.suffix == '.json':
                json.loads(content)
            if path.suffix == '.cif' and not content.lstrip().startswith(b'data_'):
                raise ValueError('Not an mmCIF file')
            temporary = path.with_suffix(path.suffix + '.part')
            temporary.write_bytes(content)
            temporary.replace(path)
            provenance = {'url': url, 'retrieved_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                          'sha256': hashlib.sha256(content).hexdigest(), 'bytes': len(content)}
            path.with_suffix(path.suffix + '.source.json').write_text(json.dumps(provenance, indent=2))
            return path
        except Exception as exc:
            last = exc
            time.sleep(attempt + 1)
    raise RuntimeError(f'{url}: {last}')


def acquire(root):
    root = Path(root)
    tasks = [(BASE + 'target/CHEMBL5023.json', root/'sources/target.json'),
             ('https://files.rcsb.org/download/1YCR.cif', root/'sources/1YCR.cif')]
    for pdb in ('5Z02', '5TRF', '8PWC'):
        tasks.append((f'https://files.rcsb.org/download/{pdb}.cif', root/f'evaluation/structures/{pdb}.cif'))
    for offset in range(0, 14000, 1000):
        tasks.append((BASE + f'activity.json?target_chembl_id=CHEMBL5023&limit=1000&offset={offset}',
                      root/f'evaluation/raw/activity_{offset:05d}.json'))
    # A real source pool, not assumed inactive. Final property matching is performed later.
    for offset in range(0, 8000, 1000):
        tasks.append((BASE + f'molecule.json?molecule_properties__full_mwt__gte=300&molecule_properties__full_mwt__lte=800&limit=1000&offset={offset}',
                      root/f'evaluation/raw/background_{offset:05d}.json'))
    for pmc in ('PMC4247201', 'PMC11612618', 'PMC2660734'):
        tasks.append((f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML',
                      root/f'evaluation/literature/{pmc}.xml'))
    for offset in range(0,15000,1000):
        tasks.append((BASE+f'molecule.json?molecule_properties__full_mwt__gte=500&molecule_properties__full_mwt__lte=620&molecule_properties__hbd__gte=2&molecule_properties__hbd__lte=5&molecule_properties__alogp__gte=2.5&molecule_properties__alogp__lte=7&limit=1000&offset={offset}',root/f'evaluation/raw/background_matched_{offset:05d}.json'))
    for offset in range(0,15000,1000):
        tasks.append((BASE+f'molecule.json?molecule_properties__full_mwt__gte=500&molecule_properties__full_mwt__lte=620&molecule_properties__hbd__gte=3&molecule_properties__hbd__lte=5&molecule_properties__rtb__lte=6&molecule_properties__psa__lte=120&limit=1000&offset={offset}',root/f'evaluation/raw/background_rigid_{offset:05d}.json'))
    for ccd in ('NUT','7HC','G7I'):
        tasks.append((f'https://files.rcsb.org/ligands/download/{ccd}.cif',root/f'evaluation/ccd/{ccd}.cif'))
    for pmc in ('PMC4247201','PMC4281096'):
        tasks.append((f'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/{pmc}/unicode',root/f'evaluation/literature/{pmc}.bioc.json'))
    for record in list(range(2250131,2250151))+[2250154]:
        tasks.append((BASE+f'compound_record/{record}.json',root/f'evaluation/records/{record}.json'))
    for url,name in [('https://www.osti.gov/servlets/purl/1345500','SAR405838_primary.pdf'),
                     ('https://ndownloader.figshare.com/files/50939479','Brigimadlin_methods.pdf'),
                     ('https://ndownloader.figshare.com/files/50939455','Brigimadlin_tableS4.pdf')]:
        tasks.append((url,root/'evaluation/literature'/name))
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
        futures = {pool.submit(fetch, url, path): (url, path) for url, path in tasks}
        for future in concurrent.futures.as_completed(futures):
            url, path = futures[future]
            try:
                future.result()
                results.append({'path': str(path.relative_to(root)), 'status': 'downloaded'})
            except Exception as exc:
                results.append({'path': str(path.relative_to(root)), 'status': 'failed', 'error': str(exc)})
            print(json.dumps(results[-1]), flush=True)
    (root/'sources/acquisition.json').write_text(json.dumps(results, indent=2))
    return results
