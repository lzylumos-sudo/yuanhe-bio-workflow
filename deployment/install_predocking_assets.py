#!/usr/bin/env python3
"""Download pinned official assets into an explicit destination; no runtime edits."""
import argparse,hashlib,json,os,shutil,ssl,sys,time
from http.client import IncompleteRead
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import urlopen
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'02_docking/src'))
from yuanhe_preflight.config import defaults

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024**2),b''):h.update(block)
    return h.hexdigest()

def transient_download_failure(exc):
    if isinstance(exc, HTTPError):
        return exc.code in (408, 500, 502, 503, 504)
    if isinstance(exc, URLError):
        exc = exc.reason
    return isinstance(exc, (TimeoutError, ConnectionError, ssl.SSLEOFError, IncompleteRead))


def download_verified(item, path):
    """Retry transient transport failures; never publish unchecked bytes."""
    temp = path.with_suffix(path.suffix + '.part')
    attempts = 3
    for attempt in range(1, attempts + 1):
        print(f"Downloading pinned asset {item['name']} from {urlsplit(item['url']).hostname} "
              f"(attempt {attempt}/{attempts})", flush=True)
        try:
            with urlopen(item['url'], timeout=120) as src, temp.open('wb') as out:
                shutil.copyfileobj(src, out)
            if digest(temp) != item['sha256']:
                raise ValueError('Official asset hash mismatch: ' + item['name'])
            os.replace(temp, path)
            return
        except Exception as exc:
            if attempt == attempts or not transient_download_failure(exc):
                raise
            print(f"Transient download failure for {item['name']} ({type(exc).__name__}); "
                  f"retrying in {attempt}s", file=sys.stderr, flush=True)
            time.sleep(attempt)
        finally:
            temp.unlink(missing_ok=True)


def install_asset(item, path, cache=None):
    """Reuse bytes only after checking the current repository's original pin."""
    if path.exists():
        if digest(path) != item['sha256']:
            raise ValueError('Existing asset differs; preserve and choose another destination: '+str(path))
        return 'existing_verified'
    cached = cache / item['filename'] if cache is not None else None
    if cached is not None and cached.is_file():
        temp = path.with_suffix(path.suffix + '.part')
        try:
            shutil.copyfile(cached, temp)
            if (temp.stat().st_size != item['bytes'] or
                    digest(temp) != item['sha256']):
                raise ValueError('Cached asset differs from pinned official asset: '+item['name'])
            os.replace(temp, path)
        finally:
            temp.unlink(missing_ok=True)
        print('Reused hash-verified official asset '+item['name'], flush=True)
        return 'verified_local_cache'
    download_verified(item, path)
    return 'official_download'


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--with-retro',action='store_true');ap.add_argument('--tox-python',default=sys.executable)
    ap.add_argument('--retro-python',default='python')
    ap.add_argument('--cache',type=Path,default=os.environ.get('YUANHE_ASSET_CACHE'),
                    help='Optional directory of official asset files; every reused file is checked against the current size and SHA-256 pins')
    args=ap.parse_args()
    if args.cache is not None and not args.cache.is_dir():
        ap.error('Asset cache must be an existing directory')
    dest=args.output.resolve();dest.mkdir(parents=True,exist_ok=True)
    specs=json.loads((ROOT/'deployment/predocking-assets.json').read_text())['assets'];installed={}
    for item in specs:
        if item['name']!='scscore' and not args.with_retro:continue
        path=dest/item['filename']
        delivery=install_asset(item, path, args.cache)
        installed[item['name']]={**item,'path':str(path),'delivery':delivery}
    cfg=defaults(ROOT);cfg['engine']['toxicity']['python']=args.tox_python
    sc=installed['scscore'];cfg['engine']['synthesis'].update(scscore_path=sc['path'],scscore_sha256=sc['sha256'])
    if args.with_retro:
        cfg['retro']['python']=args.retro_python
        for short,key in [('model','retro_model'),('templates','retro_templates'),('stock','retro_stock')]:
            item=installed[key];cfg['retro'].update({short+'_path':item['path'],short+'_version':item['sha256'],short+'_license':item['license']})
    (dest/'predocking-config.json').write_text(json.dumps(cfg,indent=2)+'\n')
    (dest/'asset-receipt.json').write_text(json.dumps(installed,indent=2)+'\n')
    print('Verified assets and predocking-config.json written to '+str(dest))
if __name__=='__main__':main()
