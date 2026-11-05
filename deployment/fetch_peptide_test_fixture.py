"""Fetch only the pinned sampling source for CPU contract tests; no weights/inference."""
import hashlib
import json
from pathlib import Path
import urllib.request
import apply_vendor_patches


def main():
    root = Path(__file__).resolve().parents[1]
    row = json.loads((root / 'deployment/vendor_patches/manifest.json').read_text())[0]
    config = json.loads((root / 'deployment/model_dependencies.json').read_text())['peptide2mol']
    url = ('https://raw.githubusercontent.com/BLUE-Flowing/Peptide2Mol/' + config['repo_commit'] + '/' + row['file'])
    path = root / '01_molecule_generation/vendor/Peptide2Mol' / row['file']
    if not path.exists():
        data = urllib.request.urlopen(url, timeout=60).read()
        if hashlib.sha256(data).hexdigest() != row['before_sha256']:
            raise ValueError('Pinned sampling source hash mismatch')
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    apply_vendor_patches.main('Peptide2Mol')
    print('Verified CPU sampling-loop fixture; no model or GPU execution')


if __name__ == '__main__':
    main()
