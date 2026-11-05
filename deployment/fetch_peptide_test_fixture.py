"""Prepare all pinned Peptide2Mol patch sources for CPU contracts, no inference."""
import hashlib
import json
from pathlib import Path
import urllib.request
import apply_vendor_patches


def main():
    root = Path(__file__).resolve().parents[1]
    rows = [row for row in json.loads((root / 'deployment/vendor_patches/manifest.json').read_text())
            if row['vendor'] == 'Peptide2Mol']
    config = json.loads((root / 'deployment/model_dependencies.json').read_text())['peptide2mol']
    for row in rows:
        path = root / '01_molecule_generation/vendor/Peptide2Mol' / row['file']
        if path.exists():
            continue  # The patcher verifies every existing source before any mutation.
        if row['file'] == 'src/eval.py':
            data = (root / 'deployment/tests/fixtures/peptide2mol/eval.py.txt').read_bytes()
        else:
            url = ('https://raw.githubusercontent.com/BLUE-Flowing/Peptide2Mol/' +
                   config['repo_commit'] + '/' + row['file'])
            data = urllib.request.urlopen(url, timeout=60).read()
        if hashlib.sha256(data).hexdigest() != row['before_sha256']:
            raise ValueError('Pinned sampling source hash mismatch: ' + row['file'])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    apply_vendor_patches.main('Peptide2Mol')
    print('Verified CPU sampling and seed fixtures; no model or GPU execution')


if __name__ == '__main__':
    main()
