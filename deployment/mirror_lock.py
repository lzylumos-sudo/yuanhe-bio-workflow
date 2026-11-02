"""Rewrite only package download locations, preserving every locked artifact/hash."""
import argparse
from pathlib import Path
from urllib.parse import urlsplit


def mirrored_lock(text, mirror):
    base = urlsplit(mirror)
    if base.scheme != 'https' or not base.netloc or base.query or base.fragment or base.username:
        raise ValueError('Use an HTTPS conda-forge mirror base URL without credentials/query')
    lines = []
    for line in text.splitlines():
        if line.startswith('https://conda.anaconda.org/conda-forge/'):
            tail = line.removeprefix('https://conda.anaconda.org/conda-forge/')
            if '#' not in tail:
                raise ValueError('Refusing a package without its pinned checksum')
            line = mirror.rstrip('/') + '/' + tail
        lines.append(line)
    return '\n'.join(lines) + '\n'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mirror', required=True)
    args = parser.parse_args()
    result = mirrored_lock(args.source.read_text(), args.mirror)
    with args.output.open('x') as output:
        output.write(result)
