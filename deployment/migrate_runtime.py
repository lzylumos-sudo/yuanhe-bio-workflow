#!/usr/bin/env python3
"""Export a compatible runtime configuration to a new file; never overwrite history."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / '00_workbench'))
from yuanhe.runtime import migrated_runtime_config


def migrate(source, destination):
    source, destination = Path(source), Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError('Input and output must be different files; preserve the original')
    raw = source.read_bytes()
    original = json.loads(raw)
    result = migrated_runtime_config(original)
    result['migration'] = {'source_schema': original.get('schema_version', 'legacy-unversioned'),
                           'source_sha256': hashlib.sha256(raw).hexdigest(),
                           'execution_validation': 'not_run'}
    encoded = json.dumps(result, indent=2, allow_nan=False) + '\n'
    # Exclusive creation also protects against concurrent writers and symlink targets.
    with destination.open('x', encoding='utf-8', newline='\n') as handle:
        handle.write(encoded)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        result = migrate(args.input, args.output)
    except (OSError, ValueError, TypeError) as exc:
        parser.exit(2, str(exc) + '\n')
    print(json.dumps({'output': str(args.output), 'schema_version': result['schema_version'],
                      'migration': result['migration']}, indent=2))
    print('Original preserved. Inspect the new file, then run deployment diagnostics before use.')


if __name__ == '__main__':
    main()
