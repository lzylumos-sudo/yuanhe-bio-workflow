"""Export legacy SBVS metadata without machine paths or nonstandard lock headers."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def explicit_export(text):
    lines = [line.strip() for line in text.splitlines()]
    urls = [line for line in lines if line.startswith('https://')]
    unexpected = [line for line in lines if line and not (
        line.startswith(('#', 'https://', 'List of packages in environment:'))
        or line == '@EXPLICIT')]
    if not urls or unexpected:
        raise ValueError('Empty or unrecognized explicit export; existing files retained')
    if any('/linux-64/' not in url and '/noarch/' not in url for url in urls):
        raise ValueError('This exporter requires a Linux x86-64 environment')
    return '# platform: linux-64\n# Conda packages only; install pip supplements separately.\n@EXPLICIT\n' + '\n'.join(urls) + '\n'


def main():
    import yaml
    root = Path(sys.argv[1]).resolve()
    prefix = os.environ.get('SBVS_ENV_PREFIX', '/opt/sbvs/env')
    mamba = os.environ.get('SBVS_MICROMAMBA', '/opt/sbvs/vendor/micromamba')
    lock = explicit_export(subprocess.check_output(
        [mamba, 'list', '-p', prefix, '--explicit'], text=True))
    spec = yaml.safe_load(subprocess.check_output(
        [mamba, 'env', 'export', '-p', prefix], text=True))
    if not isinstance(spec, dict) or not spec.get('dependencies'):
        raise ValueError('Missing exported dependencies; existing files retained')
    spec.pop('prefix', None)
    spec['name'] = 'sbvs-legacy'
    snapshot = '# Historical Linux build snapshot, not a portable minimal specification.\n' + yaml.safe_dump(spec, sort_keys=False)
    for name, content in [('conda-linux-64.lock', lock), ('environment.yml', snapshot)]:
        destination = root / 'env' / name
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=destination.parent, delete=False) as f:
            temporary = Path(f.name)
            f.write(content)
        try:
            os.replace(temporary, destination)
        finally:
            temporary.unlink(missing_ok=True)


if __name__ == '__main__':
    main()
