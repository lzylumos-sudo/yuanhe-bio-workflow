"""Reject non-portable exports without changing recorded package identities."""
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    'environment_export', ROOT / '02_docking/scripts/export_environment.py')
exporter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exporter)


def test_export_preserves_urls_and_hashes_but_removes_local_banner():
    url = 'https://conda.anaconda.org/conda-forge/linux-64/example-1.0.conda#abc123'
    output = exporter.explicit_export('List of packages in environment: "/local/env"\n\n' + url)
    records = [line for line in output.splitlines() if not line.startswith('#')]
    assert records == ['@EXPLICIT', url]
    assert '/local/env' not in output


@pytest.mark.parametrize('text', [
    '', 'Error: export failed',
    'https://conda.anaconda.org/conda-forge/osx-arm64/example.conda',
    'https://conda.anaconda.org/conda-forge/linux-64/example.conda\nunrecognized output',
])
def test_invalid_exports_fail_closed(text):
    with pytest.raises(ValueError):
        exporter.explicit_export(text)
