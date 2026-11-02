from pathlib import Path
import subprocess

from zlis import pymol_scene


def test_selected_renderer_environment_and_png_are_required(tmp_path, monkeypatch):
    seen = []
    monkeypatch.setenv('PYTHONPATH', '/wrong/abi')
    monkeypatch.setenv('ZLIS_NO_PNG', '1')
    def run(argv, **kwargs):
        seen.append((argv, kwargs))
        (tmp_path / 'scene.pse').write_bytes(b'fixture scene')
        return subprocess.CompletedProcess(argv, 0, '', '')
    monkeypatch.setattr(pymol_scene.subprocess, 'run', run)
    result = pymol_scene.build_scene(tmp_path, executable='/isolated/pymol/bin/python')
    assert not result['ok']  # A successful process without a PNG is incomplete.
    argv, kwargs = seen[0]
    assert argv[0] == '/isolated/pymol/bin/python'
    assert 'PYTHONPATH' not in kwargs['env'] and 'ZLIS_NO_PNG' not in kwargs['env']
    (tmp_path / 'binding_site_main.png').write_bytes(b'fixture image')
    assert pymol_scene.build_scene(tmp_path, executable=argv[0])['ok']


def test_missing_renderer_returns_actionable_failed_result(tmp_path):
    result = pymol_scene.build_scene(tmp_path, executable=str(tmp_path / 'missing-python'))
    assert not result['ok'] and 'missing-python' in result['stderr']
