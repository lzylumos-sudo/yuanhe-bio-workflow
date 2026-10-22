"""Administrator archives retain exact committed source and detect tampering."""
import importlib.util
import json
from pathlib import Path
import subprocess
import zipfile
import pytest
ROOT = Path(__file__).resolve().parents[2]

def module(name):
    spec=importlib.util.spec_from_file_location(name, ROOT/'deployment'/(name+'.py'))
    value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value);return value

@pytest.fixture
def repository(tmp_path):
    root=tmp_path/'source';root.mkdir()
    (root/'VERSION').write_text('0.11.0-dev.2\n')
    (root/'payload.txt').write_text('committed source\n')
    def git(*args):
        return subprocess.check_output(['git','-C',str(root),*args],stderr=subprocess.PIPE)
    git('init');git('add','.')
    git('-c','user.name=Test','-c','user.email=test@example.invalid','commit','-m','fixture')
    return root

def test_package_uses_committed_tree_and_preserves_existing_output(repository,tmp_path):
    package=module('package_release');verify=module('verify_release')
    (repository/'payload.txt').write_text('local modification')
    (repository/'private.txt').write_text('not committed')
    archive=tmp_path/'release.zip';result=package.build(archive,root=repository)
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        assert not any(n.endswith('private.txt') for n in z.namelist())
        z.extractall(tmp_path/'extracted')
    root=tmp_path/'extracted/yuanhe-bio-platform-0.11.0-dev.2'
    assert (root/'payload.txt').read_text()=='committed source\n'
    assert verify.verify(root)['source_commit']==result['git_commit']
    original=archive.read_bytes()
    with pytest.raises(FileExistsError):package.build(archive,root=repository)
    assert archive.read_bytes()==original
    (root/'payload.txt').write_text('modified after extraction')
    with pytest.raises(ValueError,match='Modified file'):verify.verify(root)

def test_verifier_rejects_manifest_path_escape(tmp_path):
    (tmp_path/'VERSION').write_text('0.11.0-dev.2\n')
    (tmp_path/'SOURCE_MANIFEST.json').write_text(json.dumps({'schema':'yuanhe.source-package/1','software_version':'0.11.0-dev.2','git_commit':'a'*40,'files':{'../escape':'0'*64}}))
    with pytest.raises(ValueError,match='unsafe file'):module('verify_release').verify(tmp_path)

def test_archive_cannot_inherit_parent_repository_identity(repository):
    spec=importlib.util.spec_from_file_location('release_version',ROOT/'00_workbench/yuanhe/version.py')
    version=importlib.util.module_from_spec(spec);spec.loader.exec_module(version)
    extracted=repository/'unpacked';extracted.mkdir()
    assert version.git_identity(extracted)=={'commit':'unknown','branch':'unknown','dirty':None}
    assert len(version.git_identity(repository)['commit'])==40
