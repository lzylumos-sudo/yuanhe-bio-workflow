"""Adapter contract tests using fake executables only — NOT model inference tests."""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

import pytest

from src.adapters import peptide2mol, reinvent4

pytestmark = pytest.mark.adapter_contract


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def executable(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    path.chmod(0o755)
    return path


def test_peptide2mol_contract_runs_official_prepare_and_inference(tmp_path, monkeypatch):
    pytest.importorskip("rdkit")
    repo = tmp_path / "Peptide2Mol"
    (repo / "notebooks").mkdir(parents=True)
    (repo / "scripts").mkdir()
    (repo / "notebooks" / "deal_with_mol_test_from_pdb_pal_5A.py").write_text(
        "import pathlib,sys\np=pathlib.Path(sys.argv[1])\n"
        "x=next(p.glob('*poc.pdb'))\n(x.with_suffix('.pt')).write_bytes(b'prepared')\n",
        encoding="utf-8",
    )
    (repo / "src").mkdir()
    (repo / "src" / "eval.py").write_text(
        "import pathlib,sys,json\n"
        "arg=lambda p: next(x.split('=',1)[1] for x in sys.argv if x.startswith(p+'='))\n"
        "out=pathlib.Path(json.loads(arg('++model.net.sample.log_dir'))+'_SDF'); out.mkdir(parents=True)\n"
        "(out/'attempt_budget.json').write_text(json.dumps({'requested_attempts':2,'started_attempts':2,'valid':1,'failed':1,'status':'finished'}))\n"
        "(out/'generated_0.sdf').write_text('generated\\n  contract\\n\\n  1  0  0  0  0  0  0  0  0  0999 V2000\\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\\nM  END\\n$$$$\\n')\n",
        encoding="utf-8",
    )
    checkpoint = tmp_path / "model.ckpt"; checkpoint.write_bytes(b"real checkpoint stand-in")
    pocket = tmp_path / "target.pdb"
    pocket.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n", encoding="utf-8")
    monkeypatch.setattr(peptide2mol.platform, "system", lambda: "Linux")

    result = peptide2mol.run(
        pocket_pdb=pocket, checkpoint=checkpoint, checkpoint_sha256=sha(checkpoint),
        output_dir=tmp_path / "run", budget=2, seed=19, executable=sys.executable, repo=repo,
    )

    assert result["status"] == "completed"
    assert result["repo_commit"] == peptide2mol.REPO_COMMIT
    assert len(result["commands"]) == 2
    assert "deal_with_mol_test_from_pdb_pal_5A.py" in result["commands"][0]["command"][1]
    assert result["commands"][1]["command"][1].endswith("src/eval.py")
    assert "++model.net.sample.num_mols=2" in result["commands"][1]["command"]
    assert "++seed=19" in result["commands"][1]["command"]
    assert result["report"] == {"requested": 2, "completed": 2, "sdf_outputs": 1, "valid": 1, "unique": 1, "failure_reasons": {"model_reconstruction_rejected": 1}}
    assert list((tmp_path / "run" / "raw_sdf").glob("*.sdf"))


def test_peptide2mol_wrong_hash_is_nonzero_and_never_executes(tmp_path):
    checkpoint = tmp_path / "model.ckpt"; checkpoint.write_bytes(b"checkpoint")
    pocket = tmp_path / "poc.pdb"; pocket.write_text("", encoding="utf-8")
    output = tmp_path / "out"
    rc = peptide2mol.main([
        "--pocket-pdb", str(pocket), "--checkpoint", str(checkpoint),
        "--checkpoint-sha256", "0" * 64, "--output-dir", str(output),
        "--budget", "1", "--seed", "1", "--executable", sys.executable,
    ])
    manifest = json.loads((output / "manifest.json").read_text())
    assert rc != 0
    assert manifest["status"] == "failed"
    assert manifest["commands"] == []
    assert "mismatch" in manifest["error"]


def test_reinvent4_contract_writes_v4_toml_and_preserves_raw_rows(tmp_path):
    pytest.importorskip("rdkit")
    repo = tmp_path / "REINVENT4"; repo.mkdir()
    checkpoint = tmp_path / "libinvent.prior"; checkpoint.write_bytes(b"real prior stand-in")
    fake = executable(tmp_path / "reinvent", """#!/usr/bin/env python3
import csv,re,sys
cfg=open(sys.argv[-1]).read()
def val(name):
    return re.search(r'^'+name+r' = "(.*)"$', cfg, re.M).group(1)
scaffold=open(val('smiles_file')).read().strip()
with open(val('output_file'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['SMILES','SMILES_state','Scaffold','R-groups','NLL'])
    w.writeheader(); w.writerow({'SMILES':'Cc1ccccc1','SMILES_state':'valid','Scaffold':scaffold,'R-groups':'C','NLL':'1.25'})
print('fake adapter contract invocation')
""")
    out = tmp_path / "run"
    result = reinvent4.run(
        scaffold_smiles="c1ccc([*:1])cc1", checkpoint=checkpoint,
        checkpoint_sha256=sha(checkpoint), output_dir=out, budget=7, seed=123,
        executable=fake, repo=repo,
    )

    config = (out / "sampling.toml").read_text()
    candidates = json.loads((out / "candidates.json").read_text())
    assert result["status"] == "completed"
    assert 'run_type = "sampling"' in config and 'version =' not in config
    assert "num_smiles = 7" in config and "unique_molecules = false" in config
    assert result["commands"][0]["command"] == [str(fake.resolve()), "--seed", "123", str((out / "sampling.toml").resolve())]
    assert candidates["candidates"][0]["SMILES"] == "Cc1ccccc1"
    assert reinvent4.scaffold_identity(candidates["candidates"][0]["Input_SMILES"]) == reinvent4.scaffold_identity("c1ccc([*:1])cc1")
    assert "[*]" in (out / "scaffold.smi").read_text()
    assert "[*:1]" not in (out / "scaffold.smi").read_text()
    assert candidates["candidates"][0]["NLL"] == "1.25"
    assert (out / "raw" / "raw_sampling.csv").is_file()


def test_reinvent4_rejects_bad_hash_before_fake_executable(tmp_path):
    repo = tmp_path / "REINVENT4"; repo.mkdir()
    checkpoint = tmp_path / "prior"; checkpoint.write_bytes(b"prior")
    marker = tmp_path / "was-run"
    fake = executable(tmp_path / "reinvent", f"#!/bin/sh\ntouch {marker}\n")
    output = tmp_path / "out"
    rc = reinvent4.main([
        "--scaffold-smiles", "c1ccc([*:1])cc1", "--checkpoint", str(checkpoint),
        "--checkpoint-sha256", "f" * 64, "--output-dir", str(output),
        "--budget", "2", "--seed", "9", "--executable", str(fake), "--repo", str(repo),
    ])
    manifest = json.loads((output / "manifest.json").read_text())
    assert rc != 0 and not marker.exists()
    assert manifest["commands"] == []


@pytest.mark.parametrize("smiles", ["c1ccccc1", "c1ccc(*)cc1", "c1cc([*:1])cc([*:1])c1", "[*:1]C[*:2]"])
def test_reinvent4_scaffold_attachment_contract_rejects_invalid(smiles):
    pytest.importorskip("rdkit")
    with pytest.raises(ValueError):
        reinvent4.validate_scaffold(smiles)
