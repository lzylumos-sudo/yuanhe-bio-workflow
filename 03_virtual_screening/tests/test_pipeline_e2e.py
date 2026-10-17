"""End-to-end pipeline on a tiny fixture cut from real GNINA output."""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
REAL_RUN = os.path.join(PROJECT, "runs", "t4")
MEMPHILOT = r"F:\MemPilot\runs\public20_engineering"
ENTRIES = ["P000002_S01_C01_seed20261015", "P000011_S01_C01_seed20261015"]

pytestmark = pytest.mark.integration


def _fixture(root: str) -> str:
    inputs = os.path.join(root, "inputs")
    os.makedirs(inputs, exist_ok=True)
    receptor = os.path.join(REAL_RUN, "reference", "1DRF.pdb")
    if not os.path.exists(receptor):
        pytest.skip("reference receptor from run t4 is missing")
    shutil.copy2(receptor, os.path.join(inputs, "receptor.pdb"))
    box = os.path.join(MEMPHILOT, "protocol_inputs", "box.json")
    if os.path.exists(box):
        shutil.copy2(box, os.path.join(inputs, "box.json"))
    for name in ENTRIES:
        src = os.path.join(MEMPHILOT, "entries", name, "poses.sdf")
        dst_dir = os.path.join(inputs, "entries", name)
        os.makedirs(dst_dir, exist_ok=True)
        shutil.copy2(src, os.path.join(dst_dir, "poses.sdf"))
    return inputs


def _run(root: str, inputs: str, run_id: str) -> subprocess.CompletedProcess:
    code = (
        "from zvss.cli import main;"
        f"raise SystemExit(main(['--root', {root!r}, 'run', '--input', {inputs!r},"
        f" '--target', {os.path.join(PROJECT, 'configs', 'targets', 'generic.yaml')!r},"
        f" '--budget', {os.path.join(PROJECT, 'configs', 'budgets', 'quick.yaml')!r},"
        f" '--run-id', {run_id!r}]))"
    )
    env = dict(os.environ)
    env.setdefault("PYTHONUNBUFFERED", "1")
    return subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env, timeout=1800)


@pytest.fixture(scope="module")
def workspace(tmp_path_factory):
    root = str(tmp_path_factory.mktemp("zvss-root"))
    return root, _fixture(root)


def test_pipeline_imports_every_pose_and_selects(workspace) -> None:
    root, inputs = workspace
    result = _run(root, inputs, "e2e1")
    assert result.returncode == 0, result.stdout + result.stderr
    summary = json.load(open(os.path.join(root, "runs", "e2e1", "stages", "import", "summary.json"), encoding="utf-8"))
    assert summary["poses_imported"] == 6
    assert summary["poses_failed"] == 0
    selection = json.load(open(os.path.join(root, "runs", "e2e1", "stages", "selection", "summary.json"), encoding="utf-8"))
    assert selection["selection_version"].startswith("SEL")
    shortlist = selection["primary"] + selection.get("preliminary", []) + selection["exploratory"]
    assert len(shortlist) >= 1
    # Formal candidates must be a subset of the qualified set: in this fixture the
    # re-dock, region and membrane evidence cannot exist, so primary is empty and
    # the visible candidates are explicitly preliminary.
    qualified = {cid for cid, verdict in selection["qualification"].items() if verdict["formal_qualified"]}
    assert {entry["compound_id"] for entry in selection["primary"]} <= qualified
    assert any(entry.get("qualification") == "preliminary" for entry in selection.get("preliminary", []))
    report = os.path.join(root, "runs", "e2e1", "reports", "report.md")
    assert os.path.exists(report)


def test_second_run_reproduces_selection_and_candidates(workspace) -> None:
    root, inputs = workspace
    result = _run(root, inputs, "e2e2")
    assert result.returncode == 0, result.stdout + result.stderr
    paths = {}
    for run_id in ("e2e1", "e2e2"):
        base = os.path.join(root, "runs", run_id)
        paths[run_id] = (
            json.load(open(os.path.join(base, "stages", "selection", "summary.json"), encoding="utf-8"))["selection_version"],
            open(os.path.join(base, "reports", "candidates.csv"), encoding="utf-8").read(),
        )
    assert paths["e2e1"][0] == paths["e2e2"][0]
    assert paths["e2e1"][1] == paths["e2e2"][1]


def test_md_package_contract(workspace) -> None:
    root, _ = workspace
    md_dir = os.path.join(root, "runs", "e2e1", "md_handoff")
    packages = [d for d in sorted(os.listdir(md_dir)) if os.path.isdir(os.path.join(md_dir, d))]
    if not packages:
        pytest.skip("no MD package produced for this fixture")
    first = os.path.join(md_dir, packages[0])
    for name in ("complex.pdb", "ligand.sdf", "manifest.json", "task.json", "README.md", "checksums.sha256"):
        assert os.path.exists(os.path.join(first, name)), name
    manifest = json.load(open(os.path.join(first, "manifest.json"), encoding="utf-8"))
    assert manifest["ligand"]["resname"] == "LIG"
    assert manifest["ligand"]["atoms"] > 0
    assert manifest["compound_id"].startswith("CMP")
    task = json.load(open(os.path.join(first, "task.json"), encoding="utf-8"))
    codes = {b["code"] for b in (task.get("execution") or {}).get("blockers", [])}
    assert "blocked_missing_membrane_capability" in codes or "blocked_platform_linux_required" in codes
    assert task["execution"]["status"] == "blocked"
    complex_pdb = open(os.path.join(first, "complex.pdb"), encoding="utf-8").read()
    assert "HETATM" in complex_pdb and " LIG " in complex_pdb
    assert "REMARK" in complex_pdb
