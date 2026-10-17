"""Single place where all paths are resolved (plan §05 RuntimePaths). No hard-coded roots."""
from __future__ import annotations
import os
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class RuntimePaths:
    home: Path          # SBVS_HOME: project root on the (external) archive volume
    project_root: Path  # SBVS_PROJECT_ROOT: where projects/<id>/ live
    scratch: Path       # SBVS_SCRATCH: local high-frequency temp (should be local ext4/APFS)

    @classmethod
    def resolve(cls, home: str | os.PathLike | None = None) -> "RuntimePaths":
        h = Path(home or os.environ.get("SBVS_HOME") or Path(__file__).resolve().parents[2]).resolve()
        pr = Path(os.environ.get("SBVS_PROJECT_ROOT") or h / "projects").resolve()
        sc = Path(os.environ.get("SBVS_SCRATCH") or Path.home() / ".sbvs2_scratch").resolve()
        return cls(h, pr, sc)

    def project(self, project_id: str) -> Path:
        if not project_id or any(c in project_id for c in "/\\") or project_id in (".", ".."):
            raise ValueError(f"unsafe project id: {project_id!r}")
        return self.project_root / project_id

    def run_dir(self, project_id: str, run_id: str) -> Path:
        if any(c in run_id for c in "/\\") or run_id in (".", ".."):
            raise ValueError(f"unsafe run id: {run_id!r}")
        return self.project(project_id) / "runs" / run_id


PROJECT_LAYOUT = ("00_original", "01_preparation", "runs", "exports")
RUN_LAYOUT = ("00_run_record", "01_logs", "02_work_products", "03_results", "04_analysis")


def init_project(paths: RuntimePaths, project_id: str) -> Path:
    p = paths.project(project_id)
    for d in PROJECT_LAYOUT:
        (p / d).mkdir(parents=True, exist_ok=True)
    return p


def init_run(paths: RuntimePaths, project_id: str, run_id: str) -> Path:
    r = paths.run_dir(project_id, run_id)
    if r.exists():
        raise FileExistsError(f"run already exists, runs are never overwritten: {r}")
    for d in RUN_LAYOUT:
        (r / d).mkdir(parents=True)
    return r


def safe_join(base: Path, rel: str) -> Path:
    """Reject path traversal for user supplied relative names."""
    target = (base / rel).resolve()
    if base.resolve() not in target.parents and target != base.resolve():
        raise ValueError(f"path escapes base directory: {rel!r}")
    return target
