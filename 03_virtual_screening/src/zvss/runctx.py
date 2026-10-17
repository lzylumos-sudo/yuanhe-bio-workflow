"""Run context: directories, config snapshot, database handle, stage bookkeeping."""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, field
from typing import Any

from . import __version__
from .config import HostConfig, ResolvedConfig, build_resolved, load_budget, load_host_config, load_stages, load_target_profile
from .db import Database
from .util import atomic_write_json, atomic_write_text, ensure_dir, now_iso, utc_stamp


def make_run_id(label: str | None = None) -> str:
    stamp = utc_stamp().replace("-", "").replace(":", "")
    suffix = f"_{label}" if label else ""
    return f"{stamp}{suffix}"


@dataclass
class RunContext:
    run_id: str
    root: str
    config: ResolvedConfig
    db: Database
    extra: dict[str, Any] = field(default_factory=dict)

    # -- paths ------------------------------------------------------------- #
    @property
    def run_dir(self) -> str:
        return os.path.join(self.root, "runs", self.run_id)

    def path(self, *parts: str) -> str:
        target = os.path.join(self.run_dir, *parts)
        ensure_dir(os.path.dirname(target)) if os.path.splitext(target)[1] else ensure_dir(target)
        return target

    def stage_dir(self, *parts: str) -> str:
        return ensure_dir(os.path.join(self.run_dir, "stages", *parts))

    def inputs_dir(self) -> str:
        return ensure_dir(os.path.join(self.run_dir, "inputs"))

    def poses_dir(self) -> str:
        return ensure_dir(os.path.join(self.run_dir, "poses"))

    def reports_dir(self) -> str:
        return ensure_dir(os.path.join(self.run_dir, "reports"))

    def md_dir(self) -> str:
        return ensure_dir(os.path.join(self.run_dir, "md_handoff"))

    def ref_dir(self) -> str:
        return ensure_dir(os.path.join(self.run_dir, "reference"))

    # -- provenance -------------------------------------------------------- #
    def write_resolved_config(self) -> str:
        path = os.path.join(self.run_dir, "resolved_config.yaml")
        return atomic_write_text(path, self.config.to_yaml())

    def write_manifest(self, payload: dict[str, Any]) -> str:
        return atomic_write_json(os.path.join(self.run_dir, "run_manifest.json"), payload)

    def log(self, message: str, *, stage: str | None = None, level: str = "info", code: str | None = None, payload: dict | None = None, run_id: str | None = None) -> None:
        self.db.log(message, run_id=run_id or self.run_id, stage=stage, level=level, code=code, payload=payload)

    def archive_input(self, path: str) -> str:
        """Copy an input into the run for provenance; keep the original read-only."""
        if not os.path.exists(path):
            return path
        target_dir = self.inputs_dir()
        base = os.path.basename(path)
        target = os.path.join(target_dir, base)
        if not os.path.exists(target):
            shutil.copy2(path, target)
        return target


def load_context(
    *,
    run_id: str | None = None,
    root: str | None = None,
    host_config: str | None = None,
    target_profile: str | None = None,
    budget_config: str | None = None,
    stage_config: str | None = None,
    overrides: dict[str, Any] | None = None,
    create: bool = True,
) -> RunContext:
    root = os.path.abspath(root or os.environ.get("ZVSS_ROOT") or os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
    host = load_host_config(host_config)
    if not target_profile:
        target_profile = os.path.join(root, "configs", "targets", "generic.yaml")
    if not budget_config:
        budget_config = os.path.join(root, "configs", "budgets", "full.yaml")
    target = load_target_profile(target_profile, structures_dir=host.maybe("structures_dir") or "")
    budget = load_budget(budget_config)
    stages = load_stages(stage_config)
    resolved = build_resolved(
        host=host, target=target, budget=budget, stages=stages,
        code_version=__version__, overrides=overrides or {},
    )
    run_id = run_id or make_run_id()
    runtime_dir = ensure_dir(os.path.join(root, "runtime"))
    db = Database(os.path.join(runtime_dir, "zvss.sqlite3"))
    db.run_id = run_id
    ctx = RunContext(run_id=run_id, root=root, config=resolved, db=db)
    if create:
        ensure_dir(ctx.run_dir)
        ctx.db.upsert_run(
            run_id,
            config_hash=resolved.digest(),
            code_version=__version__,
            target_id=target.target_id,
            status="initialised",
            budget=budget.budget_id,
            meta={"created": now_iso()},
        )
        ctx.write_resolved_config()
    return ctx
