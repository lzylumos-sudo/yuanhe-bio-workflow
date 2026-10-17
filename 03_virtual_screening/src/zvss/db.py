"""SQLite persistence: identity, statuses, assessments, cache and event log.

A single file per project (``runtime/zvss.sqlite3``) keeps the run database
portable.  Every stage writes *rows*, never half-written status flags, so an
interrupted stage is visible as ``pending``/``running`` rather than as success.
"""

from __future__ import annotations

import json
import os
import sqlite3
import threading
from contextlib import contextmanager
from typing import Any, Iterable, Iterator, Sequence

from .constants import AssessmentStatus, ExecutionStatus, SelectionStatus
from .util import atomic_write_json, now_iso

SCHEMA_VERSION = 5

SCHEMA = """
CREATE TABLE IF NOT EXISTS schema_info (key TEXT PRIMARY KEY, value TEXT);

CREATE TABLE IF NOT EXISTS run (
    run_id TEXT PRIMARY KEY,
    created TEXT NOT NULL,
    updated TEXT NOT NULL,
    config_hash TEXT,
    code_version TEXT,
    target_id TEXT,
    status TEXT NOT NULL DEFAULT 'pending',
    budget TEXT,
    meta TEXT
);

CREATE TABLE IF NOT EXISTS compound (
    compound_id TEXT PRIMARY KEY,
    canonical_smiles TEXT NOT NULL,
    inchikey TEXT,
    formula TEXT,
    heavy_atoms INTEGER,
    net_charge INTEGER,
    parent_label TEXT,
    first_seen TEXT,
    meta TEXT
);

CREATE TABLE IF NOT EXISTS chemical_state (
    state_id TEXT PRIMARY KEY,
    compound_id TEXT NOT NULL,
    canonical_smiles TEXT NOT NULL,
    formal_charge INTEGER,
    protonation_note TEXT,
    tautomer_note TEXT,
    stereo_note TEXT,
    source TEXT,
    created TEXT,
    meta TEXT,
    FOREIGN KEY (compound_id) REFERENCES compound(compound_id)
);

CREATE TABLE IF NOT EXISTS docking_job (
    job_id TEXT PRIMARY KEY,
    state_id TEXT NOT NULL,
    receptor_version TEXT,
    receptor_sha256 TEXT,
    box_json TEXT,
    box_hash TEXT,
    protocol_json TEXT,
    protocol_hash TEXT,
    seed INTEGER,
    engine TEXT,
    engine_version TEXT,
    model_id TEXT,
    scoring_function TEXT,
    source TEXT,
    created TEXT,
    meta TEXT,
    FOREIGN KEY (state_id) REFERENCES chemical_state(state_id)
);

CREATE TABLE IF NOT EXISTS pose (
    pose_id TEXT PRIMARY KEY,
    job_id TEXT NOT NULL,
    state_id TEXT NOT NULL,
    compound_id TEXT NOT NULL,
    rank_in_job INTEGER,
    record_index INTEGER,
    model_index INTEGER,
    source_file TEXT,
    source_hash TEXT,
    coords_path TEXT,
    coords_hash TEXT,
    scores TEXT,
    exec_status TEXT NOT NULL DEFAULT 'pending',
    assess_status TEXT NOT NULL DEFAULT 'not_evaluated',
    select_status TEXT NOT NULL DEFAULT 'eligible',
    reason_code TEXT,
    reason TEXT,
    created TEXT,
    meta TEXT,
    UNIQUE (job_id, record_index, model_index, coords_hash)
);

CREATE TABLE IF NOT EXISTS assessment (
    run_id TEXT NOT NULL DEFAULT '',
    scope_type TEXT NOT NULL,
    scope_id TEXT NOT NULL,
    stage TEXT NOT NULL,
    metric TEXT NOT NULL,
    status TEXT NOT NULL,
    value REAL,
    value_text TEXT,
    unit TEXT,
    reason_code TEXT,
    reason TEXT,
    evidence_path TEXT,
    config_source TEXT,
    created TEXT,
    meta TEXT,
    PRIMARY KEY (run_id, scope_type, scope_id, stage, metric)
);

CREATE TABLE IF NOT EXISTS cache (
    cache_key TEXT PRIMARY KEY,
    stage TEXT NOT NULL,
    scope_id TEXT,
    input_hash TEXT,
    config_hash TEXT,
    code_version TEXT,
    tool_versions TEXT,
    output_path TEXT,
    output_hash TEXT,
    completed INTEGER NOT NULL DEFAULT 0,
    created TEXT,
    meta TEXT
);

CREATE TABLE IF NOT EXISTS event (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id TEXT,
    ts TEXT,
    level TEXT,
    stage TEXT,
    code TEXT,
    message TEXT,
    payload TEXT
);

CREATE TABLE IF NOT EXISTS md_artifact (
    handoff_id TEXT PRIMARY KEY,
    pose_id TEXT,
    compound_id TEXT,
    selection_tier TEXT,
    package_dir TEXT,
    prep_status TEXT,
    sched_status TEXT,
    handoff_validated INTEGER DEFAULT 0,
    schema_version TEXT,
    created TEXT,
    updated TEXT,
    meta TEXT
);

CREATE INDEX IF NOT EXISTS idx_pose_state ON pose(state_id);
CREATE INDEX IF NOT EXISTS idx_pose_compound ON pose(compound_id);
CREATE INDEX IF NOT EXISTS idx_assess_scope ON assessment(scope_type, scope_id);
CREATE INDEX IF NOT EXISTS idx_event_run ON event(run_id);
"""


def _json(value: Any) -> str | None:
    if value is None:
        return None
    return json.dumps(value, ensure_ascii=False, sort_keys=True, default=str)


class Database:
    """Thin, thread-safe wrapper around a single SQLite file."""

    def __init__(self, path: str):
        self.path = path
        os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
        self._lock = threading.RLock()
        self.run_id = ""
        self._conn = sqlite3.connect(path, check_same_thread=False, timeout=30)
        self._conn.row_factory = sqlite3.Row
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("PRAGMA foreign_keys=ON")
        with self._lock:
            self._conn.executescript(SCHEMA)
            self._migrate()
            self._conn.execute(
                "INSERT OR REPLACE INTO schema_info(key, value) VALUES ('schema_version', ?)",
                (str(SCHEMA_VERSION),),
            )
            self._conn.commit()

    # -- plumbing ---------------------------------------------------------- #
    def _migrate(self) -> None:
        """Bring an older database file up to the current schema.

        Assessments used to be keyed without a run, so two runs over the same
        deterministic compound/pose ids silently overwrote each other's numbers.
        Legacy rows are preserved under an empty run id rather than dropped.
        """
        columns = {row["name"] for row in self.query("PRAGMA table_info(assessment)")}
        if not columns or "run_id" in columns:
            return
        self._conn.executescript(
            """
            ALTER TABLE assessment RENAME TO assessment_legacy;
            CREATE TABLE assessment (
                run_id TEXT NOT NULL DEFAULT '',
                scope_type TEXT NOT NULL,
                scope_id TEXT NOT NULL,
                stage TEXT NOT NULL,
                metric TEXT NOT NULL,
                status TEXT NOT NULL,
                value REAL,
                value_text TEXT,
                unit TEXT,
                reason_code TEXT,
                reason TEXT,
                evidence_path TEXT,
                config_source TEXT,
                created TEXT,
                meta TEXT,
                PRIMARY KEY (run_id, scope_type, scope_id, stage, metric)
            );
            INSERT INTO assessment(run_id, scope_type, scope_id, stage, metric, status, value, value_text,
                    unit, reason_code, reason, evidence_path, config_source, created, meta)
                SELECT '', scope_type, scope_id, stage, metric, status, value, value_text,
                       unit, reason_code, reason, evidence_path, config_source, created, meta
                FROM assessment_legacy;
            DROP TABLE assessment_legacy;
            CREATE INDEX IF NOT EXISTS idx_assess_scope ON assessment(scope_type, scope_id);
            """
        )
        self._conn.commit()

    def close(self) -> None:
        with self._lock:
            self._conn.commit()
            self._conn.close()

    def __enter__(self) -> "Database":
        return self

    def __exit__(self, *exc: Any) -> None:
        self.close()

    @contextmanager
    def tx(self) -> Iterator[sqlite3.Connection]:
        with self._lock:
            try:
                yield self._conn
                self._conn.commit()
            except Exception:
                self._conn.rollback()
                raise

    def query(self, sql: str, params: Sequence[Any] = ()) -> list[dict[str, Any]]:
        with self._lock:
            cursor = self._conn.execute(sql, tuple(params))
            return [dict(row) for row in cursor.fetchall()]

    def one(self, sql: str, params: Sequence[Any] = ()) -> dict[str, Any] | None:
        rows = self.query(sql, params)
        return rows[0] if rows else None

    def execute(self, sql: str, params: Sequence[Any] = ()) -> None:
        with self.tx() as conn:
            conn.execute(sql, tuple(params))

    # -- events ------------------------------------------------------------ #
    def log(
        self,
        message: str,
        *,
        run_id: str | None = None,
        stage: str | None = None,
        level: str = "info",
        code: str | None = None,
        payload: dict | None = None,
    ) -> None:
        self.execute(
            "INSERT INTO event(run_id, ts, level, stage, code, message, payload) VALUES (?,?,?,?,?,?,?)",
            (run_id, now_iso(), level, stage, code, message, _json(payload)),
        )

    def events(self, run_id: str | None = None, limit: int = 500) -> list[dict[str, Any]]:
        if run_id:
            return self.query(
                "SELECT * FROM event WHERE run_id=? ORDER BY id DESC LIMIT ?", (run_id, limit)
            )
        return self.query("SELECT * FROM event ORDER BY id DESC LIMIT ?", (limit,))

    # -- run --------------------------------------------------------------- #
    def upsert_run(self, run_id: str, **fields: Any) -> None:
        stamp = now_iso()
        existing = self.one("SELECT run_id FROM run WHERE run_id=?", (run_id,))
        if existing:
            keys = ", ".join(f"{key}=?" for key in fields)
            self.execute(
                f"UPDATE run SET updated=?, {keys} WHERE run_id=?",
                (stamp, *[_json(v) if isinstance(v, (dict, list)) else v for v in fields.values()], run_id),
            )
        else:
            columns = ["run_id", "created", "updated", *fields.keys()]
            values = [run_id, stamp, stamp, *[_json(v) if isinstance(v, (dict, list)) else v for v in fields.values()]]
            placeholders = ",".join("?" for _ in columns)
            self.execute(
                f"INSERT INTO run({','.join(columns)}) VALUES ({placeholders})", values
            )

    def get_run(self, run_id: str) -> dict[str, Any] | None:
        row = self.one("SELECT * FROM run WHERE run_id=?", (run_id,))
        if row and row.get("meta"):
            row["meta"] = json.loads(row["meta"])
        return row

    def list_runs(self, limit: int = 50) -> list[dict[str, Any]]:
        return self.query("SELECT * FROM run ORDER BY created DESC LIMIT ?", (limit,))

    # -- identity ---------------------------------------------------------- #
    def upsert_compound(self, record: dict[str, Any]) -> None:
        self.execute(
            """INSERT INTO compound(compound_id, canonical_smiles, inchikey, formula, heavy_atoms,
                                    net_charge, parent_label, first_seen, meta)
               VALUES (?,?,?,?,?,?,?,?,?)
               ON CONFLICT(compound_id) DO UPDATE SET
                 parent_label=COALESCE(compound.parent_label, excluded.parent_label),
                 meta=COALESCE(excluded.meta, compound.meta)""",
            (
                record["compound_id"], record["canonical_smiles"], record.get("inchikey"),
                record.get("formula"), record.get("heavy_atoms"), record.get("net_charge"),
                record.get("parent_label"), now_iso(), _json(record.get("meta")),
            ),
        )

    def upsert_state(self, record: dict[str, Any]) -> None:
        self.execute(
            """INSERT INTO chemical_state(state_id, compound_id, canonical_smiles, formal_charge,
                    protonation_note, tautomer_note, stereo_note, source, created, meta)
               VALUES (?,?,?,?,?,?,?,?,?,?)
               ON CONFLICT(state_id) DO UPDATE SET
                 source=COALESCE(excluded.source, chemical_state.source),
                 meta=COALESCE(excluded.meta, chemical_state.meta)""",
            (
                record["state_id"], record["compound_id"], record["canonical_smiles"],
                record.get("formal_charge"), record.get("protonation_note"), record.get("tautomer_note"),
                record.get("stereo_note"), record.get("source"), now_iso(), _json(record.get("meta")),
            ),
        )

    def upsert_job(self, record: dict[str, Any]) -> None:
        self.execute(
            """INSERT INTO docking_job(job_id, state_id, receptor_version, receptor_sha256, box_json,
                    box_hash, protocol_json, protocol_hash, seed, engine, engine_version, model_id,
                    scoring_function, source, created, meta)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
               ON CONFLICT(job_id) DO UPDATE SET
                 source=COALESCE(excluded.source, docking_job.source),
                 meta=COALESCE(excluded.meta, docking_job.meta)""",
            (
                record["job_id"], record["state_id"], record.get("receptor_version"),
                record.get("receptor_sha256"), _json(record.get("box")), record.get("box_hash"),
                _json(record.get("protocol")), record.get("protocol_hash"), record.get("seed"),
                record.get("engine"), record.get("engine_version"), record.get("model_id"),
                record.get("scoring_function"), record.get("source"), now_iso(), _json(record.get("meta")),
            ),
        )

    def upsert_pose(self, record: dict[str, Any]) -> str:
        existing = self.one(
            "SELECT pose_id FROM pose WHERE job_id=? AND record_index IS ? AND model_index IS ? AND coords_hash=?",
            (record["job_id"], record.get("record_index"), record.get("model_index"), record.get("coords_hash")),
        )
        pose_key = existing["pose_id"] if existing else record["pose_id"]
        self.execute(
            """INSERT INTO pose(pose_id, job_id, state_id, compound_id, rank_in_job, record_index,
                    model_index, source_file, source_hash, coords_path, coords_hash, scores,
                    exec_status, assess_status, select_status, reason_code, reason, created, meta)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
               ON CONFLICT(pose_id) DO UPDATE SET
                 scores=excluded.scores,
                 exec_status=excluded.exec_status,
                 coords_path=COALESCE(excluded.coords_path, pose.coords_path),
                 meta=COALESCE(excluded.meta, pose.meta)""",
            (
                pose_key, record["job_id"], record["state_id"], record["compound_id"],
                record.get("rank_in_job"), record.get("record_index"), record.get("model_index"),
                record.get("source_file"), record.get("source_hash"), record.get("coords_path"),
                record.get("coords_hash"), _json(record.get("scores")),
                str(record.get("exec_status", ExecutionStatus.COMPLETED)),
                str(record.get("assess_status", AssessmentStatus.NOT_EVALUATED)),
                str(record.get("select_status", SelectionStatus.ELIGIBLE)),
                record.get("reason_code"), record.get("reason"), now_iso(), _json(record.get("meta")),
            ),
        )
        return pose_key

    def set_pose_status(
        self,
        pose_id: str,
        *,
        exec_status: str | None = None,
        assess_status: str | None = None,
        select_status: str | None = None,
        reason_code: str | None = None,
        reason: str | None = None,
    ) -> None:
        fields: list[str] = []
        values: list[Any] = []
        for column, value in (
            ("exec_status", exec_status),
            ("assess_status", assess_status),
            ("select_status", select_status),
            ("reason_code", reason_code),
            ("reason", reason),
        ):
            if value is not None:
                fields.append(f"{column}=?")
                values.append(str(value))
        if not fields:
            return
        values.append(pose_id)
        self.execute(f"UPDATE pose SET {', '.join(fields)} WHERE pose_id=?", values)

    # -- assessments ------------------------------------------------------- #
    def record_assessment(
        self,
        scope_type: str,
        scope_id: str,
        stage: str,
        metric: str,
        value: Any = None,
        status: AssessmentStatus | str | None = None,
        reason_code: str | None = None,
        *,
        value_text: str | None = None,
        unit: str | None = None,
        reason: str | None = None,
        note: str | None = None,
        evidence_path: str | None = None,
        config_source: str | None = None,
        meta: dict | None = None,
    ) -> None:
        # Non-numeric values (lists, fingerprints, strings) are stored as text so a
        # structured metric can never crash the run.
        if value is not None and not isinstance(value, (int, float, bool)):
            if value_text is None:
                value_text = str(value)
            value = None
        if status is None:
            status = AssessmentStatus.NOT_EVALUATED
        if reason is None and note is not None:
            reason = note
        self.execute(
            """INSERT INTO assessment(run_id, scope_type, scope_id, stage, metric, status, value, value_text,
                    unit, reason_code, reason, evidence_path, config_source, created, meta)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
               ON CONFLICT(run_id, scope_type, scope_id, stage, metric) DO UPDATE SET
                 status=excluded.status, value=excluded.value, value_text=excluded.value_text,
                 unit=excluded.unit, reason_code=excluded.reason_code, reason=excluded.reason,
                 evidence_path=excluded.evidence_path, config_source=excluded.config_source,
                 created=excluded.created, meta=excluded.meta""",
            (
                getattr(self, "run_id", "") or "", scope_type, scope_id, stage, metric, str(status), value,
                value_text, unit, reason_code, reason, evidence_path, config_source, now_iso(), _json(meta),
            ),
        )

    def assessments(
        self,
        scope_id: str | None = None,
        scope_type: str | None = None,
        stage: str | None = None,
        run_id: str | None = None,
    ) -> list[dict[str, Any]]:
        clauses: list[str] = []
        params: list[Any] = []
        effective_run = run_id if run_id is not None else getattr(self, "run_id", "")
        if effective_run:
            clauses.append("run_id=?")
            params.append(effective_run)
        if scope_id:
            clauses.append("scope_id=?")
            params.append(scope_id)
        if scope_type:
            clauses.append("scope_type=?")
            params.append(scope_type)
        if stage:
            clauses.append("stage=?")
            params.append(stage)
        where = f" WHERE {' AND '.join(clauses)}" if clauses else ""
        rows = self.query(f"SELECT * FROM assessment{where} ORDER BY scope_id, stage, metric", params)
        for row in rows:
            if row.get("meta"):
                row["meta"] = json.loads(row["meta"])
        return rows

    def metric_map(self, scope_type: str, scope_id: str, stage: str) -> dict[str, dict[str, Any]]:
        return {row["metric"]: row for row in self.assessments(scope_id, scope_type, stage)}

    # -- cache ------------------------------------------------------------- #
    def cache_get(self, key: str) -> dict[str, Any] | None:
        row = self.one("SELECT * FROM cache WHERE cache_key=? AND completed=1", (key,))
        if not row:
            return None
        if row.get("output_path") and not os.path.exists(row["output_path"]):
            return None
        return row

    def cache_put(
        self,
        *,
        key: str,
        stage: str,
        scope_id: str | None,
        input_hash: str,
        config_hash: str,
        code_version: str,
        tool_versions: dict,
        output_path: str | None,
        output_hash: str | None,
        completed: bool,
        meta: dict | None = None,
    ) -> None:
        self.execute(
            """INSERT INTO cache(cache_key, stage, scope_id, input_hash, config_hash, code_version,
                    tool_versions, output_path, output_hash, completed, created, meta)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
               ON CONFLICT(cache_key) DO UPDATE SET
                 output_path=excluded.output_path, output_hash=excluded.output_hash,
                 completed=excluded.completed, created=excluded.created, meta=excluded.meta""",
            (
                key, stage, scope_id, input_hash, config_hash, code_version, _json(tool_versions),
                output_path, output_hash, 1 if completed else 0, now_iso(), _json(meta),
            ),
        )

    def invalidate(self, stage: str | None = None) -> int:
        if stage:
            rows = self.query("SELECT cache_key FROM cache WHERE stage=?", (stage,))
            self.execute("DELETE FROM cache WHERE stage=?", (stage,))
        else:
            rows = self.query("SELECT cache_key FROM cache")
            self.execute("DELETE FROM cache")
        return len(rows)

    # -- md artifacts ------------------------------------------------------ #
    def upsert_md_artifact(self, record: dict[str, Any]) -> None:
        stamp = now_iso()
        self.execute(
            """INSERT INTO md_artifact(handoff_id, pose_id, compound_id, selection_tier, package_dir,
                    prep_status, sched_status, handoff_validated, schema_version, created, updated, meta)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
               ON CONFLICT(handoff_id) DO UPDATE SET
                 package_dir=excluded.package_dir, prep_status=excluded.prep_status,
                 sched_status=excluded.sched_status, handoff_validated=excluded.handoff_validated,
                 updated=excluded.updated, meta=excluded.meta""",
            (
                record["handoff_id"], record.get("pose_id"), record.get("compound_id"),
                record.get("selection_tier"), record.get("package_dir"), record.get("prep_status"),
                record.get("sched_status"), int(bool(record.get("handoff_validated"))),
                record.get("schema_version"), stamp, stamp, _json(record.get("meta")),
            ),
        )

    def md_artifacts(self, run_id: str | None = None) -> list[dict[str, Any]]:
        rows = self.query("SELECT * FROM md_artifact ORDER BY created")
        for row in rows:
            if row.get("meta"):
                row["meta"] = json.loads(row["meta"])
        return rows


def open_db(path: str) -> Database:
    return Database(path)


def dump_state(db: Database, path: str) -> str:
    """Human-inspectable snapshot used by tests and hand-off records."""
    payload = {
        "schema_version": SCHEMA_VERSION,
        "runs": db.query("SELECT * FROM run"),
        "compounds": db.query("SELECT COUNT(*) AS n FROM compound"),
        "states": db.query("SELECT COUNT(*) AS n FROM chemical_state"),
        "jobs": db.query("SELECT COUNT(*) AS n FROM docking_job"),
        "poses": db.query("SELECT exec_status, assess_status, select_status, COUNT(*) AS n FROM pose GROUP BY 1,2,3"),
        "assessments": db.query(
            "SELECT scope_type, stage, status, COUNT(*) AS n FROM assessment GROUP BY 1,2,3 ORDER BY 1,2,3"
        ),
    }
    return atomic_write_json(path, payload)
