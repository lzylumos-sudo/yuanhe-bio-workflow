"""Stage orchestration with resumable, content-addressed caching."""

from __future__ import annotations

import os
import json
import traceback
from typing import Any, Callable

from . import __version__
from . import compat
from .constants import ExecutionStatus, Reason
from .ids import sha256_file, sha256_text
from .util import atomic_write_json, atomic_write_text, ensure_dir, now_iso, read_json, read_text, write_pdb

STAGE_ORDER = [
    "target",
    "import",
    "pose_qa",
    "mechanism",
    "motion",
    "redock",
    "strain",
    # membrane is a property of the target/pose environment, not a result of
    # selection, and candidate qualification reads it; it must run before
    # scoring/selection or the requirement can only ever be "not evaluated".
    "membrane",
    "scoring",
    "selection",
    "md",
    "report",
]


def live_counts(rows, selection=None, qa_finished=False):
    total={r.get('compound_id') for r in rows if r.get('compound_id')}
    counts={'imported':len(total)}
    if qa_finished:
        usable={r['compound_id'] for r in rows if r.get('compound_id') and r.get('assess_status') in {'pass','flag'}}
        pending={r['compound_id'] for r in rows if r.get('compound_id') and r.get('assess_status') not in {'pass','flag','fail'}}-usable
        counts.update(qa_usable=len(usable),qa_pending=len(pending),qa_failed=len(total-usable-pending))
    if selection is not None:
        formal={r['compound_id'] for r in selection.get('primary',[])}
        prelim={r['compound_id'] for k in ('preliminary','exploratory') for r in selection.get(k,[])}-formal
        counts.update(remaining=len(formal|prelim),formal=len(formal),preliminary=len(prelim))
    return counts



# --------------------------------------------------------------------------- #
# target preparation
# --------------------------------------------------------------------------- #


def prepare_target(ctx) -> dict[str, Any]:
    """Fetch/verify the reference structure, map residues and derive the receptor."""
    from .receptor import pdbio, sifts
    from .util import parse_pdb_text

    profile = ctx.config.target
    cache_dir = ensure_dir(os.path.join(ctx.root, "runtime", "reference_cache"))
    summary: dict[str, Any] = {
        "target_id": profile.target_id,
        "pdb_id": profile.primary.pdb_id,
        "chain": profile.primary.chain,
        "uniprot": profile.uniprot,
    }
    structure_path = ""
    fetch_record: dict[str, Any] = {}
    local = ""
    if profile.structures_dir and os.path.isdir(profile.structures_dir):
        stem = profile.primary.pdb_id
        for name in sorted(os.listdir(profile.structures_dir)):
            lowered = name.lower()
            if lowered.startswith(stem.lower()) and lowered.endswith((".cif", ".pdb", ".ent")):
                local = os.path.join(profile.structures_dir, name)
                break
    cached = ""
    for name in sorted(os.listdir(cache_dir)) if os.path.isdir(cache_dir) else []:
        if name.lower().startswith(profile.primary.pdb_id.lower()) and name.lower().endswith((".cif", ".pdb")):
            cached = os.path.join(cache_dir, name)
            break
    if local and os.path.exists(local):
        structure_path = local
        fetch_record = {"status": "local", "path": local, "sha256": sha256_file(local)}
    elif cached:
        structure_path = cached
        fetch_record = {"status": "cached", "path": cached, "sha256": sha256_file(cached)}
    else:
        fetch_record = sifts.fetch_structure(profile.primary.pdb_id, cache_dir)
        if fetch_record.get("status") != "ok":
            summary["status"] = "blocked_missing_input"
            summary["reason_code"] = str(Reason.REFERENCE_POSE_MISSING)
            summary["fetch"] = fetch_record
            return summary
        structure_path = fetch_record["path"]

    sequence = profile.native_sequence
    if not sequence and profile.uniprot:
        fasta = sifts.fetch_uniprot_sequence(profile.uniprot, cache_dir)
        sequence = fasta.get("sequence", "")
        summary["uniprot_fetch"] = {k: v for k, v in fasta.items() if k != "sequence"}
    summary["native_sequence_length"] = len(sequence)

    pdb_path = os.path.join(ctx.ref_dir(), f"{profile.primary.pdb_id}.pdb")
    if structure_path.lower().endswith((".cif", ".mmcif")):
        _cif_to_pdb(structure_path, pdb_path)
    else:
        with open(structure_path, "r", encoding="utf-8", errors="replace") as handle:
            write_pdb(pdb_path, parse_pdb_text(handle.read()).atoms)

    mapping_report = None
    if sequence:
        try:
            mapping_report = sifts.build_mapping(
                structure_path,
                pdb_id=profile.primary.pdb_id,
                chain=profile.primary.chain,
                accession=profile.uniprot,
                native_sequence=sequence,
            )
            paths = sifts.write_mapping(mapping_report, ctx.ref_dir())
            summary["mapping"] = mapping_report.summary()
            summary["mapping_paths"] = paths
            summary["mapping_provenance"] = mapping_report.provenance
            summary["mutations"] = mapping_report.mutations[:40]
            summary["mapping_warnings"] = mapping_report.warnings
        except Exception as exc:  # noqa: BLE001
            summary["mapping_error"] = f"{type(exc).__name__}: {exc}"
    else:
        summary["mapping"] = {"status": "skipped", "reason": "no native sequence available"}

    native_to_auth: dict[str, int] = {}
    if mapping_report:
        for row in mapping_report.rows:
            if row.uniprot_resseq is not None and row.auth_seq_id is not None:
                native_to_auth[str(row.uniprot_resseq)] = row.auth_seq_id

    inventory = pdbio.inventory(pdb_path)
    receptor_path = os.path.join(ctx.ref_dir(), "receptor.pdb")
    complex_reference_path = os.path.join(ctx.ref_dir(), "complex_full.pdb")
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        atomic_write_text(complex_reference_path, handle.read())

    # Restrict the receptor to residues that the verified mapping places inside
    # the declared native range, so expression fusions and tags never reach
    # docking or contact analysis. Without a mapping (generic targets) the
    # chain filter alone applies and the audit reports that honestly.
    keep_residues = None
    if mapping_report and profile.construct.native_range:
        low, high = profile.construct.native_range
        keep_residues = sorted(
            {
                (row.auth_seq_id and profile.primary.chain, int(row.auth_seq_id))
                for row in mapping_report.rows
                if row.auth_seq_id is not None
                and row.uniprot_resseq is not None
                and low <= int(row.uniprot_resseq) <= high
            }
        )

    derived = pdbio.derive_receptor(
        pdb_path,
        receptor_path,
        keep_chains=[profile.primary.chain] if profile.primary.chain else None,
        keep_residues=keep_residues,
        recipe=(
            "native_domain_only_strip_water_ions_lipids_ligands_and_other_chains"
            if keep_residues
            else "strip_water_ions_lipids_ligands_and_other_chains"
        ),
        coordinate_frame="deposited_frame",
    )
    derived["inventory"] = inventory.to_dict()

    from .receptor import audit

    receptor_audit = audit.audit_receptor(
        original_path=pdb_path,
        receptor_path=receptor_path,
        reference_complex_path=complex_reference_path,
        chain=profile.primary.chain,
        native_range=profile.construct.native_range,
        mapping_rows=[
            {
                "auth_seq_id": row.auth_seq_id,
                "uniprot_resseq": row.uniprot_resseq,
                "resname": row.resname,
                "uniprot_resname": row.uniprot_resname,
            }
            for row in (mapping_report.rows if mapping_report else [])
        ],
        mutations=profile.construct.mutations,
        c_terminal_tag=profile.construct.ctag,
        strip_log=derived.get("strip_log"),
    )
    derived["audit_paths"] = audit.write_audit(receptor_audit, os.path.join(ctx.ref_dir(), "audit"))
    derived["audit_status"] = receptor_audit["status"]
    atomic_write_json(os.path.join(ctx.stage_dir("target"), "receptor.json"), derived)
    atomic_write_json(os.path.join(ctx.stage_dir("import"), "receptor.json"), derived)

    reference = {
        "path": pdb_path,
        "source": structure_path,
        "pdb_id": profile.primary.pdb_id,
        "chain": profile.primary.chain,
        "native_to_auth": native_to_auth,
        "reference_ligand": _reference_ligand(pdb_path, profile),
        "sha256": sha256_file(pdb_path),
    }
    atomic_write_json(os.path.join(ctx.ref_dir(), "reference.json"), reference)
    summary["receptor"] = {k: v for k, v in derived.items() if k != "inventory"}
    summary["receptor_audit"] = {
        "status": receptor_audit["status"],
        "paths": derived["audit_paths"],
        "checks": receptor_audit["checks"],
        "original_atoms": receptor_audit["original"]["atoms"],
        "receptor_atoms": receptor_audit["receptor"]["atoms"],
        "receptor_residues": receptor_audit["receptor"]["residues"],
        "receptor_sha256": receptor_audit["receptor"]["sha256"],
        "segments": receptor_audit["original"]["segments"],
        "mutations": receptor_audit["mutations"],
        "reference_complex": receptor_audit["reference_complex"]["path"],
    }
    summary["inventory"] = {
        "chains": inventory.chains,
        "residues_per_chain": inventory.residues_per_chain,
        "hetero": [g.to_dict() for g in inventory.hetero][:20],
        "water_count": inventory.water_count,
        "heavy_atom_count": inventory.heavy_atom_count,
    }
    summary["status"] = "ok" if mapping_report else "partial"
    summary["fetch"] = fetch_record
    return summary


def _cif_to_pdb(cif_path: str, pdb_path: str) -> str:
    import gemmi

    structure = gemmi.read_structure(cif_path)
    structure.setup_entities()
    structure.write_pdb(pdb_path)
    return pdb_path


def _reference_ligand(pdb_path: str, profile) -> dict[str, Any]:
    """Coordinates of the first configured reference ligand, if present."""
    from .util import parse_pdb_text

    if not profile.region or not profile.region.reference_poses:
        return {"status": "not_configured"}
    structure = parse_pdb_text(read_text(pdb_path), keep_hydrogens=False)
    for pose in profile.region.reference_poses:
        atoms = [
            atom for atom in structure.atoms
            if atom.record == "HETATM"
            and atom.resname.upper() == (pose.resname or "").upper()
            and (not pose.chain or atom.chain == pose.chain)
            and (pose.native_resseq is None or True)
        ]
        if atoms:
            return {
                "status": "ok",
                "label": f"{pose.resname}:{pose.chain}",
                "source": os.path.basename(pdb_path),
                "note": pose.note,
                "atoms": [[a.x, a.y, a.z] for a in atoms if not a.is_hydrogen()],
            }
    return {"status": "not_found", "reason": "configured reference ligand is not present in the structure"}


# --------------------------------------------------------------------------- #
# individual stages
# --------------------------------------------------------------------------- #


def stage_import(ctx, inputs: dict[str, Any]) -> dict[str, Any]:
    from .ingest.discover import discover
    from .ingest.import_run import import_inputs

    roots = inputs.get("input_dirs") or []
    manifest = inputs.get("manifest")
    bundle = discover(
        roots,
        manifest=manifest,
        recursive=bool(inputs.get("recursive", True)),
        hint_files=inputs.get("hint_files"),
    )
    report = import_inputs(ctx, bundle, hints=inputs.get("hints") or {}, selector=inputs.get("selector"))
    return report.to_dict()


def stage_pose_qa(ctx) -> dict[str, Any]:
    from .qa import pose_qa

    return pose_qa.run(ctx)


def stage_mechanism(ctx) -> dict[str, Any]:
    from .mech import interactions

    return interactions.run(ctx)


def stage_motion(ctx) -> dict[str, Any]:
    from .motion import interference

    return interference.run(ctx)


def stage_redock(ctx, selection: dict[str, Any] | None = None) -> dict[str, Any]:
    from .redock import gnina

    return gnina.run_redock(ctx, _redock_plan(ctx, selection))


def _redock_plan(ctx, selection) -> list[dict[str, Any]]:
    from . import data

    budget = ctx.config.budget
    poses = data.pose_rows(ctx)
    by_state: dict[str, dict[str, Any]] = {}
    for pose in poses:
        by_state.setdefault(pose["state_id"], pose)
    plan = []
    for state_id in sorted(by_state)[: int(getattr(budget, "redock_max_parents", 200))]:
        pose = by_state[state_id]
        plan.append(
            {
                "state_id": state_id,
                "compound_id": pose.get("compound_id"),
                "pose_ids": [pose["pose_id"]],
                "box": pose.get("box") or {},
                "ligand_path": pose.get("coords_path")
                or os.path.join(ctx.poses_dir(), f"{pose['pose_id']}.sdf"),
                "seeds": list(getattr(budget, "seeds", [17, 41, 97])),
            }
        )
    return plan


def stage_strain(ctx) -> dict[str, Any]:
    from .strain import strain

    return strain.run(ctx)


def stage_scoring(ctx) -> dict[str, Any]:
    from .scoring import aggregate

    return aggregate.run(ctx)


def stage_selection(ctx, ranking: dict[str, Any]) -> dict[str, Any]:
    from .select import pipeline

    pipeline.prime_fingerprints(ctx, ranking["ranking"])
    return pipeline.run(ctx, ranking)


def stage_membrane(ctx) -> dict[str, Any]:
    from .membrane import environment

    return environment.run(ctx)


def stage_md(ctx, selection: dict[str, Any]) -> dict[str, Any]:
    from .md import adapters, package

    summary = package.build(ctx, selection)
    summary["adapters"] = adapters.build_adapters(ctx, summary)
    return summary


def stage_report(ctx, state: dict[str, Any]) -> dict[str, Any]:
    from .report import render

    return render.write_all(ctx, state)


# --------------------------------------------------------------------------- #
# runner
# --------------------------------------------------------------------------- #


_CAPABILITY_CACHE: dict[str, str] = {}


def _stage_capabilities(stage: str) -> str:
    """Fingerprint the optional tools a stage depends on.

    Without this, a run that was cached while PoseBusters was unavailable keeps
    reporting ``unsupported`` forever: the cache key only saw configuration and
    code version, not the existence of the tool.
    """
    if stage in _CAPABILITY_CACHE:
        return _CAPABILITY_CACHE[stage]
    fingerprint = ""
    if stage == "pose_qa":
        available = "unknown"
        try:
            if compat.ascii_scratch_available():
                compat.bootstrap()
                import posebusters

                available = getattr(posebusters, "__version__", "present")
        except Exception:  # noqa: BLE001 - absence is a legitimate answer
            available = "absent"
        fingerprint = f"posebusters={available}"
    _CAPABILITY_CACHE[stage] = fingerprint
    return fingerprint


def _cache_key(ctx, stage: str, extra: str = "") -> str:
    return sha256_text(
        "|".join([stage, ctx.config.digest(), __version__, extra, _stage_capabilities(stage)])
    )[:32]


def _input_digest(inputs):
    """Cache by bytes, not a mutable directory name."""
    records = []
    for directory in inputs.get('input_dirs') or []:
        for base, dirs, names in os.walk(directory):
            dirs.sort()
            for name in sorted(names):
                path = os.path.join(base, name)
                if os.path.isfile(path):
                    records.append((os.path.abspath(path), sha256_file(path)))
    manifest = inputs.get('manifest')
    if manifest and os.path.isfile(manifest):
        records.append((os.path.abspath(manifest), sha256_file(manifest)))
    return sha256_text(json.dumps({'files': records, 'inputs': inputs}, sort_keys=True, default=str))


def run_all(
    ctx,
    inputs: dict[str, Any],
    *,
    stages: list[str] | None = None,
    resume: bool = False,
    progress: Callable[[str, dict[str, Any]], None] | None = None,
) -> dict[str, Any]:
    unknown = set(stages or []) - set(STAGE_ORDER)
    if unknown:
        raise ValueError('Unknown stages: ' + ','.join(sorted(unknown)))
    wanted = [s for s in (stages or STAGE_ORDER) if s in STAGE_ORDER]
    input_digest = _input_digest(inputs)
    from pathlib import Path
    code_digest = sha256_text(json.dumps([(str(p.relative_to(Path(__file__).parent)), sha256_file(str(p))) for p in sorted(Path(__file__).parent.rglob('*.py'))]))
    state: dict[str, Any] = {"run_id": ctx.run_id, "stages": {}, "errors": {}}
    ranking: dict[str, Any] | None = read_json(os.path.join(ctx.stage_dir("scoring"), "summary.json"), {}) or None
    selection: dict[str, Any] | None = read_json(os.path.join(ctx.stage_dir("selection"), "summary.json"), {}) or None
    ctx.db.upsert_run(ctx.run_id, status="running")
    ctx.write_manifest({"run_id": ctx.run_id, "inputs": inputs, "stages": wanted, "started": now_iso()})

    for stage in wanted:
        out_dir = ctx.stage_dir(stage)
        summary_path = os.path.join(out_dir, "summary.json")
        key = _cache_key(ctx, stage, ctx.run_id + input_digest + code_digest)
        cached = ctx.db.cache_get(key) if resume else None
        if cached and os.path.isfile(summary_path) and cached.get("output_hash") == sha256_file(summary_path):
            payload = read_json(summary_path, {}) or {}
            payload["cached"] = True
            state["stages"][stage] = payload
            if stage == "scoring":
                ranking = payload
            elif stage == "selection":
                selection = payload
            if progress:
                progress(stage, {"status": "cached"})
            continue

        prior_progress=read_json(os.path.join(ctx.run_dir,"live_progress.json"),{}) or {}
        atomic_write_json(os.path.join(ctx.run_dir,"live_progress.json"),{**prior_progress,"current_stage":stage,"status":"running"})
        ctx.log(f"stage {stage} started", run_id=ctx.run_id, stage=stage)
        try:
            if stage == "target":
                payload = prepare_target(ctx)
            elif stage == "import":
                payload = stage_import(ctx, inputs)
            elif stage == "pose_qa":
                payload = stage_pose_qa(ctx)
            elif stage == "mechanism":
                payload = stage_mechanism(ctx)
            elif stage == "motion":
                payload = stage_motion(ctx)
            elif stage == "redock":
                payload = stage_redock(ctx, selection)
            elif stage == "strain":
                payload = stage_strain(ctx)
            elif stage == "scoring":
                payload = stage_scoring(ctx)
                ranking = payload
            elif stage == "selection":
                if ranking is None:
                    ranking = state["stages"].get("scoring") or {}
                payload = stage_selection(ctx, ranking)
                selection = payload
            elif stage == "membrane":
                payload = stage_membrane(ctx)
            elif stage == "md":
                payload = stage_md(ctx, selection or {"primary": [], "exploratory": []})
            elif stage == "report":
                state["stages"]["report"] = {}
                payload = stage_report(ctx, state)
            else:  # pragma: no cover
                raise ValueError(stage)
            from . import data
            atomic_write_json(os.path.join(ctx.run_dir,"live_progress.json"),{"current_stage":stage,"status":"completed","counts":live_counts(data.pose_rows(ctx),selection,STAGE_ORDER.index(stage)>=STAGE_ORDER.index('pose_qa'))})
            payload["cached"] = False
            state["stages"][stage] = payload
            atomic_write_json(summary_path, payload)
            ctx.db.cache_put(
                key=key,
                stage=stage,
                scope_id=ctx.run_id,
                input_hash=input_digest,
                config_hash=ctx.config.digest(),
                code_version=__version__,
                tool_versions={},
                output_path=summary_path,
                output_hash=sha256_file(summary_path),
                completed=True,
            )
            if progress:
                progress(stage, {"status": "done"})
        except Exception as exc:  # noqa: BLE001
            state["errors"][stage] = {"error": f"{type(exc).__name__}: {exc}", "traceback": traceback.format_exc()[-4000:]}
            ctx.log(f"stage {stage} failed: {exc}", run_id=ctx.run_id, stage=stage, level="error")
            atomic_write_json(os.path.join(out_dir, "error.json"), state["errors"][stage])
            ctx.db.upsert_run(ctx.run_id, status="failed")
            if progress:
                progress(stage, {"status": "error", "error": str(exc)})
            break

    ok = not state["errors"]
    ctx.db.upsert_run(ctx.run_id, status="completed" if ok else "failed")
    atomic_write_json(os.path.join(ctx.run_dir, "run_state.json"), state)
    return state
