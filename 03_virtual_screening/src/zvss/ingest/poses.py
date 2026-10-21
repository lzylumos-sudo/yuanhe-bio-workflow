"""Pose reading, score normalisation and score-to-pose binding.

The single most important rule in this module: a score is attached to a pose by
an *explicit key found in the data*, never by list position, file name order or
lexicographic sorting.  Concretely:

* GNINA SDF records carry their own score properties -> authoritative.
* GNINA stdout logs carry a per-mode table -> parsed and cross-checked against
  the SDF values; a mismatch is reported, not silently resolved.
* Sidecar tables (CSV/TSV) must expose a joining key that also exists in the
  pose records (``pose_file``, ``pose_id``, ``record_index`` or
  ``(parent_id, rank)``).  If no key exists, the table is treated as
  *compound-level only* and is never broadcast onto individual poses.
"""

from __future__ import annotations

import csv
import gzip
import os
import re
from dataclasses import dataclass, field
from typing import Any, Iterable

from rdkit import Chem

from ..chem import rdio
from ..constants import InputError, Reason
from ..ids import coords_digest
from ..util import AMINO_ACIDS, ION_NAMES, LIPID_NAMES, WATER_NAMES, parse_pdb_text, read_text, safe_float, safe_int

# --------------------------------------------------------------------------- #
# score field normalisation
# --------------------------------------------------------------------------- #

# canonical name -> accepted aliases (lower-cased, non-alphanumerics stripped)
SCORE_ALIASES: dict[str, tuple[str, ...]] = {
    "vina_affinity_kcal_mol": ("vinaaffinity", "vinaaffinitykcalmol"),
    "minimizedAffinity": (
        "minimizedaffinity", "minimizedaffinitykcalmol", "minimizedaffinitykcalmole",
        "affinity", "affinitykcalmol", "empiricalscorekcalmol",
        "empirical", "empiricalscore",
    ),
    "CNNscore": ("cnnscore", "cnnposescore", "cnnpose", "cnn"),
    "CNNaffinity": ("cnnaffinity", "cnnaffinitypk", "cnn_affinity", "cnnlogaffinity"),
    "CNN_VS": ("cnnvs", "cnnvsprod", "cnnvirtualscreening"),
    "minimizedRMSD": ("minimizedrmsd", "rmsd", "rmsdlb", "minimizedrmsdangstrom"),
    "intramolecularEnergy": ("intramol", "intramolenergy", "intramolecularenergy"),
    "GridScore": ("gridscore", "grid"),
}

SCORE_SEMANTICS: dict[str, dict[str, Any]] = {
    "vina_affinity_kcal_mol": {"unit":"kcal/mol", "direction":"lower_is_better", "meaning":"Official Vina empirical score; not binding free energy"},
    "minimizedAffinity": {
        "unit": "kcal/mol",
        "direction": "lower_is_better",
        "meaning": "GNINA/Vina empirical (force-field) score for this pose; NOT a binding free energy",
    },
    "CNNscore": {
        "unit": "dimensionless",
        "direction": "higher_is_better",
        "meaning": "CNN pose-quality score in [0,1]; a pose score, NOT an inhibition probability",
    },
    "CNNaffinity": {
        "unit": "pK (log scale)",
        "direction": "higher_is_better",
        "meaning": "CNN predicted affinity as a pK-like value; never convert to kcal/mol",
    },
    "CNN_VS": {
        "unit": "dimensionless",
        "direction": "higher_is_better",
        "meaning": "CNN_VS = CNNscore * CNNaffinity as reported by GNINA; a ranking aid, not evidence",
    },
    "minimizedRMSD": {
        "unit": "angstrom",
        "direction": "informational",
        "meaning": "deviation of the minimised pose from the input conformer",
    },
    "intramolecularEnergy": {
        "unit": "kcal/mol",
        "direction": "lower_is_better",
        "meaning": "internal strain inside the docking scoring function",
    },
    "GridScore": {
        "unit": "arbitrary",
        "direction": "lower_is_better",
        "meaning": "generic grid score",
    },
}


def _norm_key(text: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(text).lower())


ALIAS_LOOKUP: dict[str, str] = {}
for canonical, aliases in SCORE_ALIASES.items():
    ALIAS_LOOKUP[_norm_key(canonical)] = canonical
    for alias in aliases:
        ALIAS_LOOKUP[_norm_key(alias)] = canonical


def normalise_scores(raw: dict[str, Any]) -> tuple[dict[str, float | None], dict[str, str]]:
    """Map recognised score aliases to canonical names.

    Returns ``(scores, provenance)``; ``provenance`` records which source field
    produced each canonical value so the report can show it.
    """
    scores: dict[str, float | None] = {key: None for key in SCORE_ALIASES}
    provenance: dict[str, str] = {}
    for key, value in raw.items():
        canonical = ALIAS_LOOKUP.get(_norm_key(key))
        if not canonical:
            continue
        number = safe_float(value)
        if number is None:
            continue
        if scores.get(canonical) is None:
            scores[canonical] = number
            provenance[canonical] = str(key)
    return scores, provenance


# --------------------------------------------------------------------------- #
# GNINA log parsing
# --------------------------------------------------------------------------- #


@dataclass
class GninaModeRow:
    mode: int
    affinity: float | None
    intramol: float | None
    cnn_pose_score: float | None
    cnn_affinity: float | None
    extra: list[float] = field(default_factory=list)


@dataclass
class GninaLog:
    path: str
    version: str = ""
    build: str = ""
    argv: list[str] = field(default_factory=list)
    seed: int | None = None
    model: str = ""
    scoring: str = ""
    modes: list[GninaModeRow] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "version": self.version,
            "build": self.build,
            "argv": self.argv,
            "seed": self.seed,
            "model": self.model,
            "scoring": self.scoring,
            "modes": [row.__dict__ for row in self.modes],
            "warnings": self.warnings,
        }


MODE_ROW = re.compile(
    r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$"
)


def parse_gnina_log(text: str, path: str = "<memory>") -> GninaLog:
    log = GninaLog(path=path)
    in_table = False
    for line in text.splitlines():
        if not log.version:
            match = re.search(r"gnina v([0-9][^\s]*)", line)
            if match:
                log.version = match.group(1)
                build = re.search(r"Built ([^\.]+)\.", line)
                if build:
                    log.build = build.group(1).strip()
        if "Commandline:" in line:
            import shlex

            try:
                log.argv = shlex.split(line.split("Commandline:", 1)[1].strip())
            except ValueError:
                log.argv = line.split("Commandline:", 1)[1].split()
            for index, token in enumerate(log.argv):
                if token in {"--cnn", "-c"} and index + 1 < len(log.argv):
                    log.model = log.argv[index + 1]
                if token == "--scoring" and index + 1 < len(log.argv):
                    log.scoring = log.argv[index + 1]
                if token == "--seed" and index + 1 < len(log.argv):
                    log.seed = safe_int(log.argv[index + 1])
        if "Using random seed:" in line:
            log.seed = safe_int(line.split(":", 1)[1].strip())
        if re.search(r"\bmode\b.*affinity", line, re.IGNORECASE):
            in_table = True
            continue
        if in_table:
            if set(line.strip()) <= {"-", "+", "|", " "}:
                continue
            match = MODE_ROW.match(line)
            if match:
                values = [safe_float(match.group(i)) for i in range(2, 6)]
                log.modes.append(
                    GninaModeRow(
                        mode=int(match.group(1)),
                        affinity=values[0],
                        intramol=values[1],
                        cnn_pose_score=values[2],
                        cnn_affinity=values[3],
                    )
                )
            elif log.modes:
                in_table = False
        if "initial pose not within box" in line or "not within box" in line:
            log.warnings.append(line.strip())
        if "WARNING" in line or "ERROR" in line:
            log.warnings.append(line.strip())
    log.notes.append(f"parsed {len(log.modes)} mode rows")
    return log


def read_gnina_log(path: str) -> GninaLog:
    return parse_gnina_log(read_text(path), path)


def match_log_modes_to_records(
    log: GninaLog, records: list[dict[str, Any]], *, tolerance: float = 5e-3
) -> dict[str, Any]:
    """Pair log mode rows with SDF records by score values.

    GNINA writes modes in descending score order and the SDF follows that order,
    but relying on order is exactly the failure mode the spec forbids.  Pairing
    on the numeric tuple makes the association verifiable; unpaired rows are
    reported.
    """
    unmatched_modes = list(range(len(log.modes)))
    pairs: list[dict[str, Any]] = []
    for record_index, record in enumerate(records):
        scores = record.get("scores") or {}
        affinity = scores.get("minimizedAffinity")
        cnn = scores.get("CNNscore")
        best_index = None
        best_delta = None
        for mode_index in unmatched_modes:
            row = log.modes[mode_index]
            if affinity is None or row.affinity is None or cnn is None or row.cnn_pose_score is None:
                continue
            delta = abs(row.affinity - affinity) + abs(row.cnn_pose_score - cnn)
            if delta <= tolerance * 4 and (best_delta is None or delta < best_delta):
                best_delta = delta
                best_index = mode_index
        if best_index is not None:
            unmatched_modes.remove(best_index)
            row = log.modes[best_index]
            pairs.append(
                {
                    "record_index": record_index,
                    "mode": row.mode,
                    "value_delta": best_delta,
                    "seed": log.seed,
                }
            )
    return {
        "paired": pairs,
        "unpaired_modes": [log.modes[i].__dict__ for i in unmatched_modes],
        "unpaired_records": [
            index for index in range(len(records)) if index not in {p["record_index"] for p in pairs}
        ],
        "log_modes": len(log.modes),
        "records": len(records),
    }


# --------------------------------------------------------------------------- #
# tables
# --------------------------------------------------------------------------- #


@dataclass
class Table:
    path: str
    columns: list[str]
    rows: list[dict[str, Any]]
    delimiter: str = ","

    def canonical_columns(self) -> dict[str, str]:
        out: dict[str, str] = {}
        for column in self.columns:
            canonical = ALIAS_LOOKUP.get(_norm_key(column))
            if canonical:
                out[column] = canonical
        return out

    def to_dict(self) -> dict[str, Any]:
        return {"path": self.path, "columns": self.columns, "rows": len(self.rows)}


def read_table(path: str, *, max_rows: int = 2_000_000) -> Table:
    opener = gzip.open if path.lower().endswith(".gz") else open
    delimiter = "\t" if path.lower().endswith((".tsv", ".tsv.gz")) else ","
    with opener(path, "rt", encoding="utf-8-sig", errors="replace", newline="") as handle:
        sample = handle.read(8192)
        handle.seek(0)
        if "\t" in sample and "," not in sample.splitlines()[0]:
            delimiter = "\t"
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows: list[dict[str, Any]] = []
        for index, row in enumerate(reader):
            if index >= max_rows:
                raise InputError(f"table {path} exceeds {max_rows} rows")
            rows.append(row)
        return Table(path=path, columns=list(reader.fieldnames or []), rows=rows, delimiter=delimiter)


JOIN_KEYS = ("pose_file", "pose_id", "record_index", "model_index")


def bind_table_to_poses(table: Table, poses: list["RawPose"], *, parent_key: str = "parent_id") -> dict[str, Any]:
    """Attach table rows to poses through an explicit key.

    Returns a report; unmapped rows/poses are listed, and when no usable key
    exists the caller must keep the scores at compound level only.
    """
    lower_columns = {column.lower(): column for column in table.columns}
    usable_key = None
    for key in JOIN_KEYS:
        if key in lower_columns:
            usable_key = lower_columns[key]
            break

    report: dict[str, Any] = {
        "table": table.path,
        "rows": len(table.rows),
        "poses": len(poses),
        "join_key": usable_key,
        "mapped": 0,
        "unmapped_rows": [],
        "unmapped_poses": [],
        "level": "pose" if usable_key else "table_only",
        "reason": None,
    }
    if not usable_key:
        report["reason"] = str(Reason.SCORE_MAPPING_UNVERIFIED if hasattr(Reason, "SCORE_MAPPING_UNVERIFIED") else Reason.NEEDS_INPUT)
        report["note"] = (
            "table exposes no pose_file/pose_id/record_index column; scores stay at compound level "
            "and are NOT broadcast onto individual poses"
        )
        return report

    pose_lookup: dict[str, RawPose] = {}
    for pose in poses:
        if usable_key == "pose_file" or usable_key.endswith("file"):
            pose_lookup[os.path.basename(pose.source_file)] = pose
        elif usable_key == "pose_id":
            pose_lookup[str(pose.extras.get("pose_id", ""))] = pose
        elif usable_key == "record_index":
            pose_lookup[str(pose.record_index)] = pose
        elif usable_key == "model_index":
            pose_lookup[str(pose.model_index)] = pose

    for row_index, row in enumerate(table.rows):
        raw = row.get(usable_key, "")
        key = os.path.basename(str(raw).replace("\\", "/")) if usable_key.endswith("file") else str(raw)
        pose = pose_lookup.get(key)
        if pose is None:
            report["unmapped_rows"].append(row_index)
            continue
        scores, provenance = normalise_scores(row)
        for name, value in scores.items():
            if value is not None and pose.scores.get(name) is None:
                pose.scores[name] = value
                pose.score_provenance.setdefault(name, f"{os.path.basename(table.path)}:{usable_key}={key}")
        if row.get(parent_key):
            pose.extras.setdefault("table_parent_id", row[parent_key])
        pose.extras.setdefault("table_row_index", row_index)
        report["mapped"] += 1

    mapped_basenames = {
        os.path.basename(p.source_file) for p in poses if p.extras.get("table_row_index") is not None
    }
    report["unmapped_poses"] = [
        os.path.basename(p.source_file) for p in poses if os.path.basename(p.source_file) not in mapped_basenames
    ]
    return report


# --------------------------------------------------------------------------- #
# poses
# --------------------------------------------------------------------------- #


@dataclass
class RawPose:
    source_file: str
    source_hash: str
    record_index: int
    model_index: int | None = None
    mol: Chem.Mol | None = None
    props: dict[str, Any] = field(default_factory=dict)
    scores: dict[str, float | None] = field(default_factory=dict)
    score_provenance: dict[str, str] = field(default_factory=dict)
    parent_label: str = ""
    ligand_block: str = ""
    atom_names: list[str] = field(default_factory=list)
    resname: str = "LIG"
    chain: str = "L"
    resid: str = "1"
    needs_template: bool = False
    issues: list[str] = field(default_factory=list)
    extras: dict[str, Any] = field(default_factory=dict)

    @property
    def coords_digest(self) -> str | None:
        if self.mol is None or self.mol.GetNumConformers() == 0:
            return None
        return coords_digest(rdio.coordinates(self.mol))

    def to_dict(self) -> dict[str, Any]:
        return {
            "source_file": self.source_file,
            "record_index": self.record_index,
            "model_index": self.model_index,
            "parent_label": self.parent_label,
            "scores": self.scores,
            "score_provenance": self.score_provenance,
            "needs_template": self.needs_template,
            "issues": self.issues,
            "atoms": self.mol.GetNumAtoms() if self.mol else 0,
        }


def _mol_name(mol: Chem.Mol) -> str:
    try:
        return mol.GetProp("_Name").strip()
    except KeyError:
        return ""


def read_sdf_poses(discovered, *, table: Table | None = None, title: str | None = None) -> tuple[list[RawPose], dict[str, Any]]:
    """Read every record of an SD file as a pose candidate."""
    poses: list[RawPose] = []
    report: dict[str, Any] = {"file": discovered.path, "records": 0, "parse_failures": [], "score_bound": 0}
    record_index = 0
    for mol, props in rdio.iter_sdf(discovered.path):
        if mol is None:
            report["parse_failures"].append(record_index)
            poses.append(
                RawPose(
                    source_file=discovered.path,
                    source_hash=discovered.sha256,
                    record_index=record_index,
                    issues=["record could not be parsed by RDKit"],
                    extras={"parse_error": True},
                )
            )
            record_index += 1
            continue
        clean_props = {k: v for k, v in props.items() if not k.startswith("_")}
        scores, provenance = normalise_scores(clean_props)
        pose = RawPose(
            source_file=discovered.path,
            source_hash=discovered.sha256,
            record_index=record_index,
            model_index=None,
            mol=mol,
            props=clean_props,
            scores=scores,
            score_provenance=provenance,
            parent_label=(title or _mol_name(mol)),
            extras={"record_index": record_index},
        )
        if not any(value is not None for value in scores.values()):
            pose.issues.append("no score properties on this SDF record")
        poses.append(pose)
        record_index += 1
    report["records"] = record_index
    report["score_bound"] = sum(
        1 for pose in poses if any(v is not None for v in pose.scores.values())
    )
    if table is not None:
        report["table_binding"] = bind_table_to_poses(table, poses)
    return poses, report


def read_pdb_poses(discovered, *, selector: dict[str, str] | None = None) -> tuple[list[RawPose], dict[str, Any]]:
    """Read ligand poses out of a ligand-only or complex PDB, per MODEL."""
    text = read_text(discovered.path)
    if discovered.path.lower().endswith(".gz"):
        with gzip.open(discovered.path, "rt", encoding="utf-8", errors="replace") as handle:
            text = handle.read()
    structure = parse_pdb_text(text, keep_hydrogens=True)
    report: dict[str, Any] = {"file": discovered.path, "models": structure.models, "candidates": []}

    poses: list[RawPose] = []
    models = structure.models or [1]
    for model in models:
        extraction = extract_ligand_candidates(structure.by_model(model), selector=selector)
        report["candidates"].append(
            {"model": model, "selected": extraction.get("selected"), "candidates": extraction.get("candidates")}
        )
        selected = extraction.get("selected")
        if not selected:
            poses.append(
                RawPose(
                    source_file=discovered.path,
                    source_hash=discovered.sha256,
                    record_index=len(poses),
                    model_index=model,
                    issues=[extraction.get("reason", "no unique ligand could be identified")],
                    extras={"candidates": extraction.get("candidates")},
                )
            )
            continue
        block = extraction["block"]
        names = extraction["atom_names"]
        pose = RawPose(
            source_file=discovered.path,
            source_hash=discovered.sha256,
            record_index=len(poses),
            model_index=model,
            props={},
            parent_label=selected.get("resname", ""),
            ligand_block=block,
            atom_names=names,
            resname=selected["resname"],
            chain=selected["chain"],
            resid=str(selected["resseq"]),
            needs_template=True,
            extras={"record_index": len(poses), "model_index": model},
        )
        # A ligand-only PDB can often be parsed directly; a complex PDB cannot.
        direct = rdio.mol_from_pdb_block(block, sanitize=False, remove_hs=False)
        if direct is not None and direct.GetNumAtoms() > 0:
            pose.mol = direct
        poses.append(pose)
    report["poses"] = len(poses)
    return poses, report


def extract_ligand_candidates(structure, *, selector: dict[str, str] | None = None) -> dict[str, Any]:
    """Identify the intended small molecule among HETATM groups.

    Waters, ions, lipids and additives are excluded by name; remaining groups are
    ranked and reported.  When the choice is not unique the caller must ask for
    an explicit selector instead of taking "the first HETATM".
    """
    groups: dict[tuple[str, int, str, str], list] = {}
    for atom in structure.atoms:
        if atom.record != "HETATM":
            continue
        if atom.resname.upper() in WATER_NAMES:
            continue
        groups.setdefault((atom.chain, atom.resseq, atom.icode, atom.resname), []).append(atom)

    candidates = []
    for (chain, resseq, icode, resname), atoms in groups.items():
        upper = resname.upper()
        kind = "ligand"
        if upper in ION_NAMES:
            kind = "ion"
        elif upper in LIPID_NAMES:
            kind = "lipid"
        elif upper in AMINO_ACIDS:
            kind = "modified_residue"
        heavy = [atom for atom in atoms if not atom.is_hydrogen()]
        candidates.append(
            {
                "resname": resname,
                "chain": chain,
                "resseq": resseq,
                "icode": icode,
                "atoms": len(atoms),
                "heavy_atoms": len(heavy),
                "elements": sorted({(atom.element or atom.name[:1]).upper() for atom in heavy}),
                "kind": kind,
                "atoms_obj": atoms,
            }
        )

    ligand_like = [c for c in candidates if c["kind"] == "ligand" and c["heavy_atoms"] > 1]
    if selector:
        wanted = {
            "resname": (selector.get("resname") or "").upper(),
            "chain": selector.get("chain"),
            "resseq": safe_int(selector.get("resseq") or selector.get("resid")),
        }
        picked = [
            c
            for c in candidates
            if (not wanted["resname"] or c["resname"].upper() == wanted["resname"])
            and (wanted["chain"] in (None, "") or c["chain"] == wanted["chain"])
            and (wanted["resseq"] is None or c["resseq"] == wanted["resseq"])
        ]
        if len(picked) == 1:
            return {
                "selected": _public(picked[0]),
                "candidates": [_public(c) for c in candidates],
                "block": _ligand_block_from_atoms(picked[0]["atoms_obj"]),
                "atom_names": [atom.name for atom in picked[0]["atoms_obj"] if not atom.is_hydrogen()],
            }
        return {
            "selected": None,
            "candidates": [_public(c) for c in candidates],
            "reason": f"selector matched {len(picked)} groups; needs_input" if picked else "selector matched nothing",
        }

    if len(ligand_like) == 1:
        chosen = ligand_like[0]
        return {
            "selected": _public(chosen),
            "candidates": [_public(c) for c in candidates],
            "block": _ligand_block_from_atoms(chosen["atoms_obj"]),
            "atom_names": [atom.name for atom in chosen["atoms_obj"] if not atom.is_hydrogen()],
        }
    reason = (
        "no ligand-like HETATM group found"
        if not ligand_like
        else f"{len(ligand_like)} ligand-like groups; provide an explicit selector"
    )
    return {"selected": None, "candidates": [_public(c) for c in candidates], "reason": reason}


def _public(candidate: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in candidate.items() if key != "atoms_obj"}


def _ligand_block_from_atoms(atoms: list) -> str:
    lines = ["HEADER    ZVSS LIGAND EXTRACT", "MODEL        1"]
    for atom in atoms:
        lines.append(atom.to_line())
    lines.append("ENDMDL")
    lines.append("END")
    return "\n".join(lines) + "\n"


def collect_poses(discovered_files: Iterable, *, table: Table | None = None, selector: dict | None = None,
                  titles: dict[str, str] | None = None) -> dict[str, Any]:
    """Turn discovered pose files into pose candidates plus a provenance report.

    A directory that exposes both ``poses.sdf`` (all modes) and per-mode
    ``pose_00N.sdf`` files is read once, from ``poses.sdf``; otherwise the same
    physical pose would be imported twice.
    """
    files = list(discovered_files)
    directories_with_combined = {
        os.path.dirname(item.path)
        for item in files
        if os.path.basename(item.path).lower().startswith(("poses.sdf", "poses.sdf.gz"))
    }
    poses: list[RawPose] = []
    reports: list[dict[str, Any]] = []
    skipped: list[str] = []
    for item in files:
        base = os.path.basename(item.path).lower()
        if os.path.dirname(item.path) in directories_with_combined and not base.startswith("poses.sdf"):
            skipped.append(item.path)
            continue
        if item.kind == "pose_sdf":
            file_poses, report = read_sdf_poses(item, table=table if len(poses) == 0 else None,
                                                 title=(titles or {}).get(os.path.abspath(item.path)) or (titles or {}).get(os.path.basename(item.path)))
            poses.extend(file_poses)
            reports.append(report)
        elif item.kind == "pose_pdb":
            file_poses, report = read_pdb_poses(item, selector=selector)
            poses.extend(file_poses)
            reports.append(report)
    for index, pose in enumerate(poses):
        pose.extras.setdefault("global_index", index)
    return {"poses": poses, "reports": reports, "skipped_duplicates": skipped}
