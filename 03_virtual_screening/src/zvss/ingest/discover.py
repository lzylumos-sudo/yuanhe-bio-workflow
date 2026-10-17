"""Input discovery and manifest handling.

Requirements this implements:
* a directory, a list of directories, explicit files, or an explicit manifest;
* arbitrary file names (Chinese, spaces, parentheses, non-consecutive numbers);
* GNINA native SDF/SDF.gz with several molecules and several poses per molecule;
* small-molecule PDB and protein-ligand complex PDB, including multi-MODEL;
* score CSV/TSV and GNINA logs as *sidecars*, never as a replacement for
  per-pose identity.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any, Iterable

import yaml

from ..constants import InputError
from ..ids import sha256_file
from ..util import read_text, safe_float

SDF_SUFFIXES = (".sdf", ".sdf.gz", ".sd", ".mol2")
PDB_SUFFIXES = (".pdb", ".ent", ".pdb.gz")
TABLE_SUFFIXES = (".csv", ".tsv")
LOG_SUFFIXES = (".log", ".txt", ".out")
MANIFEST_SUFFIXES = (".yaml", ".yml", ".json")

GNINA_LOG_HINTS = ("gnina", "Commandline", "CNN pose score", "minimizedAffinity")

#: SD files that hold a *pre-docking* structure rather than a pose.
INPUT_SDF_STEMS = {
    "input", "ligand", "ligand_template", "template", "ligands", "reference", "ref_ligand", "crystal_ligand",
}


@dataclass
class DiscoveredFile:
    path: str
    kind: str          # pose_sdf | pose_pdb | table | log | receptor | template | unknown
    size: int
    sha256: str
    note: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "kind": self.kind,
            "size": self.size,
            "sha256": self.sha256,
            "note": self.note,
        }


@dataclass
class InputBundle:
    root: str
    files: list[DiscoveredFile] = field(default_factory=list)
    receptor_file: str | None = None
    box_file: str | None = None
    manifest_path: str | None = None
    hints: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def by_kind(self, kind: str) -> list[DiscoveredFile]:
        return [item for item in self.files if item.kind == kind]

    def to_dict(self) -> dict[str, Any]:
        return {
            "root": self.root,
            "receptor_file": self.receptor_file,
            "box_file": self.box_file,
            "manifest_path": self.manifest_path,
            "hints": self.hints,
            "warnings": self.warnings,
            "files": [item.to_dict() for item in self.files],
        }


def _classify(path: str) -> str:
    name = os.path.basename(path).lower()
    stem = os.path.splitext(name)[0].replace(".gz", "")
    if name.startswith("._"):
        return "metadata"
    if name.endswith(SDF_SUFFIXES):
        # Pre-docking inputs and templates are NOT poses: counting them would
        # multiply every hit and break score-to-pose binding.
        if stem in INPUT_SDF_STEMS or stem.startswith(("ligand", "template", "ref_")):
            return "template"
        return "pose_sdf"
    if name.endswith(PDB_SUFFIXES):
        if "receptor" in name or "protein" in name:
            return "receptor"
        if stem.startswith("complex"):
            return "complex"
        return "pose_pdb"
    if name.endswith(TABLE_SUFFIXES):
        return "table"
    if name.endswith(LOG_SUFFIXES):
        return "log"
    return "unknown"


def _looks_like_gnina_log(path: str) -> bool:
    try:
        head = read_text(path)[:4000]
    except OSError:
        return False
    hits = sum(1 for hint in GNINA_LOG_HINTS if hint in head)
    return hits >= 2


def iter_candidate_files(root: str, *, recursive: bool = True, max_files: int = 20000) -> Iterable[str]:
    if os.path.isfile(root):
        yield root
        return
    if not os.path.isdir(root):
        raise InputError(f"input root does not exist: {root}")
    count = 0
    if recursive:
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
            for filename in sorted(filenames):
                if filename.startswith("._") or filename.startswith("."):
                    continue
                if filename in {"Thumbs.db", "desktop.ini"}:
                    continue
                count += 1
                if count > max_files:
                    raise InputError(
                        f"more than {max_files} candidate files under {root}; pass an explicit manifest"
                    )
                yield os.path.join(dirpath, filename)
    else:
        for filename in sorted(os.listdir(root)):
            full = os.path.join(root, filename)
            if os.path.isfile(full) and not filename.startswith("."):
                yield full


def load_manifest(path: str) -> dict[str, Any]:
    text = read_text(path)
    if path.lower().endswith((".yaml", ".yml")):
        payload = yaml.safe_load(text) or {}
    else:
        import json

        payload = json.loads(text)
    if not isinstance(payload, dict):
        raise InputError(f"manifest {path} must be a mapping")
    return payload


def discover(
    roots: list[str],
    *,
    manifest: str | None = None,
    recursive: bool = True,
    hint_files: dict[str, str] | None = None,
) -> InputBundle:
    """Collect input files and classify them without reading chemistry yet."""
    if not roots and not manifest:
        raise InputError("give at least one input root or a manifest")
    bundle = InputBundle(root=os.path.abspath(roots[0]) if roots else os.path.dirname(os.path.abspath(manifest or ".")))
    seen: set[str] = set()

    hints: dict[str, Any] = {}
    if manifest:
        hints = load_manifest(manifest)
        bundle.manifest_path = os.path.abspath(manifest)
        manifest_dir = os.path.dirname(os.path.abspath(manifest))
        for entry in hints.get("inputs", []) or []:
            candidate = entry if isinstance(entry, str) else entry.get("path")
            if not candidate:
                continue
            candidate = os.path.expandvars(os.path.expanduser(str(candidate)))
            if not os.path.isabs(candidate):
                candidate = os.path.join(manifest_dir, candidate)
            if not os.path.exists(candidate):
                bundle.warnings.append(f"manifest input missing: {candidate}")
                continue
            for path in iter_candidate_files(candidate, recursive=recursive):
                _add_file(bundle, path, seen)
        receptor = hints.get("receptor_file") or hints.get("receptor")
        if receptor:
            receptor = str(receptor)
            if not os.path.isabs(receptor):
                receptor = os.path.join(manifest_dir, receptor)
            if os.path.exists(receptor):
                bundle.receptor_file = os.path.normpath(receptor)
            else:
                bundle.warnings.append(f"manifest receptor missing: {receptor}")
        box = hints.get("box_file") or hints.get("box")
        if isinstance(box, str):
            if not os.path.isabs(box):
                box = os.path.join(manifest_dir, box)
            if os.path.exists(box):
                bundle.box_file = os.path.normpath(box)
            else:
                bundle.warnings.append(f"manifest box missing: {box}")

    for root in roots:
        for path in iter_candidate_files(root, recursive=recursive):
            _add_file(bundle, path, seen)

    for logical, path in (hint_files or {}).items():
        if path and os.path.exists(path):
            bundle.hints[logical] = os.path.abspath(path)

    # Promote a log that is actually a GNINA log; drop macOS resource forks.
    for item in bundle.files:
        if item.kind == "log" and _looks_like_gnina_log(item.path):
            item.note = "gnina_log"
        if item.kind == "metadata":
            bundle.warnings.append(f"ignored macOS metadata file: {item.path}")

    if bundle.receptor_file is None:
        receptors = bundle.by_kind("receptor")
        if len(receptors) == 1:
            bundle.receptor_file = receptors[0].path
        elif len(receptors) > 1:
            bundle.warnings.append(
                "several files look like receptors; pass an explicit receptor_file "
                f"(candidates: {[os.path.basename(r.path) for r in receptors][:6]})"
            )

    bundle.hints.update({k: v for k, v in hints.items() if k not in {"inputs"}})
    return bundle


def _add_file(bundle: InputBundle, path: str, seen: set[str]) -> None:
    real = os.path.abspath(path)
    if real in seen:
        return
    seen.add(real)
    kind = _classify(real)
    if kind == "metadata":
        return
    try:
        size = os.path.getsize(real)
        digest = sha256_file(real)
    except OSError as exc:
        bundle.warnings.append(f"cannot read {real}: {exc}")
        return
    bundle.files.append(DiscoveredFile(path=real, kind=kind, size=size, sha256=digest))


def detect_box(source: dict[str, Any] | None) -> dict[str, Any] | None:
    """Normalise a GNINA/Vina box record from config or a JSON sidecar."""
    if not source:
        return None
    if "center" in source and "size" in source:
        return {
            "center": [float(x) for x in source["center"]],
            "size": [float(x) for x in source["size"]],
            "units": source.get("units", "angstrom"),
            "method": source.get("method", "declared"),
            "source": source.get("source", ""),
        }
    if all(key in source for key in ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")):
        return {
            "center": [float(source["center_x"]), float(source["center_y"]), float(source["center_z"])],
            "size": [float(source["size_x"]), float(source["size_y"]), float(source["size_z"])],
            "units": source.get("units", "angstrom"),
            "method": source.get("method", "declared"),
            "source": source.get("source", ""),
        }
    return None


def box_from_log(text: str) -> dict[str, Any] | None:
    """Extract center/size from a GNINA command line, when present."""
    import shlex

    for line in text.splitlines():
        if "Commandline:" not in line:
            continue
        try:
            tokens = shlex.split(line.split("Commandline:", 1)[1].strip())
        except ValueError:
            continue
        values: dict[str, float] = {}
        for index, token in enumerate(tokens):
            key = token.lstrip("-")
            if key in {
                "center_x", "center_y", "center_z", "size_x", "size_y", "size_z",
            } and index + 1 < len(tokens):
                value = safe_float(tokens[index + 1])
                if value is not None:
                    values[key] = value
        if len(values) == 6:
            return {
                "center": [values["center_x"], values["center_y"], values["center_z"]],
                "size": [values["size_x"], values["size_y"], values["size_z"]],
                "units": "angstrom",
                "method": "parsed_from_gnina_commandline",
                "source": "gnina_log",
            }
    return None
