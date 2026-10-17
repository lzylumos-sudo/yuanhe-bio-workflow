"""Residue mapping between a construct, the natural sequence and Uniprot.

Primary source: the PDBe *updated* mmCIF, which carries ``_pdbx_sifts_xref_db``.
Fallback: ``_struct_ref_seq`` + ``_pdbx_poly_seq_scheme`` with an explicit
sequence verification, then a global sequence alignment against the natural
sequence.  Every path is recorded in the output, and an empty/failed SIFTS
response is *never* replaced by a guessed constant offset.
"""

from __future__ import annotations

import os
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from typing import Any, Iterable

from ..ids import sha256_file
from ..util import atomic_write_json, atomic_write_text, ensure_dir, read_text

PDBE_UPDATED = "https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}_updated.cif"
PDBE_ORIGINAL = "https://files.rcsb.org/download/{pdb_id}.cif"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M",
    "SEC": "U", "PYL": "O",
}


# --------------------------------------------------------------------------- #
# fetching
# --------------------------------------------------------------------------- #


def fetch(url: str, target: str, *, retries: int = 3, timeout: int = 60) -> dict[str, Any]:
    """Download with retries; returns a provenance record, never raises on 404."""
    ensure_dir(os.path.dirname(os.path.abspath(target)))
    last_error = ""
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "zvss/0.9"})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                status = getattr(response, "status", 200)
                data = response.read()
            if status != 200 or len(data) < 64:
                last_error = f"HTTP {status} with {len(data)} bytes"
                continue
            atomic_write_text(target, data.decode("utf-8", errors="replace"))
            return {
                "url": url,
                "path": target,
                "status": "ok",
                "bytes": len(data),
                "sha256": sha256_file(target),
                "attempts": attempt,
            }
        except urllib.error.HTTPError as exc:
            last_error = f"HTTPError {exc.code}"
            if exc.code in (403, 404):
                break
        except Exception as exc:  # network, TLS, timeout
            last_error = f"{type(exc).__name__}: {exc}"
    return {"url": url, "path": None, "status": "failed", "error": last_error, "attempts": retries}


def fetch_structure(pdb_id: str, cache_dir: str) -> dict[str, Any]:
    """Prefer the SIFTS-annotated file, fall back to the deposited mmCIF."""
    pdb_id = pdb_id.lower()
    updated = os.path.join(cache_dir, f"{pdb_id}_updated.cif")
    record = fetch(PDBE_UPDATED.format(pdb_id=pdb_id), updated)
    if record["status"] == "ok" and _has_sifts(updated):
        record["kind"] = "pdbe_updated_with_sifts"
        return record
    original = os.path.join(cache_dir, f"{pdb_id}.cif")
    fallback = fetch(PDBE_ORIGINAL.format(pdb_id=pdb_id), original)
    if fallback["status"] == "ok":
        fallback["kind"] = "rcsb_original_without_sifts"
        fallback["note"] = (
            "PDBe updated mmCIF unavailable or without _pdbx_sifts_xref_db; "
            "mapping falls back to _struct_ref_seq plus sequence verification"
        )
        fallback["sifts_attempt"] = record
        return fallback
    original_pdb = os.path.join(cache_dir, f"{pdb_id}.pdb")
    third = fetch(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb", original_pdb)
    if third["status"] == "ok":
        third["kind"] = "rcsb_pdb_no_mapping_tables"
        third["note"] = "only coordinates available; mapping must come from sequence alignment"
        return third
    return {"status": "failed", "sifts_attempt": record, "rcsb_attempt": fallback, "pdb_attempt": third}


def fetch_uniprot_sequence(accession: str, cache_dir: str) -> dict[str, Any]:
    target = os.path.join(cache_dir, f"{accession}.fasta")
    if os.path.exists(target) and os.path.getsize(target) > 32:
        return {"status": "cached", "path": target, "sequence": _parse_fasta(read_text(target)),
                "sha256": sha256_file(target)}
    record = fetch(UNIPROT_FASTA.format(acc=accession), target)
    if record["status"] != "ok":
        return {**record, "sequence": ""}
    return {**record, "sequence": _parse_fasta(read_text(target))}


def _parse_fasta(text: str) -> str:
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))


def _has_sifts(path: str) -> bool:
    try:
        return "_pdbx_sifts_xref_db." in read_text(path)[:2_000_000]
    except OSError:
        return False


# --------------------------------------------------------------------------- #
# CIF tables
# --------------------------------------------------------------------------- #

try:  # gemmi is a hard dependency, but keep import errors explicit
    import gemmi
except Exception as _exc:  # pragma: no cover
    gemmi = None
    _GEMMI_ERROR = _exc


def _block(path: str):
    if gemmi is None:  # pragma: no cover
        raise RuntimeError(f"gemmi is required for mmCIF parsing: {_GEMMI_ERROR}")
    document = gemmi.cif.read(path)
    return document.sole_block()


def _table(block, category: str) -> list[dict[str, str]]:
    try:
        table = block.find_mmcif_category(category)
    except Exception:
        return []
    if table is None:
        return []
    tags = [str(tag) for tag in table.tags]
    short = [tag.split(".", 1)[-1] for tag in tags]
    rows: list[dict[str, str]] = []
    for row in table:
        rows.append({short[i]: (row[i] or "") for i in range(len(short))})
    return rows


# --------------------------------------------------------------------------- #
# mapping
# --------------------------------------------------------------------------- #


@dataclass
class ResidueMapRow:
    pdb_id: str
    chain: str
    label_asym_id: str
    label_seq_id: int | None
    auth_seq_id: int | None
    icode: str
    resname: str
    uniprot_resseq: int | None
    uniprot_resname: str
    construct_match: bool
    source: str

    def to_dict(self) -> dict[str, Any]:
        return self.__dict__.copy()


@dataclass
class MappingReport:
    pdb_id: str
    chain: str
    accession: str
    native_length: int
    source: str
    rows: list[ResidueMapRow] = field(default_factory=list)
    mutations: list[dict[str, Any]] = field(default_factory=list)
    mismatches: list[dict[str, Any]] = field(default_factory=list)
    unmapped_construct: list[dict[str, Any]] = field(default_factory=list)
    missing_native: list[int] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "pdb_id": self.pdb_id,
            "chain": self.chain,
            "accession": self.accession,
            "native_length": self.native_length,
            "source": self.source,
            "rows": [row.to_dict() for row in self.rows],
            "mutations": self.mutations,
            "mismatches": self.mismatches,
            "unmapped_construct": self.unmapped_construct,
            "missing_native": self.missing_native,
            "warnings": self.warnings,
            "provenance": self.provenance,
        }

    def summary(self) -> dict[str, Any]:
        mapped = [row for row in self.rows if row.uniprot_resseq is not None]
        return {
            "pdb_id": self.pdb_id,
            "chain": self.chain,
            "source": self.source,
            "mapped_residues": len(mapped),
            "construct_residues": len(self.rows),
            "mutations": len(self.mutations),
            "mismatches": len(self.mismatches),
            "missing_native": len(self.missing_native),
            "native_range": [min((r.uniprot_resseq for r in mapped), default=None),
                             max((r.uniprot_resseq for r in mapped), default=None)],
        }


def build_mapping(
    cif_path: str,
    *,
    pdb_id: str,
    chain: str,
    accession: str,
    native_sequence: str,
) -> MappingReport:
    block = _block(cif_path)
    scheme = _table(block, "_pdbx_poly_seq_scheme.")
    xref = _table(block, "_pdbx_sifts_xref_db.")
    ref_seq = _table(block, "_struct_ref_seq.")
    ref = _table(block, "_struct_ref.")
    diffs = _table(block, "_struct_ref_seq_dif.")

    # construct residues in label order for the requested chain
    construct: list[dict[str, Any]] = []
    for row in scheme:
        asym = row.get("asym_id") or row.get("pdb_strand_id") or ""
        if asym != chain:
            continue
        construct.append(
            {
                "label_asym_id": asym,
                "label_seq_id": _int(row.get("seq_id")),
                "auth_seq_id": _int(row.get("pdb_seq_num")) if row.get("pdb_seq_num") else _int(row.get("auth_seq_num")),
                "icode": (row.get("pdb_ins_code") or "").strip(),
                "resname": (row.get("mon_id") or "").strip().upper(),
            }
        )
    construct.sort(key=lambda item: (item["label_seq_id"] if item["label_seq_id"] is not None else 10**9))

    sifts_pairs: dict[tuple[str, int], tuple[int, str]] = {}
    db_names: set[str] = set()
    for row in xref:
        db_name = _val(row, "xref_db_name").upper()
        if db_name:
            db_names.add(db_name)
        acc = _val(row, "unp_acc", "xref_db_acc").upper()
        # Two SIFTS-in-mmcIF layouts exist: the current one uses unp_acc/unp_num,
        # older files use xref_db_name/xref_db_acc/db_num.  Both are accepted,
        # and '?' placeholders are treated as absent rather than as values.
        # The accession match alone decides membership: xref_db_name is an
        # auxiliary label and rejecting rows on it silently dropped 90% of the
        # real 9NU4 mapping.
        if acc != accession.upper():
            continue
        asym = _val(row, "asym_id")
        seq_id = _int(_val(row, "seq_id"))
        db_num = _int(_val(row, "unp_num", "db_num"))
        if asym and seq_id is not None and db_num is not None:
            sifts_pairs[(asym, seq_id)] = (db_num, _val(row, "mon_id").upper())
    report = MappingReport(
        pdb_id=pdb_id,
        chain=chain,
        accession=accession,
        native_length=len(native_sequence),
        source="sifts_pdbx_sifts_xref_db" if sifts_pairs else "struct_ref_seq",
    )

    if sifts_pairs:
        for item in construct:
            key = (item["label_asym_id"], item["label_seq_id"])
            pair = sifts_pairs.get(key)
            native_pos = pair[0] if pair else None
            report.rows.append(_make_row(pdb_id, chain, item, native_pos, native_sequence, "sifts"))
        report.provenance["sifts_rows"] = len(sifts_pairs)
        report.provenance["sifts_db_names"] = sorted(db_names)
    else:
        offset = _offset_from_struct_ref_seq(ref_seq, ref, chain, accession)
        report.provenance["struct_ref_seq"] = ref_seq
        if offset is not None:
            report.source = "struct_ref_seq_offset_verified_by_sequence"
            for item in construct:
                native_pos = None
                if item["auth_seq_id"] is not None:
                    candidate = item["auth_seq_id"] + offset
                    if 1 <= candidate <= len(native_sequence):
                        native_pos = candidate
                report.rows.append(_make_row(pdb_id, chain, item, native_pos, native_sequence, "struct_ref_seq"))
        else:
            report.source = "global_sequence_alignment"
            report.warnings.append(
                "no SIFTS mapping and no usable _struct_ref_seq range; used a global alignment "
                "between the observed chain sequence and the natural sequence"
            )
            alignment = align_sequences(
                "".join(THREE_TO_ONE.get(item["resname"], "X") for item in construct), native_sequence
            )
            for index, item in enumerate(construct):
                native_pos = alignment["map"].get(index)
                report.rows.append(_make_row(pdb_id, chain, item, native_pos, native_sequence, "sequence_alignment"))
            report.provenance["alignment"] = {
                "identity": alignment["identity"],
                "gaps": alignment["gaps"],
                "score": alignment["score"],
            }

    _annotate(report, diffs, native_sequence)
    return report


def _val(row: dict[str, str], *keys: str) -> str:
    """First present, non-placeholder value among ``keys``."""
    for key in keys:
        value = (row.get(key) or "").strip()
        if value and value not in {"?", "."}:
            return value
    return ""


def _int(value: Any) -> int | None:
    try:
        text = str(value).strip()
        if not text or text in {"?", "."}:
            return None
        return int(float(text))
    except (TypeError, ValueError):
        return None


def _make_row(pdb_id, chain, item, native_pos, native_sequence, source) -> ResidueMapRow:
    native_resname = ""
    match = False
    if native_pos is not None and 1 <= native_pos <= len(native_sequence):
        native_resname = THREE_TO_ONE.get(native_sequence[native_pos - 1], native_sequence[native_pos - 1])
        match = THREE_TO_ONE.get(item["resname"], "X") == native_resname
    return ResidueMapRow(
        pdb_id=pdb_id,
        chain=chain,
        label_asym_id=item["label_asym_id"],
        label_seq_id=item["label_seq_id"],
        auth_seq_id=item["auth_seq_id"],
        icode=item["icode"],
        resname=item["resname"],
        uniprot_resseq=native_pos,
        uniprot_resname=native_resname,
        construct_match=bool(match),
        source=source,
    )


def _offset_from_struct_ref_seq(ref_seq, ref, chain, accession) -> int | None:
    acc_by_id = {}
    for row in ref:
        if (row.get("db_code") or "").upper() == accession.upper():
            acc_by_id[row.get("id") or row.get("align_id") or ""] = row
    for row in ref_seq:
        if (row.get("pdbx_strand_id") or "") != chain:
            continue
        ref_id = row.get("ref_id") or row.get("align_id") or ""
        if ref_id not in acc_by_id:
            continue
        auth_beg = _int(row.get("pdbx_auth_seq_align_beg"))
        db_beg = _int(row.get("db_align_beg"))
        if auth_beg is not None and db_beg is not None:
            return db_beg - auth_beg
    return None


def _annotate(report: MappingReport, diffs: list[dict[str, Any]], native_sequence: str) -> None:
    for row in report.rows:
        if row.uniprot_resseq is None:
            report.unmapped_construct.append(
                {"label_seq_id": row.label_seq_id, "auth_seq_id": row.auth_seq_id,
                 "resname": row.resname, "icode": row.icode}
            )
        elif not row.construct_match:
            report.mismatches.append(
                {
                    "label_seq_id": row.label_seq_id,
                    "auth_seq_id": row.auth_seq_id,
                    "resname": row.resname,
                    "uniprot_resseq": row.uniprot_resseq,
                    "uniprot_resname": row.uniprot_resname,
                    "kind": "sequence_difference",
                }
            )
    # engineered mutations / tags declared in the deposition
    for diff in diffs:
        strand = diff.get("pdbx_pdb_strand_id") or diff.get("pdbx_strand_id") or ""
        if strand != report.chain:
            continue
        seq_num = _int(diff.get("pdbx_auth_seq_num") or diff.get("seq_num") or diff.get("pdb_seq_num"))
        entry = {
            "auth_seq_id": seq_num,
            "label_seq_id": _int(diff.get("seq_num")),
            "db_mon_id": (diff.get("db_mon_id") or "").strip(),
            "pdb_mon_id": (diff.get("mon_id") or diff.get("pdb_mon_id") or "").strip(),
            "details": (diff.get("details") or "").strip().strip("'\""),
            "source": "struct_ref_seq_dif",
        }
        if entry not in report.mutations:
            report.mutations.append(entry)
    mapped = sorted({row.uniprot_resseq for row in report.rows if row.uniprot_resseq})
    if mapped and native_sequence:
        present = set(mapped)
        report.missing_native = [i for i in range(1, len(native_sequence) + 1) if i not in present]


# --------------------------------------------------------------------------- #
# alignment (deterministic Needleman-Wunsch)
# --------------------------------------------------------------------------- #


def align_sequences(query: str, reference: str, *, match: int = 2, mismatch: int = -1, gap: int = -2) -> dict[str, Any]:
    """Global alignment; returns query-index -> reference-index map (1-based)."""
    n, m = len(query), len(reference)
    if n == 0 or m == 0:
        return {"map": {}, "identity": 0.0, "gaps": max(n, m), "score": 0}
    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = gap * i
        trace[i][0] = 1
    for j in range(1, m + 1):
        score[0][j] = gap * j
        trace[0][j] = 2
    for i in range(1, n + 1):
        qi = query[i - 1]
        for j in range(1, m + 1):
            diagonal = score[i - 1][j - 1] + (match if qi == reference[j - 1] else mismatch)
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap
            best = max(diagonal, up, left)
            score[i][j] = best
            trace[i][j] = 0 if best == diagonal else (1 if best == up else 2)
    mapping: dict[int, int] = {}
    i, j = n, m
    identical = 0
    gaps = 0
    while i > 0 or j > 0:
        step = trace[i][j] if (i > 0 and j > 0) else (1 if i > 0 else 2)
        if step == 0:
            mapping[i - 1] = j
            if query[i - 1] == reference[j - 1]:
                identical += 1
            i -= 1
            j -= 1
        elif step == 1:
            gaps += 1
            i -= 1
        else:
            gaps += 1
            j -= 1
    return {
        "map": mapping,
        "identity": identical / max(n, 1),
        "gaps": gaps,
        "score": score[n][m],
    }


def write_mapping(report: MappingReport, out_dir: str) -> dict[str, str]:
    ensure_dir(out_dir)
    rows_path = os.path.join(out_dir, f"residue_map_{report.pdb_id}_{report.chain}.csv")
    header = "pdb,chain,label_asym_id,label_seq_id,auth_seq_id,icode,resname,uniprot_resseq,uniprot_resname,construct_match,source\n"
    lines = [header]
    for row in report.rows:
        lines.append(
            f"{row.pdb_id},{row.chain},{row.label_asym_id},{row.label_seq_id or ''},{row.auth_seq_id if row.auth_seq_id is not None else ''},"
            f"{row.icode},{row.resname},{row.uniprot_resseq or ''},{row.uniprot_resname},{row.construct_match},{row.source}\n"
        )
    atomic_write_text(rows_path, "".join(lines))
    json_path = os.path.join(out_dir, f"mapping_{report.pdb_id}_{report.chain}.json")
    atomic_write_json(json_path, report.to_dict())
    return {"csv": rows_path, "json": json_path}


def require_chain_coverage(report: MappingReport, residues: Iterable[int]) -> dict[str, Any]:
    """Check that configured native residues actually exist in the mapping."""
    mapped = {row.uniprot_resseq: row for row in report.rows if row.uniprot_resseq is not None}
    missing = [r for r in residues if r not in mapped]
    return {
        "requested": sorted(residues),
        "missing": missing,
        "found": {r: mapped[r].auth_seq_id for r in residues if r in mapped},
    }


