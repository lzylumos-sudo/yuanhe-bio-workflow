"""Self-contained chemistry in PDB REMARK records; no external services.

This container transports an already verified SDF; it does not infer chemistry.
Consumers must still validate its chemical graph and atom mapping against PDB.
The implementation is copied into each independent Nexus project.
"""
from __future__ import annotations
import base64
import gzip
import hashlib
import json
import re
import zlib

PREFIX = "REMARK 950 NEXUS1 "
MAX_PAYLOAD = 8 * 1024 * 1024
MAX_CHUNKS = 100000

class EmbeddedChemistryError(ValueError):
    pass

def _validate_payload(payload):
    if not isinstance(payload, dict) or payload.get("schema") != "nexus-complex-1":
        raise EmbeddedChemistryError("Unsupported embedded Nexus chemistry schema")
    if not isinstance(payload.get("manifest"), dict):
        raise EmbeddedChemistryError("Embedded manifest must be an object")
    sdf = payload.get("ligand_sdf")
    if sdf is not None and not isinstance(sdf, str):
        raise EmbeddedChemistryError("Embedded ligand SDF must be text or null")
    expected = hashlib.sha256(sdf.encode("utf-8")).hexdigest() if sdf is not None else None
    if payload.get("ligand_sdf_sha256") != expected:
        raise EmbeddedChemistryError("Embedded SDF hash mismatch")
    if payload["manifest"].get("ligand") and sdf is None:
        raise EmbeddedChemistryError("A ligand manifest requires its exact embedded SDF")
    return payload

def write_embedded(pdb_text: str, manifest: dict, sdf_text: str | None = None) -> str:
    """Prepend <=80-column REMARK records; leave all coordinate records untouched."""
    if not isinstance(pdb_text, str) or not isinstance(manifest, dict):
        raise EmbeddedChemistryError("PDB text and a manifest object are required")
    payload = _validate_payload({
        "schema": "nexus-complex-1", "manifest": manifest,
        "ligand_sdf": sdf_text,
        "ligand_sdf_sha256": hashlib.sha256(sdf_text.encode("utf-8")).hexdigest() if sdf_text is not None else None,
    })
    raw = json.dumps(payload, ensure_ascii=True, sort_keys=True, separators=(",", ":")).encode("utf-8")
    if len(raw) > MAX_PAYLOAD:
        raise EmbeddedChemistryError("Embedded chemistry exceeds the 8 MiB limit")
    encoded = base64.b64encode(gzip.compress(raw, mtime=0)).decode("ascii")
    chunks = [encoded[i:i + 40] for i in range(0, len(encoded), 40)]
    if len(chunks) > MAX_CHUNKS:
        raise EmbeddedChemistryError("Too many embedded chemistry records")
    digest = hashlib.sha256(raw).hexdigest()
    records = [PREFIX + "FORMAT GZIP_JSON_V1", PREFIX + "SHA256 1 " + digest[:32], PREFIX + "SHA256 2 " + digest[32:]]
    records += [PREFIX + f"DATA {i:06d} {len(chunks):06d} {chunk}" for i, chunk in enumerate(chunks, 1)]
    assert all(len(line) <= 80 for line in records)
    # Only replace our own namespace, preserving every other PDB line verbatim.
    clean = "".join(line for line in pdb_text.splitlines(keepends=True) if not line.startswith(PREFIX))
    return "\n".join(records) + "\n" + clean

def read_embedded(pdb_text: str) -> dict | None:
    records = [line.rstrip("\r\n") for line in pdb_text.splitlines() if line.startswith(PREFIX)]
    if not records:
        return None
    if records.count(PREFIX + "FORMAT GZIP_JSON_V1") != 1:
        raise EmbeddedChemistryError("Missing, duplicate, or unsupported Nexus metadata header")
    hashes, data = {}, {}
    total = None
    for line in records:
        if line == PREFIX + "FORMAT GZIP_JSON_V1":
            continue
        suffix = line[len(PREFIX):]
        match = re.fullmatch(r"SHA256 ([12]) ([0-9a-f]{32})", suffix)
        if match:
            key = int(match.group(1))
            if key in hashes:
                raise EmbeddedChemistryError("Duplicate embedded hash record")
            hashes[key] = match.group(2)
            continue
        match = re.fullmatch(r"DATA (\d{6}) (\d{6}) ([A-Za-z0-9+/=]{1,40})", suffix)
        if not match:
            raise EmbeddedChemistryError("Malformed embedded Nexus record")
        index, count = int(match.group(1)), int(match.group(2))
        if not 1 <= count <= MAX_CHUNKS or not 1 <= index <= count or index in data:
            raise EmbeddedChemistryError("Invalid or duplicate embedded chunk index")
        if total is not None and total != count:
            raise EmbeddedChemistryError("Inconsistent embedded chunk counts")
        total = count
        data[index] = match.group(3)
    if set(hashes) != {1, 2} or total is None or len(data) != total:
        raise EmbeddedChemistryError("Incomplete embedded chemistry records")
    try:
        compressed = base64.b64decode("".join(data[i] for i in range(1, total + 1)), validate=True)
        decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)
        raw = decoder.decompress(compressed, MAX_PAYLOAD + 1)
        if len(raw) > MAX_PAYLOAD or decoder.unconsumed_tail or not decoder.eof or decoder.unused_data:
            raise EmbeddedChemistryError("Oversized, truncated, or concatenated embedded payload")
        if hashlib.sha256(raw).hexdigest() != hashes[1] + hashes[2]:
            raise EmbeddedChemistryError("Embedded chemistry payload hash mismatch")
        payload = json.loads(raw.decode("utf-8"))
    except (ValueError, zlib.error, UnicodeError) as exc:
        if isinstance(exc, EmbeddedChemistryError):
            raise
        raise EmbeddedChemistryError("Embedded chemistry decoding failed") from exc
    return _validate_payload(payload)
