"""Ligand CSV ingest (plan §12). Every source row stays traceable; bad rows are isolated, never dropped."""
from __future__ import annotations
import csv, io
from dataclasses import dataclass, field
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolHash
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")
ALLOWED = {"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"}   # intersection of prep/docking/FF support (conservative)
KNOWN_COUNTERIONS = {"[Na+]", "[K+]", "[Cl-]", "[Br-]", "[I-]", "O", "[Li+]", "[Ca+2]", "[Mg+2]"}


@dataclass
class Failure:
    source_row: int; compound_id: str; stage: str; reason_code: str; retryable: bool; original_smiles: str; detail: str = ""


@dataclass
class IngestResult:
    parents: list = field(default_factory=list)       # one per unique stereo-aware structure
    source_map: list = field(default_factory=list)    # every source row -> parent_id (or failure)
    failures: list = field(default_factory=list)
    flags: list = field(default_factory=list)         # review flags (salts, unknown stereo, ...)


def read_csv(data: bytes, id_col: str | None = None, smiles_col: str | None = None):
    text = data.decode("utf-8-sig")
    rdr = csv.DictReader(io.StringIO(text, newline=""))
    cols = rdr.fieldnames or []
    low = {c.lower().strip(): c for c in cols}
    smiles_col = smiles_col or low.get("smiles") or low.get("smi")
    if smiles_col is None or smiles_col not in cols:
        raise ValueError(f"SMILES column not found; columns={cols}. Provide explicit mapping.")
    id_col = id_col or low.get("compound_id") or low.get("id")
    return list(rdr), smiles_col, id_col


def _salt_policy(mol):
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) == 1:
        return mol, None
    keep, drop = [], []
    for f in frags:
        (drop if Chem.MolToSmiles(f) in KNOWN_COUNTERIONS else keep).append(f)
    if len(keep) == 1:                      # explicit known counter-ions only; no "largest fragment" guess
        return keep[0], f"counterions_removed:{sorted(Chem.MolToSmiles(f) for f in drop)}"
    return None, "MULTI_FRAGMENT_NEEDS_POLICY"


def ingest_ligands(data: bytes, id_col=None, smiles_col=None) -> IngestResult:
    rows, scol, icol = read_csv(data, id_col, smiles_col)
    res, seen_ids, by_key = IngestResult(), set(), {}
    for i, row in enumerate(rows, start=2):               # line 1 is header
        cid = (row.get(icol) or "").strip() if icol else ""
        cid = cid or f"ROW{i:06d}"
        smi = (row.get(scol) or "").strip()

        def fail(stage, code, detail="", retry=False):
            res.failures.append(Failure(i, cid, stage, code, retry, smi, detail))
            res.source_map.append({"source_row": i, "compound_id": cid, "parent_id": None, "status": "failed", "reason_code": code})
        if cid in seen_ids:
            fail("validate", "DUPLICATE_ID", "compound_id repeated"); continue
        seen_ids.add(cid)
        if not smi:
            fail("validate", "EMPTY_SMILES"); continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None or mol.GetNumAtoms() == 0:
            fail("validate", "INVALID_SMILES"); continue
        mol2, note = _salt_policy(mol)
        if mol2 is None:
            fail("fragments", "SALT_POLICY_REQUIRED", note); continue
        bad = sorted({a.GetSymbol() for a in mol2.GetAtoms()} - ALLOWED)   # checked on the kept parent only
        if bad:
            fail("validate", "UNSUPPORTED_ELEMENT", ",".join(bad)); continue
        if note:
            res.flags.append({"source_row": i, "compound_id": cid, "flag": "SALT_STRIPPED", "detail": note})
        std = rdMolStandardize.Cleanup(mol2)
        iso = Chem.MolToSmiles(std)                                   # stereo-aware canonical
        key = iso                                                      # different known stereoisomers never merge
        unspec = len(Chem.FindMolChiralCenters(std, includeUnassigned=True, useLegacyImplementation=False)) \
            - len(Chem.FindMolChiralCenters(std, includeUnassigned=False, useLegacyImplementation=False))
        if key not in by_key:
            pid = f"P{len(by_key) + 1:06d}"
            by_key[key] = {"parent_id": pid, "isomeric_smiles": iso, "first_row": i, "source_rows": [], "unspecified_stereocenters": unspec}
            res.parents.append(by_key[key])
            if unspec:
                res.flags.append({"source_row": i, "compound_id": cid, "flag": "UNKNOWN_STEREO", "detail": f"{unspec} unspecified center(s): review or enumerate"})
        by_key[key]["source_rows"].append({"source_row": i, "compound_id": cid, "original_smiles": smi, "name": row.get("name", "")})
        res.source_map.append({"source_row": i, "compound_id": cid, "parent_id": by_key[key]["parent_id"], "status": "ok", "reason_code": ""})
    return res
