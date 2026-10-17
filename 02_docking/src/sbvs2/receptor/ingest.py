"""Protein ingest + structure identity (plan §08, §09). Pure python, no third-party deps.

Never guesses: missing residues are only *reported* and classified; unknown sequence -> MISSING_SEQUENCE.
"""
from __future__ import annotations
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from ..hashing import sha256_bytes

WATER = {"HOH", "WAT", "DOD", "H2O"}
IONS = {"NA", "CL", "K", "MG", "CA", "ZN", "MN", "FE", "FE2", "CU", "CU1", "NI", "CO", "CD", "HG", "BR", "IOD", "F"}
METALS = {"ZN", "MN", "FE", "FE2", "CU", "CU1", "NI", "CO", "CD", "HG", "MG", "CA"}
BUFFERS = {"SO4", "PO4", "GOL", "EDO", "PEG", "ACT", "DMS", "TRS", "MES", "EPE", "FMT", "BME", "PGE", "PG4"}
SUGARS = {"NAG", "MAN", "BMA", "FUC", "GAL", "GLC", "NDG"}
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
       "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
       "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C"}
SIDECHAIN_HEAVY = {  # expected heavy atoms (incl. backbone N,CA,C,O)
    "ALA": 5, "ARG": 11, "ASN": 8, "ASP": 8, "CYS": 6, "GLN": 9, "GLU": 9, "GLY": 4, "HIS": 10, "ILE": 8,
    "LEU": 8, "LYS": 9, "MET": 8, "PHE": 11, "PRO": 7, "SER": 6, "THR": 7, "TRP": 14, "TYR": 12, "VAL": 7}


@dataclass
class Atom:
    record: str; serial: int; name: str; altloc: str; resname: str; chain: str; resseq: int; icode: str
    x: float; y: float; z: float; occ: float; bfac: float; element: str; model: int; line: str

    @property
    def rkey(self):
        return (self.chain, self.resseq, self.icode)


def _elem(line: str, name: str) -> str:
    e = line[76:78].strip() if len(line) >= 78 else ""
    if e:
        return e.upper()
    n = name.strip().lstrip("0123456789")
    return (n[:1] if n else "X").upper()


def parse_pdb(text: str) -> tuple[list[Atom], dict]:
    atoms, meta, model = [], {"seqres": defaultdict(list), "remark465": [], "ssbond": [], "link": [], "models": 0,
                              "title": ""}, 1
    for ln in text.splitlines():
        rec = ln[:6].strip()
        if rec == "MODEL":
            model = int(ln[10:14]); meta["models"] += 1
        elif rec in ("ATOM", "HETATM"):
            try:
                a = Atom(rec, int(ln[6:11]), ln[12:16].strip(), ln[16].strip(), ln[17:20].strip(), ln[21].strip() or "A",
                         int(ln[22:26]), ln[26].strip(), float(ln[30:38]), float(ln[38:46]), float(ln[46:54]),
                         float(ln[54:60] or 1.0), float(ln[60:66] or 0.0), _elem(ln, ln[12:16]), model, ln)
            except (ValueError, IndexError) as e:
                raise ValueError(f"unparseable coordinate line: {ln[:40]!r}: {e}")
            if not all(map(math.isfinite, (a.x, a.y, a.z))):
                raise ValueError("non-finite coordinate")
            atoms.append(a)
        elif rec == "SEQRES":
            meta["seqres"][ln[11].strip() or "A"].extend(ln[19:].split())
        elif rec == "REMARK" and ln[7:10] == "465":
            toks = ln[10:].split()
            if len(toks) >= 3 and toks[0] in AA3 and toks[-1].lstrip("-").isdigit():
                meta["remark465"].append({"resname": toks[0], "chain": toks[1], "resseq": int(toks[2])}) if len(toks) == 3 else None
        elif rec == "SSBOND":
            meta["ssbond"].append(ln.rstrip())
        elif rec == "LINK":
            meta["link"].append(ln.rstrip())
        elif rec == "TITLE":
            meta["title"] += ln[10:].strip() + " "
    meta["models"] = max(meta["models"], 1)
    meta["seqres"] = dict(meta["seqres"])
    return atoms, meta


def classify_hetero(resname: str) -> str:
    if resname in WATER: return "water"
    if resname in METALS: return "metal"
    if resname in IONS: return "ion"
    if resname in BUFFERS: return "buffer"
    if resname in SUGARS: return "glycan"
    return "ligand_or_cofactor"


@dataclass
class IngestReport:
    sha256: str
    n_models: int
    chains: dict
    hetero: dict
    altlocs: list
    missing_heavy_atoms: list
    gaps: list
    missing_from_seqres: list
    issues: list = field(default_factory=list)   # {code, severity, detail}
    plan: list = field(default_factory=list)     # per-case graded decisions


def ingest_pdb(path: str | Path, pocket_center: tuple | None = None) -> IngestReport:
    raw = Path(path).read_bytes()
    atoms, meta = parse_pdb(raw.decode("utf-8", errors="replace"))
    if not atoms:
        raise ValueError("no ATOM/HETATM records")
    issues, plan = [], []
    if meta["models"] > 1:
        issues.append({"code": "MULTI_MODEL", "severity": "review",
                       "detail": f"{meta['models']} models; models are never merged, pick one explicitly"})
    first_model = atoms[0].model
    atoms = [a for a in atoms if a.model == first_model]

    prot = [a for a in atoms if a.record == "ATOM" and a.resname in AA3]
    # altlocs: report only, caller chooses a *whole-residue* occupancy-consistent set
    alt = defaultdict(set)
    for a in atoms:
        if a.altloc:
            alt[a.rkey + (a.resname,)].add(a.altloc)
    altlocs = [{"chain": k[0], "resseq": k[1], "icode": k[2], "resname": k[3], "altlocs": sorted(v)} for k, v in alt.items()]

    residues = defaultdict(list)
    for a in prot:
        residues[a.rkey].append(a)
    chains = defaultdict(lambda: {"residues": [], "sequence_from_atoms": ""})
    for (c, n, i), ats in residues.items():
        chains[c]["residues"].append((n, i, ats[0].resname))
    for c, d in chains.items():
        d["residues"].sort()
        d["length"] = len(d["residues"])
        d["sequence_from_atoms"] = "".join(AA3.get(r[2], "X") for r in d["residues"])
        d["seqres_length"] = len(meta["seqres"].get(c, []))
        d["seqres_sequence"] = "".join(AA3.get(r, "X") for r in meta["seqres"].get(c, []))
        d["residues"] = f"{d['residues'][0][0]}-{d['residues'][-1][0]}"

    # missing heavy atoms (only residue types we know the expected count for)
    missing_atoms = []
    for (c, n, i), ats in residues.items():
        names = {a.name for a in ats if a.element != "H" and not a.altloc or a.altloc in ("", "A")}
        rn = ats[0].resname
        if rn in SIDECHAIN_HEAVY and len(names) < SIDECHAIN_HEAVY[rn]:
            missing_atoms.append({"chain": c, "resseq": n, "resname": rn, "have": len(names), "expect": SIDECHAIN_HEAVY[rn]})
    # chain breaks: CA-CA > 4.2 A between sequence neighbours or residue number gap
    gaps = []
    for c in chains:
        cas = sorted([(k[1], k[2], [a for a in v if a.name == "CA"]) for k, v in residues.items() if k[0] == c], key=lambda t: (t[0], t[1]))
        for p, q in zip(cas, cas[1:]):
            if not p[2] or not q[2]:
                continue
            d = math.dist((p[2][0].x, p[2][0].y, p[2][0].z), (q[2][0].x, q[2][0].y, q[2][0].z))
            if q[0] - p[0] > 1 or d > 4.2:
                gaps.append({"chain": c, "after_resseq": p[0], "before_resseq": q[0], "n_missing_by_numbering": q[0] - p[0] - 1,
                             "ca_ca": round(d, 2)})
    missing_seq = []
    for c, d in chains.items():
        if d["seqres_length"] and d["seqres_length"] > d["length"]:
            missing_seq.append({"chain": c, "seqres_len": d["seqres_length"], "modeled": d["length"]})

    # hetero
    het = defaultdict(list)
    for a in atoms:
        if a.record == "HETATM":
            het[(a.resname, a.chain, a.resseq, a.icode)].append(a)
    hetero = defaultdict(list)
    for (rn, c, n, i), ats in het.items():
        hetero[classify_hetero(rn)].append({"resname": rn, "chain": c, "resseq": n, "n_atoms": len(ats)})

    # graded plan (plan §09)
    near_pocket = lambda g: False
    if pocket_center is not None:
        pc = pocket_center
        def near_pocket(g):
            for k, v in residues.items():
                if k[0] == g["chain"] and g["after_resseq"] <= k[1] <= g["before_resseq"]:
                    for a in v:
                        if math.dist((a.x, a.y, a.z), pc) <= 6.0 + 10:
                            return True
            return False
    for m in missing_atoms:
        plan.append({"case": "missing_sidechain_atoms", "target": f"{m['chain']}:{m['resname']}{m['resseq']}",
                     "action": "PDBFixer/template completion", "auto": True})
    for g in gaps:
        n = g["n_missing_by_numbering"]
        has_seq = g["chain"] in meta["seqres"]
        if not has_seq:
            code, auto, act = "MISSING_SEQUENCE", False, "pause: request FASTA/source; never guess residues"
            issues.append({"code": code, "severity": "needs_input", "detail": f"gap {g['chain']}:{g['after_resseq']}-{g['before_resseq']} but no SEQRES"})
        elif n >= 6 or near_pocket(g):
            code, auto, act = "LOOP_REVIEW", False, "MODELLER-class loop modelling -> candidate set -> human choice"
            issues.append({"code": code, "severity": "review", "detail": f"gap of {n} residues {'near pocket' if near_pocket(g) else '>=6'}"})
        elif n >= 1:
            code, auto, act = "SHORT_LOOP", True, "short-loop pipeline, 10 candidates, only if away from pocket and QC passes"
        else:
            code, auto, act = "CHAIN_BREAK", False, "break with no numbering gap (possible cross-chain/modified residue): pause"
            issues.append({"code": code, "severity": "review", "detail": f"CA-CA {g['ca_ca']} A at {g['chain']}:{g['after_resseq']}"})
        plan.append({"case": code, "target": f"{g['chain']}:{g['after_resseq']}-{g['before_resseq']}", "action": act, "auto": auto})
    if altlocs:
        issues.append({"code": "ALTLOC", "severity": "review", "detail": f"{len(altlocs)} residues with alternate locations; choose a whole-residue consistent set"})
    for c in ("metal", "ligand_or_cofactor", "glycan"):
        if hetero.get(c):
            issues.append({"code": f"HETERO_{c.upper()}", "severity": "review",
                           "detail": f"{len(hetero[c])} {c} component(s) need explicit keep/remove decision (no blanket HETATM delete)"})
    if hetero.get("water"):
        plan.append({"case": "waters", "target": f"{len(hetero['water'])} waters", "action": "default delete bulk; list pocket/bridging waters as keep-candidates", "auto": False})
    return IngestReport(sha256_bytes(raw), meta["models"], dict(chains), dict(hetero), altlocs, missing_atoms, gaps, missing_seq, issues, plan)
