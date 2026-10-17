"""Persist prepared ligand states as files and bind approvals to their exact content (plan §14-15).

Per state:  states/<state_id>/{state.sdf (all kept conformers + Gasteiger charges as atom props), state.json, review.svg}
The approval hash set is {sdf, svg, state_json} so a changed structure, conformer, charge or image invalidates the approval.
"""
from __future__ import annotations
import json
from pathlib import Path
from rdkit import Chem
from ..hashing import sha256_file, atomic_write_json, atomic_write_bytes
from .review import render_review_svg


def _conf_mol_list(rec):
    mh, keep = rec["_mol"], rec["_keep"]
    out = []
    for k, cid in enumerate(keep):
        m = Chem.Mol(mh, confId=cid) if False else Chem.Mol(mh)
        conf = mh.GetConformer(cid)
        m.RemoveAllConformers(); m.AddConformer(Chem.Conformer(conf), assignId=True)
        q = rec["charges"]["values"]
        m.SetProp("atom.dprop.gasteiger_charge", " ".join(f"{float(x):.6f}" for x in q))   # SDF atom-property list
        m.SetProp("_Name", f'{rec["state_id"]}_c{k}'); m.SetProp("state_id", rec["state_id"]); m.SetProp("state_hash", rec["state_hash"])
        m.SetProp("charge_method", rec["charges"]["method"]); m.SetIntProp("conformer_index", k)
        out.append(m)
    return out


def save_state(rec: dict, parent: dict, out_root: Path) -> dict:
    """Write one state. Returns {'dir','hashes':{sdf,svg,state_json}}. Refuses to overwrite a different existing state."""
    d = Path(out_root) / rec["state_id"]
    meta = {k: v for k, v in rec.items() if not k.startswith("_") and k != "charges"}
    meta["charges"] = {"method": rec["charges"]["method"], "check": rec["charges"]["check"]} if "charges" in rec else None
    meta["parent_smiles"] = parent.get("smiles") or parent.get("canonical_smiles")
    if d.exists():
        old = json.loads((d / "state.json").read_text())
        if old.get("state_hash") != rec["state_hash"]:
            raise FileExistsError(f"STATE_DIR_CONFLICT {d}: different content; use a new output root, never overwrite")
        return {"dir": str(d), "hashes": {"sdf": sha256_file(d / "state.sdf"), "svg": sha256_file(d / "review.svg"), "state_json": sha256_file(d / "state.json")}, "reused": True}
    d.mkdir(parents=True)
    if "_mol" in rec:
        w = Chem.SDWriter(str(d / "state.sdf"))
        for m in _conf_mol_list(rec):
            w.write(m)
        w.close()
    else:
        (d / "state.sdf").write_text("")
    atomic_write_bytes(d / "review.svg", render_review_svg(rec, meta["parent_smiles"]).encode("utf-8"))
    atomic_write_json(d / "state.json", meta)
    return {"dir": str(d), "hashes": {"sdf": sha256_file(d / "state.sdf"), "svg": sha256_file(d / "review.svg"), "state_json": sha256_file(d / "state.json")}, "reused": False}


def load_state_sdf(d: Path):
    return [m for m in Chem.SDMolSupplier(str(Path(d) / "state.sdf"), removeHs=False) if m is not None]
