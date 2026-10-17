"""Excel / CSV / JSON exports built from the single interactions.json event list."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

TYPE_LABEL = {
    "hydrogen_bond": "Hydrogen bond", "ionic": "Salt bridge / ionic", "metal_coordination": "Metal coordination",
    "halogen_bond": "Halogen bond", "pi_pi": "pi-pi", "cation_pi": "Cation-pi",
    "water_bridge": "Water bridge (candidate)", "hydrophobic": "Hydrophobic", "vdw_contact": "van der Waals contact",
    "steric_clash": "Steric clash (warning)",
}


def _names(atoms):
    return ", ".join(a["name"] for a in atoms)


def events_to_frame(result: dict) -> pd.DataFrame:
    rows = []
    for e in result["events"]:
        row = {
            "event_id": e["event_id"],
            "type": TYPE_LABEL.get(e["type"], e["type"]),
            "chain": e["residue"]["chain"], "residue": e["residue"]["resname"],
            "resseq": e["residue"]["resseq"], "icode": e["residue"]["icode"],
            "protein_part": e.get("protein_part", ""),
            "protein_atoms": _names(e["protein_atoms"]),
            "ligand_atoms": _names(e["ligand_atoms"]),
            "protein_atom_uid": "; ".join(a["uid"] for a in e["protein_atoms"]),
            "ligand_atom_uid": "; ".join(a["uid"] for a in e["ligand_atoms"]),
            "distance_A": e.get("distance"),
            "donor": "", "acceptor": "", "donor_side": e.get("donor_side", ""),
            "H_A_distance_A": e.get("h_a_distance"), "angle_DHA_deg": e.get("angle_DHA"),
            "plane_angle_deg": e.get("plane_angle"), "other_angle_deg": e.get("angle"),
            "geometry": e.get("geometry", ""), "water": e.get("water", {}).get("label", ""),
            "confidence": e.get("confidence", ""), "flags": ";".join(e.get("flags", [])),
            "engine": e.get("engine", ""), "rule": e.get("rule", ""),
        }
        if e["type"] == "hydrogen_bond":
            row["donor"] = f'{e["donor"]["uid"]}'
            row["acceptor"] = f'{e["acceptor"]["uid"]}'
        rows.append(row)
    return pd.DataFrame(rows)


def summary_frame(result: dict) -> pd.DataFrame:
    df = events_to_frame(result)
    if df.empty:
        return pd.DataFrame(columns=["residue_label", "n_events", "types"])
    df["residue_label"] = df["chain"] + ":" + df["residue"] + df["resseq"].astype(str) + df["icode"]
    g = df.groupby("residue_label").agg(
        n_events=("event_id", "count"), types=("type", lambda s: ", ".join(sorted(set(s)))),
        min_distance_A=("distance_A", "min"),
    ).reset_index()
    return g.sort_values("n_events", ascending=False)


def write_tables(result: dict, out_dir: Path) -> dict:
    out_dir = Path(out_dir)
    df = events_to_frame(result)
    sm = summary_frame(result)
    csv = out_dir / "interactions.csv"
    xlsx = out_dir / "interactions.xlsx"
    df.to_csv(csv, index=False, encoding="utf-8-sig")
    chem = result["chemistry"]
    meta = pd.DataFrame([
        ("ZLIS version", result["zlis_version"]), ("Input", result["input_pdb"]),
        ("Ligand", result["ligand"]["label"]), ("Model", result["model"]),
        ("Ligand SMILES", chem["smiles"]), ("Ligand charge", chem["formal_charge"]),
        ("Chemistry source", chem["source"]), ("Chemistry status", chem["status"]),
        ("Notes", " | ".join(chem["notes"] + result["protein_prep"].get("notes", []) + result["warnings"])),
        ("Waters in input", result["n_waters"]),
        ("Disclaimer", "Static-structure geometric analysis. Not a free-energy or force calculation."),
    ], columns=["field", "value"])
    with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
        df.to_excel(xw, sheet_name="Interactions", index=False)
        sm.to_excel(xw, sheet_name="Residue summary", index=False)
        meta.to_excel(xw, sheet_name="Run info", index=False)
        for name, d in (("Interactions", df), ("Residue summary", sm), ("Run info", meta)):
            ws = xw.sheets[name]
            ws.freeze_panes = "A2"
            for col in ws.columns:
                width = max(len(str(c.value)) if c.value is not None else 0 for c in col[:200])
                ws.column_dimensions[col[0].column_letter].width = min(max(10, width + 2), 60)
    return {"csv": str(csv), "xlsx": str(xlsx)}
