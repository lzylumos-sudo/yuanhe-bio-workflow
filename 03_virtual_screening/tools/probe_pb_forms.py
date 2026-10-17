"""Find the PoseBusters input form that works in this environment."""

from __future__ import annotations

import os
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PROJECT, "src"))

from zvss import compat  # noqa: E402

compat.bootstrap()
import pandas as pd  # noqa: E402
from rdkit import Chem  # noqa: E402
from posebusters import PoseBusters  # noqa: E402

LIGAND = r"F:\MemPilot\runs\public20_engineering\entries\P000005_S01_C01_seed20261015\input.sdf"
RECEPTOR = r"F:\MemPilot\runs\public20_engineering\protocol_inputs\receptor_prepared.pdb"
mol = Chem.MolFromMolBlock(open(LIGAND, encoding="utf-8").read().split("$$$$")[0], removeHs=False)
print("atoms:", mol.GetNumAtoms(), "conformers:", mol.GetNumConformers())
buster = PoseBusters(config="dock")

def attempt(name, func):
    try:
        out = func()
        row = out.iloc[0]
        items = {k: bool(v) for k, v in row.to_dict().items() if isinstance(v, bool)}
        print(f"{name}: OK items={len(items)} failed={[k for k, v in items.items() if not v]}")
    except Exception as exc:  # noqa: BLE001
        print(f"{name}: {type(exc).__name__}: {exc}")

attempt("positional", lambda: buster.bust(mol, None, RECEPTOR))
attempt("frame_receptor_path", lambda: buster.bust(pd.DataFrame([{"mol_pred": mol, "mol_true": None, "mol_cond": RECEPTOR}]), None, None))
attempt("frame_no_cond", lambda: buster.bust(pd.DataFrame([{"mol_pred": mol}]), None, None))
attempt("frame_mol_true_none_cond_none", lambda: buster.bust(pd.DataFrame([{"mol_pred": mol, "mol_true": None, "mol_cond": None}]), None, None))
attempt("frame_only_pred_key", lambda: buster.bust(pd.DataFrame([{"mol_pred": mol, "mol_cond": RECEPTOR}]), None, None))
