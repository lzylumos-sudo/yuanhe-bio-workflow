"""Probe PoseBusters directly on one real pose and one deliberately broken pose.

Run with the project interpreter.  Writes a JSON record next to this file so the
exception text, the failing path and the per-check verdicts stay auditable.
"""

from __future__ import annotations

import json
import os
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
SRC = os.path.join(PROJECT, "src")
sys.path.insert(0, SRC)

from zvss import compat  # noqa: E402


def main() -> int:
    report: dict[str, object] = {}
    print("scratch_root:", compat.scratch_root())
    print("ascii_scratch_available:", compat.ascii_scratch_available())
    try:
        state = compat.bootstrap()
        report["bootstrap"] = state
    except Exception:
        report["bootstrap"] = traceback.format_exc()

    from rdkit import Chem
    from rdkit import RDConfig

    print("RDDataDir:", RDConfig.RDDataDir)
    print("RDDataDir ascii:", compat._is_ascii(RDConfig.RDDataDir))
    print("BaseFeatures exists:", os.path.exists(os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")))

    try:
        import posebusters

        report["posebusters_version"] = getattr(posebusters, "__version__", "unknown")
    except Exception:
        report["posebusters_import_error"] = traceback.format_exc()
        print(json.dumps(report, indent=1))
        return 2

    ligand = os.path.join(r"F:\MemPilot\runs\public20_engineering\entries",
                          "P000005_S01_C01_seed20261015", "input.sdf")
    receptor = os.path.join(r"F:\MemPilot\runs\public20_engineering", "protocol_inputs", "receptor_prepared.pdb")
    mol = Chem.MolFromMolFile(ligand, removeHs=False) if os.path.exists(ligand) else None
    if mol is None:
        mol = Chem.MolFromMolBlock(open(ligand, encoding="utf-8").read(), removeHs=False)
    print("ligand atoms:", mol.GetNumAtoms() if mol else None)

    from posebusters import PoseBusters

    buster = PoseBusters(config="dock")
    real = buster.bust(mol, None, receptor)
    report["real_pose"] = {"failed": [c for c in real.columns if real[c].iloc[0] is False or real[c].iloc[0] == 0]}
    print(real.T.to_string())

    broken = Chem.Mol(mol)
    conf = broken.GetConformer()
    for i in range(broken.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (pos.x + 8.0, pos.y, pos.z))
    bad = buster.bust(broken, None, receptor)
    report["broken_pose"] = {"failed": [c for c in bad.columns if bad[c].iloc[0] is False or bad[c].iloc[0] == 0]}
    print(bad.T.to_string())

    out = os.path.join(HERE, "probe_posebusters.json")
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=1, default=str)
    print("written:", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
