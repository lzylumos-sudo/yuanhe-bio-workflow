"""Acceptance run: PoseBusters on a real pose and on deliberately invalid poses.

This is the executable evidence behind the PoseBusters section of the completion
matrix.  It uses the same pinned PoseBusters build, the same receptor and the
same input form as the pipeline (positional ``bust(mol_pred, mol_true, mol_cond)``)
and writes both a JSON record and a Markdown table.
"""

from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PROJECT, "src"))

from zvss import compat  # noqa: E402

compat.bootstrap()

from rdkit import Chem  # noqa: E402

ENTRY = r"F:\MemPilot\runs\public20_engineering\entries\P000005_S01_C01_seed20261015"
LIGAND = os.path.join(ENTRY, "input.sdf")
RECEPTOR = r"F:\MemPilot\runs\public20_engineering\protocol_inputs\receptor_prepared.pdb"
OUT_DIR = os.path.join(PROJECT, "docs", "acceptance")


def load_ligand():
    text = open(LIGAND, encoding="utf-8").read().split("$$$$")[0]
    return Chem.MolFromMolBlock(text, removeHs=False)


def translated(mol, dx: float):
    out = Chem.Mol(mol)
    conf = out.GetConformer()
    for index in range(out.GetNumAtoms()):
        position = conf.GetAtomPosition(index)
        conf.SetAtomPosition(index, (position.x + dx, position.y, position.z))
    return out


def stretched_bond(mol):
    out = Chem.Mol(mol)
    conf = out.GetConformer()
    bond = out.GetBondWithIdx(0)
    begin = bond.GetBeginAtomIdx()
    end = bond.GetEndAtomIdx()
    start = conf.GetAtomPosition(begin)
    stop = conf.GetAtomPosition(end)
    vector = (stop.x - start.x, stop.y - start.y, stop.z - start.z)
    length = max((vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2) ** 0.5, 1e-6)
    scale = 2.2 / length
    conf.SetAtomPosition(
        end, (start.x + vector[0] * scale, start.y + vector[1] * scale, start.z + vector[2] * scale)
    )
    return out


def overlapped_atoms(mol):
    out = Chem.Mol(mol)
    conf = out.GetConformer()
    first = conf.GetAtomPosition(0)
    for index in range(1, out.GetNumAtoms()):
        conf.SetAtomPosition(index, (first.x + 0.4, first.y + 0.4, first.z + 0.4))
    return out


def main() -> int:
    from posebusters import PoseBusters

    buster = PoseBusters(config="dock")
    base = load_ligand()
    cases = {
        "real_pose": base,
        "translated_8A_into_protein": translated(base, 8.0),
        "stretched_bond_2.2A": stretched_bond(base),
        "all_atoms_overlapped": overlapped_atoms(base),
    }
    results = {}
    for name, mol in cases.items():
        table = buster.bust(mol, None, RECEPTOR)
        row = table.iloc[0].to_dict()
        items = {key: bool(value) for key, value in row.items() if isinstance(value, bool)}
        failed = [key for key, value in items.items() if not value]
        results[name] = {
            "atoms": mol.GetNumAtoms(),
            "checks": len(items),
            "passed": len(items) - len(failed),
            "failed": failed,
            "items": items,
        }
    os.makedirs(OUT_DIR, exist_ok=True)
    payload = {
        "posebusters_version": __import__("posebusters").__version__,
        "scratch_root": compat.scratch_root(),
        "ligand": LIGAND,
        "receptor": RECEPTOR,
        "results": results,
    }
    with open(os.path.join(OUT_DIR, "posebusters_run.json"), "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=1)

    lines = [
        "# PoseBusters 真实执行记录",
        "",
        f"- PoseBusters 版本: `{payload['posebusters_version']}`",
        f"- ASCII 运行时临时目录: `{payload['scratch_root']}`",
        f"- 真实位姿来源: `{LIGAND}`",
        f"- 受体: `{RECEPTOR}`",
        "",
        "| 用例 | 原子数 | 检查项 | 通过 | 失败项 |",
        "| --- | --- | --- | --- | --- |",
    ]
    for name, record in results.items():
        lines.append(
            f"| {name} | {record['atoms']} | {record['checks']} | {record['passed']} | "
            f"{', '.join(record['failed']) or '无'} |"
        )
    lines += [
        "",
        "## 判定",
        "",
        f"- 真实位姿失败项数: {len(results['real_pose']['failed'])}（应为 0）",
        f"- 平移进蛋白位姿失败项: {', '.join(results['translated_8A_into_protein']['failed'])}",
        f"- 拉长键位姿失败项: {', '.join(results['stretched_bond_2.2A']['failed'])}",
        f"- 原子重叠位姿失败项: {', '.join(results['all_atoms_overlapped']['failed'])}",
        "",
        "失败项逐条来自 PoseBusters 的布尔检查列，未经任何后处理或裁剪。",
        "",
    ]
    with open(os.path.join(OUT_DIR, "POSEBUSTERS_RUN.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
    print(json.dumps({name: record["failed"] for name, record in results.items()}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
