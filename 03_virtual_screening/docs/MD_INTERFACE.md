# MD 交付接口

## 1. 任务包目录

```
runs/<run_id>/md_handoff/<role>_<order>_<compound_id>/
├─ complex.pdb        # 受体 + 配体，配体为 HETATM，resname LIG，chain L，resid 1
├─ ligand.sdf         # 与 complex.pdb 中配体逐原子一致的 3D 结构
├─ manifest.json      # 机器契约：身份、原子映射、来源
├─ task.json          # 任务描述、环境要求、阻塞项、必检项
├─ README.md          # 人工交接说明（含如何喂给 04/05）
├─ adapter04.json     # 到 05_molecular_dynamics 的适配记录
├─ adapter05.json     # 到 04_interaction_analysis 的适配记录
└─ checksums.sha256   # 上述文件校验和
```

`complex.pdb` 头部写 `REMARK 900–903`（run_id / compound_id / state_id / pose_id），供下游追溯。

## 2. `manifest.json` 契约

```json
{
  "format": "zvss-md-manifest-1",
  "compound_id": "CMP...", "state_id": "STA...", "pose_id": "POS...",
  "ligand": {"resname": "LIG", "chain": "L", "resid": "1",
             "atoms": 9, "smiles": "...", "atom_mapping": {"C1": 1, "...": 0}},
  "receptor": {"path": "receptor.pdb", "version": "REC...", "atoms": 4598},
  "checksums": {"complex.pdb": "sha256:...", "ligand.sdf": "sha256:..."}
}
```

`atom_mapping` 把 PDB 原子名映射到 `ligand.sdf` 中的 1-based 原子序号，下游无需再猜配体原子顺序。

## 3. 与 `05_molecular_dynamics`（OpenMM）的接口

已有项目契约（只读参考，未修改）：`bash run_md.sh --input <complex.pdb>`，输入为 `complex.pdb` + 相邻 `ligand.sdf`/`manifest.json`（或 PDB 内嵌 `REMARK 950`）。

ZVSS 适配器（`md/adapters.py::adapter04`）：

- 校验 `complex.pdb` / `ligand.sdf` / `manifest.json` 三者存在且相互一致；
- 生成可直接执行的命令模板，写入 `adapter04.json`；
- 记录两条硬阻塞：`blocked_platform_linux_required`（该项目在非 Linux 直接退出）、`blocked_missing_membrane_capability`（该项目明确无膜/非水体系协议）；
- `executed: false`，`verified_offline: true`（只验证了文件契约，没有假装跑过 MD）。

## 4. 与 `04_interaction_analysis`（ZLIS/ProLIF）的接口

已有项目契约：`bash run.sh analyze <pdb> [--ligand A:BEN1] [--template lig.sdf]`。

ZVSS 适配器（`md/adapters.py::adapter05`）：

- 由 `manifest.json` 生成选择器（例如 `L:LIG1`）与完整命令；
- 记录输出清单（interactions.csv/json、2D SVG、scene.pse、ligand.sdf）；
- `executed: false`：不替下游伪造相互作用结果。

## 5. 阻塞码（写进 `task.json.execution.blockers`）

| code | 含义 |
| --- | --- |
| `blocked_missing_membrane_capability` | 现有 MD 运行器没有膜/非水协议，膜体系未构建 |
| `blocked_platform_linux_required` | 现有 MD 运行器仅 Linux 可用 |

## 6. 交接流程

1. `zvss run ...` 产出任务包与 `adapter04.json` / `adapter05.json`；
2. 在具备 Linux + OpenMM 的机器上，把包目录拷到 04 项目可访问位置，按 `adapter04.json.command_readable` 执行；
3. 需要相互作用分析时按 `adapter05.json.command_readable` 执行；
4. 任何一步没有实际执行，都保持 `executed: false` 与阻塞项，不做“已完成”表述。
