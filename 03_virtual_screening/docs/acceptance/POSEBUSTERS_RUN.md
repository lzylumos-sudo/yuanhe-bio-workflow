# PoseBusters 真实执行记录

- PoseBusters 版本: `0.6.5`
- ASCII 运行时临时目录: `C:\Users\user\AppData\Local\Temp\dsh-iKCsDM\zvss-scratch`
- 真实位姿来源: `F:\MemPilot\runs\public20_engineering\entries\P000005_S01_C01_seed20261015\input.sdf`
- 受体: `F:\MemPilot\runs\public20_engineering\protocol_inputs\receptor_prepared.pdb`

| 用例 | 原子数 | 检查项 | 通过 | 失败项 |
| --- | --- | --- | --- | --- |
| real_pose | 12 | 22 | 22 | 无 |
| translated_8A_into_protein | 12 | 22 | 20 | minimum_distance_to_protein, volume_overlap_with_protein |
| stretched_bond_2.2A | 12 | 22 | 20 | bond_lengths, internal_energy |
| all_atoms_overlapped | 12 | 21 | 18 | bond_lengths, bond_angles, protein-ligand_maximum_distance |

## 判定

- 真实位姿失败项数: 0（应为 0）
- 平移进蛋白位姿失败项: minimum_distance_to_protein, volume_overlap_with_protein
- 拉长键位姿失败项: bond_lengths, internal_energy
- 原子重叠位姿失败项: bond_lengths, bond_angles, protein-ligand_maximum_distance

失败项逐条来自 PoseBusters 的布尔检查列，未经任何后处理或裁剪。
