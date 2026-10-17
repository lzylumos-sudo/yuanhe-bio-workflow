# Zhiyu Ligand Interaction Studio (ZLIS) 0.1.0

蛋白-小分子相互作用分析：原子级表格、可编辑 2D 图、PyMOL 场景、真实原子 PDB。
一套核心代码（`src/zlis`），Mac 与 Windows 只是两个启动入口。

## 启动
- macOS：双击 `Start_Mac.command`，或 `bash run.sh web`；命令行 `bash run.sh analyze inputs/3PTB.pdb --out results`
- Windows：`Start_Windows.cmd` / `run.ps1`（**未测试**，未在真实 Windows 上运行过）
- 首次安装环境：`bash env/setup_mac.sh`（conda 环境 `zlis` 必须建在内置磁盘，不能建在 exFAT 移动硬盘上）
- 路径含空格或中文已测试（见 `tests/test_pipeline.py::test_moved_project_dir_with_spaces_and_chinese`）。

## 命令
- `ligands <pdb>` 列出候选配体
- `analyze <pdb> [--ligand A:BEN1] [--template lig.sdf|SMILES] [--charge N] [--model N] [--offline] [--plip] [--no-png] [--layout layout.json]`
- `web [--port 8765]` 本地网页，仅监听 127.0.0.1，需要启动时打印的随机令牌；无账号体系，不要暴露到网络。

## 输出（每个配体一个文件夹）
`interactions.csv/.xlsx/.json`，`interaction_2d_{publication,detailed,audit}.svg/.png/.pdf`，`interaction_2d_editor.html`（拖动标签，保存 `layout.json`，再用 `--layout` 重绘），`interaction_2d_audit.html`，
`complex_analysis.pdb`，`ligand.sdf`，`scene.pse`，`scene.pml`，`rebuild_scene.py`，`scene_config.json`，`binding_site_main.png`，可选 `plip_crosscheck.json`。

## 设计约定
- 相互作用只由 ProLIF（`count=True`）加 ZLIS 自带几何（金属配位、水桥候选、位阻冲突）判定；PyMOL 与 2D 图只读 `interactions.json`，不重新判定。
- `atom_uid = M{model}:{chain}:{resname}{resseq}{icode}:{atomname}`，全流程一致；加上去的氢带 `:H{n}` 后缀。
- 配体坐标永远来自输入 PDB。SDF/SMILES 模板只提供键级与电荷。重原子最大位移 0.0 Å（每次运行写入 `protein_prep`）。
- 氢键箭头方向是供体指向受体，仅表示角色，不表示电子流。盐桥、疏水等没有箭头。
- 水桥只是候选（两段 O…极性原子距离都在 3.5 Å 内），因结晶水通常无氢，不指定供受体；报告距离取较长一段。无水输入不会产生水桥，也不会往 PDB 里添水。
- 不计算也不编造能量。
- 配体来源于 RCSB 化学组分字典默认质子化时标 `needs_confirmation`；`configs/ligand_templates.json` 中的 BEN 模板（苯甲脒阳离子 +1）是人工设定的假设。
- 蛋白氢由 OpenMM Modeller 在 pH 7.0 下添加；His/Asp/Glu/Lys 状态为模型假设，写入 notes。缺 OXT 的链端用临时占位原子加氢后删除，末端羧基氢与电荷未建模。

## 已知限制
见 `docs/TEST_REPORT.md`。
