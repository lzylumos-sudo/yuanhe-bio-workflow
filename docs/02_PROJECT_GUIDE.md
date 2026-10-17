# 项目架构与代码阅读指南

## 项目定位与贡献

项目组合现有生成模型、对接和模拟工具，贡献集中在阶段组织、输入输出契约、科学输入审核、失败恢复、结构与化学检查、结果可追溯性和研究界面。Peptide2Mol、REINVENT4、GNINA、RDKit、OpenMM、GROMACS 等不是本项目原创模型或引擎。

## 推荐阅读顺序

1. 根 README：问题、流程、目录和已知边界。
2. `00_workbench/webui/workflow.py`：阶段状态、上游哈希验证、设置冻结、重试归档。
3. `00_workbench/pipeline_bridge.py`：真实子进程适配与各阶段产物检查。
4. `01_molecule_generation/src/`：结构准备、核心审核、化学身份、模型来源记录。
5. `02_docking/src/sbvs2/docking/` 与 `03_virtual_screening/src/zvss/`：对接记录、位姿与分数绑定、筛选。
6. `04_interaction_analysis/src/zlis/`、`05_molecular_dynamics/modules/`：独立下游工具。
7. 各模块 tests 与本次验证记录：判断测试覆盖的具体范围。

## 实现状态与验证范围

当前实现包括统一界面、阶段编排、执行门禁与工程回归测试。真实模型采样、目标 GNINA 批量筛选、膜体系准备/模拟与跨机运行，必须有匹配版本的输入、命令、环境、原始结果及日志后才列为已验收。

工程测试用于验证软件行为；模型性能、候选分子质量与科学结论需要独立的计算和实验依据。数据来源、第三方依赖与授权说明见根目录 THIRD_PARTY_NOTICES.md。
