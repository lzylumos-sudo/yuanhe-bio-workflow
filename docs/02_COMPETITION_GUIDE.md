# 比赛代码阅读指南

## 项目贡献应如何表述

项目组合现有生成模型、对接和模拟工具，贡献集中在阶段组织、输入输出契约、科学输入审核、失败恢复、结构与化学检查、结果可追溯性和研究界面。Peptide2Mol、REINVENT4、GNINA、RDKit、OpenMM、GROMACS 等不是本项目原创模型或引擎。

## 推荐阅读顺序

1. 根 README：问题、流程、目录和已知边界。
2. `00_workbench/webui/workflow.py`：阶段状态、上游哈希验证、设置冻结、重试归档。
3. `00_workbench/pipeline_bridge.py`：真实子进程适配与各阶段产物检查。
4. `01_molecule_generation/src/`：结构准备、核心审核、化学身份、模型来源记录。
5. `02_docking/src/sbvs2/docking/` 与 `03_virtual_screening/src/zvss/`：对接记录、位姿与分数绑定、筛选。
6. `04_interaction_analysis/src/zlis/`、`05_molecular_dynamics/modules/`：独立下游工具。
7. 各模块 tests 与本次验证记录：判断测试覆盖的具体范围。

## 可展示与未完成

可展示统一界面、代码结构、阶段门禁和已通过的工程测试。真实模型采样、目标 GNINA 批量筛选、膜体系准备/模拟与跨机运行，必须有匹配版本的输入、命令、环境、原始结果及日志后才列为已验收。

没有在本次整理中补造模型结果或科学结论。未给出具体比赛章程，因此不声称已满足其开源许可证、公开时间、数据授权、团队署名或 AI 辅助开发披露要求。
