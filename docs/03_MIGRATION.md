# 目录与命名迁移

主流程使用两位数字加英文功能名；Python 包名、CLI 参数、科学产物文件名和第三方模型名保持原接口。改名同步修改桥接路径、启动脚本、结构选择与生成模块配置，避免只改外层文件夹导致断链。

| 原目录/入口 | 新目录/入口 |
|---|---|
| 元核生物_UI | 00_workbench |
| 元核生物_对接前分子生成/v2 | 01_molecule_generation |
| 01_GNINA_SBVS | 02_docking |
| 元核生物_虚拟筛选方案/06_Virtual_Screening_Studio | 03_virtual_screening |
| 05_Ligand_Interaction_Studio | 04_interaction_analysis |
| 04_Universal_MD_Simulation | 05_molecular_dynamics |
| 02_GNINA_SBVS_v2 | archive/docking_prototype |
| 根目录 PDB | data/structures/，文件名统一小写 |
| START_LINUX.sh | start_linux.sh |
| START_PROJECT.cmd / .ps1 | start_windows.cmd / .ps1 |

`FILE_MAP.csv` 保存复制文件与原始 SHA-256；改名后修改的源码哈希自然不同。旧根目录的开发提示词、临时诊断脚本、论文全文、重复 v1 生成副本、历史结果不纳入发布。工作台仍保留历史页面依赖的 pepwedge 包，但历史生成入口不是主流程模型后端。

本次整理在独立副本完成，原外置盘内容保持原状。旧运行记录的绝对路径及哈希不批量重写，避免伪造可恢复性；需要继续旧运行时使用原项目。
