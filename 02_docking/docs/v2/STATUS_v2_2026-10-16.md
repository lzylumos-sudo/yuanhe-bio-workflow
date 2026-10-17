# 实现状态（诚实记录，2026-10-16）

判定口径沿用方案 §01：未装显示“未安装”，无靶点验证显示“需验证”，未测试不计入通过。

| 方案部分 | 状态 | 证据 |
|---|---|---|
| M0 旧项目基线（只读哈希） | 已完成 | reports/baseline/legacy_01_manifest.json（2259 文件） |
| RuntimePaths / 无硬编码路径 / 中文空格路径 | 已实现并测试 | tests/test_core.py::test_paths_chinese_space |
| 原子写入 + checksums + 篡改检测 | 已实现并测试 | test_checksums_detect_tamper |
| 状态机（审批不可跳过） | 已实现并测试 | test_state_machine |
| 审批账本、哈希绑定、入队硬门禁 | 已实现并测试 | test_approval_invalidated_by_hash_change, test_queue_gate |
| PDB 摄取：断链/缺原子/altloc/HETATM 分级 | 已实现；真实 1DRF 与旧记录一致（GLN12、FOL、249 水） | reports/smoke/1DRF_ingest.json |
| 缺失片段 **建模** 本身 | **未实现**：缺失残基只报告、不建模（长/短 loop、MODELLER 均未做） | — |
| 受体：补缺失原子、去水/杂原子（须显式确认）、PDBFixer 加氢、OpenMM 受约束最小化、漂移 QC | 已实现并在真实 1DRF 上测试：补 GLN12，1511 H，骨架漂移 0.27 Å。**不做 pKa 预测**，His/Asp/Glu 互变需人工审核；保留辅因子/金属明确拒绝 | tests/test_prepare_dock.py, projects/20261016_d888a546/01_preparation |
| 配体 CSV：逐行溯源、失败隔离、盐、立体 | 已实现并测试 | test_csv_traceability_and_failures |
| 状态生成 | 仅 **规则回退**（method=rdkit_rules_v0，酸/脂肪胺）；MolScrub 未安装，未与其对比 | reports/smoke/ligands_100_states.json |
| 3D 构象 + Gasteiger 电荷校验 | 已实现；100 母体 → 132 状态，0 电荷失败 | 同上 |
| 2D 审核 SVG + 批量规则不得越过异常 | 已实现并测试（无 Web 界面） | test_unknown_stereo_enumerated... |
| 盒子、协议哈希、argv 数组命令 | 已实现并测试 | test_box, test_protocol_hash_and_argv |
| 排序/聚合/EF/稳定性（不跨量纲加权） | 已实现并测试 | test_ranking_*, test_ef_and_na |
| 引擎/GPU 检测 | 已实现；本机 arm64 → 诚实报 unsupported_arch | `./run.sh doctor` |
| **Vina 实际对接** | 已执行（vina 1.2.7 arm64 python 绑定，Meeko 0.8.0 PDBQT）。每 job 独立目录、原始输出+哈希、失败隔离、目录不复用 | projects/20261016_d888a546/runs/…/02_work_products |
| 中断恢复 | 已实测一次：会话/硬盘中断后恢复，复核协议+输入哈希，11 个已验证 job 复用，1 个不完整 job 隔离到 `.quarantine_*`（未覆盖），其余重跑 | state_log.jsonl（RUNNING→INTERRUPTED→RUNNING） |
| **GNINA 对接 / CNN 重打分 / GPU** | **未执行**：无 arm64 版 gnina，本机无 NVIDIA GPU；适配器只生成命令（argv），未对真实二进制验证 | `./run.sh doctor` |
| Web 界面、SSH 隧道、WSL/Ubuntu 安装器 | **未实现** | — |
| PoseBusters、PDBQT 导出、报告生成、共价/柔性/金属/MD | **未实现** | — |
| FOL 重对接（Vina，3 构象×2 状态×3 seed=18 job，起始 3D 来自 SMILES） | 已做：**top1 RMSD 1.88–3.10 Å，仅 1/18 ≤2 Å**；任一 pose ≤2 Å 为 6/18，最小 0.93 Å。属 Vina 单靶点结果，**不是**富集验证，也不能代表 GNINA | run `redock_results.json`, `summary.json` |
| 靶点富集验证（EF/AUC） | 未做，无活性标签 → “尚无目标富集验证” | — |

## 已知注意
- 本机实测环境：Python 3.11 + RDKit 2026.03.4（借用另一项目的解释器，未安装到本项目内）。
- 酸/胺规则是粗略 pKa 窗口，**不能替代 MolScrub/Epik 类工具**，2D 审核即为此设。
- 盘为 exFAT：无符号链接、无 POSIX 权限；高频运行数据应放计算主机本地盘。
