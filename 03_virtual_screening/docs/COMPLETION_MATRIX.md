# 功能完成矩阵

图例：**DONE** = 已实现；**VERIFIED** = 已用真实数据（非示例）验证；**BLOCKED** = 已实现但被本机外部条件阻塞；**NOTDONE** = 未实现。

## 本轮验收（七项要求逐条状态）

| # | 要求 | 状态 | 证据 / 阻塞条件 |
| --- | --- | --- | --- |
| 1 | MurJ 受体审计 | **DONE + VERIFIED** | `runs/murj2/reference/audit/`：沉积 4824 原子（A 链 4598 = 天然域 3–500 共 3785 + N 端融合 −101–2 共 813；B 链 29 残基 Sgl 肽）；受体裁剪后 **3785 原子 / 498 残基 / auth 3–500**；删除 104 残基见 `deletions.csv`；原生复合体留档 `complex_full.pdb`；WT 侧链恢复 `not_evaluated`（未做 rotamer 重建，保留沉积突变 K5F/S12I/M13A）；编号往返/仅天然域/完整覆盖/参考留档四项 pass；原子数是测量值，不作为建模后的固定目标 |
| 2 | 候选准入审计 | **DONE + VERIFIED** | `select/qualify.py` + `reports/qualification.csv`；正式主候选与展示性初步候选分离；`runs/gq5`：正式主候选 **0**、展示性初步 2、探索 1，未满足项逐条列出（无局部区域、无复对接、无膜来源）；`md` 对未合格候选不下发任务包（`blocked_no_formally_qualified_candidates`） |
| 3 | PoseBusters 修复 | **DONE + VERIFIED** | 根因：`tempfile.gettempdir()` 在 TEMP 不可写时回退到 CJK 项目路径 → RDKit `BaseFeatures.fdef` 窄字符打开失败；修复 `compat.scratch_root()`（探测平台 TEMP 及会话子目录，CLI/Web 共用）；真实位姿 22/22 通过，平移进蛋白失败 `minimum_distance_to_protein`+`volume_overlap_with_protein`，拉长键失败 `bond_lengths`+`internal_energy`，原子重叠失败 3 项（`docs/acceptance/POSEBUSTERS_RUN.md`）；异常全量落盘；缓存键加入能力指纹，`runs/gq2` 旧 `unsupported` 记录已重算 |
| 4 | 真实 GNINA 执行环境 | **BLOCKED（CUDA 运行时）** | 唯一 Linux 入口 WSL2 `Ubuntu-22.04`、GPU `RTX 5070`；按上游 vendor lock 下载并校验 GNINA 1.3.3（sha256 `3340c1f4…`）；补齐 cuDNN 9.8.0.87 + `libcudart.so.12` + `libcublas.so.12` 后仍缺 `libcusparse.so.12`（`docs/acceptance/GNINA_ENV.md`）。适配器已实现多 seed、固定受体坐标系 `CalcRMS` 聚类、planned/attempted/succeeded/supporting 四类计数，但无真实运行 → 未验证 |
| 5 | 七类指标真实执行证据 | **部分 VERIFIED** | 第 1 类 QA（含 PoseBusters 22 项）已验证；第 2/3/5/6 类真实执行（`runs/gq5`）；第 4 类被 CUDA 阻塞；第 6 类仅有几何代理（`proxy_only`）→ `not_evaluated`；异常应变判 `not_evaluated` + `strain_not_converged`，不再记为 0 应变通过 |
| 6 | MurJ 试筛（约 20 分子） | **BLOCKED（缺局部对接盒）** | 受体已冻结（3785 原子/498 残基）；MurJ 无可核验局部对接盒来源（9NU4 无共结晶配体，B 链 29 残基 Sgl 肽不足以定义位点），故不生成任何 MurJ 对接或位姿级指标；候选为空、解释保留、规则未放宽。缺失输入：经验证位点定义（共结晶配体 / 定点突变位点 / 上游冻结盒） |
| 7 | MD 包与规格等价性 | **DONE + VERIFIED（交接包层面）** | `tools/validate_md_handoff.py`：`runs/handoff1` 2/2 包通过全部 8 项（文件、checksum、原子映射、复合体组成、坐标与 `ligand.sdf` 一致 <0.001 Å、来源、protocol+units、阻塞项显式）；`runs/t6` 历史包仅缺 `protocol/units`（本轮已修复 `md/package.py`，历史包保留为缺项证据）；膜引擎预检仍 `blocked_missing_membrane_capability` + `blocked_platform_linux_required`，未在 06 另建 MD 引擎 |

| # | 规格要求 | 状态 | 证据 / 说明 |
| --- | --- | --- | --- |
| 1 | 消费既有 GNINA 批量结果（SDF/CSV/log） | DONE + VERIFIED | `runs/t4`：313 文件、60/60 位姿导入、0 失败；盒子从 GNINA 日志 `argv` 解析 |
| 2 | 输入发现与分类、`._*`/未知文件忽略 | DONE + VERIFIED | `ingest/discover.py`；`input.sdf`/`ligand.sdf` 不再被当作位姿，`poses.sdf` 与 `pose_00N.sdf` 去重 |
| 3 | 分数与位姿绑定（禁止按顺序猜） | DONE + VERIFIED | `ingest/poses.py`；日志 mode 与 SDF 记录按分数值配对，失败则 `table_only` + reason |
| 4 | 分数单位语义（CNNaffinity 是 pK） | DONE + VERIFIED | `SCORE_SEMANTICS` + 单测 `test_score_semantics.py`；报告不把 pK 折算成 kcal/mol |
| 5 | 四层身份 + 确定性 ID | DONE + VERIFIED | `ids.py` + `test_ids.py`；两轮运行产生相同 `selection_version` |
| 6 | 位姿质量校验（PoseBusters + 自实现） | DONE + VERIFIED（本轮修复） | PoseBusters 0.6.5 已真实执行：60/60 位姿各 22 项检查、受体参与蛋白相关项；`runs/gq5`；ASCII 临时目录可移植（CLI/Web 共用） |
| 7 | MurJ 局部机制特征（残基映射→接触） | DONE + VERIFIED | `runs/murj1`：498 残基映射、TM2/TM7 区域解析；`runs/t4`：60 位姿接触 Jaccard（无参考配体时记录 `not_evaluated`） |
| 8 | 构象干扰 / 位姿一致性聚类 | DONE + VERIFIED | `motion/interference.py`，单链聚类阈值来自配置，t4 中每条记录都有 `pose_cluster_id` |
| 9 | 多 seed 再对接 | DONE but BLOCKED | `redock/gnina.py` 能力探测 + `blocked_external_tool`；本机无 `gnina`，成功分支未运行 |
| 10 | 分数整合 + 证据分层 + 多目标排序 | DONE + VERIFIED | t4：B=18、C=2、D=0；缺失分量在 `components_missing` 与报告中列出 |
| 11 | 膜环境评估 | DONE（仅代理） | 跨膜螺旋主轴几何代理 `proxy_only`；不声称知道真实膜内取向 |
| 12 | 分子性质 / PAINS / 骨架 / 多样性 | DONE + VERIFIED | `strain/strain.py` + `chem/props.py`，t4 全部 20 个化合物有性质与骨架分组 |
| 13 | 配体应变 | DONE + VERIFIED | MMFF94（UFF 回退）；非物理单点能量会 `not_evaluated` 而不是伪造数值 |
| 14 | 预算约束（≤200 再对接 / ≤50 深度 / ≤12 MD / ≤2 位姿） | DONE + VERIFIED | 来自 `budgets/*.yaml`；t4 产出 3 个任务包，未超限 |
| 15 | 候选选择（层内多目标 + 多样性 + 骨架上限） | DONE + VERIFIED | t4：2 primary + 1 exploratory，`SEL6830370C0C42E694` |
| 16 | 排除报告与理由码 | DONE + VERIFIED | `reports/exclusions.csv`；`Reason` 枚举覆盖化学/几何/分数/预算/证据不足 |
| 17 | 报告（人读 + 机读） | DONE + VERIFIED | `report.md` / `candidates.csv` / `ranking.csv` / `exclusions.csv` / `report.json` |
| 18 | MD 任务包（complex.pdb/ligand.sdf/manifest/task/README/checksums） | DONE + VERIFIED | t4/t5 共 3 个包；`checksums.sha256` 可校验；`task.json` 带真实阻塞码 |
| 19 | 与 `05_molecular_dynamics`（OpenMM）适配 | DONE but BLOCKED | `adapter04.json` 记录命令模板与两条阻塞（Linux-only、无膜协议）；`executed: false` |
| 20 | 与 `04_interaction_analysis`（ZLIS/ProLIF）适配 | DONE（未执行） | `adapter05.json` 记录 `run.sh analyze` 命令与选择器；不在下游环境执行 |
| 21 | SIFTS 编号映射（PDB↔UniProt） | DONE + VERIFIED | 真实 9NU4：`sifts_pdbx_sifts_xref_db`，498 映射、3 处差异、缺失区正确 |
| 22 | 突变 / 标签 / 融合区识别 | DONE + VERIFIED | K5F/S12I/M13A 与 C 端 `SAWSHPQFEK`（auth 501–510）由确证法与 SIFTS diff 一致得出 |
| 23 | CLI（doctor/run/serve/runs/report） | DONE + VERIFIED | 全部在真实数据上跑过 |
| 24 | 本机 Web UI（127.0.0.1 只读） | DONE + VERIFIED | 8770 端口端点冒烟通过（health/runs/首页/运行页/报告页/404） |
| 25 | 缓存与断点续跑 | DONE + VERIFIED | 缓存键 `(stage, input_hash, config_hash, code_version)`；`--only` 重跑不丢其它阶段 |
| 26 | 测试集 | DONE + VERIFIED | 29 passed（26 单元 + 3 端到端） |
| 27 | 文档 | DONE | README / ARCHITECTURE / MD_INTERFACE / DEPENDENCIES / TEST_REPORT / 本矩阵 |
| 28 | 膜体系构建 + 膜 MD | NOTDONE + BLOCKED | 本机与既有 04 项目均无膜协议；只做代理评估并如实标 `blocked_missing_membrane_capability` |
| 29 | 真实 MurJ 对接位姿筛选 | NOTDONE（无输入） | 工作区不存在 MurJ 对接位姿；位姿级结论为 `not_evaluated`，不做替代性结论 |
| 30 | 在 Linux 上实跑 OpenMM MD | NOTDONE + BLOCKED | 需 Linux + OpenMM 环境；本工作区为 Windows |

## 结论口径

- 可以对外说的：**“已用真实 GNINA 输出（60 位姿 / 20 化合物）跑通全流程，产出可解释候选、排除报告与 3 个 MD 任务包，且两次独立运行结果逐字节一致；MurJ 侧完成真实结构制备与编号映射验证。”**
- 不可以说的：**“完成了膜体系 MD”“完成了 MurJ 分子筛选”“再对接验证通过”**——这三件事在本环境被外部条件阻塞或缺少输入，报告中均以阻塞码/`not_evaluated` 呈现。
