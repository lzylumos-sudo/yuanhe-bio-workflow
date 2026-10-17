# 架构

## 1. 分层

```
CLI (cli.py)  ──►  stages.run_all(ctx)  ──► 各阶段模块 ──► SQLite + runs/<id>/stages/<stage>/summary.json
Web (api.py)  ──►  只读：DB + 已生成产物/报告
```

- **输入层** `ingest/`：`discover.py` 只做文件发现与分类（位姿 / 表格 / GNINA 日志 / 受体 / 模板 / 复合物 / 未知），`poses.py` 负责 SDF/PDB 位姿解析、分数规范化与“分数↔位姿”绑定，`import_run.py` 负责身份生成与入库。
- **结构层** `receptor/`：`pdbio.py` 做受体清空（水/离子/脂/配体/其它链）与几何检查；`sifts.py` 处理 PDB↔UniProt 编号（SIFTS → `_struct_ref_seq` 偏移 → 全局比对三级回退）。
- **计算层**：`qa/pose_qa.py`、`mech/interactions.py`、`motion/interference.py`、`redock/gnina.py`、`strain/strain.py`、`membrane/environment.py`、`scoring/aggregate.py`、`select/pipeline.py`。
- **交付层**：`md/package.py`（任务包）、`md/adapters.py`（04/05 适配器）、`report/render.py`（人类可读 + 机器可读报告）。

每个阶段写 `summary.json`（结构化、可回溯）并可写 `error.json`；`report` 阶段会合并磁盘上的所有阶段摘要，因此**单独重跑某阶段不会让报告丢失其它阶段的事实**。

## 2. 四层身份

| 标识 | 生成方式 | 为什么 |
| --- | --- | --- |
| `compound_id` (CMP) | canonical isomeric SMILES 的 sha256 前 16 位 | 同一分子在不同名字/文件里必须同 id |
| `chemical_state_id` (STA) | SMILES + 形式电荷 + 质子化说明 | 同一分子不同电荷态不是同一个可比较实体 |
| `docking_job_id` (JOB) | 受体版本 + 盒子 + 协议摘要 | 不同协议的结果不可混算 |
| `pose_id` (POS) | job + 记录序号 + 坐标摘要 + rank | 位姿级唯一，且与文件顺序无关 |

坐标摘要按 1e-3 Å 取整，避免浮点格式噪声产生新位姿。

## 3. 三种状态族（互不混用）

- `execution_status`：`pending/running/completed/failed/cached/skipped`
- `assessment_status`：`pass/flag/fail/needs_input/not_evaluated/unsupported`
- `selection_status`：`eligible/excluded/exploratory/selected/not_selected`

规则：**没有真正算过的检查一律 `not_evaluated` / `needs_input` / `unsupported`，永不写 `pass`。**

## 4. 七个筛选计算（规格 §）与实现位置

| # | 计算 | 位置 | 关键诚实点 |
| --- | --- | --- | --- |
| 1 | 化学身份 + 位姿质量 | `chem/`, `qa/` | PoseBusters 不可用时写 `unsupported` 并说明原因；盒子包含性、蛋白内部冲突、键长离群为自实现确定性检查 |
| 2 | MurJ 局部机制特征 | `mech/interactions.py` | 参考配体缺失（MurJ 的 Sgl 是蛋白/肽）时接触 Jaccard = `not_evaluated`，不编造 |
| 3 | 构象干扰 / 位姿一致性 | `motion/interference.py` | 单链聚类阈值来自配置，聚类按位姿 id 排序，结果稳定 |
| 4 | 多 seed 再对接 + 位姿聚类 | `redock/gnina.py` | 无 GNINA 时 `blocked_external_tool`，不产生“再对接证据” |
| 5 | 分数整合 | `scoring/aggregate.py` | 权重只在存在的分量上归一；缺失分量在报告与 `components_missing` 中列出 |
| 6 | 膜环境 + 分子性质 | `membrane/`, `strain/` | 膜只给跨膜螺旋主轴几何代理（`proxy_only`），明确“未构建膜体系” |
| 7 | 应变 + 多样性 | `strain/` | 力场不可参数化或能量非物理时 `not_evaluated`，绝不当作 0 应变 |

## 5. 证据分层与排序

先分层（A/B/C/D），再在层内用多目标分排序：

- 层级由“必须存在的证据”决定（如 A 需要局部机制证据 + 可靠分数 + 通过质检），阈值与要求来自 `budgets/*.yaml` 的 `tier_requirements`。
- 层内分数 = 各分量在候选集合内的百分位加权和；缺分量的权重被重新归一，并记录 `components_missing`。
- 平局用 `compound_id` 兜底，保证顺序确定。
- Pareto 前沿按 `pareto_objectives` 在可用分量上计算，用于报告“哪些候选是权衡上的非支配解”。

**同一 parent 的多个位姿不会当成多个命中**：排序在 parent 级别做，每个 parent 只取一个代表位姿（`representative_pose`：最优 affinity → 最优 CNN → pose_id 兜底），其余位姿只贡献“一致性/替代位姿”信息。

## 6. 预算与选择

`budgets/*.yaml` 约束：`redock_max_parents ≤ 200`、`deep_max_parents ≤ 50`、`md_max_parents ≤ 12`（primary + exploratory）、`poses_per_parent_max = 2`。选择阶段用贪心 max-min 多样性在 Tier 合格集合上挑选，同时满足骨架上限、相似度上限、膜告警负担上限；每个排除都带 reason code 与细节。

## 7. 确定性与可复现

- 所有遍历按路径/ID 排序；不依赖 `hash()`（字符串哈希随机化会破坏复现）。
- 阈值全部来自 YAML，并在 `resolved_config.yaml` 中回显来源。
- 缓存键 = `(stage, input_hash, config_hash, code_version)`；`--resume` 只复用完全匹配的缓存。
- 同一配置 + 同一输入 + 同一代码版本 ⇒ 相同 `selection_version` 与相同的候选/排名文件（已在 `runs/t4` 与 `runs/t5` 上逐字节验证）。

## 8. Web 层

FastAPI + Jinja2，仅绑定 `127.0.0.1`，只读：运行列表、阶段状态、候选、排名、排除、原始 `report.md`。没有写接口，页面不会触发计算或修改数据库。
