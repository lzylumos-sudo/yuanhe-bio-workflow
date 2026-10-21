# 双 RTX PRO 5000 Blackwell 研究配置

本次交付仅完成代码修改和静态验收。没有运行单元测试、集成测试、模型推理、Vina、GROMACS、checkpoint 实测、GUI 验收或性能测试。不得据此声称双卡环境已验证或科学结果有效。

## 统一配置与默认值

`00_workbench/configs/research_dual_blackwell.json` 是正式定义；`00_workbench/research/config.py` 是 UI/API/CLI 共用校验器。独立 `00_workbench/configs/smoke.json` 保持 MD 为一条 500 ns 生产任务，smoke 不会自动启动 MD。

|字段|默认|
|---|---:|
|Peptide2Mol 原始尝试|5000|
|审核核心上限|200|
|REINVENT4 全任务尝试|200000|
|生成唯一母体上限|50000|
|外部母体上限|50000|
|初筛母体上限|100000|
|精筛母体上限|2000|
|详细审核上限|50|
|MD 候选上限|8|
|每候选生产轨迹|1 × 500 ns|
|MD GPU 并发|2|
|单 MD 作业阶段超时上限|60 天（执行保护上限，不是 ETA）|
|湿实验目标|4–8|
|Vina worker / 线程 / 总预算|24 / 2 / 48|
|初筛构象 / 状态 / 种子 / exhaustiveness / modes|1 / 1 / 1 / 8 / 5|
|精筛构象 / 状态 / 种子 / exhaustiveness / modes|2 / 2 / 2 / 32 / 10|

全部为上限或目标。精筛实际仅有一个受体或化学状态时按一个计算，保留缺口；不得复制结构凑数。额外状态必须通过 `ligand_states` 提供审核、文件哈希和独立前置筛选证据。

## 实际调用链

工作台创建 → `workflow.create` → `research.config.resolve` → 冻结配置/哈希/环境/资产清单。

Peptide2Mol → 核心人工审核 → REINVENT4 全任务配额 → `research.qualification.prepare_risk` 标准化母体 → 共享 `yuanhe_preflight` 风险筛查 → 供应证据/本地 AiZynthFinder 路线审核 → `research.docking.run` 两级官方 Vina → 各受体独立 Y-MAVS → 结构审核 → 既有膜构建 → `research.md_queue` → 既有独立 GROMACS/分析入口 → 湿实验人工审核。

默认规则为已有前置模块的可追溯规则子集，不是完整毒性预测。没有供应目录、模型、模板、原料库存或合格路线时保持 pending。SA/QED 不产生可获取标签。没有合格分子时不会启动 Vina。

MD 的 `build` 包装器仅绑定真实对接受体；`build_prepared`、力场、电荷、膜模型、MDP、run、checkpoint_probe、analyze 保持基线。现有独立 20 ps checkpoint 自检属于后端恢复检查，不是新增 MD 筛选阶段。

## 复用、恢复与证据

配置不可原地改变。`/api/research/run/<id>/derive` 只接受 Top N、MD 数量、湿实验数量与选择模式变化。生成及前置证据先验证再复用；Vina 复用键含输入 SDF、化学状态、受体、搜索盒、种子、参数、软件与适配器哈希。所有缓存输出重新核验。只为缓存未命中的组合调用 Vina。

已完成 MD 通过原任务文件及哈希引用复用，不复制/改写原始轨迹。未完成 MD 在同任务内通过原 checkpoint 恢复。生成批次完成后可按哈希复用；若中断时原始尝试量无法确定，保留预算预留并阻断重采样，需人工审计，不能无限重试消耗新预算。

GPU 用 UUID 独占锁；每个 MD 子进程仅看见一张卡。作业不会为了优先级强杀别的计算。资源不足缩减并发或阻断。ETA 无实测时为 null。MD 日志中的实际生产时间、Performance 和 checkpoint 状态显示在研究面板。

## 双卡机器上的验收命令（本次未执行）

在已配置的化学环境中从仓库根目录运行：

```bash
python -m pytest 00_workbench/tests/test_research_dual_blackwell.py 02_docking/tests/preflight -q
python 00_workbench/run_research.py preview --config 00_workbench/configs/smoke.json
python deployment/preflight_doctor.py
python 00_workbench/run_webui.py --help
```

任务创建命令接收整个工作台配置：`python 00_workbench/run_research.py create --config <完整配置.json>`。仅预览不创建任务；创建也不会执行计算。开始阶段需要显式 `start --run-id <ID> --stage <阶段>`。请先填真实执行环境、已审核核心、受体、搜索盒与获取证据。

## 尚未完成或需目标环境验收

- 本次没有任何运行验收结果或实测速率；Blackwell CUDA、模型、GROMACS/RESP 和资源识别必须在目标机验证。
- 初筛采用 Vina 分数与真实骨架分组的多样性轮选；精筛后保留 Y-MAVS。初筛尚未计算全部候选的构象应变/关键接触多目标排序。
- REINVENT4 当前按核心均分总配额、依次执行；未实现根据有效率自适应分配或双 GPU 核心并行。
- CPU worker 按配额限流，未实现基于实测吞吐量的在线自适应调参；队列不实现跨节点集群或自动抢占。
- 不完整生成批次不自动续采；保留预算与原始输出等待审计。失败 Vina 缓存也保留，需检查后显式处理。
- 各阶段仍保留部分候选元数据列表；姿势文件按作业保存，但尚未把所有历史 JSON/CSV 读取改为数据库分页。
- 交付 ZIP 包含结构、审核与轨迹原件哈希清单；大体积轨迹保持原目录，需要按清单一并交付。
- 化学状态需要显式导入并审核，不自动预测额外质子化/互变异构状态。
- 已提供新测试代码；没有执行，也没有将既有另一任务的测试结果计为本次整合版本通过。

这是一版接入原项目的研究流程升级，仍须上述目标环境验收和后续能力补齐，不能视为所有 P1/P2 均已完成。
