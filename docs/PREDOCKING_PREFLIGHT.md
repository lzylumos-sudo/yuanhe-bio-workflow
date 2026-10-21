# 对接前风险与获取路径筛选

执行顺序：REINVENT4 → 结构标准化与精确去重 → `risk` → `availability` → AutoDock Vina → 原有筛选／MD。Peptide2Mol 种子在 risk 前并入输入，不能在对接入口追加。两阶段实际执行，不以 UI 标签或通过布尔值代替资格证明。默认未提供供应目录、模型和库存时，不会产生获取路径通过标签。

## 配置与入口

- 规则：`02_docking/src/yuanhe_preflight/configs/risk_rules.v1.json`。
- 配置模板：`02_docking/src/yuanhe_preflight/configs/predocking.example.json`。
- 工作流配置：任务 `config.json` 的 `predocking`；阶段设置保存在 `reviews/risk.json`、`reviews/availability.json`，执行时再冻结为阶段 `review.json` 和 `config.snapshot.json`。有效配置逐阶段合并后传给真实模块。
- 依赖规格：`deployment/preflight-environment.yml`，与 MD 环境分开。AiZynthFinder 适配器针对 4.4.1；实际资产未随代码提供，也不会自动下载。
- 许可清单：`deployment/preflight_licenses.json`；上游许可证及接口核对哈希在 `docs/preflight_sources/`。

从项目根目录运行 UI：

```sh
python 00_workbench/run_webui.py --no-browser
```

在“风险子结构过滤／可购买／可合成性筛选 → Stage settings”编辑 `predocking` JSON。`rules_path`、`cache_dir`、目录路径、模型／模板／库存路径和逆合成 Python 建议使用绝对路径。可配置规则动作、人工例外、目录有效期、搜索时间、最大深度、迭代数、工作进程和批大小。配置不包含 MD 参数。

外部供应分子库：risk 阶段设置 `imported_candidates` 为 SDF 或带表头 CSV/TSV/SMILES 文件。实际对接前仍需完成原有受体准备和审查。外部库允许跳过生成阶段，但不能跳过 risk／availability。

独立 CLI 只执行前置筛选，不运行生成、对接或 MD：

```sh
PYTHONPATH=02_docking/src python -m yuanhe_preflight.cli \
  --input /absolute/library.sdf \
  --config /absolute/predocking.json \
  --out /absolute/new_preflight_run
python deployment/preflight_doctor.py --config /absolute/predocking.json
```

CSV 必须包含 `smiles`，可含 `id`、`core`、`seed_id` 等来源字段。`.smi/.smiles/.tsv` 使用带表头的制表分隔格式。没有供货字段的普通 SMILES 列表仅是候选输入，不是采购证据。未指定的立体中心不会被解释为某一确定对映体；生成构象造成 SMILES/SDF 身份不一致时进入复核。

## 风险规则

默认是经过逐条正反例验证的小型精选规则集，非完整毒理预测器；不启用 `ALL`。

| 类别 | 默认规则／动作 |
|---|---|
| 毒性风险 | 项目定义 N-亚硝胺模式，默认排除；FDA 来源支持风险背景，不代表 FDA 提供此 SMARTS |
| 非特异性反应性 | Brenk 酰卤、异氰酸酯默认排除；过氧化物与项目改编的 NIH 环氧基警报默认复核 |
| 实验干扰 | PAINS catechol_A(92)，默认保留警告；不归类为毒性结论 |
| 其他药物化学警报 | Brenk 硝基，默认警告 |

每条规则包含 ID、SMARTS、类别、严重程度、动作、来源、版本、正例和反例。启动时检查 SMARTS 和 FilterCatalog，逐个执行正反例；无效规则导致任务失败，不静默忽略。匹配索引对应标准化 SMILES 的零基重原子编号，同时保存到原 SDF 原子编号的映射。SVG 显示索引及高亮。

标准化只做 sanitization、去除原子映射号、规范化异构 SMILES；不自动脱盐、中和、归并互变异构体或去除同位素。重复结构保留原始行与全部来源。排除是排除完整分子，原始输入和排除 SDF 仍保存。未命中时只显示“未命中当前规则”。

`risk_policies` 形如 `{"BRENK_NITRO":"review"}`，动作可为 `exclude/review/warn`。人工例外通过 `risk_exceptions` 记录 `rule_id`、规范化 `structure_smiles`、`rules_sha256`、`reviewer`、`reason`、`reviewed_at`；例外不能抹去原警报。risk 完成后需要补充复核，可点“修订前置证据”，先归档旧阶段，再编辑配置、重新运行前置阶段；不重新生成分子。已有对接或 MD 结果的任务不就地追改证据。

## 真实采购证据

`catalogs`（成品）和 `stock_catalogs`（逆合成末端原料）都是描述符数组，例如：

```json
{"path":"/absolute/vendor.csv","version":"vendor-2026-10-21","source":"supplier-authorized-export","license":"contract/reference","trusted":true}
```

支持 SDF 属性、CSV 或带表头 TSV/SMILES。每条有效供货记录需要：

- `smiles`、`vendor`、`sku`、`supply_type`（`in_stock` 或 `made_to_order`）。
- `catalog_date`、`checked_at`（ISO 日期），可选 `valid_until`。
- `evidence_kind` 为 `stock_confirmation`、`dated_quote` 或 `orderable_confirmation`。
- `evidence_file` 为相对目录文件或绝对文件路径；`evidence_sha256` 为该真实报价／库存确认文件的 SHA256。

目录与查询日期均受 `evidence_max_age_days` 限制（默认 30 天）。系统核对文件存在、内容哈希、结构和日期，不替代对供应商文件真实性的人工核实。单独的产品页 URL、CAS、PubChem、相似分子或旧目录都不构成通过证据。相同骨架、另一对映体、不同同位素和不同盐型不会误匹配。可选 `parent_smiles` 只记录盐型／母体映射，不据此授予目标盐型可购买资格。

供货待确认与过期状态单独显示；现货和按需合成可订购产生不同标签。正式采购前仍需再次确认。

## 本地逆合成与路线完整性

采购优先。无有效采购证据才启动 AiZynthFinder 子进程，实际加载 `model_path`、`templates_path`、`stock_path`。每个资产需要版本和独立许可字段。仅允许本地 template-based 模型和本地 inchiset 库，不接受任意远程模型配置。`external_interfaces` 当前必须为空：本版本没有第三方接口客户端，也不上传结构。

原料库存用于真实搜索；`stock_catalogs` 用于独立核验每一个末端原料的精确供货证据。仅有 AiZynthFinder `in_stock=true` 不够。

路线核对包括根节点精确结构、分子／反应树完整性、无循环、每步模板来源、RDKit 重放逆反应得到的精确前体组合、全部末端原料供货记录。保存完整树、步数、最大深度、模板和元数据、软件版本、原始搜索、资产 SHA256。对含立体化学潜在问题的路线采取保守复核；后端声明的未解决选择性／缺步骤不自动通过。SA／SCScore 不参与资格判断。

`route_reviews` 可包含 `route_sha256`、`target_smiles`、`reviewer`、`reason`、`reviewed_at`、`identity_stereo_selectivity_confirmed:true`。人工审核不能补造缺失的步骤或原料证据。可用 `search_artifact` 链接已有搜索结果，使已购买目标同时展示已有合成路线标签；不会因此强制重做搜索。

通过标签为“可买｜现货证据”“可买｜按需合成可订购”“可合成｜计算路线支持，待化学确认”“可合成｜路线已人工审核”，可以并存。计算路线可进入对接，但不会改变最终湿实验／实际委托合成所需的获取路径审核。

`model_missing`、`templates_missing`、`stock_missing`、`dependency_missing`、`search_timeout`、`no_route_found`、`route_needs_review` 分别记录。未找到路线不等于无法合成。进程超时预算是搜索时间加模型启动宽限。搜索并发上限为 workers，每批最多 batch_size；每个结果独立落盘。跨任务缓存键包含精确结构、规则哈希、模型／模板／原料库哈希与版本、目录版本与哈希、搜索参数、工作进程代码哈希。文件锁和原子落盘保证中断后只复用完整搜索结果；超时／错误不缓存为成功。复用路线时重新核验供货时效，过期证据不会永久放行。

## 对接门槛、身份与导出

`vina_official.run_vina_job` 首先调用 `require_ligand`。正式工作流、调度器、新任务及缓存恢复、远程 batch CLI 和二次对接均执行同一验证；无证据的历史 PDBQT-only Python 入口拒绝执行。旧 `sbvs` CLI/API 需配置 `qualification_manifest` 指向新的 `availability/qualified.json`，枚举后每个确切化学状态必须具有自己的证据。

验证不是读取 `decision=true`：重新核对结构身份、原 SDF 哈希、筛选记录哈希、当前规则匹配与未处理阻断项、目录／原料库版本、有效供货或完整路线、原始搜索和模型资产。改变结构、修改证据或证据到期都会拒绝。旧历史结果不因此获得新资格标签。

每条合格记录包含 compound_id、state_id、来源核心、生成批次、原生成 ID 和证据哈希。沿用原 zvss 的 compound/state ID 配方；pose_id 仍由原机制生成。证据进入 docking job、pose sidecar、筛选数据库、candidates.csv/json、`predocking_handoff.json` 和最终 `predocking_evidence.json`／报告。MD 源码和数值参数没有改动。

主要阶段产物：

- `risk/raw_input.json`、`records.json`、`rules.snapshot.json`、`inputs/*.sdf`、`passed/excluded/review.sdf`、分子高亮 SVG。
- `availability/records.json`、`qualified.json`、`identity_evidence.json`、`evidence/*.json`、`ligands/*.sdf` 及 `.qualification.json` 旁文件。
- `source_evidence/` 保存供应目录、供货文件和原始搜索快照；所有阶段有 CSV 和计数摘要。

UI 原有“Export reproduction records”会包含这些阶段文件及哈希清单。轻量导出仍沿用原有 20 MB 截断策略，遗漏文件列在 EXPORT_MANIFEST；需要完整目录／原始证据时选择 full 导出。模型权重不打包，但保留哈希和版本。跨主机需重新配置本地资产路径并运行前置验证，不能编辑通过布尔值复用旧资格。数量不足按实际合格数量继续；零合格时不启动 Vina；不自动追加生成。

## 验收范围

新增测试在 `02_docking/tests/preflight/test_preflight.py`。其中供应商、报价和路线均明确为模拟单元夹具；RDKit 匹配、模板重放、文件哈希和门槛判断实际执行。真实模型搜索、真实库存／报价和生产目标环境验收需要提供相应资产后执行。本次未启动正式筛选批次或 MD。
