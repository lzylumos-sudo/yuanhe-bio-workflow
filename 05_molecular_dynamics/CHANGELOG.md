# 2026-10-16 — v2 通用输入和任务管理

- `run_md.sh`、`scripts/run_workflow.py`：统一输入、决定、准备、建系、OpenMM、分析与验收。旧入口和计算模块保留。
- `modules/automation/input_resolver.py`：指定输入优先、限制目录发现、任意名称、无人值守选择阻塞。
- `pdb_inventory.py`：固定列无损原子清单、稳定 ID、组分/连接/缺失/质子化变体识别、逐次规范化与双向映射。
- `ligand_chemistry.py`：本地模板核验、有超时的 RDKit 推断、原始电荷位点/氢/手性/坐标核验、精确候选确认及缓存。
- `backend.py`、`forcefields.py`：保留 HIS 质子化证据、显式末端处理、参数与配置验证、原模拟后端适配、真实 checkpoint 与轨迹验收。
- `analysis_mapping.py`：原始原子/残基至最终拓扑与分析映射，重建氢不冒充原原子。
- `state.py`、`workflow.py`：指纹、运行锁、科学决定、分阶段日志、完整性检查、恢复和状态。
- `setup.sh`、`scripts/check_environment.py`：安装不再只检查 Python 存在；验证库/工具/极小 CPU 与 checkpoint。显式包锁补正 `@EXPLICIT` 格式，保留原 338 个包 URL。
- `scripts/validate_v2_ubuntu.py`：真实 Ubuntu 短分子链路、中断恢复和幂等验收；本轮未执行 native MD。
- 新增自动化回归测试、原样 R2 夹具、中文操作说明和短测配置。

已实际复现：R2 空配体字段、重复名称、HIE 超出原输入协议。新版保留 3094 原子，配体 92 原子/98 边/+2 电荷位点，候选需确认；HXT 未决定不改写。还保留四个由现有环氢识别的 HID。包锁旧首行格式错误已修复。相关测试记录见 `V2_VALIDATION.json`。

静态已接通、尚待 Ubuntu 验收：新准备适配器、强制 AmberTools AM1-BCC、建系、完整模拟、实际中断恢复、最终拓扑映射与分析。原 v1 历史 PASS 不自动升级成新版 PASS；Mac 缺少 Linux 运行能力不是目标流程本身失败。
