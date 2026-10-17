# 整理后验证记录

日期：2026-10-17；本机 macOS / Apple Silicon。以下是本次整理副本的实际检查，不是原目录历史数字的转抄。

| 检查范围 | 本次结果 |
|---|---|
| 全部项目 Python 源文件 AST 解析 | 302 个通过（按 UTF-8 BOM 兼容读取） |
| 分子生成 tests | 57 passed |
| 主对接 tests/v2 | 29 passed, 2 skipped |
| 工作台 test_workflow_regressions | 15 passed |
| 虚拟筛选非 integration 测试 | 24 passed, 3 skipped, 3 deselected |
| ZLIS tests | 40 passed, 1 warning |
| OpenMM inventory / chemistry / workflow / mapping / backend 工程契约 | 76 passed |
| 独立结构 prepare 入口 | 成功；使用迁移后的集中结构路径 |
| 本机 HTTP 首页、结构列表、新建运行 | 首页 200；4 个结构；新运行所有阶段 pending |
| 启动 shell 语法 | 通过 |
| 常见令牌/私钥/凭据赋值模式扫描 | 无匹配；不等同于全面安全审计 |

合计 241 项测试通过，5 项跳过，3 项集成测试未选。ZLIS 有一条 NumPy 二进制兼容性 warning，未影响本次测试通过，但目标环境仍应按依赖锁验证。

测试 Python 为已有 zlis 环境 Python 3.11，RDKit 2026.09.1；筛选及 OpenMM 工程测试临时增加 PyYAML 6.0.3 和 gemmi 0.7.5。未修改既有环境。各模块在各自目录独立运行，以避免同名包混用。

没有运行：Peptide2Mol/REINVENT4 正式推理、目标 GNINA 批量对接、GROMACS 膜模拟、OpenMM 生产模拟、Windows 启动、第二台机器迁移。未把这些列为通过。模型适配器测试包含子进程替身，工程通过不构成科学结论。

原始审查说明保留于 `history/integration_audit_2026-10-17.txt`，其中历史文件数、结果位置与测试范围不代表本次发布副本。
