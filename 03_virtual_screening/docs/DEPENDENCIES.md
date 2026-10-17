# 依赖

运行时（`pyproject.toml` 声明，`pip freeze` 实测版本）：

| 包 | 实测版本 | 用途 |
| --- | --- | --- |
| rdkit | 2026.9.1 | 化学身份、构象、力场、指纹、骨架 |
| numpy | 2.5.4 | 几何与统计 |
| scipy | 1.18.1 | 距离/统计辅助 |
| pandas | 3.0.6 | 表格读取与报告 |
| gemmi | 0.7.5 | mmCIF→PDB、结构解析 |
| MDAnalysis | 2.10.0 | 结构读取辅助（与下游 MD 生态一致） |
| prolif | 2.2.3 | 与 05 项目一致的相互作用引擎（可选） |
| posebusters | 0.6.5 | 位姿校验（可选，需 ASCII 临时目录） |
| fastapi | 0.143.3 | 本机 Web |
| uvicorn | 0.54.0 | ASGI 服务器 |
| pydantic | 2.14.1 | 配置/接口校验 |
| PyYAML | 6.0.3 | 配置解析 |
| Jinja2 | 3.1.6 | 页面模板 |
| matplotlib | 3.11.2 | 报告图（可选） |
| requests | 2.34.2 | 参考结构抓取（实际走 urllib 回退） |
| psutil | 7.2.2 | 主机能力探测 |

开发：`pytest 9.1.1`。

## 许可注意

- RDKit（BSD-3-Clause）、Gemmi（MPL-2.0）、MDAnalysis（GPL-2.0+）、ProLIF（Apache-2.0）、PoseBusters（MIT）、FastAPI（MIT）、Uvicorn（BSD-3-Clause）、PyYAML（MIT）。
- **MDAnalysis 为 GPL-2.0+**：只在本地进程内作为可选读取器使用，不随本项目再分发其代码；若将来以服务形式对外分发，需要重新审视许可证义务。
- GNINA 本体不在本环境内，也没有被再分发；ZVSS 只读取其**输出文件**。

## 安装注意

- 本机 `pip` 在 PEP 517 隔离构建下会因临时目录权限（`WinError 5`）失败，因此使用 `pip install --no-build-isolation -e .`。
- RDKit 的 C++ 层无法打开非 ASCII 路径，`zvss.compat` 会把 `rdkit/Data` 镜像到 ASCII 临时目录并设置 `RDConfig.RDDataDir`；若没有可用 ASCII 目录，相关可选检查标记为 `unsupported` 而不是失败。
