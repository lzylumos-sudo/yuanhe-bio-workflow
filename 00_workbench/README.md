# 统一研究工作台

从仓库根目录运行 `python3 00_workbench/run_webui.py --no-browser`。

`webui/workflow.py` 管理主流程；`pipeline_bridge.py` 调用真实阶段。`webui/static/index.html` 是主页面。历史页面及 pepwedge 仅保留旧接口，历史计算数据不打包。

请以根 README 和 `docs/01_SETUP.md` 为准。`README.historical.md` 是整理前文档，不代表当前版本验收。
