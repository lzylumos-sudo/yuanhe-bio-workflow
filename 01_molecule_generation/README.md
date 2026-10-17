# 分子生成与准备

主实现位于 `src/`，独立 CLI 是 `murj_pipeline.py`。从本目录执行：

```bash
python murj_pipeline.py --run-id demo_prepare prepare
python murj_pipeline.py --run-id demo_prepare status
PYTHONPATH=. python -m pytest tests -q
```

准备读取 `../data/structures/9nu4.pdb`。模型与权重须恢复到 `vendor/`，见根目录 `docs/01_SETUP.md`。准备阶段不会调用生成模型；真实生成入口不会回退为旧规则枚举。

`README.historical.md` 保留旧开发记录，其中环境位置与状态描述可能过时。当前入口和验证范围以根目录文档为准。
