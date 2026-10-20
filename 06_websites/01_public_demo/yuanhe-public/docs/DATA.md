# 公开数据与后续接入

本版只有下列公开样例：

| 文件 | 来源 / 依据 | 含义 |
|---|---|---|
| public/data/9NU4.pdb | https://files.rcsb.org/download/9NU4.pdb | 原始公开坐标，2026-10-18 独立下载 |
| public/data/sample.json | scripts/prepare_data.py 从上项复算 | SHA-256、来源、方法、链摘要及最小接触距离 |
| public/structure.svg | 同一脚本从原始 Cα 坐标投影 | 仅用于轻量首页预览，非任意装饰蛋白 |

PDB 数据依据 https://www.wwpdb.org/about/usage-policies 使用 CC0。结构引用：https://doi.org/10.2210/pdb9NU4/pdb；原论文：https://doi.org/10.1038/s41586-026-10163-w。

方法：第一模型，ATOM 重原子，altloc 空或 A；受体 A:3–500；B 链可见片段；每个受体残基对所有 B 重原子取最小欧氏距离，保留 ≤8 Å，四位小数存储。当前 PDB 未见重复 altloc；换结构时必须重新审核选模与 altloc 策略。原始 A 链记录包括构建体额外序列（编号可为负），不要将全链残基数等同于本演示选区。

```sh
python3 scripts/prepare_data.py
pnpm test
pnpm build
```

接触是坐标派生量，不是氢键判定、自由能、活性或药效证明。原文件包含工程突变，B 只有 2–30 可见残基，不把未解析片段补成已观测结构。

## 后续内部样例审批记录（必须先得到用户明确确认）

每个拟公开文件应单独登记：样例 ID、靶点、原始输出与关联分子、运行程序/版本、参数/种子、SHA-256、来源许可、审核人、明确批准的公开范围与日期、排除的敏感字段。没有批准记录的文件不复制到 public/。

建议扩展 sample schema 为具有 `status: approved | unavailable`、`provenance`、`license`、`approvedAt` 的分阶段对象。候选需稳定 ID、SMILES、SDF、真实分数及单位，MD 需轨迹与时间单位和分析方法。当前 results 字段均为 null，页面未实现任意导入内部结果的隐式通道。接入新数据时同时编写展示组件和一致性测试，再发布新版本。
