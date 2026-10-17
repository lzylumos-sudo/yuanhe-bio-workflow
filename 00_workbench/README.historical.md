# WEDGE-Design：从噬菌体肽到 MurJ 小分子骨架

> MurJ（lipid II flippase）是革兰氏阴性菌细胞壁合成的必需转运蛋白。
> 三种独立进化、序列无同源性的噬菌体裂解蛋白（Sgl^M / Sgl^PP7 / Sgl^CJ3）
> 趋同地把它锁在**向外开放（periplasm-open）**构象。
> 本项目把这套"肽的抑制逻辑"翻译成可枚举、可解释、可合成的小分子**骨架**。

本仓库是该设计算法的**完整可运行原型**：输入三个冷冻电镜结构，输出排名化的分子骨架、3D 位姿与全部中间证据链。

---

## 快速开始

```bash
pip install numpy matplotlib

python run_evidence.py --out evidence    # 证据链：30 秒，把每个数字从 PDB 算一遍
python run_pipeline.py --out results     # 设计：50 秒
python run_filters.py --out results      # 结构衍生过滤器：538 -> 32（可选，81 秒）
python tests/test_smoke.py               # 自检（结构/先验/证据链/组装/枚举/口袋/端到端）
```

**先跑 `run_evidence.py`。** 它会把 24 条文献先验逐条对坐标核对，
输出 `evidence/report.md` 与 7 个数据文件。若 PASS 数不是 22，说明结构或先验被改动过，
此时任何下游设计结论都不可信。

约 **80 秒**（证据链 30 s + 设计 50 s，单核 CPU）后得到：

```
evidence/
  report.md                  # 【自动生成】证据链报告：输入结构、界面规模、重要残基、锚点几何、先验核对
  00_census.csv              # 链/残基/原子数、缺失区间、构建体残基（已剔除）+ 由此得到的实际靶点范围
  01_interface_murj.csv      # 每个界面残基：接触、dSASA、收敛度、CA 坐标
  02_interface_peptide.csv   # 肽侧每个残基：接触、埋藏面积
  03_anchors.csv             # 盐桥 + 锚点两两距离（逐结构）
  04_priors_check.csv        # 先验逐条核对：PASS / FAIL / FLAG / MISMATCH
  05_atom_budget.json        # 界面规模 + 肽 vs 小分子原子预算
  06_priors_resolved.json    # 先验文件解析后的形态（含每条核对规则的实际参数）
results/
  interface_hotspots.csv     # 52 个界面残基的遗传学锚定热点打分
  pharmacophore.json         # 12 个跨结构共识药效团特征 + 轴对距离
  pocket.json                # 5,440 Å³ 腔体、4 个亚口袋、埋藏度、化学组成
  axes.json                  # 15 条设计轴（14 条有候选）的枚举统计
  scaffolds_ranked.csv       # Top 200 骨架：SMILES + 性质 + 打分分解 + 接触残基
  diverse_scaffolds.sdf      # 分层去冗余后的 3D 位姿
  top_scaffold_poses.sdf     # Top 25 的 3D 位姿
  summary.json               # 一页纸汇总（漏斗、热点、特征、最佳骨架）
figures/
  fig0_flow.png              # 算法流程
  fig1_hotspots.png          # 热点排序
  fig2_pocket_map.png        # 口袋体素 + 三肽路径 + 药效团
  fig3_top_scaffold.png      # 最佳骨架在口袋中的位置
  fig4_funnel.png            # 筛选漏斗
```

## 可视化控制台（零依赖 Web UI）

```bash
python run_webui.py                             # 自动打开 http://127.0.0.1:8777/
python run_webui.py --port 9000 --no-browser    # 换端口 / 不自动开浏览器
python run_webui.py --host 0.0.0.0              # 局域网访问（无认证，勿直接暴露公网）
```

九个页面把 `results/`、`evidence/`、`docs/` 里的文件直接变成可读界面：
概览漏斗、骨架库（200 行可排序、可搜索、可导出 CSV）、界面热点（52 残基 + 逐残基详情）、
共识药效团与设计轴、口袋与四个亚口袋、**3D 结构浏览器**（三肽叠合 + 25 个位姿 +
体素腔体 + 药效团/热点球，纯 Canvas 2D，无 WebGL）、证据链（报告 / 24 条先验逐条核对 /
普查 / 界面测量 / 锚点）、运行管线（在浏览器里点按钮重跑并实时看日志）、文档。

- 只用 Python 标准库 + 原生 ES 模块：**不需要 npm、不需要联网、不需要 RDKit 或 WebGL**，
  在生成结果的同一环境里直接起。
- 服务端只读：页面不修改 `results/`、`evidence/`；唯一的写操作是“运行管线”按钮，
  且只能执行仓库里四个入口脚本（参数逐项校验，同时最多 2 个任务）。
- 自检：`python webui/test_api.py`（55 项 API/静态文件/越权断言）、
  `node webui/test_pages.mjs`（15 项页面渲染与路由断言，离线）。
- 用法、API 参考、排错与扩展见 **[webui/README.md](webui/README.md)**。

## 算法一句话

> 药效团对给出**一条设计轴**；骨架是沿轴排列的刚性单元链
> `头帽–[桥]–环–[桥]–…–尾帽`，每个单元有固定跨距，于是"能否同时够到两个锚点"
> 变成一个**一维长度约束**，用区间查找精确求解——确定性、可穷举、可复现。

详见 [docs/01_分子骨架设计算法.md](docs/01_分子骨架设计算法.md)；
过程与被否掉的方案见 [docs/04_工作日志.md](docs/04_工作日志.md)。

## 结果速览

| 指标 | 数值 |
|---|---|
| 结构输入 | 9NU4 / 9NU5 / 9NU8（3.6–3.7 Å） |
| MurJ | 5–500 共 496 残基（另有 106/7/14 个融合构建体残基已被剔除） |
| 界面残基 | 22 / 36 / 40 个；三种肽共同接触 14 个 |
| 肽的接触面 | 肽建模 29/29/40 残基，进入界面 18/19/28 个；接触重原子 64/111/121 |
| 小分子骨架 | 20–27 个重原子（中位数 24），约为肽接触原子的 1/5 |
| 埋藏腔体 | 5,440 Å³；最深亚口袋埋藏度 0.94 |
| 先验核对 | 24 条：PASS 22、FLAG 1（SglM D18 距 K46 6.67 Å）、WEAK 1 |
| 药效团 | 12 个共识特征；4 条主设计轴 |
| 枚举骨架 | 20,307 |
| 位阻 + 类药 + 肽模拟三重达标 | **538** |
| 最佳骨架 | `OC(=O)c1ccc(c2ccc(C(=O)Nc3cnccc3)cc2)cc1`（MW 318，cLogP 3.2，TPSA 79） |

最佳骨架与三种肽的空间重合度为 **1.00**，同时接触 K368/R53（阴离子头）与
L247/I248/T251/Y266（TM7 底板）——即"用刚性联苯把 TM2/TM7 的闭合面撑住"。

## 目录结构

```
pepwedge/
  core.py           PDB 解析、Shrake–Rupley SASA、接触/盐桥、热点打分、Kabsch
  pocket.py         体素化、射线埋藏度、亚口袋划分
  pharmacophore.py  药效团提取、跨结构共识聚类、设计轴对
  design.py         单元库、长度约束枚举、SMILES 与 3D 同步生成、性质包络
  pose.py           刚体放置、接触型打分
  priors.py         先验文件加载与校验
  evidence.py       全部结构测量（不依赖任何文献数字）
  verify.py         先验 × 坐标 的逐条核对
  pipeline.py       端到端编排与导出
  figures.py        图表
priors/
  murj_priors.json          文献先验（带原文引用 + 机器可执行的核对规则）
docs/
  01_分子骨架设计算法.md     方法学主文档
  02_下游接口指南.md         生成模型 / 对接 / MD 的对接协议
  03_证据链与可复现性.md     每个数字怎么算出来的、两个已修 bug、如何扩展
  04_工作日志.md             按时间顺序的过程记录：错误、被否掉的方案、数字快照
run_evidence.py      一键证据链（先跑这个）
run_pipeline.py      一键端到端
run_webui.py         一键 Web 控制台（可视化全部结果，也可在浏览器里重跑）
make_flow_figure.py  架构图
tests/test_smoke.py  冒烟测试（无 pytest 依赖）
webui/               控制台：server.py（标准库 HTTP）+ data.py（只读数据层）
                     + jobs.py（白名单任务）+ static/（原生 ES 模块前端）
```

## 设计原则

1. **遗传学锚定**：药效团点来自丙氨酸扫描 + 自发耐药突变 + 三肽趋同接触，不从形状拟合。
2. **机制可量化**：用"与三肽的空间重合度"把"像不像肽"变成 0–1 的数。
3. **约束而非采样**：骨架枚举确定、可复现、可归因。
4. **透膜性内生化**：性质包络是枚举期硬约束，且针对"胞质周靶点"重新调整（见主文档 §1.2）。
5. **测量与声称分离**：结构数字由 `evidence.py` 现算，文献断言放在 `priors/*.json`
   并逐条核对；复现不了的会以 `FLAG`/`FAIL` 显式列出，而不是被静默接受。
6. **接受局限**：位姿是刚体放置而非对接，最终必须由下游 MD / MM-PBSA 精修。

## 依赖

`numpy`、`matplotlib`。无 GPU、无 RDKit、无 PyMOL、无 BioPython。

> 说明：cLogP/TPSA 为碎片加和估计（与 RDKit 有 ±0.5 偏差），生产环境建议替换为 RDKit 的 Crippen 与 TPSA 实现。

## 引用

- Li, Y. E. et al. Convergent MurJ flippase inhibition by phage lysis proteins. *Nature* **652**, 1274–1280 (2026). DOI: 10.1038/s41586-026-10163-w
- 结构：PDB 9NU4 / 9NU5 / 9NU8
