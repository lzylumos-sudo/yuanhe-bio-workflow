# WEDGE-Design 控制台（webui/）

把 `results/`、`evidence/`、`docs/` 里已经算好的文件变成一个可以点、可以看、可以重跑的界面：
漏斗、骨架库、界面热点、药效团、口袋、**3D 结构浏览器**、证据链、运行管线、文档。

设计约束只有一条：**零外部依赖**。
Python 侧只用标准库（`http.server`），前端是原生 ES 模块 + Canvas 2D，没有构建步骤、
不需要 npm、不联网、不需要 RDKit 或 WebGL。这样它能在生成结果的同一台机器上、
同一种环境里直接跑起来——控制台本身就是可复现性的一部分。

---

## 启动

```bash
python run_webui.py                      # http://127.0.0.1:8777/ 并自动打开浏览器
python run_webui.py --port 9000          # 换端口（被占用会自动顺延）
python run_webui.py --no-browser         # 不自动开浏览器
python run_webui.py --host 0.0.0.0       # 局域网访问（见下方“安全”）
```

首次访问 3D 页面时，口袋体素要从 PDB 现算（约 15 秒，单核），结果缓存到
`webui/cache/pocket_<hash>.json`；之后打开是即时的。缓存键包含 PDB 与
`pepwedge/pocket.py` 的参数，改了参数会自动失效。

**先跑过 `run_evidence.py` 与 `run_pipeline.py`**，否则页面会提示数据缺失
（这本身是有用的信息：控制台读的是文件，不猜）。

---

## 九个页面

| 页面 | 读的文件 | 能做什么 |
|---|---|---|
| 概览 | `results/summary.json` 等 | 四个 KPI、筛选漏斗（20,307 → 538 → 32）、Top-10 骨架、热点条形图、药效团类型分布、亚口袋体积、产物文件清单、五张图 |
| 骨架库 | `scaffolds_ranked.csv` + `structure_filters.csv`（按 SMILES 连接） | 200 行可排序/搜索/按轴过滤；勾选“仅三重达标”“仅过滤存活”；SMILES 正则过滤；导出 CSV；点行看打分分解（位阻/盐桥/氢键/疏水/肽重合/埋藏/性质惩罚/刚性）、F1–F4 判定、接触残基、单元装配配方；一键跳到 3D 看它的位姿 |
| 界面热点 | `interface_hotspots.csv` + `evidence/01_interface_murj.csv` | 52 个残基的遗传学锚定打分、三种肽的共同接触、逐残基的接触原子数 / dmin / dSASA / 收敛度 |
| 药效团 | `pharmacophore.json` + `axes.json` | 12 个共识特征（类型、权重、趋同度、容差、来源肽）、特征对的距离—权重散点与 9–18 Å 设计窗口、15 条设计轴的枚举统计（`A7-H12` 无候选也保留） |
| 口袋 | `pocket.json` + `pharmacophore.json` | 5,440 Å³、4 个亚口袋的体积/包络/埋藏度/化学组成（正电、负电、极性、芳香、疏水）、亚口袋与药效团的对应 |
| 3D 结构 | `9NU4/9NU5/9NU8.pdb`、`top_scaffold_poses.sdf`、`pocket.json` | 三种肽经 Kabsch 叠合到同一坐标系；MurJ Cα 轨迹按亚口袋着色；25 个位姿可切换；体素腔体按化学或亚口袋着色；12 个药效团球、24 个热点球；悬停拾取残基/特征 |
| 证据链 | `evidence/report.md`、`04_priors_check.csv`、`05_atom_budget.json`、`01/02/03_*.csv` | 报告正文 + 目录；24 条先验逐条核对（点开看原文引用、核对规则、观测值、出处）；结构普查与原子预算；界面测量（MurJ 侧 / 肽侧）；锚点几何与盐桥 |
| 运行管线 | 通过 API 触发 | 四个白名单脚本一键重跑，实时滚动日志、显示退出码与耗时、可取消、可导出日志；任务完成后自动刷新数据缓存 |
| 文档 | `docs/*.md`、`README.md` | 方法学主文档、下游对接协议、证据链说明、工作日志，按仓库原文渲染 |

导航是 hash 路由（`#/scaffolds?axis=A1-H3`），可以直接把链接发给别人。
右上角 `☾` 切换深/浅色（存在 localStorage）。

---

## 3D 结构浏览器的操作

| 操作 | 效果 |
|---|---|
| 拖拽 | 旋转（在相机坐标系里做，任意角度都不会翻转） |
| 右键拖拽 / Shift+拖拽 | 平移 |
| 滚轮 | 缩放（以视图中心为基准） |
| 双击 | 回到初始视角 |
| 悬停 | 显示 Cα 的残基名或药效团编号 |
| 单击 | 右侧显示该残基的亚口袋、dSASA、证据类型（或特征的类型/权重/接触残基） |

左侧显示开关分四层：**受体**（Cα 轨迹 / 侧链球棍 / 界面 13 Å 内的重原子）、
**配体**（25 个位姿任选，原子半径放大以突出）、**肽**（Sgl^M / Sgl^PP7 / Sgl^CJ3
各自可开关）、**口袋**（体素云：关闭 / 按化学 / 按亚口袋；药效团球；热点球 / 标签）。

渲染是画家算法 + 球/棍/点三种图元，在 `Canvas 2D` 上手写投影与深度排序。
没有 WebGL：某些远程桌面或虚拟机上 WebGL 不可用，而这块画布在任何浏览器里都能画。

坐标系：所有结构都在 **Sgl^CJ3（9NU8）** 的坐标系里。切换“参考结构”只是重新
用该复合物做一次 Kabsch 叠合，位姿与药效团跟着一起变换，所以口袋不会错位。

---

## 运行管线

只能启动仓库里的四个入口脚本，参数逐项校验（路径白名单正则、整数/浮点范围）：

| 任务 | 脚本 | 参数 |
|---|---|---|
| 证据链重建 | `run_evidence.py` | `--out`（默认 `evidence`） |
| 设计管线 | `run_pipeline.py` | `--out`（默认 `results`） |
| 结构衍生过滤器 | `run_filters.py` | `--out`、`--targets`（0–5000）、`--max-spread`（0.5–50） |
| 冒烟测试 | `tests/test_smoke.py` | 无 |

- 同时最多 2 个任务；同类任务不允许并发（避免两个进程写同一个 `--out`）。
- 输出目录必须是相对路径且不含 `..`；不接受任意命令、不接受 shell。
- 子进程的 stdout/stderr 合并后按行流式回传，`PYTHONIOENCODING=utf-8` 保证中文不乱码。
- 取消 = `SIGTERM`；任务状态为 `cancelled`，不会显示成失败。
- 任务完成后前端 `invalidate()` 掉内存缓存，下次进页面读到的是新文件。

---

## API

全部返回 JSON（`Content-Type: application/json; charset=utf-8`），错误统一为
`{"error": "...", "status": 4xx}`。

| 方法 | 路径 | 说明 |
|---|---|---|
| GET | `/api/health` | 版本、时间、任务数（探活） |
| GET | `/api/meta` | 工作目录、三个结构、产物文件 mtime、可用性、figure 清单 |
| GET | `/api/overview` | 概览页一次性载荷（漏斗、最佳骨架、热点前 12、口袋、证据状态） |
| GET | `/api/hotspots` | 52 个界面热点（按打分降序）+ 共同接触数 |
| GET | `/api/pharmacophore` | 12 个特征、特征对、设计轴、口袋 |
| GET | `/api/axes` | 15 条设计轴 |
| GET | `/api/pocket` | 腔体体积、亚口袋、化学组成、埋藏度 |
| GET | `/api/pocket/points?step=0.8` | 体素点云（`status` = `computing` 时返回 **202**，算完返回 200） |
| GET | `/api/scaffolds` | 200 个骨架（含 F1–F4 过滤结果与位姿下标） |
| GET | `/api/filters` | 538 行结构过滤明细 + 汇总 + 校准 |
| GET | `/api/poses` | 25 个 Top 位姿 + 24 个去冗余位姿（原子/键/SMILES） |
| GET | `/api/scene?ref=SglCJ3&ligand=0&radius=13` | 3D 场景：受体、三条肽、配体、药效团、热点、亚口袋 |
| GET | `/api/evidence/report` | 报告 Markdown + 证据状态 + 原子预算 |
| GET | `/api/evidence/priors` | 24 条先验核对（含 `claim`/`quote`/`locator`/`observed`） |
| GET | `/api/evidence/census` | 结构普查 + 原子预算 |
| GET | `/api/evidence/interface` | MurJ 侧 / 肽侧界面测量 + 锚点 + 热点 |
| GET | `/api/docs` | 文档清单（目录、标题、大小、mtime） |
| GET | `/api/docs/<name>` | 单篇 Markdown 原文 |
| GET | `/api/jobs` | 最近 40 个任务（不含日志行与进程句柄）+ 任务类型 |
| GET | `/api/jobs/<id>?since=N` | 任务状态 + 从第 N 行起的日志（增量轮询） |
| POST | `/api/jobs` | `{"kind": "pipeline", "params": {"out": "results"}}` → 201 `{"id": ...}` |
| POST | `/api/jobs/<id>/cancel` | 取消运行中的任务 |
| GET | `/files/<dir>/<path>` | 下载 `results/`、`evidence/`、`figures/`、`docs/`、`priors/` 下的文件；`?download=1` 加附件头 |

前端静态文件在 `/`、`/css/*`、`/js/*`。缺失的静态资源返回 404（不是 SPA 首页），
无扩展名的未知路径回落到 `index.html`，因为路由在 hash 里。

---

## 代码结构

```
run_webui.py            入口：--host/--port/--no-browser/--quiet
webui/
  server.py             路由表 + 静态文件 + 端口顺延；~450 行，只用标准库
  data.py               只读数据层：CSV/JSON/SDF/PDB 解析、按 mtime 记忆化、
                        Kabsch 叠合、体素点云缓存（后台线程 + 锁）
  jobs.py               白名单任务运行器：参数校验、并发上限、行缓冲、取消
  static/
    index.html          外壳：侧栏 + 顶栏 + #view
    css/app.css         设计系统（深/浅色、栅格、表格、图表、抽屉、控制台）
    js/api.js           fetch 封装 + 内存缓存 + invalidate
    js/ui.js            hyperscript、格式化、DataTable（排序/分页/选中）、
                        条形图/漏斗/图例、Markdown 渲染、抽屉、CSV 导出
    js/viewer.js        Canvas 2D 3D 渲染器（投影、深度排序、拾取、交互）
    js/app.js           hash 路由、导航、主题、任务轮询
    js/pages/*.js       九个页面，各自 mount(root, ctx) 并返回清理函数
  test_api.py           API 端到端断言（可加 --job 跑一次真实任务）
  test_pages.mjs        前端渲染/路由断言（Node + 自制 DOM 垫片，离线）
  capture_fixtures.py   把 API 响应快照到 webui/testdata/（test_pages 用）
  check_js.mjs          语法检查所有 ES 模块（不需要浏览器）
```

数据层只读：`data.py` 不做任何写入，除非是 `webui/cache/` 下的体素缓存。

---

## 自检

```bash
python webui/test_api.py          # 55 项：每个接口的形状与语义、静态文件、越权与路径穿越
python webui/test_api.py --job    # 再加 4 项：真实启动一个任务并等它完成
node webui/test_pages.mjs         # 15 项：九个页面渲染、表格排序/过滤、抽屉、
                                  #      Markdown、路由切换、主题切换、3D 投影
node --experimental-vm-modules webui/check_js.mjs   # 13 个模块的语法
```

`test_pages.mjs` 用一个 ~600 行的 DOM/Canvas 垫片把页面代码跑在 Node 里，
数据来自 `webui/testdata/*.json`（用 `python webui/capture_fixtures.py` 重新生成）。
它不能验证像素，但能抓住“字段名写错、导入漏了、渲染时抛异常、路由没挂上”这类问题——
在没有浏览器可用的环境里，这是能做到的最强验证。

`tests/test_smoke.py` 仍然是算法本身的回归测试；控制台的测试只覆盖界面与 API。

---

## 排错

| 现象 | 原因 / 处理 |
|---|---|
| 页面顶部显示“数据不完整” | `results/` 或 `evidence/` 缺文件；跑对应的 `run_*.py` |
| 3D 页面腔体一直“正在计算” | 首次体素化约 15 s；看终端是否有异常；删 `webui/cache/` 重试 |
| 端口被占用 | `--port` 会自动顺延到空端口，启动横幅会打印实际地址 |
| 中文变成方块/乱码 | Windows 控制台编码问题；页面本身是 UTF-8，子任务进程已强制 `PYTHONIOENCODING=utf-8` |
| 改了 `pepwedge/` 但界面没变 | 重跑脚本；界面读的是文件，不缓存业务数据（只缓存 JSON 载荷，刷新即失效） |
| 任务按钮点了没反应 | 检查终端日志；同类任务并发会被 409 拒绝；`jobs.py` 里 `MAX_LINES` 只保留最近 6000 行 |
| 想清空任务历史 | 重启服务（历史在内存里） |

---

## 安全

- 默认只绑定 `127.0.0.1`，没有登录，**不要直接暴露到公网**；需要局域网访问时
  至少放到反向代理后面并加认证。
- `POST /api/jobs` 只能跑四个白名单脚本，路径参数走正则白名单，不拼接 shell。
- `/files/` 只允许 `results|evidence|figures|docs|priors`，`..` 与绝对路径被拒绝
  （已用测试覆盖）。

## 已知限制

- 位姿是刚体放置（rigid placement），不是对接；控制台展示的是几何可行性，不是结合自由能。
- 3D 视图没有卡通/表面渲染，也没有真正的遮挡剔除（体素云按深度排序近似）。
- 3D 画布始终是深色（像显微镜视野），浅色主题只作用于界面外壳。
- 没有做移动端专门优化：<900 px 时侧栏折叠为抽屉，表格横向滚动。
- 图（`figures/*.png`）由 matplotlib 生成，控制台只负责展示，不重新绘制。
