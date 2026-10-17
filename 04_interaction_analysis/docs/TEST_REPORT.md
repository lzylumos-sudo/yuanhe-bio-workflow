# ZLIS 测试与安装记录

## 环境（macOS 26.5.1 arm64，conda 环境 `/opt/miniconda3/envs/zlis`）
Python 3.11.17，RDKit 2026.9.1，ProLIF 2.2.3，MDAnalysis 2.10.0，Biopython 1.88，OpenMM 8.6.1，pandas 3.0.6，openpyxl 3.1.5，numpy 2.4.6，gemmi 0.7.5，cairosvg 2.7.1，pytest 8.3.4，PyMOL 3.1.0（pymol-open-source）。
可选 PLIP 3.0.0 在独立环境（路径写在 `env/plip_python.txt`）。

## 自动测试：40 项全部通过（从移动硬盘上的副本运行，原始输出见 `pytest_record_macos.txt`）
- `tests/test_thresholds.py`（16）：氢键（配体供体、蛋白供体、3.49/3.55 Å 距离边界、125°/135° 角度边界）、盐桥（两种方向、4.49/4.6 Å、同号电荷负例）、疏水、π-π、阳离子-π（含中性胺负例）、卤键、金属配位、水桥（含无水负例）、event_id 唯一。
- `tests/test_pipeline.py`（17）：同一残基两条氢键都列出、供体侧与主链/侧链、重原子坐标不变、atom_uid 贯通、PDB 仅含真实原子、SDF 原子顺序与坐标不同时坐标不变、错误模板被报告、未知配体标 needs_confirmation、缺氢补氢并记录、插入码、多 MODEL 不混合、多链同 resseq、含空格与中文路径的复制目录可运行、PSE 重新打开对象仍在、PyMOL 脚本不含相互作用判定逻辑。
- `tests/test_diagram.py`（7）：三种视图标签不越界不重叠、氢键箭头方向、非氢键无箭头、layout.json 往返、SVG/PNG/PDF/HTML 生成。

## 真实演示
| 体系 | 配体 | 结果 |
|---|---|---|
| 3PTB（胰蛋白酶） | A:BEN1（模板库，苯甲脒阳离子） | 氢键 3、盐桥 4、水桥候选 4、疏水 2、vdW 5 |
| 1HSG（HIV 蛋白酶） | B:MK1902（RCSB 化学组分默认质子化，needs_confirmation） | 氢键 4、水桥候选 13、疏水 7、vdW 9 |

3PTB 氢键：ASP189 N2–OD2 2.87 Å ∠147.0°，ASP189 N1–OD1 2.92 Å ∠151.3°，GLY219 N1…O 2.82 Å ∠144.6°（主链羰基）。
1HSG 氢键：A:ASP25（配体供体）、B:ARG8（蛋白供体）、B:GLY27（配体供体，主链）、B:ASP29（蛋白供体，主链 N）。结果在 `results/`。

## PLIP 交叉检查（独立规则，差异不自动调和）
- 3PTB：PLIP 找到 2 条（ASP189），未报 GLY219 的氢键；ZLIS 报 3 条。PLIP 未报盐桥，ZLIS 报 ASP189。
- 1HSG：PLIP 6 条，ZLIS 4 条；PLIP 报 B:ASP25，ZLIS 报 A:ASP25（链归属不同，未查明原因）。两条引擎的差异需要人工核对，尚未裁决谁对。

## 已知问题与未完成
- **Windows 未测试**：`Start_Windows.cmd`、`run.ps1`、`env/setup_windows.ps1` 只写完，没有在 Windows 上运行。
- 2D 图与 PyMOL PNG 我无法目视检查（所用模型无图像输入），只做了程序化几何检查（越界、重叠、箭头方向）。请在发布前人工看一遍。
- 2D 编辑器的浏览器拖动与保存按钮没有做浏览器自动化测试，只测了 `layout.json` 往返。
- 离子型相互作用与氢键在同一原子对上会同时出现（盐桥带氢键），属于两种几何判据各自成立，不是重复计数。
- 配体质子化与 His 状态（均为 HID）是模型假设，会影响氢键与盐桥；PLIP 与 ProLIF 的差异可能部分来自此处。
- 水桥只是候选，不分供受体。
- ProLIF 的 π-π 判据对 3.8 Å 面对面叠合会同时给出疏水，已原样保留。
- 驱动器为 exFAT，macOS 会生成 `._*` 文件，已清理但访问时可能再次出现。
- 未实现：多构象（altloc）选择界面（默认取占有率最高者）、批量多文件运行界面。
