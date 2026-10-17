# 测试报告

环境：Windows（中文路径工作区）、Python 3.13.5、独立 venv `C:\Users\user\Desktop\元核生物\.venvs\zvss`。
命令：`python -m pytest tests -q`（在 `06_Virtual_Screening_Studio` 下执行）。

## 1. 结果

```
29 passed
```

- 单元测试 26 项：`pytest -m "not integration"`
- 端到端测试 3 项：`pytest -m integration`（从真实 GNINA 输出裁剪 2 个 entry 组成 fixture，跑完整 12 阶段流水线）

## 2. 覆盖内容

| 文件 | 断言要点 |
| --- | --- |
| `tests/test_ids.py` | CMP/STA/POS/SEL 可复现；坐标摘要抗 1e-6 噪声；`selection_version` 与成员顺序无关；缓存键对每个分量敏感 |
| `tests/test_config.py` | generic 无局部区域；ecmurj 声明 P0AF16 / native_range 3–500 / 三处突变 / TM2+TM7；缺必需键与超范围序列必须报 `ConfigError`；预算上下界 |
| `tests/test_chem.py` | 规范 SMILES 等价；cLogP 明确标注 “NOT logD” 且 `logd_available=false`；Murcko 骨架；模板回收保身份；对称等价 RMSD 为 0、真实平移后 > 0.1 Å |
| `tests/test_score_semantics.py` | CNNaffinity 语义为 pK 且不可折算 kcal/mol；分数缺失必须为 `None` 而不是 0；真实 GNINA 日志可解析出 argv 与 mode 表 |
| `tests/test_sifts_mapping.py` | 真实 9NU4 mmCIF + P0AF16 序列：498 残基映射、614 construct 残基、3 处序列差异恰为 auth 5/12/13、缺失 13 个 native 位点、auth 3↔native 3、auth 500↔native 500 |
| `tests/test_pipeline_e2e.py` | fixture 恰好导入 6 个位姿且 0 失败；两次独立运行 `selection_version` 相同、`candidates.csv` 相同；MD 包文件齐全、`ligand.resname == LIG`、`complex.pdb` 含 `HETATM`/`LIG`/`REMARK`、`task.json` 带真实阻塞码且状态为 `blocked` |

## 3. 真实数据运行记录（比测试更强的证据）

| run | 输入 | 结果 |
| --- | --- | --- |
| `runs/t4` | `F:\MemPilot\runs\public20_engineering`（GNINA 1.3.3，20 parents × 3 modes） | 313 文件 → 60/60 位姿导入、0 失败；质检 56 pass / 4 flag；60 位姿接触分析；tier B=18 C=2 D=0；2 primary + 1 exploratory；3 个 MD 任务包 |
| `runs/t5` | 同上，独立重跑 | `selection_version` 与 t4 相同（`SEL6830370C0C42E694`），`candidates.csv`、`ranking.csv` 逐字节一致 |
| `runs/murj1` | 真实 9NU4 + UniProt P0AF16 | SIFTS 映射 498、序列差异 3（=K5F/S12I/M13A）、缺失 13（1–2、501–511）、受体 4598 原子、膜 `proxy_only`、MD `blocked` |

Web 冒烟（本地 8770 端口，测后关闭）：`/api/health`、`/api/runs`、`/`、`/run/t4`、`/report/t4`、`/api/run/t4` 均 200，不存在的 run 返回 404。

## 4. 已知未覆盖

- 多 seed 再对接路径：本机无 `gnina`，只有“不可用 → 阻塞”的分支被真实执行过，成功分支未运行。
- 膜体系构建与膜 MD：能力缺失，仅警告/阻塞路径覆盖。
- PoseBusters 成功路径：当前会话临时目录为中文路径，RDKit 数据文件无法加载，因此报告为 `unsupported`；把 `ZVSS_SCRATCH` 指到 ASCII 目录后可导入，但未在此环境完成整轮带 PoseBusters 的运行。
- Web 端未做浏览器自动化测试，仅做 HTTP 端点冒烟。
