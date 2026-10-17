# ZVSS · Zhiyu Virtual Screening Studio

本地虚拟筛选工作台：吃进**已有的 GNINA 批量对接结果**，做自动化的七类筛选计算，产出**可解释候选 + 排除理由 + MD 任务包**，并把任务包适配到已有的 `05_molecular_dynamics`（OpenMM）与 `04_interaction_analysis`（ZLIS/ProLIF）项目。

本目录是独立新项目，不改动 F 盘任何既有项目；所有跨项目衔接都通过显式适配器与文件契约完成。

## 目录

```
06_Virtual_Screening_Studio/
├─ src/zvss/                 # 包源码（ingest / receptor / qa / mech / motion / redock /
│                            #   scoring / strain / select / membrane / md / report / api）
├─ configs/
│  ├─ host.yaml              # 主机能力：工具路径、执行器、开关
│  ├─ targets/ecmurj.yaml    # MurJ（9NU4 / P0AF16），含工程突变、标签、TM2/TM7 区域
│  ├─ targets/generic.yaml   # 通用靶点（1DRF），无局部区域
│  └─ budgets/{full,quick}.yaml
├─ tests/                    # pytest：单元 + 端到端（integration 标记）
├─ runs/<run_id>/            # 每次运行的完整证据：配置、日志、阶段产物、报告、MD 包
├─ runtime/                  # SQLite 库、参考结构缓存
└─ docs/                     # ARCHITECTURE / MD_INTERFACE / TEST_REPORT / COMPLETION_MATRIX
```

## 安装

```powershell
# 已创建的独立虚拟环境（Python 3.13）
# C:\Users\user\Desktop\元核生物\.venvs\zvss
$py = "C:\Users\user\Desktop\元核生物\.venvs\zvss\Scripts\python.exe"

# 本机 pip 在 PEP 517 隔离构建下会因临时目录权限失败，因此禁用构建隔离：
& $py -m pip install --no-build-isolation -e "C:\Users\user\Desktop\元核生物\06_Virtual_Screening_Studio"
& $py -m pip install pytest            # 仅测试需要
```

依赖版本见 `docs/DEPENDENCIES.md`；RDKit / PoseBusters 需要 ASCII 临时目录，见下文“已知限制”。

## 启动

```powershell
$zvss = "C:\Users\user\Desktop\元核生物\.venvs\zvss\Scripts\zvss.exe"
$root = "C:\Users\user\Desktop\元核生物\06_Virtual_Screening_Studio"

# 1) 环境自检
& $zvss --root $root doctor --with-tools

# 2) 一次完整筛选（消费既有 GNINA 结果目录）
& $zvss --root $root run `
    --input "F:\MemPilot\runs\public20_engineering" `
    --target "$root\configs\targets\generic.yaml" `
    --budget "$root\configs\budgets\quick.yaml"

# 3) 本机只读 Web 界面（127.0.0.1）
& $zvss --root $root serve --port 8770   # 浏览器打开 http://127.0.0.1:8770

# 4) 重渲染报告 / 查看历史运行
& $zvss --root $root report --run-id <run_id>
& $zvss --root $root runs
```

注意：`--root` 是全局参数，必须写在子命令之前。

## 真实运行记录（本工作区，非示例数据）

输入：`F:\MemPilot\runs\public20_engineering`（GNINA 1.3.3，1DRF chain A，20 parents × 3 modes = 60 poses）

| run | 靶点 | 结果 |
| --- | --- | --- |
| `runs/t4` | generic（1DRF） | 313 个文件 → 60 位姿全部导入；质检 56 pass / 4 flag；60 位姿接触分析；tier B=18 C=2；2 primary + 1 exploratory；3 个 MD 任务包 |
| `runs/t5` | generic（1DRF） | 独立第二次运行：`selection_version` 与 `t4` **完全相同**，`candidates.csv` / `ranking.csv` 逐字节相同 |
| `runs/murj1` | ecmurj（9NU4 / P0AF16） | 真实结构制备：SIFTS 映射 498 残基、3 处序列差异即工程突变、缺失 13 个 native 位点（1–2、501–511）、受体去水去离子去配体后 4598 原子、膜几何代理 `proxy_only`、MD 明确 `blocked` |

`t4` 的候选（`SEL6830370C0C42E694`）：

```
primary,1,CMP7047E1960DBD3C23,...,B,0.7555
primary,2,CMPCBBA17C44FCCD2B3,...,B,0.5695
exploratory,1,CMP75C7A4B0B41E4390,...,B,0.6725
```

## 能力边界（诚实声明）

- **已用真实数据验证**：GNINA 结果解析与分数绑定、化学身份与位姿几何、PoseBusters 之外的确定性 QC、接触分析、位姿一致性聚类、性质/骨架/多样性、证据分层与多目标排序、候选与排除、MD 包结构、SIFTS 编号映射、CLI 与 Web。
- **本机不可用（已如实标注，不当作通过）**：多 seed 再对接（无可用 `gnina` → `blocked_external_tool`）、膜体系构建与膜 MD（`blocked_missing_membrane_capability`）、Linux-only MD 执行（`blocked_platform_linux_required`）。
- **需要 ASCII 临时目录**：RDKit 数据文件与 PoseBusters 无法从含中文的路径加载；`ZVSS_SCRATCH` 指向 ASCII 目录时 PoseBusters 可用，否则报告写 `unsupported`。
- **MurJ 结论边界**：本工作区没有任何 MurJ 对接位姿，因此 MurJ 的位姿级筛选是 `not_evaluated`；已验证的是受体/编号/区域/参考结构这条结构证据链，结论仅为“可进入 MD 的候选优先级 + 假设支持性证据”。

更细的内容见 `docs/ARCHITECTURE.md`、`docs/MD_INTERFACE.md`、`docs/TEST_REPORT.md`、`docs/COMPLETION_MATRIX.md`。
