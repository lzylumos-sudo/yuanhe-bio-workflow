# Nexus MD v2 中文使用说明

## 安装与首次验收

在旁边的 x86-64 Ubuntu 电脑，或 Windows 的 WSL2 Ubuntu 中，打开终端，进入移动硬盘实际挂载位置中的 `Nexus_Structural_Biology/05_molecular_dynamics`。Linux 电脑不需要登录 GPT。

```bash
bash setup.sh
source runtime.sh
python scripts/check_environment.py --native --output installation_evidence/my-environment-check.json
python scripts/validate_v2_ubuntu.py
```

安装优先校验并恢复本模块随附的 Linux 环境包。环境解压到 Linux 用户缓存目录；源码、配置与结果随硬盘保存。没有环境包时按固定包 URL 锁联网重建。锁已修正为 `@EXPLICIT`，但旧清单没有逐包校验和；随附环境包有完整 SHA-256。已有环境会检查实际导入、工具和极小 CPU/checkpoint 运行，不能仅凭 Python 文件存在就算成功。

环境损坏时，`bash setup.sh --repair` 将原环境移到带时间戳的备份位置再重建，不删除模拟结果。不要移动解压后的 Conda 目录；搬迁后在新机器重新运行安装脚本。`NEXUS_ENV_PREFIX` 可指定 Linux 文件系统中的安装位置。项目路径不固定为任何盘符或 Mac 挂载点。

最后一条命令使用公开 1FKJ 小体系，固定 1 ps NVT + 1 ps NPT + 10 ps 生产，包含真实中断、二进制 checkpoint 载入与恢复、分析和完成任务重复续跑检查。默认超时 1200 秒，可用 `--max-seconds 3600` 调整墙钟超时，模拟时长仍不变。没有中断成功会报告 NOT_RUN，不能当成续跑通过。本轮在 Mac 仅检查了入口的平台保护，尚未执行这项 Ubuntu 分子验收。

## 日常命令

```bash
# 交互向导；只扫描本项目 inputs/，多个候选时让你选择
bash run_md.sh

# 输入任意文件名和路径，.PDB 大写也支持
bash run_md.sh --input "/实际路径/我的复合物 (R2).pdb"

# 使用保存的设置；配置路径可以在项目外
bash run_md.sh --input "/实际路径/复合物.pdb" --config config/v2_short.yaml

# 只检查/整理输入，不参数化，不模拟
bash run_md.sh --input "/实际路径/复合物.pdb" --inspect-only

# 结构准备、参数化、溶剂化和体系构建，不推进模拟
bash run_md.sh --input "/实际路径/复合物.pdb" --prepare-only

# 强制工程短测：1 次、CPU 2 线程、1+1 ps 平衡、10 ps 生产
bash run_md.sh --input "/实际路径/复合物.pdb" --smoke-test

# 无人值守遇到信息不足时保存 pending.json 并退出，不等待键盘
bash run_md.sh --input "/实际路径/复合物.pdb" --config config/v2_short.yaml --non-interactive

# 用实际打印的任务 ID 恢复
bash run_md.sh --resume "runs/<任务ID>" --decisions "/实际路径/decisions.yaml"
```

默认新入口：**1 次 100 ps 生产、100 ps NVT + 100 ps NPT、CPU 2 线程**。每次显示实际设置；`--smoke-test` 明确覆盖为更短测试。正式时长、重复数、温度、压力、盐浓度、输出间隔和 CPU/GPU 在 YAML 中指定。旧 `run.sh` 仍保留原默认长任务，不能与新入口混用。

指定 `--input` 时只使用该文件。不指定时仅在明确的输入目录查找，不能把 `runs/` 或 `results/` 当新输入。每次新建唯一任务目录；相同文件名不共用结果，已有 `--run-id` 拒绝覆盖。旧任务使用 `run.sh --run-id ID --resume`，见 `LEGACY_README.md`。

## 你需要准备什么

通常只需一份完整复合物 PDB，以及希望的模拟设置。不要求改名、拆蛋白/配体或统一手工准备 SDF、ITP/TOP。流程优先读取前面项目的内嵌 SDF/manifest 或旁置模板；模板必须通过元素、连接、映射、电荷与立体化学检查。只有 PDB 时，自动尝试根据 CONECT、显式氢、电荷和坐标生成候选，全部在本地进行。

仍需你决定的是化学和研究对象本身：

- 蛋白是否就是实际目标构建体；只有 SEQRES 或连续编号不算完整序列验证。
- 仅由 PDB 推断的配体候选是否确实为目标分子。核对键级、质子化、原始电荷位点和手性后确认候选哈希；也可补充准确 SMILES，不必统一制作新配体文件。
- 多 MODEL、altLoc、保留哪些组分、特殊末端、输入水是否可被新溶剂化替代。
- 蛋白默认加氢 pH、力场整套预设、模拟时长与资源。pH 不是位点 pKa 预测。

结构名称兼容问题由流程处理；这些科学决定记录在 `decisions.json`。缺少决定就暂停，原始输入和已完成工作保留。以下只是决定文件的**格式示意**，ID 与哈希必须来自实际任务，不能对所有输入照抄批准：

```yaml
construct_confirmed: true  # 只有核对构建体后才设为 true
termini:
  r0188: standard_charged  # 仅在你已选择该残基的常规带电末端时
chemistry:
  c002:
    accept_candidate_sha256: "填写已核对候选的完整哈希"
    # 或 smiles: "包含正确电荷和立体化学的 SMILES"
# model: 1
# altloc: {r0010: A}
# keep_components: [c001, c002]
# water_policy: replace_bulk
```

也可在 YAML 配置的 `decisions:` 下写同样内容。蛋白原有重原子在准备时核对保留位置；简单缺失侧链可按模板补齐并记录。HID/HIE/HIP，以及 HIS 已有环氢证据保留。HXT 不会直接改名 OXT，只有明确选择常规末端后才删除 HXT、按模板补 OXT并记录变化。大段主链缺失或异常连接不自动补造。

## 支持范围

`bash run_md.sh --list-presets` 显示完整预设状态。现在只启用 AMBER ff14SB + GAFF 2.11/AM1-BCC + amber14 TIP3P 及相应 Na+/Cl− 路线。ff19SB/OPC、CHARMM36m/CGenFF、膜或非水溶剂尚未完成接入验证，不能随意混搭。模拟参数可手动设定；当前不能任意切换力场和溶剂。

| 输入 | 状态 |
|---|---|
| 普通非共价蛋白—单小分子 | 输入层及原 OpenMM 链路已接通；v2 仍须 Ubuntu 完整验收 |
| 蛋白—蛋白，多链蛋白 | 可识别并接原后端；小分子 MM/GBSA 分区不适用 |
| 多蛋白链 + 一个分子胶 | 可识别三元体系；可在 manifest/YAML 给链分配分析角色 |
| 多个配体/同配体多个拷贝 | 输入层全部识别保留；当前模拟适配器只支持单配体，明确阻塞 |
| 金属中心、共价配体、辅因子、复杂糖基化、膜、主链缺失 | 识别并阻塞/分流；没有冒充普通体系支持 |

## 输出与续跑

`runs/<任务ID>/` 保存原始输入和 SHA-256、完整原子清单、每次 `normalized/attempt_*`、化学候选/确认、配置、决定、逐阶段日志、体系和拓扑、轨迹/checkpoint、编号映射与验收。

`READY` 只表示输入检查通过；`PREPARED` 表示准备和建系完成；`COMPLETE` 必须经过实际 MD、分析、DCD/CSV/最终 checkpoint 验收。WAITING_FOR_DECISION、DEPENDENCY_BLOCKED、FAILED、INTERRUPTED 分别表示等待决定、依赖阻塞、失败和中断。退出码 0 为完成本次模式，2 失败，3 待决定，4 依赖/平台阻塞，5 已有运行锁，130 用户中断。

续跑检查输入、配置、决定、代码、依赖及产物指纹。改变质子化、力场或上游结构后建立新任务，不复用旧拓扑。生产恢复使用真实二进制 checkpoint，保留速度和运行状态；未提交的部分段可能从上个检查点重算。跨硬件 checkpoint 兼容性需实测，失败时不偷偷改成最后一帧重新赋速度。

断电留下 `.running` 锁时，先查看 `owner.json`，在对应主机确认进程已经停止，再移除该任务的锁目录；不要删 checkpoint 或 progress。运行前检查至少 2 GiB 可用空间，这不是长任务容量估算。

分析复用原有蛋白 Cα 对齐 RMSD、蛋白对齐后的配体重原子 RMSD、Cα RMSF、接触/氢键/距离/SASA，处理周期边界和分子完整性。`mappings/` 追踪原始原子至最终拓扑，`*_original.csv` 补回原始残基标识；重建氢明确记为新原子。蛋白 RMSD 的选择为 Cα，不声称是 N/CA/C/O 全主链 RMSD。MM/GBSA 无熵，三元体系默认全部蛋白对小分子分区，不等于协同性。

短测试不证明充分平衡、收敛、结合更强或降解活性。

## 开发测试与回滚

```bash
source runtime.sh
python -B -m unittest discover -s tests -p 'test_*v2.py' -v
python -B -m unittest discover -s tests -p test_embedded_chemistry.py -v
python scripts/validate_v2_ubuntu.py
```

真实输入回归为 `tests/fixtures/complex_R2.pdb`。测试里模拟候选批准仅验证软件机制，不构成科学确认。当前已执行/未执行项见 `V2_VALIDATION.json`。

部署前原件在 `.backups/` 的带时间戳目录，附部署清单和 `rollback.py`。回滚前停止相关任务，在备份目录运行 `python3 rollback.py --apply`。它先检查现文件是否仍为本次部署版本，遇到后续改动会拒绝覆盖；原轨迹和新 runs 不删除。
