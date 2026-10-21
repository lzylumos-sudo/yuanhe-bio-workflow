# Linux / WSL2 显卡电脑：公开最小验收

这套入口重建独立环境，运行 **1CRN + 乙醇 → 真实 Vina 对接 → 筛选与导出 → HF/RESP + GAFF2 建模 → 极短 CUDA MD → 分析与检查点审计**。

这是命令行计算链路验收，不是整个平台全部模块的验收。生成模型、网页操作、膜体系、长时间 MD、MM/PBSA 和 MM/GBSA 不在本次范围内。乙醇是公开工程对照，不能据此评价亲和力或稳定性；筛选正式合格数可以为零，MD 使用明确标记的工程对照，不修改正式筛选门槛。

## 目标电脑准备

- Linux x86_64，或 Windows 下的 WSL2 Linux；本次实测 Ubuntu 22.04 / RTX 5080。
- 已安装 NVIDIA 驱动，Linux/WSL 内能识别 GPU。锁包含 CUDA 13.4 用户态库，目标驱动必须兼容；脚本实际创建 CUDA 上下文检查，失败即停止，不静默回退 CPU。本包不安装驱动或更改系统环境。
- Git、Bash、micromamba，以及可访问 conda-forge、PyPI、GitHub、RCSB 和 Molekula 的网络。micromamba 按其[官方安装说明](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)准备即可，不要求初始化 shell。
- 完整科学依赖环境需要数 GB 存储；“最小”指计算量，不代表下载体积小。

## 安装和运行

以下命令均在 Linux / WSL 终端执行。使用仓库中包含本目录的版本；在仓库根目录执行后两条命令。

```bash
git clone --branch deployment/public-minimal-gpu-20261021 https://github.com/lzylumos-sudo/yuanhe-bio-workflow.git
cd yuanhe-bio-workflow

# 若 micromamba 不在 PATH，先设置为它的实际可执行文件路径：
# export MICROMAMBA=/absolute/path/to/micromamba
bash deployment/public_minimal/install.sh "$HOME/yuanhe-public-env"
bash deployment/public_minimal/run_public.sh "$HOME/yuanhe-public-env" "$HOME/yuanhe-public-run-001"
```

环境路径和结果路径必须是新的绝对路径，脚本拒绝覆盖。安装失败后保留现场，修复网络或依赖原因，再选择一个新环境路径重试；运行失败后检查日志，再选择新的结果路径。没有自动断点续跑，也不会自动安装、激活或修改原有研究环境。所有缓存、量化临时文件、公开来源快照和结果与原研究项目分开。

可设置 `CUDA_VISIBLE_DEVICES=0` 选择一张卡。该验收只使用单 GPU；不包含多 GPU 扩展验证。Vina 为 CPU 单线程，MD 为 OpenMM CUDA。安装中存在 CPU PyTorch 依赖，但没有运行生成模型。

## 固定版本与来源

- `linux-64.lock`：300 个 Conda 包的确切构建 URL 和校验值；不是重新求解的宽松环境文件。
- `pip-requirements.txt`：16 个 pip 包的固定版本/提交。RESP 固定上游提交和源码 SHA-256；安装后强制 `pip check`。
- 主要版本：Python 3.11.17、NumPy 1.26.4、SciPy 1.15.2、Biopython 1.85、OpenMM 8.6.1、AmberTools 26.0、Psi4 1.11、RDKit 2026.03.1、Meeko 0.8.0、Vina 1.2.7。
- 公开 1CRN 和 Vina 下载执行 SHA-256 检查。供应商库存页面每次重新获取，核对乙醇身份和 InStock 字段，留存抓取时间与哈希。页面变更、无法访问或无库存时失败，不伪造供应证据；未购买任何商品。
- 不附带第三方二进制、供应商网页、模型权重或私有分子。各依赖沿用上游许可证，参阅[研究依赖许可说明](../../docs/16_RESEARCH_DEPENDENCY_LICENSES.md)。

## 最小计算参数与通过条件

| 步骤 | 设置 / 结果要求 |
| --- | --- |
| 对接 | 1 个分子、1 个姿势、CPU 1 线程、exhaustiveness 1、seed 17 |
| 对接区域 | 1CRN 表面 10 Å 立方框，仅作工程检查，不是经验证的结合口袋 |
| 参数化 | 真正 HF/6-31G* 优化和两阶段 RESP，GAFF2，单线程 |
| MD | ff14SB / TIP3P；最小化上限 100 步；NVT 0.1 ps；NPT 0.1 ps；生产 0.2 ps（200 步，1 fs） |
| 分析 | 4 个生产帧，RMSD / RMSF / 接触 / 氢键 / SASA 等，原生审计通过且检查点可重载 |

成功时结果目录生成 `COMPLETED.json`，状态为 `PUBLIC_ENGINEERING_PATH_PASSED`。同时检查 `pip-check.txt`、`gpu.json`、`md_audit.json`、各阶段日志和 `artifact_hashes.json`。失败记录为 `FAILED.json`，不得把部分文件存在当成运行完成。

通过后，说明这台电脑能运行这条公开计算链路。正式研究仍需使用自己的经审查蛋白、口袋、化学状态和科学参数，在独立研究任务中执行；不要把本目录的 0.2 ps 参数用于正式结论。原始科研默认配置没有被这套验收修改。
