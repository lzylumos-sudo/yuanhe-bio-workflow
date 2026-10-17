# MurJ v2 — Peptide2Mol + REINVENT4 可审计重建

本目录是与旧 `pepwedge/` 结果隔离的重建实现。旧版结果仅作历史基线；旧 SDF 已被 RDKit 证实不可回读，因此不得作为 v2 化学输入。

## 范围

- 主结构：项目根目录的 `9nu4.pdb`（与用户指定 `/Users/lizhiyu/Desktop/9nu4.pdb` SHA-256 完全一致）。
- 主生成器：官方 Peptide2Mol，固定提交和 checkpoint 哈希。
- 扩展器：官方 REINVENT4 / LibInvent，固定提交和 prior 哈希。
- 交付目标：对接前 SMILES CSV、旁表、manifest、指标与失败记录。
- 明确不运行：GNINA/Vina 正式对接、分子动力学、自由能计算。
- 禁止：模型失败时回退到旧手工枚举并冒充模型输出。

## 固定的官方依赖

| 组件 | 固定版本 |
|---|---|
| Peptide2Mol | commit `2c5af9628aeeab5ed2dff5dc5326ce2463316291` |
| PMT checkpoint | SHA-256 `7dc82c4f895efea83764ce25c2fae5787f9f25b7935590e5119bde7199540129` |
| REINVENT4 | commit `3cf206ea5b385ac955dc44bbae6063b55ee31e0a` |
| LibInvent prior | Zenodo version DOI `10.5281/zenodo.20701824`;本地 SHA-256 由下载后 manifest 记录 |

## 目录

- `inputs/raw/`：原始输入引用/清单，不改写根目录 PDB。
- `inputs/prepared/`：稳定准备产物。
- `configs/project.json`：项目参数与证据边界。
- `src/structure.py`、`hypotheses.py`：结构准备与局部假设。
- `src/adapters/`：真实模型子进程适配器。
- `src/chemistry.py`、`cores.py`、`synthesis.py`：身份、核心及路线证据。
- `src/geometry.py`、`diversity.py`、`export.py`：构象、硬门槛、多样性及导出。
- `runs/<run_id>/`：原始日志与分阶段产物；成功阶段才有 `stage_success.json`。
- `exports/<run_id>/`：通过正式门槛的交付目录；官方 demo、mock 和旧版分子不会进入。

## 环境

各模型独立环境，使用文件契约通信：

```bash
# 化学/QC 测试使用现有 RDKit 环境
/opt/miniconda3/envs/zlis/bin/python -m unittest discover -s v2/tests -v

# Peptide2Mol（Python 3.9；官方主要验证 CUDA/NVIDIA）
v2/envs/peptide2mol/bin/python ...

# REINVENT4（Python 3.11；本机 macOS/MPS 或 CPU）
v2/envs/reinvent4/bin/reinvent --help
```

本机为 Apple Silicon，无 NVIDIA/CUDA。允许完成结构、化学、适配器、LibInvent CPU/MPS 及所有非 GNINA/非 MD 检查；Peptide2Mol 若无法在 MPS/CPU 正确推理，将明确报告为平台受限，不生成替代候选。

## 入口

```bash
python v2/murj_pipeline.py --run-id murj_9nu4_pilot prepare
python v2/murj_pipeline.py --run-id murj_9nu4_pilot hypotheses
python v2/murj_pipeline.py --run-id murj_9nu4_pilot generate-seeds --budget 20
python v2/murj_pipeline.py --run-id murj_9nu4_pilot status
```

模型阶段会校验 checkpoint 哈希、官方仓库提交、输入和命令，并保存 stdout/stderr。错误哈希、缺权重、非法输入或模型错误均非零退出。

## 正式导出门槛

主表 `murj_predocking_library.csv` 只允许同时满足以下条件的真实候选：

1. 真实 Peptide2Mol 种子与真实 LibInvent 扩展事件可追溯；
2. RDKit 可解析、单组分、无 dummy/query atom，标准化身份唯一；
3. 路线状态为完整支持（`supported_forward` 或经规则允许的完整 `supported_retro`）；
4. 至少一个状态/构象满足一个局部假设并通过严重碰撞硬门槛；
5. 所有外键、SDF 化学图、形式电荷、键级、立体化学和坐标通过全量回读。

数量不足时交付真实合格子集与缺口原因，不复制分子、不把构象当新化合物。

## 环境位置与当前状态

- 项目位于外置盘 exFAT。conda 环境含绝对路径和符号链接，留在 `~/Desktop/01_molecule_generation/envs`，可用 `MURJ_ENV_ROOT` 覆盖。
- 测试：`PYTHONPATH=v2 <envs/reinvent4>/bin/python -m pytest v2/tests -q`（structure、cores、adapters、chemistry 已通过；geometry 测试未写）。
- Peptide2Mol 在 macOS CPU 上可启动采样（约 2.2 s/步 x 1000 步），但本次未完成采样验证；未跑 REINVENT4 采样；未做 GNINA/MD。
- 经 Hydra 传含中文的路径必须加引号：`"ckpt_path='...'"`。适配器尚未包含该修复。
- v1 旧代码的已知缺陷（SDF 写出、锚点计数、解析器）尚未修复。
