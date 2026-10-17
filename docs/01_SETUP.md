# 安装与配置

## 1. 网页预览

在 Python 3.11+ 运行 `python3 00_workbench/run_webui.py --no-browser`，打开 `http://127.0.0.1:8777`。默认仅监听本机。首页读取集中结构目录；历史页面没有随仓库打包旧结果，不能把空历史页面当作本次计算失败。

## 2. 化学和科学依赖

不同模块/模型使用独立环境，不应强行安装到同一个环境。

- 生成模块工程测试：Python 3.11、pytest、RDKit、NumPy。
- 对接环境：见 `02_docking/env/sbvs2.environment.yml` 与对应平台锁文件；实际 GNINA 可执行程序另行安装。
- 筛选环境：依赖列表见 `03_virtual_screening/pyproject.toml`，可在独立环境内执行 `pip install -e ./03_virtual_screening`。
- ZLIS：见 `04_interaction_analysis/env/` 和该模块 README。
- OpenMM：x86-64 Ubuntu/WSL2，进入 `05_molecular_dynamics` 执行 `bash setup.sh`；仓库保留 `env/linux-64-explicit.lock`，不包含离线环境包。新版入口为 `bash run_md.sh`。
- Peptide2Mol、REINVENT4：分别按官方固定版本环境安装，来源与目标路径见 `deployment/model_dependencies.json`。不将模型运行环境装进仓库。

## 3. 恢复本地第三方模型

已有原始项目时，在仓库根目录运行：

```bash
python3 deployment/restore_models.py --source-root '/原项目路径/元核生物'
```

工具只复制第三方模型目录到忽略的 `01_molecule_generation/vendor/`，排除 `.git`、`.env` 和缓存；不覆盖现有目录；核对两个权重的 SHA-256。不安装依赖、不运行模型。缺少 Git 元数据时，适配器仅记录源码指纹，不声称验证了官方提交。

无原项目时，可按 `deployment/model_dependencies.json` 从官方获取固定提交与权重；官方来源不保证重现旧本地补丁，真实推理仍须独立验证。模型适配器拒绝缺失权重、错误哈希或空有效输出。

## 4. 本机配置

复制 `deployment/runtime.example.json` 为 `deployment/runtime.json`，填写各独立环境/程序的实际路径，后者已被 Git 忽略。然后运行：

```bash
python3 deployment/doctor.py --config deployment/runtime.json
```

doctor 是只读环境检查，检查通过也不等于完整推理或科学验证。网页中创建新运行时选择 `data/structures/9nu4.pdb`。必须明确链、残基范围与科学输入；核心选择、GNINA 受体/盒及已准备 GROMACS TPR 在对应阶段审核。

## 5. 工程测试

使用具备相应依赖的 Python，分别在各模块目录运行 pytest；不要从仓库根目录一次混收集多个同名 `src`/测试包。

```bash
cd 01_molecule_generation
PYTHONPATH=. python -m pytest tests -q
# 回到根目录后：
cd ../00_workbench
PYTHONPATH=. python -m pytest tests/test_workflow_regressions.py -q
```

当前验证范围和依赖差异见 `04_VALIDATION.md`。Windows 启动、Linux 安装和 GPU 真实计算需在目标机器再验收。
