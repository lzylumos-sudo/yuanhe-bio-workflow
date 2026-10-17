# GNINA 自动化筛选平台 v2（archive/docking_prototype）

依据《GNINA自动化筛选平台_完整实施方案》（docs/ 内附原 PDF）构建。
旧项目 `../02_docking` **只读保留**，其 2259 个文件的 SHA-256 清单见 `reports/baseline/legacy_01_manifest.json`（不含 26G runtime/）。

## 快速开始
```bash
export SBVS_PYTHON=$HOME/.sbvs2_env/bin/python   # RDKit/OpenMM/PDBFixer/Vina/Meeko，见 env/
./run.sh doctor
./run.sh new-project 我的项目
./run.sh ingest-receptor some.pdb --out receptor_report.json
./run.sh ingest-ligands ligands.csv --out ligands_report.json
./run.sh estimate --states 132 --conformers 3 --receptors 1 --seeds 3
cd source && python -m pytest -q tests
```

## 目录
- `source/sbvs2/` 代码：`paths` `hashing` `states` `approvals`（硬门禁）`receptor/ingest` `ligands/{csv_io,prepare,review}` `docking/{box,protocol,engines}` `analysis/{qc,ranking}` `project` `cli`
- `source/tests/` 单元测试
- `projects/<id>/{00_original,01_preparation,runs,exports}`；run 内 `00_run_record 01_logs 02_work_products 03_results 04_analysis`，run 永不覆盖
- `reports/baseline|smoke` 基线与真实数据冒烟证据
- `configs/default.yaml` [P] 类参数起点

## 真实重对接示例
`python scripts/redock_fol.py` 新建项目并跑 1DRF/FOL；`python scripts/redock_fol.py <project_id> <run_id>` 恢复中断的 run（复核哈希，不覆盖）。
