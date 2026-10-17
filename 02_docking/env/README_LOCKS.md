# 锁定环境 (v2)
- `sbvs2.environment.yml`: 人工维护的版本约束（已验证组合：Python 3.12.15 / RDKit 2026.03.1 / OpenMM 8.6.1 / PDBFixer 1.12 / Vina 1.2.7）。
- `sbvs2.linux-64.explicit.lock`: Ubuntu 22.04 x86_64 上实际解算结果（含 md5）。URL 已规范到 conda.anaconda.org（原始解算经过本机镜像，同名文件同 md5）。
- `sbvs2.osx-arm64.explicit.lock`: Mac 开发机实际环境。
- `sbvs2.pip.lock`: pip 部分。**Meeko 必须用 pip 0.8.0**；conda 解出的 meeko 0.5.0/0.7.1 在 rdkit 2026.x 下因 `rdkit.six` 缺失而无法导入。
- 重建: `micromamba create -p <env> --file sbvs2.linux-64.explicit.lock && <env>/bin/pip install --no-deps -r sbvs2.pip.lock`
- 旧 WSL 环境的锁 (`conda-linux-64.lock`, `requirements-pip.lock`) 原样保留，属于另一环境(Python 3.11)，**未在本次 Mac/Ubuntu 上验证**。
