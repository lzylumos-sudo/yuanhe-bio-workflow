# Environment locks

- `sbvs2.environment.yml`: maintained constraints; historical tested combination Python 3.12.15, RDKit 2026.03.1, OpenMM 8.6.1, PDBFixer 1.12 and Vina 1.2.7.
- `sbvs2.linux-64.explicit.lock`: actual Ubuntu 22.04 x86-64 solve with MD5; URLs normalized to conda.anaconda.org while retaining package identity/checksums.
- `sbvs2.osx-arm64.explicit.lock`: historical Mac development environment.
- `sbvs2.pip.lock`: pip dependencies. Meeko 0.8.0 is required here; older Conda Meeko versions failed against RDKit 2026.x because rdkit.six was absent.

Rebuild with `micromamba create -p ENV --file sbvs2.linux-64.explicit.lock`, then `ENV/bin/pip install --no-deps -r sbvs2.pip.lock` after confirming the Conda dependencies are satisfied.

Legacy `conda-linux-64.lock` and `requirements-pip.lock` describe a separate Python 3.11 WSL environment. Their preservation does not certify them on the current host.
