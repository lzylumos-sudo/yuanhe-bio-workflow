# Environment specifications and reconstruction

Keep environment definitions in Git; install the environment itself on each compute host. Installed Conda prefixes and Python virtual environments contain platform-specific binaries and often absolute paths. Copying their directories is not a general migration method. The repository excludes installed environments; no `site-packages` or `conda-meta` tree is tracked.

## What each file guarantees

| File | Purpose and limits |
|---|---|
| `02_docking/env/sbvs2.environment.yml` | V2 direct-dependency specification. A solver input, not a complete lock or proof of cross-platform compatibility. |
| `02_docking/env/sbvs2.linux-64.explicit.lock` and `sbvs2.osx-arm64.explicit.lock` | Platform-specific Conda package URLs and recorded hashes. Select exactly the target OS/architecture. Add `sbvs2.pip.lock` separately. |
| `02_docking/env/environment.yml` | Historical legacy Linux full-build snapshot. It is distinct from V2; do not combine the two. The local prefix has been removed. |
| `02_docking/env/conda-linux-64.lock` and `requirements-pip.lock` | Legacy Linux Conda URLs plus separate pip supplements. Conda URL hashes were not recorded; pip versions are pinned without artifact hashes. |
| `03_virtual_screening/pyproject.toml` | Standard Python package dependency ranges; intentionally not an exact environment lock. |
| `04_interaction_analysis/env/environment.yml` | ZLIS solver specification with several unpinned dependencies. A fresh solve can select different versions. See the module README for isolated PyMOL rendering. |
| `05_molecular_dynamics/environment.yml` and `env/linux-64-explicit.lock` | Standalone OpenMM specification and Linux URL snapshot. The URL snapshot has no recorded artifact hashes. `requirements.txt` is only a pip subset, not the complete scientific stack. |
| `deployment/membrane-environment.yml` | Separate membrane environment candidate. A successful target solve must produce its own lock and validation evidence. |
| `deployment/validation_5090/pip-freeze.txt` | Recorded Python package versions used by engineering CI, with additional explicit CI dependencies. It does not contain GNINA, CUDA drivers, model weights or the complete GPU workflow. |
| `deployment/validation_5070/peptide-inference-requirements.txt` | Candidate Peptide2Mol dependency subset; not a validated complete PyTorch/CUDA/operator environment. |
| `deployment/model_dependencies.json` | Model source/weight identities. It is not an environment lock; REINVENT4 and Peptide2Mol still require their own compatible installations. |

## Reconstruct a recorded docking environment

Example for **Linux x86-64 V2** from the repository root, using a new environment name:

```bash
conda create -n yuanhe-docking-v2 --file 02_docking/env/sbvs2.linux-64.explicit.lock
conda run -n yuanhe-docking-v2 python -m pip install -r 02_docking/env/sbvs2.pip.lock
conda run -n yuanhe-docking-v2 python -m pip check
cd 02_docking
conda run -n yuanhe-docking-v2 python -m pip install --no-deps -e .
conda run -n yuanhe-docking-v2 python -m pytest tests -q
```

For Apple Silicon use the `osx-arm64` lock instead. These are reconstruction instructions, not a claim that both environments were freshly recreated during this review. Native engine executables, system libraries, models and target GPU checks remain separate. A solver YAML can be used with `conda env create -n NEW_NAME -f FILE`, but its output must be validated before replacing a known working environment.

## Maintain lightweight, useful records

Maintain three layers: a readable direct-dependency specification, an OS/architecture-specific resolved package record, and a successful validation report. Conda and pip records describe different package sets; neither alone necessarily reconstructs a mixed environment. Capture the Python version, GPU driver, actual engine versions and model hashes with validation results. Keep local runtime interpreter paths in ignored `deployment/runtime.json`.

From an activated, validated Conda environment, export a portable starting specification and a same-platform package record:

```bash
conda env export --from-history > environment.from-history.yml
conda list --explicit --md5 > conda-explicit.txt
python -m pip freeze > pip-freeze.txt
python -m pip check
```

Review exports before committing: remove a local `prefix`, give the environment a meaningful name, and exclude credentials/local editable paths. A full pip freeze is an inventory and can duplicate Conda-managed packages; do not blindly install it over a Conda scientific stack. Maintain a separate list of pip-only additions. `--md5` preserves Conda artifact identity; hash-pinned pip artifacts require an additional platform-specific locking/download workflow. Do not invent missing hashes or claim current URL snapshots are fully hash-locked.

## When an old installation can be removed

First recreate a replacement in a separate location, verify imports and package consistency, then run a bounded real workflow for the engines actually used. Update and test runtime paths, retain package/engine/model records, and verify that no active run or other project depends on the old interpreter. Only then remove an obsolete installation. Model weights, patched source, research data and checkpoints are not disposable environments. If offline recovery matters, retain a verified package cache or platform-specific environment archive as well; a text specification still needs downloadable packages.

This repository review did not inspect or remove disk environments and did not perform a fresh complete GPU environment reconstruction.

## References

- [Conda environment management](https://docs.conda.io/projects/conda/en/stable/user-guide/tasks/manage-environments.html)
- [Conda portable history exports](https://docs.conda.io/projects/conda/en/stable/commands/env/export.html)
- [pip secure and hash-checked installs](https://pip.pypa.io/en/stable/topics/secure-installs/)
