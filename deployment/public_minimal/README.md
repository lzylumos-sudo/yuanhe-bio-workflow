# Linux / WSL2 GPU deployment: public minimal acceptance

This entrypoint recreates an independent environment and runs **1CRN + ethanol -> real Vina docking -> screening/export -> HF/RESP + GAFF2 parameterization -> short CUDA MD -> analysis and native checkpoint audit**.

This is command-line engineering acceptance, not acceptance of every platform module. Generation models, browser interactions, membrane systems, production-duration MD, MM/PBSA and MM/GBSA are outside this run. Ethanol is a public engineering control, not evidence of affinity or stability. Formal screening may qualify zero candidates: the explicit engineering MD control does not override scientific screening gates.

## Target computer

- Linux x86_64, or Linux under Windows WSL2. Locally validated on Ubuntu 22.04 / RTX 5080.
- An installed NVIDIA driver with the GPU accessible inside Linux/WSL. The lock includes CUDA 13.4 user-space libraries; the target driver must be compatible. The runner creates a CUDA context and stops on failure rather than silently using CPU MD. It does not install drivers or change system environments.
- Git, Bash, micromamba, and network access to conda-forge, PyPI, GitHub, RCSB and Molekula. Follow the [official micromamba installation instructions](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html); shell initialization is not required.
- Several GB of storage for scientific dependencies. Minimal refers to computation, not download size.

## Install and run

Run these commands in a Linux / WSL terminal:

```bash
git clone https://github.com/lzylumos-sudo/yuanhe-bio-workflow.git
cd yuanhe-bio-workflow

# If micromamba is not on PATH, set its actual executable path first:
# export MICROMAMBA=/absolute/path/to/micromamba
bash deployment/public_minimal/install.sh "$HOME/yuanhe-public-env"
bash deployment/public_minimal/run_public.sh "$HOME/yuanhe-public-env" "$HOME/yuanhe-public-run-001"
```

Installation requires a new absolute environment path; the first run requires a new absolute output path. Installation failures retain diagnostic files: resolve the cause and retry with a new environment path. Existing research environments are not installed into, activated or modified. Caches, quantum-chemistry scratch files, public source snapshots and results remain separate from research projects.

## Resume an interrupted run

For a run created by this version, use the same checkout, environment, GPU selection and output path:

```bash
bash deployment/public_minimal/run_public.sh --resume "$HOME/yuanhe-public-env" "$HOME/yuanhe-public-run-001"
```

- `RUN_STATE.json` commits completed stages and SHA-256 artifact manifests atomically. Resume verifies the code/configuration snapshot, interpreter, installed package versions, Conda metadata and GPU identity before reusing work. Changed or missing committed inputs fail closed.
- Completed stages are verified and skipped. Failed download, preparation, docking, screening, build or analysis attempts are archived under `.recovery/interrupted-*`, then retried from the last committed stage. Upstream scientific artifacts must still match their recorded hashes. Recovery never upgrades old runs without `RUN_STATE.json`.
- Interrupted MD retains its replica workspace and delegates continuation to the existing OpenMM backend. That backend verifies the system, committed trajectory segments and checkpoint before continuing; incompatible checkpoints fail without reseeding. An uncommitted segment is replayed from the last committed checkpoint. Before the first checkpoint, MD starts from the unchanged built system.
- SIGINT/SIGTERM stops the active phase process group and records interruption. A per-run OS lock prevents concurrent writers; the phase subprocess inherits the lock so a killed parent cannot immediately launch a competing resume. SIGKILL cannot record a final status; the next invocation uses the persisted active-stage record.
- Failed logs and metadata are retained. A successful finish archives stale `FAILED.json`; a rejected resume does not overwrite failure evidence. Snapshot backups require additional disk space and are intended for this tiny public control, not long production trajectories.
- Resume reuses the original supplier snapshot and its actual fetch date; it does not refresh inventory or certify current stock. Create a new run for fresh availability evidence.
- Dependency/CUDA preflight happens before a run record is created. If that preflight fails and no `RUN_STATE.json` exists, fix the environment and use a new output path. Resume does not install dependencies or repair an environment.

Recovery regression tests use injected failures and real OS subprocess/signal/locking checks. They do not establish CUDA checkpoint continuation on a target machine. The historical `VALIDATION.json` predates this recovery change; repeat interrupted CUDA acceptance on Linux/WSL before claiming that new path is GPU-validated.

Use `CUDA_VISIBLE_DEVICES=0` to select one GPU if needed. This is a single-GPU acceptance run, not multi-GPU scaling validation. Vina uses one CPU thread; MD uses OpenMM CUDA. The environment includes CPU PyTorch dependencies but does not run generation models.

## Pinned dependencies and sources

- `linux-64.lock`: exact build URLs and checksums for 300 Conda packages, without a fresh dependency solve.
- `pip-requirements.txt`: 16 pip packages pinned by version or commit. RESP uses an exact upstream commit and source SHA-256. Installation requires `pip check` to pass.
- Main versions: Python 3.11.17, NumPy 1.26.4, SciPy 1.15.2, Biopython 1.85, OpenMM 8.6.1, AmberTools 26.0, Psi4 1.11, RDKit 2026.03.1, Meeko 0.8.0, Vina 1.2.7.
- Downloads of public 1CRN and Vina are checked against fixed SHA-256 values. Each run fetches a fresh supplier page, verifies ethanol identity and InStock data, and records fetch time and hash. Unavailable or changed pages, or missing stock evidence, fail closed. No stock evidence is fabricated and no purchase is made.
- Third-party binaries, supplier pages, model weights and private molecules are not bundled in Git. Dependencies retain their upstream licenses; see [research dependency licenses](../../docs/16_RESEARCH_DEPENDENCY_LICENSES.md).

## Minimal settings and success criteria

| Stage | Setting / requirement |
| --- | --- |
| Docking | One molecule, one pose, one CPU thread, exhaustiveness 1, seed 17 |
| Region | A 10 angstrom cube at the 1CRN surface; an engineering region, not a validated binding pocket |
| Parameterization | Real HF/6-31G* optimization and two-stage RESP, GAFF2, one thread |
| MD | ff14SB / TIP3P; at most 100 minimization iterations; NVT 0.1 ps; NPT 0.1 ps; production 0.2 ps (200 steps at 1 fs) |
| Analysis | Four production frames; RMSD, RMSF, contacts, hydrogen bonds and SASA; passing native audit with checkpoint reload |

Success creates `COMPLETED.json` with status `PUBLIC_ENGINEERING_PATH_PASSED`. Inspect `pip-check.txt`, `gpu.json`, `md_audit.json`, phase logs and `artifact_hashes.json`. Failures create `FAILED.json`; partial artifacts alone do not demonstrate completion. See [recorded local validation](VALIDATION.json).

Passing validates this public computational path on the target machine. Formal research still requires reviewed proteins, pockets, chemical states and scientific settings in a separate research task. Do not use the 0.2 ps settings to support scientific conclusions. The original research default configuration is not modified by this acceptance run.
