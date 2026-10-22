# Linux / WSL2 GPU deployment: public minimal acceptance

This entrypoint recreates an independent environment and runs **1CRN + ethanol -> real Vina docking -> screening/export -> HF/RESP + GAFF2 parameterization -> short CUDA MD -> analysis and native checkpoint audit**.

This is command-line engineering acceptance, not acceptance of every platform module. Generation models, browser interactions, membrane systems, production-duration MD, MM/PBSA and MM/GBSA are outside this run. Ethanol is a public engineering control, not evidence of affinity or stability. Formal screening may qualify zero candidates: the explicit engineering MD control does not override scientific screening gates.

## Target computer

- Linux x86_64, or Linux under Windows WSL2. Locally validated on Ubuntu 22.04 / RTX 5080.
- An installed NVIDIA driver with the GPU accessible inside Linux/WSL. The lock includes CUDA 13.4 user-space libraries; the target driver must be compatible. The runner creates a CUDA context and stops on failure rather than silently using CPU MD. It does not install drivers or change system environments.
- Git, Bash, micromamba, and network access to conda-forge, PyPI, GitHub, RCSB. Follow the [official micromamba installation instructions](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html); shell initialization is not required.
- Several GB of storage for scientific dependencies. Minimal refers to computation, not download size.

## Install and run

Run these commands in a Linux / WSL terminal:

```bash
git clone https://github.com/lzylumos-sudo/yuanhe-bio-workflow.git
cd yuanhe-bio-workflow

# If micromamba is not on PATH, set its actual executable path first:
# export MICROMAMBA=/absolute/path/to/micromamba
bash deployment/public_minimal/install.sh "$HOME/yuanhe-public-env"
# Install real ADMET/SCScore/AiZynth assets as described in ../../02_docking/docs/PREDOCKING_V2.md
export YUANHE_PREDOCKING_CONFIG=/absolute/path/predocking-assets/predocking-config.json
bash deployment/public_minimal/run_public.sh "$HOME/yuanhe-public-env" "$HOME/yuanhe-public-run-001"
```

The environment and output paths must be new absolute paths; existing paths are rejected. Installation failures retain diagnostic files: resolve the cause and retry with a new environment path. For run failures, inspect the logs and retry with a new output path. Automatic resume is not implemented. Existing research environments are not installed into, activated or modified. Caches, quantum-chemistry scratch files, public source snapshots and results remain separate from research projects.

Use `CUDA_VISIBLE_DEVICES=0` to select one GPU if needed. This is a single-GPU acceptance run, not multi-GPU scaling validation. Vina uses one CPU thread; MD uses OpenMM CUDA. The environment includes CPU PyTorch dependencies but does not run generation models.

## Pinned dependencies and sources

- `linux-64.lock`: exact build URLs and checksums for 300 Conda packages, without a fresh dependency solve.
- `pip-requirements.txt`: 16 pip packages pinned by version or commit. RESP uses an exact upstream commit and source SHA-256. Installation requires `pip check` to pass.
- Main versions: Python 3.11.17, NumPy 1.26.4, SciPy 1.15.2, Biopython 1.85, OpenMM 8.6.1, AmberTools 26.0, Psi4 1.11, RDKit 2026.03.1, Meeko 0.8.0, Vina 1.2.7.
- Downloads of public 1CRN and Vina are checked against fixed SHA-256 values. The supplied model configuration runs the same toxicity, SA/SCScore and computed-route gate as production. Missing models or incomplete routes stop the run.
- Third-party binaries, model weights and private molecules are not bundled in Git. Dependencies retain their upstream licenses; see [research dependency licenses](../../docs/16_RESEARCH_DEPENDENCY_LICENSES.md).
- The prior recorded GPU validation predates the synthesis-only policy. This updated full CUDA path requires new target-host acceptance; CPU route/Vina acceptance does not establish CUDA MD validation.

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
