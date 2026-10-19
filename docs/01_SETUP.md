# Installation and configuration

The recorded compute environment is Linux x86-64 with an RTX 5090 D; final release acceptance remains open. Windows/WSL2 with RTX 5070 is a planned deployment target, not a validated platform. See [the support matrix](12_SUPPORT_MATRIX.md) for the exact boundary. Keep environments on the WSL Linux filesystem, not an exFAT external disk. The external disk carries the project, structures, model files and run records. A new host still needs compatible installed dependencies.

## Preview

With Python 3.11+, run `python3 00_workbench/run_webui.py --no-browser` from the repository root and open `http://127.0.0.1:8777`. Starting the interface does not start calculations. Historical result datasets are not bundled; an empty historical page is expected.

On Windows, use `start_windows.cmd` and select the WSL distribution. The launcher reads the chemistry interpreter from `deployment/runtime.json`. Daily preparation, review, execution, recovery and export are available in the workbench.

## Separate environments

See [environment specifications and reconstruction](14_ENVIRONMENTS.md) for the purpose, platform and limitations of each specification and lock. Installed environments stay on the compute host; the repository stores reconstruction metadata.

- Chemistry and screening: Python, RDKit, NumPy, PyYAML, gemmi and the screening package dependencies. See `02_docking/env/sbvs2.environment.yml` and `03_virtual_screening/pyproject.toml`.
- Peptide2Mol: its own Python/PyTorch/CUDA and geometric-operator environment. RTX 5070 compatibility must be tested with the actual operators, not inferred from GPU detection.
- REINVENT4: a separate environment with the pinned source and LibInvent prior.
- GNINA: a separately installed executable and compatible CNN runtime. Record actual CNN model, device and version.
- Membrane preparation and analysis: `deployment/membrane-environment.yml`, installed through the interface after the base environment is configured. The first successful solve exports an exact package lock and pip freeze; the YAML alone is not a solved lock.
- GROMACS: Linux build selected by its actual executable path. Record whether execution really used CPU or GPU.
- Optional ZLIS and standalone OpenMM: follow their module READMEs. Their environments and evidence do not validate GROMACS.

See `deployment/validation_5070/INSTALL_AND_VALIDATE.txt` for installation commands and target-machine experiments, and `USER_GUIDE.txt` for the UI workflow.

## Restore models

From the repository root:

```bash
python3 deployment/restore_models.py --source-root /path/to/source-project
```

The restorer accepts this repository layout or a source directory containing the pinned model files. It copies only model directories into ignored `01_molecule_generation/vendor/`, excludes Git metadata and secrets, refuses overwrite and checks both weight hashes. It does not install dependencies or run inference. Without upstream Git metadata, adapters record a source fingerprint rather than claiming verified upstream commit identity.

Alternatively obtain the pinned source and weights listed in `deployment/model_dependencies.json`. Upstream sources may not reproduce historical local patches; validate actual inference separately. Missing weights, wrong hashes and empty valid outputs must fail explicitly.

## Runtime settings

Copy `deployment/runtime.example.json` to ignored `deployment/runtime.json`, then enter actual paths in Environment and diagnostics. Save before creating a new run. Existing runs retain frozen settings.

`python3 deployment/doctor.py --config deployment/runtime.json` is an optional read-only check. A passing dependency check is not scientific acceptance. Select `data/structures/9nu4.pdb`, review chains and residue ranges, then review cores, receptor, docking box and membrane system at the relevant stages.

## Engineering checks

Run each module's tests in its own directory to avoid collisions between packages named `src`:

```bash
cd 01_molecule_generation
PYTHONPATH=. python -m pytest tests -q
cd ../00_workbench
PYTHONPATH=. python -m pytest tests/test_workflow_regressions.py tests/test_studio_integration.py -q
```

See `04_VALIDATION.md` for the actual acceptance boundary. Windows launch, environment solving, model inference, GNINA CNN and membrane dynamics require target-machine validation.
