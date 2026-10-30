# Installation and configuration

**Administrators deploying the complete browser workflow:** start with [DEPLOY.md](../DEPLOY.md) and the release's versioned [site configuration](../deployment/site.example.json). Existing campus installations use the [upgrade procedure](../deployment/UPGRADE_RUNTIME.md) and [campus configuration](../deployment/site.campus-upgrade.json). The [full workbench reference](../deployment/WORKBENCH_DEPLOYMENT.md) expands individual installation steps.

The recorded compute environment is Linux x86-64 with an RTX 5090 D; final release acceptance remains open. Windows/WSL2 with RTX 5070 is a planned deployment target, not a validated platform. See [the support matrix](12_SUPPORT_MATRIX.md) for the exact boundary. Keep environments on the WSL Linux filesystem, not an exFAT external disk. The external disk carries the project, structures, model files and run records. A new host still needs compatible installed dependencies.

## Preview

After configuring the chemistry environment, run `bash start_linux.sh` from the repository root and open `http://127.0.0.1:8777`. Starting the interface does not start calculations. Historical result datasets are not bundled; an empty historical page is expected.

On Windows, use `start_windows.cmd` and select the WSL distribution. The launcher reads the chemistry interpreter from `deployment/runtime.json`. Daily preparation, review, execution, recovery and export are available in the workbench.

## Separate environments

See [environment specifications and reconstruction](14_ENVIRONMENTS.md) for the purpose, platform and limitations of each specification and lock. Installed environments stay on the compute host; the repository stores reconstruction metadata.

- Chemistry and screening: Python, RDKit, NumPy, PyYAML, gemmi and the screening package dependencies. See `02_docking/env/sbvs2.environment.yml` and `03_virtual_screening/pyproject.toml`.
- Peptide2Mol: its own Python/PyTorch/CUDA and geometric-operator environment. RTX 5070 compatibility must be tested with the actual operators, not inferred from GPU detection.
- REINVENT4: a separate environment with the pinned source and LibInvent prior.
- Docking: official AutoDock Vina 1.2.7 executable, CPU threads and seeds; Meeko 0.8.0 preserves chemical mapping. No CNN runtime is required for new jobs.
- Research QM/MD: `deployment/research-environment.yml` pins Psi4, RESP and AmberTools; see [migration and research methods](15_VINA_RESP_RESEARCH.md).
- Membrane preparation and analysis: `deployment/membrane-environment.yml`, installed through the interface after the base environment is configured. The first successful solve exports an exact package lock and pip freeze; the YAML alone is not a solved lock.
- GROMACS: Linux build selected by its actual executable path. Record whether execution really used CPU or GPU.
- Optional ZLIS and standalone OpenMM: follow their module READMEs. Their environments and evidence do not validate GROMACS.

The files under `deployment/validation_5070/` record historical installation experiments; use `DEPLOY.md` for current installation. See `USER_GUIDE.txt` for the UI workflow.

## Restore models

From the repository root:

```bash
python3 deployment/restore_models.py --source-root /path/to/source-project
```

The restorer accepts this repository layout or a model bundle. It reuses only the two hash-verified weights and reconstructs each vendor directory from the exact upstream commit in `deployment/model_dependencies.json`. Existing vendor Git objects can supply that commit offline; otherwise Git fetches the pinned public source. Old working-tree edits, installed environments, wrappers and configuration are not copied. New source metadata is retained, existing destinations are refused, and no dependencies or inference are run. Apply this release's recorded vendor patches before installation, as described in `DEPLOY.md`.

Alternatively obtain the pinned source and weights using `deployment/fetch_generation_models.py` as described in `DEPLOY.md`. Actual inference still needs separate validation. Missing weights, wrong hashes and empty valid outputs must remain explicit.

## Runtime settings

For current installations, `deployment/setup_workbench.py --config /path/to/site.json --execute` generates ignored `deployment/runtime.json` from the selected release and the new environments. `site.example.json` is the administrator input; `runtime.example.json` is a reference for manually managed environments, not a replacement for the generated runtime. Keep schema fields and all actual interpreter paths together. Save any deliberate environment changes before creating a new run; existing runs retain frozen settings. Never copy an old run configuration over a new installation.

`python3 deployment/doctor.py --config deployment/runtime.json` is an optional read-only check. A passing dependency check is not scientific acceptance. Select `data/structures/9nu4.pdb`, review chains and residue ranges, then review cores, receptor, docking box and membrane system at the relevant stages.

## Engineering checks

Run each module's tests in its own directory to avoid collisions between packages named `src`:

```bash
cd 01_molecule_generation
PYTHONPATH=. python -m pytest tests -q
cd ../00_workbench
PYTHONPATH=. python -m pytest tests/test_workflow_regressions.py tests/test_studio_integration.py -q
```

See `04_VALIDATION.md` for the actual acceptance boundary. Windows launch, environment solving, model inference, Vina docking, HF/RESP and membrane dynamics require target-machine validation.
