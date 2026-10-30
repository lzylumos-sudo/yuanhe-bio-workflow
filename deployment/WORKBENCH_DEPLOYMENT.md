# Administrator deployment: complete research workbench

Start with the [self-contained installation entrypoint](../DEPLOY.md). It covers first installation without an old project, configuration, preflight, startup and upgrades. This reference expands individual installation steps. Use a complete clean Git checkout or the verified named administrator ZIP. The public minimal CLI test, historical 5070 notes, and `Dockerfile.research` alone do **not** install or validate the full workbench.

The administrator performs installation and service management. Researchers operate only the supplied UI URL; they do not need SSH or a local compute installation. A code update requires the administrator to deploy the updated revision and restart the service. Merely refreshing the browser cannot change server dependencies or its running Python.

## 1. Prerequisites and deployment contract

- Linux x86_64, natively or under WSL2, with Bash, Git, Python 3, curl, patch, compiler/build tools and CA certificates. Environments and run outputs must be on a writable Linux filesystem.
- Administrator-assigned NVIDIA GPU(s), a compatible host driver, and GPU access inside WSL/container. Never install a Linux display driver inside WSL. The shared runtime lock contains CUDA 13.4 user libraries; Peptide2Mol uses PyTorch cu128. Validate both runtimes on the actual driver. Do not infer compatibility from the GPU name.
- Install [micromamba from its official source](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html); set `MICROMAMBA` if it is not on PATH. The installer does not change drivers, firewall, system Python, unrelated Conda environments or scheduler allocations.
- Network access for approved sources: GitHub/codeload, conda-forge, PyPI, download.pytorch.org, data.pyg.org and the upstream model distribution sites. For restricted networks the administrator must provide equivalent verified mirrors/caches and preserve hashes, versions and licenses. Failure to download must stop installation.
- Storage for environments, source, checkpoints, molecules and trajectories **plus** `research.disk_reserve_gb` (default 50 GiB). A 50 GiB reserve is not the total storage estimate. Measure trajectory growth on the target system before scheduling two 500 ns jobs; configure `estimated_md_gb_per_candidate` accordingly. Preserve checkpoint/output space throughout execution.
- Keep the service on the administrator-approved private network. The built-in service has no user authentication and is not a public multi-tenant endpoint. Reuse the site's existing network/access controls.

Record the Git commit, OS, driver, GPU UUID/model/VRAM, service account, project path, environment prefix and output filesystem. Avoid changing the checkout under an active job. Preserve `deployment/runtime.json`, licensed assets and `00_workbench/workflow_runs/` when redeploying; do not overwrite frozen run records. Changed source can invalidate an existing task's source snapshot: create a new task or use an applicable verified continuation workflow.

## 2. Dependency map

| Component | Required runtime / assets | Check before delivery |
| --- | --- | --- |
| HTTP UI and drawing | RDKit, NumPy, RDKit drawing in the **service Python itself** | Startup import check and `ui_interpreter` diagnostic |
| Chemistry, risk, docking adapter | RDKit FilterCatalog, NumPy/SciPy/pandas, PyYAML, gemmi, Meeko 0.8.0, filelock | Imports, real risk-rule construction/examples, Vina version |
| Screening | MDAnalysis, ProLIF, PoseBusters, SciPy/pandas, FastAPI and package dependencies | Probe the actual `screen_python` |
| Receptor preparation | PDBFixer, OpenMM and Amber14/TIP3P data | Imports and force-field loading in `prep_python` |
| Peptide2Mol | Separate Python 3.10, PyTorch 2.7.0/cu128, matching PyG extensions, pinned requirements, patched source and exact checkpoint | CUDA/PyG operations, source patch and weight SHA-256, then real small inference |
| REINVENT4 | Separate Python 3.11, pinned source/uv.lock, patched LibInvent and exact prior | CLI/import checks, patch and prior SHA-256, then real sampling |
| Synthesis screening | AiZynthFinder 4.4.1 + model/templates/planning building-block pool | Asset/dependency inventory; per-molecule complete-route verification remains separate |
| Membrane build / HF-RESP | OpenMM CUDA, OpenFF, openmmforcefields, ParmEd, AmberTools, Psi4, RESP, RDKit | Imports, CUDA Context, force-field data, `antechamber`, `parmchk2`, `tleap` |
| Production MD | CUDA-enabled GROMACS with compatible driver/toolkit | `gmx --version` must report CUDA; actual small membrane mdrun and checkpoint recovery |
| Final trajectory and PB analysis | MDAnalysis, NumPy, matplotlib, ParmEd, AmberTools `MMPBSA.py`, `sander`, `cpptraj` | Tools in the MD environment PATH; real analysis of the acceptance trajectory |

The new installer uses the repository's exact Linux Conda lock and pinned pip supplement for a shared chemistry/preparation/MD runtime, then installs the local screening package without re-solving dependencies. The lock was reconstructed for the public docking-to-OpenMM acceptance path; its use here does not claim full browser/model/membrane acceptance on a new machine. Model environments remain separate. Optional standalone state enumeration uses molscrub 0.2.2; optional ZLIS/PLIP/PyMOL workflows retain their own module installation instructions and are outside the eight-stage UI pipeline.

## 3. Get the project and model assets

```bash
git clone https://github.com/lzylumos-sudo/yuanhe-bio-workflow.git
cd yuanhe-bio-workflow
git rev-parse HEAD
```

Use clean `main` for the latest merged fixes, recording the exact commit before installation. To reproduce a published release instead, select its immutable tag from `deployment/current_release.json`; it can predate main. For the named ZIP, follow `DEPLOY.md` and `verify_release.py`; do not mistake a source commit hash for the product version. Prefer the versioned site configuration and `setup_workbench.py` in `DEPLOY.md`. For the manual commands below, choose a new dedicated installation prefix owned by the service account, such as `/srv/yuanhe/runtime-NEW`; replace `YUANHE_PREFIX` once, consistently.

GitHub deliberately does not contain model weights, installed environments or private scientific data. For first installation, `python3 deployment/fetch_generation_models.py --destination /absolute/NEW-model-bundle --execute` obtains the exact public sources and weights, retains their licenses and verifies hashes. Review upstream terms linked in `DEPLOY.md`. Use that bundle with `setup_workbench.py`, which restores it automatically. For the manual steps below, restore the downloaded bundle or an existing compatible bundle:

```bash
python3 deployment/restore_models.py --source-root /absolute/model-bundle
```

The restorer copies only verified weights. It reconstructs vendor source from the pinned Git commits, using cached committed objects when available and fetching public source otherwise. Old working-tree modifications, environments and wrappers are not reused. Current patches are applied in the new source directory. For manual restoration, check out the exact upstream revisions into `01_molecule_generation/vendor/Peptide2Mol` and `01_molecule_generation/vendor/REINVENT4`, then place the two weights at the manifest paths. Retain the upstream LICENSE files and REINVENT4 `uv.lock`. Do not substitute a similarly named model or edit expected hashes to accept a different weight.

Run `python3 deployment/apply_vendor_patches.py` after restoring both sources and **before** installing REINVENT4. It verifies the known before/after hashes and refuses unknown revisions. The patches preserve bounded sampling and chemical stereochemistry. The installer repeats this idempotent check. Missing assets are reported by diagnostics; they are not silently replaced with historical molecules.

## 4. Install the shared runtime and model environments

```bash
export YUANHE_PREFIX=/srv/yuanhe/runtime-NEW
# If needed: export MICROMAMBA=/absolute/path/to/micromamba
bash deployment/install_workbench.sh runtime "$YUANHE_PREFIX"
bash deployment/install_workbench.sh peptide "$YUANHE_PREFIX"
bash deployment/install_workbench.sh reinvent "$YUANHE_PREFIX"
bash deployment/install_workbench.sh synthesis "$YUANHE_PREFIX"
bash deployment/install_workbench.sh vina "$YUANHE_PREFIX"
```

Run one step at a time and stop on any nonzero exit. Runtime/model installation refuses existing target environments to avoid accidental upgrades. After a failed installation preserve the log, resolve its cause and choose a new prefix; do not treat a partially populated directory as success. No step runs a target experiment. Logs and package freezes are in `$YUANHE_PREFIX/evidence/`.

The shared runtime contains CPU PyTorch for unrelated dependency resolution; **do not point `peptide_python` to it**. Peptide2Mol uses its dedicated GPU environment. Match PyTorch/CUDA/PyG builds as described in the [official PyG installation guide](https://pytorch-geometric.readthedocs.io/en/latest/install/installation.html). The repository pins an older supported combination for the recorded model; do not replace it with latest packages without testing the model.

The REINVENT4 environment is installed from the patched source with its frozen lock, followed by [pinned runtime supplements](reinvent-runtime-requirements.txt). The recorded dependency patches select official PyTorch 2.12.0+cu130 and torchvision 0.27.0+cu130 with their locked CUDA dependencies. The upstream cu126 build lacks Blackwell sm_120 kernels and failed the original installation check on an RTX 5080. Both dependency files have verified before/after hashes; unrelated package versions and sources remain unchanged. The pinned upstream lock omits SciPy even though the CLI imports it. The installer adds this supplement with `--no-deps`, then checks the CLI and dependency consistency. Readiness still requires real CUDA operations on the target GPU; installing the new build alone is not acceptance. Supply a compatible NVIDIA driver. Do not silently substitute CPU execution or change packages in a failed environment; preserve the evidence and deploy the corrected source into a new prefix.

## 5. Install or select GPU GROMACS

Set an actual absolute path to a CUDA-enabled GROMACS binary supplied by the administrator. A distribution CPU package is not sufficient for this GPU workflow. Follow the [official GROMACS build guide](https://manual.gromacs.org/current/install-guide/index.html); validate compiler, CUDA toolkit and GPU architecture together.

For a portable source build use `bash deployment/build_gromacs.sh /absolute/new-prefix /absolute/cuda 4`. It verifies the official archive, builds CUDA support and runs upstream checks. It requires CMake >= 3.28, make, curl and compatible C/C++ compilers. No drivers or site services are changed.

For the repository's recorded `/opt/yuanhe` layout, the optional helper is `deployment/validation_5070/build_gromacs_gpu.sh`. It requires the registered project host, CUDA Toolkit 12.8, GCC/G++ 12 and the base tools installation. It installs GROMACS 2026.4 into a separate prefix and runs its regression target. This helper is not required when the administrator already provides a compatible build, and it must not be run unmodified against an unrelated cluster installation.

```bash
export YUANHE_GMX=/absolute/path/to/gromacs/bin/gmx
"$YUANHE_GMX" --version
bash deployment/install_workbench.sh configure "$YUANHE_PREFIX"
```

`configure` writes actual interpreter/executable paths to ignored `deployment/runtime.json`, backing up any old file under the installation evidence directory. The shared runtime serves `chem_python`, `screen_python`, `prep_python` and `md_python`; model paths are separate. Configuration uses current synthesis-only defaults and the asset receipt from the synthesis step. Incompatible old research configuration is rejected without overwriting it. Configure is not a readiness pass.

## 6. Configure synthesis models before generation

The `synthesis` installer step installs both isolated environments and pinned screening assets before `configure`. For custom or existing environments use [the current screening guide](../02_docking/docs/PREDOCKING_V2.md) to install isolated ADMET-AI and AiZynthFinder environments, official SCScore weights, expansion templates and the planning building-block pool. Import the generated config in Pre-Docking Screening. New defaults have no supplier/catalog fields or purchase bypass. Verify a small real search before releasing a host configuration. Missing models cannot qualify candidates.


## 7. Readiness and service startup

```bash
bash deployment/install_workbench.sh check "$YUANHE_PREFIX"
# All required checks must pass before the service is handed to researchers.
bash start_linux.sh
# On the site's already-approved private interface, supply its binding explicitly:
# bash start_linux.sh --host ADMIN_APPROVED_PRIVATE_BIND_ADDRESS --port 8777
```

The launcher and direct `run_webui.py` entrypoint select the configured chemistry Python and check the service's own RDKit/NumPy/drawing imports. `YUANHE_UI_PYTHON` is an explicit administrator override. A wrong interpreter now produces a startup error with the deployment guide, instead of an apparently healthy webpage that fails on task creation. Changing environment paths in the UI does not replace a running service interpreter: a restart is required.

A service manager should use the repository as working directory, the approved service account, `bash /absolute/project/start_linux.sh` as the entrypoint, absolute runtime paths, and the site's existing logging/access controls. Keep it running independently of a terminal connection. Do not auto-update the checkout or restart the compute host during a trajectory.

The UI now provides **Environment and diagnostics -> Check saved environment before creating a run**, plus a selected-run check and the latest report. It runs independently of task creation and writes `installation_evidence/readiness_*.json`. The CLI shares the same checks:

```bash
python3 deployment/doctor.py --config deployment/runtime.json --out installation_evidence/readiness.json
```

Use `--ui-python /actual/service/python` if the service uses an explicit override. CLI diagnosis can run under standard Python even if scientific dependencies are missing; probes run in the configured environments. Missing tools, imports, weights, patches, current evidence, CUDA support and disk reserve are distinct failures. The first preparation stage repeats the full readiness check before generation, saving `deployment_readiness.json` in the run. The check cannot guarantee future disk availability, numerical convergence, service uptime or route completeness; those remain execution-time checks.

## 8. Required acceptance before the administrator delivers the URL

1. Run the dependency inventory with the same account, paths and GPU allocation as the service. Retain the report and exact package locks/freezes.
2. Create a new task through the browser; inspect frozen settings and actual input structure. Verify diagnostics also work with no task selected. No unhandled RDKit import error is acceptable.
3. In a separate explicitly small acceptance task, run real structure preparation, Peptide2Mol inference, reviewed core expansion, risk/synthesis gates, Vina, screening, reviewed membrane build, equilibration, short production and analysis/export. Keep failed/zero-qualified outcomes explicit. The public minimal CLI test is supplementary and does not replace this UI/membrane test.
4. Validate checkpoint interruption/recovery and actual GPU GROMACS execution. Check final analysis can find AmberTools after the service starts without an activated shell environment. A successful import or `gmx --version` alone does not establish these properties.
5. Preserve input/configuration/source hashes, logs, manifests, model provenance, environment evidence, candidate counts, trajectories/checkpoints and exports. Only then start the requested research task.

For the requested quarter-size research batch: Peptide2Mol 1,250 attempts; up to 50 cores; REINVENT4 50,000 total attempts; generated/external unique caps 12,500 each; initial cap 25,000; refine cap 500; review cap 13 (integer rounding); **MD limit 2, one 500 ns trajectory each, total 1,000 ns**. Keep seeds, force fields, timestep, temperature, quality criteria and the 400-500 ns analysis window unchanged. These are requested limits, not promised qualified counts. Hardware concurrency remains subject to administrator allocation and measured capacity. The repository default remains 8 candidates x 500 ns; set the quarter counts in the new task or in explicit administrator research defaults.

## 9. Failure handling and redeployment

| Symptom | Administrator action |
| --- | --- |
| `No module named rdkit` at task creation or drawing | Confirm actual service Python; deploy this fix, configure `chem_python`, restart through the checked launcher |
| Diagnostics ask for a task before checking | Deploy the independent diagnostics endpoint/UI; no placeholder research run is required |
| UI path is correct but process still uses old environment | Restart the service; saving a JSON path cannot replace the current interpreter |
| Model file/patch check fails | Restore pinned licensed assets and apply known patch; never change expected hashes to hide mismatch |
| `torch_scatter`/CUDA kernel fails | Match PyTorch/CUDA/PyG wheels and validate actual operators on the host |
| Synthesis gate blocks all candidates | Inspect model paths, route budget, terminal blocks and chemistry review outcomes |
| `MMPBSA.py`, `antechamber` or `tleap` not found | Check the complete MD environment; the worker now prepends its bin directory and sets AMBERHOME when applicable |
| GROMACS is CPU-only | Configure a tested CUDA build and confirm actual GPU mdrun |
| Installation interrupted | Keep logs/partial environment for diagnosis; recreate into a new prefix rather than upgrading arbitrary packages in place |
| Research task interrupted | Inspect task state/checkpoint and use verified recovery; do not overwrite outputs or launch duplicates |

Repository regression tests verify code contracts. This deployment revision must still pass the actual administrator-machine acceptance above; no 500 ns calculation was performed by publishing these changes.
