# Administrator deployment entrypoint

Send this page to the administrator: [current deployment instructions](https://github.com/lzylumos-sudo/yuanhe-bio-workflow/blob/main/DEPLOY.md).
It covers generation, toxicity and synthesis screening, docking, ZVSS and membrane MD. A source download alone does not install these dependencies.

## 1. Select the current source and preserve the old installation

Use the **current administrator build** download in the [README](README.md). Download its named ZIP and `SHA256SUMS.txt`, verify and extract into a **new** directory:

```bash
sha256sum -c SHA256SUMS.txt
# Enter the extracted yuanhe-bio-platform-<version> directory.
python3 deployment/verify_release.py
```

Verification reports the exact source commit and checks every packaged file. Compare the extracted `VERSION` with the current README. Do not reuse an old ZIP because it was previously called current. GitHub's automatic source ZIP lacks `SOURCE_MANIFEST.json`; use the named administrator asset or a clean Git checkout.

For Git deployment, clone into a new destination and check out the immutable release tag listed in [current_release.json](deployment/current_release.json). Record `git rev-parse HEAD`. Preserve the old installation, `runtime.json`, licensed models, task directories and checkpoints. Do not copy old frozen research configuration over the new defaults: the old purchase/catalog schema is incompatible with the current synthesis-only schema. Historical tasks remain in their original deployment.

## 2. Confirm four site-specific inputs

The administrator supplies these actual values:

1. A **new** writable installation prefix on Linux x86_64/WSL2.
2. The old project or licensed bundle containing both `01_molecule_generation/vendor/Peptide2Mol` and `vendor/REINVENT4`, with their exact checkpoints, prior, source and LICENSE files. Restoration verifies pinned weight hashes and refuses overwrites. If no bundle is available, obtain the pinned upstream sources/weights listed in [model_dependencies.json](deployment/model_dependencies.json) first.
3. A CUDA-enabled GROMACS executable. If missing, use the portable build below.
4. The already-approved private bind address and port. Retain existing access controls; this research service is not an authenticated public endpoint.

Install Git, Bash, Python 3.11+, curl, patch, build tools, CA certificates and [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html). Set `MICROMAMBA` if absent from PATH. Verify `nvidia-smi` in the service account's Linux/WSL environment. The shared runtime uses CUDA 13.4 user libraries; Peptide2Mol uses cu128. Supply a driver supporting the selected runtimes. Scripts never change drivers, firewall rules or system Python. Allocate space for environments, model assets and trajectories in addition to the configured 50 GiB reserve.

### If CUDA GROMACS is missing

Follow [GROMACS prerequisites](https://manual.gromacs.org/2026.4/install-guide/index.html): CMake >= 3.28, CUDA toolkit >= 12.1 and a compatible host compiler. For the recorded Blackwell setup use CUDA 12.8 and GCC/G++ 12. Do not install a Linux display driver inside WSL. Then:

```bash
YUANHE_CC=gcc-12 YUANHE_CXX=g++-12 bash deployment/build_gromacs.sh \
  /srv/yuanhe/gromacs-NEW /usr/local/cuda-12.8 4
```

This downloads official GROMACS 2026.4 with a fixed SHA-256, builds CUDA support and FFTW into the new prefix, runs upstream checks and retains logs. Use `/srv/yuanhe/gromacs-NEW/install/bin/gmx` below. Unlike the historical helper, it has no machine-specific project registration or fixed `/opt` prefix.

## 3. Plan, then install the complete stack

Run from the **new clean source directory**. Replace all example paths and the bind address:

```bash
python3 deployment/deploy_workbench.py \
  --prefix /srv/yuanhe/runtime-NEW \
  --model-source /absolute/old-project \
  --gmx /absolute/gromacs/bin/gmx \
  --host APPROVED_PRIVATE_BIND_ADDRESS --port 8777
```

The default is a read-only plan. Inspect the source SHA, paths and storage allocation, then repeat the command with `--execute`. It performs:

1. Copy and hash-check licensed model assets without changing the source installation.
2. Install the locked chemistry/preparation/MD runtime and ZVSS.
3. Install pinned Peptide2Mol with matching CUDA/PyG and patched REINVENT4 from its frozen lock.
4. Install isolated ADMET-AI/AiZynthFinder; download and verify SCScore, expansion model, templates and planning building-block pool.
5. Install official Vina 1.2.7, retaining its receipt outside tracked source.
6. Atomically bind current synthesis-only defaults to actual environments and model paths.
7. Run full deployment readiness with the service configuration.

No scientific stage, service restart or 500 ns trajectory starts here. Progress/failures are in `PREFIX/evidence/deployment.json` and named step logs. A failed step stops the sequence and does not create a service handoff script. Preserve failures; resolve the cause and use a new prefix and clean source. Experienced administrators may use [individual installation steps](deployment/WORKBENCH_DEPLOYMENT.md); configure preserves valid existing policies and rejects incompatible old schemas without overwriting them.

## 4. Switch the service and prove which code is serving

Only after readiness passes, the installer writes `PREFIX/evidence/start-workbench.sh`. Check jobs at OS level as well as in the UI: an empty run list does not prove the GPU is idle. During an approved idle window, point the existing service manager at:

```bash
bash /srv/yuanhe/runtime-NEW/evidence/start-workbench.sh
```

Run as the installation owner with existing private network controls. The launcher checks RDKit/NumPy/drawing in the configured service Python. Saving a different interpreter path in a running UI does not replace that process.

From the new source directory, supply the SHA from package verification or `git rev-parse HEAD`:

```bash
python3 deployment/verify_service.py --url http://APPROVED_PRIVATE_ADDRESS:8777 \
  --expected-commit FULL_SOURCE_COMMIT --out installation_evidence/service.json
```

The read-only check rejects a wrong version, same-version/wrong-commit deployment, dirty Git source, unverified ZIP source or missing independent diagnostics. System & Versions reports verified package identity even without `.git`. Service identity is separate from readiness and scientific execution.

## 5. Deliver the URL with real acceptance evidence

Open **Environment and diagnostics -> Check saved environment before creating a run** and save the report. Create a separate small run and record its ID. Execute preparation, Peptide2Mol, core review, REINVENT4, structure/toxicity/complete-route screening, receptor/box review, Vina and ZVSS. Keep zero-qualified outcomes explicit; never relax gates to force downstream execution.

MD requires reviewed candidate and membrane inputs. A short engineering trajectory, explicitly configured as a separate custom test, is not the formal 500 ns protocol. Verify actual CUDA mdrun, checkpoint recovery and analysis. Do not label skipped stages passed. Preserve force fields, gates, default protocols, historical records and requested research durations.

Return: source SHA, actual service Python/start command, `deployment.json`, readiness/service reports, run ID, stage logs, model hashes and artifact/checkpoint records. [Detailed troubleshooting](deployment/WORKBENCH_DEPLOYMENT.md) covers missing imports, GPU builds, model assets and interrupted runs.

## Release and deployment are separate

[Verified administrator release](.github/workflows/release.yml) packages an explicit main-branch SHA after both push CI workflows pass, verifies the extracted archive and publishes a new tag/ZIP/checksum. Existing tags are never replaced. It does not deploy or restart the internal workstation. A merge or green CI does not establish a server update.
