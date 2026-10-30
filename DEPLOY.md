# Administrator deployment entrypoint

This page is sufficient for a new administrator with access to a suitable machine. No old project directory, private model transfer, separate patch or conversation history is required. It covers generation, toxicity and synthesis screening, docking, ZVSS and membrane MD. Commands run inside Linux as the installation owner. A source download alone does not install the dependencies. After startup, researchers follow the [UI workflow guide](docs/18_UI_RESEARCH_WORKFLOW.md).

## 1. Select the current source and preserve the old installation

For the latest merged fixes, clone **main** into a new directory and record its exact commit:

```bash
git clone --branch main https://github.com/lzylumos-sudo/yuanhe-bio-workflow.git
cd yuanhe-bio-workflow
git rev-parse HEAD
```

Alternatively, reproduce a published snapshot using the **current administrator build** in [README](README.md). That immutable release can predate fixes on main. Download its named ZIP and `SHA256SUMS.txt`, verify and extract into a **new** directory:

```bash
sha256sum -c SHA256SUMS.txt
# Enter the extracted yuanhe-bio-platform-<version> directory.
python3 deployment/verify_release.py
```

Verification reports the exact source commit and checks every packaged file. Compare the extracted `VERSION` with the current README. Do not reuse an old ZIP because it was previously called current. GitHub's automatic source ZIP lacks `SOURCE_MANIFEST.json`; use the named administrator asset or a clean Git checkout.

For a reproducible tagged checkout instead of latest main, select the immutable tag in [current_release.json](deployment/current_release.json). Never run `git pull` or edit tracked settings inside an active deployment. For upgrades, preserve the old installation, `runtime.json`, model bundle, tasks and checkpoints. Do not copy old frozen research configuration over new defaults: the old purchase/catalog schema is incompatible with the current synthesis-only schema. Historical tasks remain in their original deployment.

## 2. Confirm four site-specific inputs

The administrator supplies these actual values:

1. A **new** writable installation prefix on Linux x86_64/WSL2.
2. A model bundle containing both `vendor/Peptide2Mol` and `vendor/REINVENT4`. First installations obtain it using the public download helper below. An old project is optional: restoration reuses verified weights, reconstructs pinned source from cached Git objects or the public upstream, and refuses overwrites. Old vendor working-tree code is not copied.
3. A CUDA-enabled GROMACS executable. If missing, use the portable build below.
4. A private bind address and unused port. The example defaults to `127.0.0.1:8777` for access on the machine. For remote users choose the site's approved private address and access controls; this service has no public authentication gateway.

Install Git, Bash, Python 3.11+, curl, patch, build tools, CA certificates and [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html). Set `MICROMAMBA` if absent from PATH. Verify `nvidia-smi` in the service account's Linux/WSL environment. The shared runtime uses CUDA 13.4 user libraries; Peptide2Mol uses cu128. Supply a driver supporting the selected runtimes. Scripts never change drivers, firewall rules or system Python. Allocate space for environments, model assets and trajectories in addition to the configured 50 GiB reserve.

On Ubuntu 24.04, install base tools with `sudo apt-get update` followed by `sudo apt-get install -y git python3 python3-venv curl patch build-essential ca-certificates tar bzip2`. Ubuntu 22.04's default Python 3.10 is insufficient: use an administrator-provided 3.11+ interpreter without replacing system Python. Place micromamba on PATH or set `export MICROMAMBA=/absolute/path/to/micromamba`; shell activation is unnecessary. Allow network access to GitHub/release assets, Zenodo, conda-forge, PyPI, PyTorch, PyG and the repository's model URLs. RCSB is used for online structure retrieval. No paid model account or API key is required for this installation path.

The `/srv/yuanhe` examples assume a parent directory already owned by the service account. If unavailable, use a writable parent such as `$HOME/yuanhe` and replace all example paths consistently. Native Windows PowerShell is not the deployment environment; run inside Linux/WSL2.

### First installation: fetch the public models

Review the upstream [Peptide2Mol terms](https://github.com/BLUE-Flowing/Peptide2Mol), [REINVENT4 terms](https://github.com/MolecularAI/REINVENT4), [pinned prior record](https://doi.org/10.5281/zenodo.20701824) and [dependency licenses](docs/16_RESEARCH_DEPENDENCY_LICENSES.md). Run from this checkout:

```bash
python3 deployment/fetch_generation_models.py --destination /srv/yuanhe/models-NEW
# Inspect the plan, then download:
python3 deployment/fetch_generation_models.py --destination /srv/yuanhe/models-NEW --execute
```

The helper fetches exact upstream source commits and two official weights (about 162 MB of weights), retains LICENSE files and verifies size/SHA-256 against [model_dependencies.json](deployment/model_dependencies.json). Success is recorded as `status: verified` in `models-NEW/model-bundle.json`. Set `model_source` to that bundle directory. It can be reused without re-downloading. No inference or dependency installation runs here. Failures preserve the receipt and partial files: fix the cause and use a new destination; never change pinned hashes. The later synthesis step downloads its own pinned assets. A generation bundle alone is not a complete air-gapped installer.

### If CUDA GROMACS is missing

Follow [GROMACS prerequisites](https://manual.gromacs.org/2026.4/install-guide/index.html): CMake >= 3.28, CUDA toolkit >= 12.1 and a compatible host compiler. For the recorded Blackwell setup use CUDA 12.8 and GCC/G++ 12. Do not install a Linux display driver inside WSL. Then:

```bash
YUANHE_CC=gcc-12 YUANHE_CXX=g++-12 bash deployment/build_gromacs.sh \
  /srv/yuanhe/gromacs-NEW /usr/local/cuda-12.8 4
```

This downloads official GROMACS 2026.4 with a fixed SHA-256, builds CUDA support and FFTW into the new prefix, runs upstream checks and retains logs. Use `/srv/yuanhe/gromacs-NEW/install/bin/gmx` below. Unlike the historical helper, it has no machine-specific project registration or fixed `/opt` prefix.

## 3. Plan, then install the complete stack

All research defaults, model pins and runtime locks are already in GitHub. Only five machine settings are needed. Copy [site.example.json](deployment/site.example.json) outside the clean source tree:

```bash
cp deployment/site.example.json /srv/yuanhe/site.json
```

Edit the five machine settings, preserving the schema and intended software version:

```json
{
  "schema_version": "yuanhe.site/2",
  "software_version": "0.11.0-dev.5",
  "prefix": "/srv/yuanhe/runtime-0.11.0-dev.5",
  "model_source": "/srv/yuanhe/models-NEW",
  "gmx": "/srv/yuanhe/gromacs-NEW/install/bin/gmx",
  "host": "127.0.0.1",
  "port": 8777
}
```

Legacy five-field site files remain readable. Versioned files reject a different source version before installation. Use the template shipped with the selected release; never copy an old research-task snapshot into site settings. The recorded campus installation has a [separate-port upgrade template](deployment/site.campus-upgrade.json) and [upgrade procedure](deployment/UPGRADE_RUNTIME.md). Review those machine-specific paths before use.

`prefix` must be new/nonexistent with a writable parent. `model_source` is the verified bundle above; `gmx` is the CUDA executable. Keep loopback access or select a private bind address and unused port. Run:

```bash
python3 deployment/setup_workbench.py --config /srv/yuanhe/site.json
# When every preflight check passes:
python3 deployment/setup_workbench.py --config /srv/yuanhe/site.json --execute
```

Preflight collects platform/Python/tool/driver/GROMACS/model/clean-source/destination checks and prints the plan. Without `--execute` it writes nothing. Any failed check blocks installation. It does not guarantee enough GPU memory for every scientific job or probe every remote package server.

For automation without a configuration file, the equivalent lower-level command is:

```bash
python3 deployment/deploy_workbench.py \
  --prefix /srv/yuanhe/runtime-NEW \
  --model-source /srv/yuanhe/models-NEW \
  --gmx /absolute/gromacs/bin/gmx \
  --host 127.0.0.1 --port 8777
```

The default is a read-only plan. Inspect the source SHA, paths and storage allocation, then repeat the command with `--execute`. It performs:

1. Reconstruct pinned vendor sources and copy hash-verified model weights without changing the source installation; apply the selected release's patches before model installation.
2. Install the locked chemistry/preparation/MD runtime and ZVSS.
3. Install pinned Peptide2Mol with matching CUDA/PyG and patched REINVENT4 from its frozen lock.
4. Install isolated ADMET-AI/AiZynthFinder; download and verify SCScore, expansion model, templates and planning building-block pool.
5. Install official Vina 1.2.7, retaining its receipt outside tracked source.
6. Atomically bind current synthesis-only defaults to actual environments and model paths.
7. Run full deployment readiness with the service configuration.

No scientific stage, service restart or 500 ns trajectory starts here. Progress/failures are in `PREFIX/evidence/deployment.json` and named step logs. A failed step stops the sequence and does not create a service handoff script. Preserve failures; resolve the cause and use a new prefix and clean source. Experienced administrators may use [individual installation steps](deployment/WORKBENCH_DEPLOYMENT.md); configure preserves valid existing policies and rejects incompatible old schemas without overwriting them.

## 4. Switch the service and prove which code is serving

Only after readiness passes (`all_runtime_checks_pass: true` in `readiness.json`, `ready_for_service_switch` in `deployment.json`), the installer writes `PREFIX/evidence/start-workbench.sh`. On a first installation, run it in a Linux terminal as the installation owner and keep the terminal open for initial verification. For unattended use, configure the site's service manager to execute the same absolute command as that owner. No existing service-unit name is assumed. For upgrades, check jobs at OS level and in the UI, then switch the existing service command during an idle window:

```bash
bash /srv/yuanhe/runtime-NEW/evidence/start-workbench.sh
```

Run as the installation owner with existing private network controls. The launcher checks RDKit/NumPy/drawing in the configured service Python. Saving a different interpreter path in a running UI does not replace that process.

Open `http://127.0.0.1:8777/` on the machine, or the configured private address from an allowed client. Do not launch a second service on the same port. For rollback, stop new jobs first and restore the original service command with its original checkout/runtime. Never run two versions against the same task directory.

From the new source directory, supply the SHA from package verification or `git rev-parse HEAD`:

```bash
python3 deployment/verify_service.py --url http://APPROVED_PRIVATE_ADDRESS:8777 \
  --expected-commit FULL_SOURCE_COMMIT --out installation_evidence/service.json
```

The read-only check rejects a wrong version, same-version/wrong-commit deployment, dirty Git source, unverified ZIP source or missing independent diagnostics. System & Versions reports verified package identity even without `.git`. Service identity is separate from readiness and scientific execution.

## 5. Deliver the URL with real acceptance evidence

Open **Environment and diagnostics -> Check saved environment before creating a run** and save the report. Create a separate small run and record its ID. Execute preparation, Peptide2Mol, core review, REINVENT4, structure/toxicity/complete-route screening, receptor/box review, Vina and ZVSS. Keep zero-qualified outcomes explicit; never relax gates to force downstream execution.

MD requires reviewed candidate and membrane inputs. A short engineering trajectory, explicitly configured as a separate custom test, is not the formal 500 ns protocol. Verify actual CUDA mdrun, checkpoint recovery and analysis. Completion requires `Finished mdrun`, final step/time, valid final checkpoint, required outputs/analysis and the workflow completion record. A process or trajectory file alone is insufficient. Do not label skipped stages passed. Preserve force fields, gates, default protocols, historical records and requested research durations.

Return: source SHA, actual service Python/start command, `deployment.json`, readiness/service reports, run ID, stage logs, model hashes and artifact/checkpoint records. [Detailed troubleshooting](deployment/WORKBENCH_DEPLOYMENT.md) covers missing imports, GPU builds, model assets and interrupted runs.

## Troubleshooting without a separate handoff

| Symptom | Action |
| --- | --- |
| Download/hash mismatch | Inspect `model-bundle.json`, official URL and proxy connectivity. Keep failed files and retry in a new directory; never change hashes. |
| Preflight failure | Resolve every `ok: false` item, then repeat. Machine paths/tools belong in site configuration or PATH, not research defaults. |
| Install failure | Inspect `PREFIX/evidence/deployment.json` and its named step log. Keep evidence, fix the cause, retry with a new clean checkout and prefix. Reuse the verified bundle. |
| CUDA/GROMACS error | Verify the driver as the service owner and `gmx --version`; use the portable build instructions above. |
| Old code still serving | Compare full source commit, actual startup command/process and port. A GitHub merge does not restart the server. |
| Large result/editor failure | Confirm the evidence UI repair is deployed; follow [its behavior and recovery reference](deployment/UI_EVIDENCE_FIX.txt). |
| No qualifying candidate | Inspect gate evidence. Missing routes remain blocked; use another authorized small generation batch without lowering thresholds. |
| Interrupted task/revision | Preserve checkpoints/archives and follow the recovery reference. Do not delete pending markers to bypass checks. |

## Release and deployment are separate

[Verified administrator release](.github/workflows/release.yml) requires engineering, predocking-model and evidence-UI push checks for an explicit main-branch SHA, verifies the extracted archive and publishes a new tag/ZIP/checksum. Existing tags are never replaced. It does not deploy or restart the internal workstation. A merge or green CI does not establish a server update.
