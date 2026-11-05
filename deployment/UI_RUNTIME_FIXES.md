# UI runtime corrections after the November 2 deployment report

Scope: source commit `baf29bcfb72d4c2d2bb2173b735b5dadfa4719f3` and its pinned deployment contracts. These are engineering fixes, not evidence of molecular activity, prediction quality or completed production MD.

## Corrections

- Load official SA_Score from either PyPI or RDKit's configured Conda Contrib directory; require its fragment data. Readiness calculates a reference score. The import failure does not establish that RDKit 2026 removed Contrib; see the [RDKit maintainer's explanation](https://greglandrum.github.io/rdkit-blog/posts/2023-12-01-using_sascore_and_npscore.html).
- Require the hash-recorded Peptide2Mol patch before preparation/inference. The base revision already shipped and installed a bounded-sampling patch; pristine upstream code is not the deployed contract. The revised patch counts every rejected attempt, handles the final partial batch, emits progress and leaves a failed audit on model/device errors. Never fabricate attempt_budget.json or suppress its validation.
- Redirect losses.txt to the run output. Run diagnostic subprocesses, including interpreter exit hooks, inside disposable directories so Psi4 cannot leave timer.dat in the checkout.
- Install a separate pinned PyMOL renderer, bind pymol_python, test actual scene/PNG creation in readiness, and render both artifacts through the UI module adapter. Partial tables/SVG remain available with explicit failure if rendering fails.
- Refuse case-colliding vendor checkouts on insensitive filesystems; verify the restored bytes against every tracked Git blob. A clean Git status alone cannot prove a lossless NTFS checkout.
- Support explicit HTTPS mirrors, hash-preserving Conda lock rewriting, effective pip index selection, and verified download retry/resume. No third-party proxy is chosen automatically.

The original scientific baseline inventory remains intact. Two reasoned implementation revisions are recorded in `docs/architecture/scientific_baseline_revisions.json`; all other protected files retain their recorded hashes. Weights, diffusion steps, qualification thresholds and MD protocols are unchanged. Budget means attempts including rejections; valid molecule counts may be smaller or zero.

## Deployment and recovery

Use the existing administrator upgrade procedure with a new checkout/prefix, preserving running jobs and their source snapshots. Fresh full installations now include the renderer step. For a deliberately maintained existing installation:

```bash
python3 deployment/apply_vendor_patches.py --vendor Peptide2Mol
bash deployment/install_workbench.sh renderer /absolute/installation-prefix
```

Known previous patch hashes upgrade idempotently; unknown modifications are rejected. Set `pymol_python` in runtime.json to `/absolute/installation-prefix/envs/pymol/bin/python`; new installations configure it automatically. The host may also need OpenGL loader libraries; readiness reports actual rendering failures. These commands do not change drivers or existing frozen job settings.

Run doctor and then a small UI acceptance task through the actual screening gates. Check SA_Score, interaction_renderer, model patch identities, OpenMM CUDA and GROMACS separately. Zero qualified molecules is legitimate, not permission to bypass screening.

Prefer WSL ext4 for environments and vendor source. On NTFS, enable case sensitivity on the **new destination before checkout**, using Windows administrator approval. The installer does not change filesystem flags.

Set `YUANHE_CONDA_MIRROR` to an approved HTTPS conda-forge mirror base and `PIP_INDEX_URL` to an approved Python index. The generated `.download.lock` preserves original package filenames/checksums and leaves the repository lock unchanged. Standard HTTPS proxy settings remain available.

Resume a bundle previously created by the acquisition tool with:

```bash
python3 deployment/fetch_generation_models.py --destination /absolute/model-bundle --resume --execute
```

Optional `--peptide2mol-url` and `--libinvent-url` select HTTPS mirrors for the same pinned bytes. Resume revalidates source/weight identities; `.part` files use HTTP Range or restart safely if Range is ignored. Mismatched final hashes are never accepted. Preserve a suspect partial file for diagnosis and remove it explicitly before restarting that file. Arbitrary manually assembled directories are not recorded resumable bundles.

The legacy membrane installer accepts an existing `MICROMAMBA` executable and otherwise gives installation guidance, avoiding reliance on one bootstrap API endpoint.

## Acceptance boundaries

Tests execute the actual patched sampling method with synthetic outputs, including empty output, origin/separation/decode/reconstruction failure, disconnected molecules, valid output, partial batches and device failure. They do not run model inference. Other regressions cover real installed SA scoring, Conda layout, diagnostic exit hooks/timeouts, renderer selection, download ranges and source identity.

`CUDA_ERROR_UNSUPPORTED_PTX_VERSION` requires resolving the host driver/PTX mismatch. On WSL, update the Windows NVIDIA driver, restart WSL, and repeat the actual OpenMM CUDA test; do not install a Linux display driver. Comparing the nvidia-smi banner to a library version is insufficient because [CUDA minor-version compatibility has PTX limitations](https://docs.nvidia.com/deploy/cuda-compatibility/minor-version-compatibility.html).

GROMACS compilation/regression tests, a small real GPU mdrun, complete weight/stock downloads and real model/UI acceptance must be recorded on the target host. Repository tests or a successful build do not substitute for those results.
