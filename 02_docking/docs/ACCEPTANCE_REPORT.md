# Implementation and acceptance report

Date:2026-10-13(Asia/Shanghai). All passed have local files,/engine/and test evidence;not_tested is not counted as passing.

## Actual deployment

Windows Project:`F:\02_docking`. Dedicated WSL: `GNINA-SBVS-1e0ea3f4`. VHDX: `F:\02_docking\runtime\wsl\ext4.vhdx`. External disk format: exFAT, Imported and started WSL2; Linux Dependencies and SQLite stored in ext4 . Original Ubuntu remains separate.

Start:`START_WINDOWS.cmd`; URL:`http://localhost:5902/`. Stop:`STOP_WINDOWS.cmd`. Before migration:`PREPARE_MIGRATION.ps1`; Register on new computer:`REGISTER_ON_THIS_PC.ps1`.

## Observed versions and devices

Python 3.11.17; RDKit 2026.9.1; Meeko 0.8.0; Vina 1.2.7; GNINA 1.3.3; Molscrub 0.2.2. Complete build/sources listed in env/conda-linux-64.lock, requirements-pip.lock, package_provenance.json, vendor.lock.json.

name, driver_version, memory.total [MiB]
NVIDIA GeForce RTX 5080, 616.56, 16303 MiB

GNINA Fixed model:crossdock_default2018_1_3. GPU Explicit device 0; Actual GPU Device trace Record 768 successful runs /dev/dxg ioctl, with finite CNN outputs.

## Three fresh default-page CCO

| Run ID | seconds | Vina affinity kcal/mol | screenshots |
|---|---:|---:|---|
|20261013_110845_44ac03aa6f|1.973|-2.546|reports/browser/cco_1.png|
|20261013_110848_ab3ffdd357|1.891|-2.546|reports/browser/cco_2.png|
|20261013_110850_6221a66fb6|1.902|-2.546|reports/browser/cco_3.png|

Each run has an independent execution log and raw Vina output; identical seeded values came from actual recomputation. The first page had no prefilled results, and local assets worked with browser network access blocked. 3D remained functional.

## Batch,GNINA and scientific checks

20-molecule classroom run:`runs/20261013_104206_d78a5f9f32`, 20/20 succeeded with 180 poses. The report is stored in run/report.html, with two distinct pose images in run/images.

GNINA GPU: `runs/20261013_105343_c036895713`, 3.120 seconds;CPU CNN: `runs/20261013_105347_17a3e16bc6`, 4.309 seconds. Both have docking SDF, empirical score, CNN pose score, CNN affinity pK.

FOL Actual redocking:`runs/20261013_104230_a2c0eaef20`; Top1 RMSD 1.8458 Å, Top3 best pose and all pose best results were 1.8458 Å. using the same receptor frame and symmetry-aware/atom mapping, without independent ligand alignment. Reference: original FOL A:187, 32 heavy atoms; do not treat score_only as redocking.

1DRF associated Gln12 and other missing heavy atoms caused the initial Meeko failure; the final repair explicitly used PDBFixer to add atoms, preserving coordinate provenance and repair records, without deleting residues or rebuilding missing loop. Initial GNINA missing CUDA user-space libraries and invalid model aliases remain documented in logs/Failed run and were subsequently repaired and retested.

Custom research receptors,mmCIF plans, unknown-chain errors,Molscrub state enumeration,TopK score_only / redock, multiple seed, minimize, cancellation,/interruption,/resume, chemical identity and manual export gates have corresponding reports evidence.

## Itemized status

| Project | Status | Description / Evidence |
|---|---|---|
|A_environment_CPU|passed|Initial raw receptor attempt failed; explicit missing-atom repair and real CPU docking subsequently passed|
|B_default_page|passed|['reports/browser_acceptance.json', 'reports/browser/initial.png', 'reports/browser/cco_1.png', 'reports/browser/cco_2.png', 'reports/browser/cco_3.png']|
|C_invalid_inputs|passed|['reports/unit_tests.log', 'reports/recovery_acceptance.json', 'reports/research_acceptance.json']|
|D_class20|passed|['reports/class_report_images.json']|
|E_GNINA_GPU|passed|['reports/gpu_device_trace.json', 'reports/gnina_gpu_smoke.json']|
|E_GNINA_no_gpu_branch|passed|Explicit --no_gpu tested on GPU-equipped host; separate physical GPU-less computer not tested|
|F_chemistry_pose_RMSD|passed|['reports/artifact_acceptance.json', 'reports/FOL_redocking.json', 'reports/unit_tests.log']|
|G_recovery|passed|['reports/recovery_acceptance.json', 'reports/corruption_acceptance.json', 'reports/stop_guard.json']|
|H_spaces_Chinese_paths|passed|['reports/path_acceptance.json', 'reports/windows_path_test.log']|
|H_other_PC_full_VHDX_migration|not_tested|No second computer or cloned-VHDX registration test; no claim of full cross-machine validation|
|I_offline|passed|Real CCO in isolated Linux network namespace, browser external requests blocked|
|J_result_consistency|passed|['reports/browser_acceptance.json', 'reports/artifact_acceptance.json', 'reports/research_acceptance.json']|
|K_MD_handoff|passed|Automated test fixture explicitly labelled; not a real scientific choice by user; no MD executed|
|state_enumeration|passed|['reports/enumeration_acceptance.json']|
|target_enrichment_validation|not_tested|No labelled benchmark provided|
|teacher_platform_remote_access|blocked|Initial sandbox request failed; bounded host request returned HTTP 401 authentication required. No case content viewed; no application dependency|
|migration_local|passed|['runtime/local_machine.json', 'reports/register_acceptance.log', 'reports/windows_path_test.log']|
|runtime_export|passed|Dedicated runtime exported after project-only safe stop. Restore import not tested.|

## Untested and unsupported scope

Not accepted on another physical computer, a host without GPU or a copied new VHDX registration;Windows ARM, WSL1 is not claimed supported. Offline export verifies snapshot creation and hashes, not successful restoration. The instructor-host request returned HTTP 401, without authentication or case access; remote-proxy acceptance was not performed.

Covalent docking, general metal coordination,/complex cofactor parameterization, flexible receptors, long missing loop reconstruction and MOL2 chemical reliability are unsupported. No activity labels or enrichment evaluation; no promise of cross-hardware bitwise identity.NVIDIA The local process list for WSL did not provide matching PID, therefore GPU evidence uses actual /dev/dxg trace and finite CNN outputs.

MD handoff Acceptance files are explicitly labeled automated test samples, not user-selected candidates; no MD, force-field parameterization or free-energy calculation was performed.
