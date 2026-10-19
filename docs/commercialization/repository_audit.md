# Repository audit

Audit date: 2026-10-19. Baseline: `841652ecd65e05aee272283126cb4c5de274a455` on `main`, Release 0.4.0-rc1. This report audits the public source and its recorded evidence; it is not a new GPU experiment or final acceptance certificate.

## Scope and method

The live GitHub metadata reports Public, default branch main. A full clone contains 15 reachable commits, main and codex/rtx5090-workflow-20261018 remote branches, and no fetched tags. The working branch is docs/commercialization-audit-20261019. The initial tree has 690 distributable files. Read source, environment manifests, CI, module tests, release records, license files, all reachable historical blobs and historical license changes before editing. Baseline English/source and repository-integrity checks both passed. Final checks are in [verification](verification.md).

The connector reported read-only access, but the existing local Git credential was separately verified to have repository write/admin access. Repository metadata edits are separate from a branch/PR and are verified separately. During the audit, main advanced through d39b7fe and 563953b (dependency/CI and docking-registry fixes). The documentation branch was fast-forwarded to 563953b without modifying those changes. The initial historical scan counts below refer to 841652e; final verification records the updated scope.

## Maturity scale

A = implemented with recorded real execution, bounded to the stated case; B = implemented but insufficiently validated; C = partial; D = planned; E = historical/problematic. A does not imply biological validation, broad platform support or production readiness. Evidence below is repository-recorded unless explicitly identified as rerun in verification.md.

| Module | Class | Code evidence | Evidence and limits |
|---|---|---|---|
| Local workbench, review, cancellation and provenance | A / B at product level | `00_workbench/webui/workflow.py`: create, review, start, cancel, verify_source; `webui/server.py`; `studio_worker.py` | Real browser execution and disconnect/reconnect recorded in docs/10_INTEGRATED_ACCEPTANCE.md; single-workstation orchestration is not multi-tenant scheduling |
| Structure/pocket preparation | A, exact MurJ input | `01_molecule_generation/src/structure.py`; `00_workbench/scientific/target_profile.py` | Recorded 3,785 receptor heavy atoms, 258 pocket atoms; chain, construct, protonation and box remain scientific review decisions |
| Peptide2Mol integration | A, bounded inference | `01_molecule_generation/src/adapters/peptide2mol.py`: REPO_COMMIT, platform_readiness | 16 attempts, 1 accepted in recorded run; released adapter requires Linux. No new model or demonstrated peptide-function retention |
| Core review and REINVENT4 expansion | A, bounded sampling | `01_molecule_generation/src/cores.py`, `src/adapters/reinvent4.py` | Recorded 16 samples, 8 unique molecules; core selection is manual; no proprietary learned optimizer |
| GNINA docking | A, bounded execution | `02_docking/src`; `00_workbench/scientific/validation.py`: docking_control | Nine molecules/81 poses for MurJ; 1DRF three-seed control has top-ranked recovery 2/3 at 2 A, best-pose 3/3. No activity enrichment evidence |
| Screening and MD handoff | A, bounded execution | `03_virtual_screening/src`; tests/test_handoff_identity.py, test_score_semantics.py, test_redock_coordinates.py | Corrected archived reassessment: 2 formal, 5 preliminary; original failed assessment retained. Formal means protocol qualification, not validated binding |
| Peptide-contact evidence | A, geometric comparison | `00_workbench/scientific/peptide_contacts.py`: compare | Seven shortlisted molecules each recover 1/22 residues; same-pool contact sorting does not improve ordering. Negative result retained |
| Integrated membrane/GROMACS | C for final release acceptance | `00_workbench/scientific/membrane.py`: build_prepared, run, checkpoint_probe, analyze; `scientific/validation.py`: summarize_replicas | Earlier bounded run recorded; two 5 ns repeats complete at last observation, third in equilibration, aggregate verdict and final UI export unobserved. Do not infer completion from elapsed time |
| ZLIS interaction analysis | B | `04_interaction_analysis/src`, tests, docs/pytest_record_macos.txt | Recorded macOS checks; optional environment, not a replacement for integrated scientific acceptance |
| Standalone OpenMM | B | `05_molecular_dynamics`, V2_VALIDATION.json, tests | Distinct optional implementation with historical records; not the integrated GROMACS engine |
| Synthesis feasibility | C | `01_molecule_generation/src/synthesis.py`: audit_route, forward_reconstruct, stock_is_stale | Route/identity contracts and stock snapshot checks exist; no validated automated synthesis planning or laboratory synthesis demonstrated |
| Packaging, installer, doctor | B | `deployment/package_component.py`, package_md.py, doctor.py, restore_models.py | Clean-environment record exists; model restore needs local vendor sources/weights; not one-command reproducibility from public Git alone |
| Independent commercial API, tenants, RBAC, enterprise queue | D | Current local server and compute lease are only reusable starting points | No enterprise acceptance evidence; local HTTP endpoints do not establish a secure commercial API |
| Proprietary generation, calibrated uncertainty, active learning loop | D | Existing adapters and filters are building blocks | No demonstrated original learned model, validated AL acquisition or autonomous design loop |
| Docking prototype | E | `archive/docking_prototype` | Historical reference; preserve, do not promote as current path |

## Engineering evidence versus scientific evidence

CI in `.github/workflows/engineering.yml` defines separate generation, workbench, screening, docking, interaction, standalone MD and packaging suites plus frontend syntax and source-integrity checks. Excluded integration tests, missing dependencies and synthetic UI fixtures must remain labelled. Historical counts in docs/04_VALIDATION.md overlap later suites and must not be added. docs/10_INTEGRATED_ACCEPTANCE.md records workbench 76 passed/4 skipped and screening 55 passed at its checkpoint; these are not reruns at this audit commit. Raw GPU runs, weights and complete trajectories are not in the public clone, so independent reconstruction of all historical claims is unavailable here.

## Reproducibility and portability findings

- `webui/workflow.py:defaults` contains /opt/yuanhe prep/MD interpreter defaults. `deployment/runtime.example.json` contains installation-specific absolute executables. These are configurable, not universal paths. Validate configuration with doctor on each host.
- `deployment/membrane-environment.yml` pins several engine versions but leaves RDKit/numpy unpinned and requests CUDA 13.1; the recorded GPU environment is not a portable solved environment for every OS.
- `deployment/model_dependencies.json` pins model source commits and weight hashes. `restore_models.py` requires a locally supplied vendor tree; a checksum proves identity, not availability or commercial permission.
- Linux x86-64 GPU is the demonstrated scientific path. macOS is suitable for preview and supported checks; use a Linux worker for the released model adapter. Windows uses WSL2; RTX 5070 target acceptance is still open. See docs/12_SUPPORT_MATRIX.md.
- The UI binds loopback by default. Do not treat the current local service as an Internet-facing enterprise deployment. No server exposure change is proposed.

## Public history and sensitive information

The first reachable commit is adfdfcb (2026-10-17), not proof of the first-ever distribution. All 15 reachable commits and 843 unique text blobs were scanned for private-key headers, common GitHub/OpenAI/AWS credential formats and sensitive key/weight/runtime-image filenames. No matches in those narrow checks. This is not a security clearance: no entropy scanner, GitHub secret alerts, removed remote refs, forks, release attachments, external hosting or unpublished files were covered. No credentials are reproduced in this report.

History already discloses adapters, screening rules, orchestration, MD logic, example structures, acceptance counts and selected candidate identifiers. A filename/content scan cannot establish ownership of scientific data. No customer dataset was identified as such in the reviewed public source; owner confirmation is still required. Deleting current files would not recover secrecy. See [license audit](license_audit.md) for historical MIT scope. No sensitive-file removal, history rewrite, license replacement or code migration was performed.

## Demonstrations and deployments

The audited tree directly provides a local Python browser workbench and Linux/WSL launchers, plus optional module interfaces. No production public website deployment manifest or verified public demo URL was found in the current README/docs/deployment search. Earlier external ModelScope/static-site work must be checked against its own deployed tree; it is not evidence that this repository's full GPU workflow is hosted. Public demo uptime, content identity, mainland accessibility and backend isolation were not tested in this audit. Keep any future demo read-only with approved public inputs until separately accepted.

## Conclusion

An integrated research candidate with meaningful bounded computational records and reusable engineering safeguards. It is not a validated general-purpose AIDD product, a proven peptide mimic generator, a demonstrated enterprise service or a commercially cleared dependency bundle. Priorities are reproducible acceptance, defensible benchmarks and a rights inventory before algorithm expansion.
