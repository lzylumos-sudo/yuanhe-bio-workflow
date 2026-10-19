# Repository review and upstream practices

Review date: 2026-10-19. Base revision: `18c7a6e5b83f19d56d9381ee06a842e2921d1b4a`. Scope: distributed source, repository documentation, portable tests, entrypoints and independent source archives. No remote scientific job was launched or recovered during this review.

## Primary-source references

| Official project | Observed practice | Applied here |
|---|---|---|
| [GNINA](https://github.com/gnina/gnina) | Installation guidance, concrete docking commands, examples/tests and separate upstream license texts. | Distinguish UI preview from engine installation; restore missing recorded license texts; check portable docking tests and exported files. |
| [REINVENT4](https://github.com/MolecularAI/REINVENT4) | Explicit platform qualification, environment/dependency instructions, changelog and contribution guide. | Add a support matrix, clarify Linux versus unverified WSL execution, and add contribution/change records. |
| [OpenMM](https://github.com/openmm/openmm) | Separate examples/tests, user/developer help and component-specific license information. | Separate actual engine acceptance tools from test collection; extend standalone MD checks and retain license boundaries in exports. |

These are engineering references, not claims of scientific equivalence, upstream endorsement or competitive performance. Upstream installation instructions do not automatically validate our integrated environment. No scientific threshold was weakened to make a check pass.

## Confirmed defects and corrections

| Finding | Correction / verification |
|---|---|
| Target installation guide retained event positioning. | Replace it with research validation and release terminology. |
| API and console reported 0.3.0-rc1 while README identified 0.4.0-rc1. | Align the API and startup message; add a repository consistency check. |
| Setup described an unverified WSL target without foregrounding the recorded Linux environment. | Add the support matrix and evidence links. |
| Screening module overview still described every membrane orientation as a proxy. | Document the exact-construct checked PPM3 path and preserve unsupported-input limits. |
| Five license texts named in the existing inventory were absent. | Retrieve each from the already recorded official tag and verify its existing SHA-256; CI now checks presence and hashes. |
| Independent exports omitted detailed notices and referenced license texts. | Share license inclusion across all four exporters; verify archive hashes, notices, CLI imports and overwrite refusal after extraction. |
| Optional PLIP interpreter configuration contained a developer's absolute path. | Remove the tracked path; retain an ignored local override and document the environment variable. |
| Test collection executed dedicated runtime or completed-run validation scripts. | Move those tools to explicit `scripts/validate_*.py` entrypoints; make integrity regressions ordinary isolated tests. |
| Frontend syntax check omitted nested flow modules. | Recursively check all 18 ES modules and run the check in CI. |
| Clean Linux docking collection failed because Meeko was absent from package metadata and the combined test environment. | Declare Meeko and its required chemistry dependencies in the docking package; install the recorded Meeko version in CI. |
| Clean Linux exposed concurrent SQLite journal initialization and a spawn test launching the same worker twice. | Serialize registry initialization, preserve existing WAL mode, close the parent log handle, reap the actual spawned child, and add a six-process initialization regression. |
| Linux ZLIS scene tests lacked PyMOL; the available renderer wheel requires NumPy 2 while chemistry is locked to 1.26. | Support an explicit separate PyMOL interpreter and install its isolated CI environment; retain actual scene creation/reopening assertions. Run one ligand comparison serially to avoid unnecessary fork warnings. |
| CI used deprecated action runtimes and a moving Ubuntu label. | Update official actions and pin Ubuntu 24.04 to keep the acceptance environment explicit. |
| CI covered only three main suites. | Add docking, ZLIS, standalone MD and standalone archive checks, with explicit dependencies and preserved JUnit results. |
| Historical execution hashes differed from the renamed Git history. | Preserve original evidence hashes and supply a file-tree-equivalent old/new commit mapping. |
| A fresh local UI displayed a hardcoded remote/GPU environment label. | Replace the implied connection status with a configuration prompt; verify preview and environment dialog in the browser. |
| A pipeline type annotation referred to an unimported type. | Import `Receptor`; no scoring or molecular behavior changes. |
| The relocation regression promised Unicode coverage but used only ASCII paths. | Exercise both relocated source and output directories with spaces and Chinese characters; retain real analysis and scene export assertions. |

A static warning about `sim` in nested OpenMM functions was inspected: `sim` is assigned before those functions are called and deleted only after completion. It was not reproduced as a runtime defect, so the MD implementation was left unchanged. Unused-import/style observations are not presented as scientific failures.

## Validation boundary

Clean Ubuntu 24.04 validation at `757c0f2` completed with **362 passed, 4 skipped** and no failures; [the preserved CI run](https://github.com/lzylumos-sudo/yuanhe-bio-workflow/actions/runs/39208990129) contains all seven JUnit reports. The skips cover optional native docking dependencies, MolScrub, an opt-in native Ubuntu MD check, and undistributed historical outputs. The subsequent Unicode-path regression also passed locally; consult the current revision's CI for its Linux result.

See the current commit's GitHub Actions run for fresh Linux JUnit reports. Local macOS checks: **361 passed, 5 skipped**, plus 18 JavaScript modules parsed and five upstream license hashes verified. Exact outcomes and skip reasons are in [the audit result](history/repository_review_2026-10-19.json). A fresh browser preview and its environment dialog opened without observed console errors; no scientific job was started. Optional native checks and absent historical output must be listed as skipped, not passed.

This review does not close the existing A1-A5 acceptance checklist. Final clean deployment, a second complete supported target, final remote repeat results, raw-data export recovery and real-UI acceptance remain subject to their existing evidence requirements. No experimental activity claim is established. Credential-pattern scanning and local HTTP regression checks are limited checks, not a full security or dependency-vulnerability audit.

`docs/FILE_MAP.csv` is the original migration inventory, not a manifest of current files. It retains original hashes even when a later change renames or edits a file. The removed local PLIP setting and relocated acceptance scripts are recorded by this maintenance commit, not by rewriting that historical inventory.
