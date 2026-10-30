# Changelog

## 0.11.0-dev.5 - 2026-10-30

- Correct the REINVENT4 CUDA build for Blackwell using hash-verified dependency patches and a frozen cu130 lock; retain unrelated package versions and scientific defaults.
- Check dependency metadata separately from installed Python modules, and update administrator templates to the new version and a fresh runtime prefix.
- Target installation and UI acceptance remain pending; see docs/releases/0.11.0-dev.5.md.

## 0.11.0-dev.4 - 2026-10-30

- Unify runtime resolution across registry probes, diagnostics and workers; retain schema-checked compatibility with legacy environment settings.
- Correct REINVENT4 interpreter resolution, official Vina version parsing and MD executable/environment probes. Keep engineering readiness separate from scientific validation.
- Return recoverable resource-busy responses without relaxing locking or provenance; isolate delayed task, table and log responses, clear resolved UI errors and show zero valid model outputs explicitly.
- Generate frontend version guards and versioned administrator site settings from VERSION. Verify served frontend hashes against the expected source during deployment.
- Provide a separate-port campus upgrade configuration and migration instructions. Historical tasks, scientific rules, model pins and force fields remain unchanged.
- Reuse only verified model weights from old installations; reconstruct vendor code from pinned upstream commits before applying current patches, so old working-tree edits and wrappers cannot enter a new runtime.

## 0.11.0-dev.3 — 2026-10-26

- Integrate synthesis environment/asset installation and current configuration defaults. Preserve old installations and fail before service handoff on dependency errors.
- Add a portable SHA-256-pinned CUDA GROMACS build, ordered administrator installer, service verifier and package identity in the UI/API.
- Keep Vina receipts outside tracked source; publish verified packages only after successful push CI.
- Preserve scientific protocols and historical runs. Target-host model/GPU/MD acceptance remains required.

## 0.11.0-dev.2 - integrated administrator build

- Integrated the existing architecture/version-governance branch with the deployment-readiness fixes from main.
- Preserved independent no-task diagnostics, checked UI startup, evidence-path resolution and MD runtime tools while retaining governed configuration dispatch and run provenance.
- Added one generated current-download entrypoint, a separate historical-version index, immutable versioned administrator packages, source manifests and SHA-256 checksums.
- Software metadata reads root VERSION; scientific protocols, MD defaults and historical results are unchanged. Target-host full workflow and 500 ns validation remain required.

## 0.11.0-dev.1 - historical architecture-governance branch

- Unified platform identity, registry-backed UI/API/CLI, independent scientific protocol IDs and generated package metadata.
- Added frozen execution configuration, new-run manifests/artifact records and compatibility facades.
- Preserved scientific defaults, algorithms, original protocols, licenses and historical records; added strict drift and lightweight regression checks.
- No stable release, target-GPU full-chain acceptance, historical-data migration or public deployment. See docs/migration/architecture_refactor.md.

## 2026-10-21 — Vina / RESP research source migration

- Official Vina 1.2.7 CPU jobs, mapped Meeko SDF export, explicit affinity semantics and engine-separated history/caches.
- Preserve Y-MAVS and local UI improvements; remove CNN dependencies from new scoring and execution paths.
- HF/6-31G* ESP, two-stage RESP and GAFF2 charge-preserving topology preparation.
- One 500 ns production trajectory, explicit GPU selection/checkpoint recovery, 400–500 ns RMSD statistics and actual PB backend.
- Update configuration, installation/container recipes, reports, public demo source and license notes.
- Code-only delivery by request. Task-generated acceptance results removed; target-host testing and production remain deferred. See docs/17_CODE_DELIVERY_20261021.md.

## Unreleased — repository maintenance, 2026-10-19

- Clarify research use cases, independent module entrypoints, software citation and community contribution routes.

- Normalize Linux Conda explicit locks without changing package URLs, remove the legacy YAML prefix, and prevent malformed future exports; document environment reconstruction limits.
- Align UI/API version reporting and deployment guidance with the 0.4.0-rc1 record.
- Remove remaining event-specific positioning from the target-machine guide.
- Separate target-only acceptance scripts from automatic test collection.
- Restore five missing upstream license texts at their recorded hashes and include license notices in independent exports.
- Remove a developer-specific optional PLIP interpreter path and the hardcoded remote/GPU label in a fresh UI.
- Serialize docking registry initialization and fix duplicate worker spawning in its regression test.
- Support an isolated PyMOL renderer and avoid multiprocessing for a single interaction comparison.
- Declare missing docking preparation dependencies and pin the CI Ubuntu version.
- Expand CI to docking, ZLIS, standalone MD contracts, package integrity and recursive JavaScript syntax checks.
- Add the support matrix, contribution guidance and repository integrity checks.

These are repository and engineering corrections. They do not complete final GPU execution, MD repeat acceptance, raw export recovery or experimental validation.

## 0.4.0-rc1 — release candidate

See [integrated execution and limitations](docs/10_INTEGRATED_ACCEPTANCE.md) and [remaining acceptance](docs/11_RELEASE_MAINTENANCE.md). This identifies the candidate source; it does not imply a completed stable release or a published GitHub Release.

## Historical evidence

The 0.3 GPU and membrane MD records remain in `docs/06_RTX5090_ACCEPTANCE.md` and `docs/08_MEMBRANE_MD_ACCEPTANCE.md`. Preserve their original conditions and scope; do not add historical test totals to new totals.
