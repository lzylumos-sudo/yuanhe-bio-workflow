# Architecture, naming and version governance migration

Baseline: main commit `e4906f2c76c711522654516b2ad15111707b1484`; no tags; clean cloned checkout. Work occurs on `refactor/architecture-governance`. Other local workspaces and the unavailable external drive were not edited.

## Implemented

- Recursive tracked-file hash/AST inventory, call graph and A-E duplication classifications before code edits.
- Six functional module identities, thirteen backend/history records, verified source-symbol mappings and runtime dependency reporting, used by real workbench/API/CLI.
- Single platform VERSION, generated standalone metadata, dynamic UI/API identity and independent protocol/model/data/run identities.
- Repository-native `yuanhe version/doctor/modules/protocols/validate/run/status`; unchanged old scripts and API routes.
- Public docking and generation facades that call the same original functions.
- Shared configuration hashing/layer provenance, reuse of existing validators, frozen stage configuration verified in the actual worker, explicit custom protocol identity for scientific overrides.
- New workbench and native soluble MD manifests, artifact registry, separate computation/review status, model execution evidence and actual per-stage environments. Historical records are not backfilled.
- Separate soluble and membrane protocol records, exact 500 ns / 2 fs step check, system protocol and parameter hashes; no force-field, charge or integration recipe changes.
- Registry-backed UI status, per-run original identity view, System & Versions, generated static-demo version footer.
- CI drift gates, governance regressions, preserved existing test suites and standalone archive checks.

## Protected behavior

127 protected configuration/scientific/data/history files remain byte-identical. Generation model adapters, state/preparation code, Vina scoring, ZVSS/Y-MAVS algorithms, ZLIS, standalone OpenMM/RESP/analysis code and formal default configurations are unchanged. The membrane MDP function remains text-identical; additions are provenance and consistency checks. No dependency lock, weight, license, historical data or original result was changed. No physical directories were moved or removed.

## Rollback

Switch new tasks to the original commit in a separate checkout/environment. Do not overwrite an active checkout or copy new checkpoints into an old environment. Old task directories remain authoritative for their original software/configuration. Preserve new task directories as evidence. Revert governance commits only through ordinary Git history; no force reset, database downgrade or data rewrite is required.

## Acceptance boundaries and remaining work

The implemented change is the safe integrated-workbench governance migration. It is **not full acceptance of every item in the broader platform program**. Independent docking/screening schedulers and native schemas remain preserved rather than fully consolidated; a universal target-configuration/schema layer across all standalone tools, automated legacy migration/export adapters, exhaustive error-class normalization and per-capability target-host validation evidence are not complete. Existing analysis definitions are preserved/documented rather than unified into one metric implementation.

No target-GPU model inference, real Vina execution, RESP calculation, membrane construction/conversion, checkpoint restart or production MD was performed in this refactor. Formal full-chain execution cannot be claimed from contract tests. No public service was deployed, no running research task interrupted, no database migrated, no Release created and no automatic merge requested. These limits are substantive and must remain visible in review.

See [validation](../validation/refactor_validation_report.md) and [compatibility](compatibility_matrix.md). File-level changes are generated in `changed_files.json` after implementation.
