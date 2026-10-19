# Verification and change inventory

Audit started at 841652e and incorporated upstream changes through 563953b by fast-forward. Branch: docs/commercialization-audit-20261019. No core algorithm, existing test, launcher, license, third-party notice or historical document is changed by this documentation patch.

## Checks

Local verification on macOS arm64 / Python 3.14 in an isolated audit environment (not the pinned Python 3.12 Linux CI runtime):

| Check | Result |
|---|---|
| English/source syntax and license digest | 701 distributable files; 0 errors |
| Repository links/version/license-copy integrity | 701 files; 0 errors |
| Frontend syntax | 18/18 JavaScript modules parsed |
| Generation | 63 passed |
| Workbench | 80 passed, 1 skipped after the single loopback test was rerun with socket permission |
| Screening, excluding integration marker | 55 passed |
| Docking, latest 563953b code | 41 passed, 2 skipped after installing missing Meeko 0.8.0 |
| Interaction analysis | 37 passed, 3 failed: missing Cairo native library (PNG output) and missing PyMOL (two scene checks) |
| Standalone MD contracts | 79 passed, 1 skipped |
| Independent packaging | 4 passed |
| Dependency consistency | pip check: no broken requirements |
| Patch whitespace | git diff --check passed |
| Updated reachable-history scan | 17 commits, 853 text blobs; no matches for the stated four credential-format families |

Unique resolved outcomes: **359 passed, 4 skipped, 3 failed**. Initial workbench result was 79 passed/1 failed/1 skipped due to sandbox loopback binding; the one failing test then passed and is counted once. Initial docking collection failed without Meeko; the complete rerun supersedes it. Do not add rerun counts to previous totals. Skips cover optional/native or absent historical fixtures and are not passes. Raw per-test JUnit/log evidence and the environment freeze are delivered with the local audit bundle, not uploaded as unreviewed run data.

The three optional-renderer failures remain visible. Application code in those modules is byte-identical to upstream; no test or scientific threshold was changed. Current tests exercise an additional local runtime, not cross-platform engine equivalence. No Peptide2Mol/REINVENT training, GPU inference, production MD, remote job recovery or wet-lab experiment was launched.

The original scan covered 15 commits/843 blobs at 841652e. The updated scan covers 17/853 through 563953b; neither is a comprehensive secret/security audit. Existing licenses, notices, launchers, scientific source and prior documentation remain byte-identical to the latest base. A final GitHub PR/CI and rendered README review are separate delivery checks; local results do not imply remote CI success.

## Patch scope

README.md plus docs/commercialization/README.md, repository_audit.md, license_audit.md, historical_mit_scope.txt, public_private_boundary.md, technical_strategy.md, closed_loop_design.md, scientific_validation.md, roadmap_18_months.md, github_metadata.json and this file. No existing technical documentation is deleted. New reports cite only already-public repository evidence and generic research plans.
