# Contributing

Use an isolated branch and one coherent change per commit. Start from the current remote branch; do not merge an older history after the 2026-10-19 commit-title rewrite. Preserve any local work before updating a checkout.

1. Record the observed problem, minimal reproduction, affected revision and environment in an issue. Remove private research inputs, tokens and local credentials from public reports.
2. Make the smallest complete correction. Preserve original input/output provenance, failed attempts and scientific thresholds. New scientific assumptions require explicit documentation and review.
3. Run the affected module's tests and repository checks in [the support matrix](docs/12_SUPPORT_MATRIX.md). Add a regression when correcting behavior. Unit tests, actual engine runs and experimental validation are separate evidence.
4. Update usage and limitations when behavior changes, and add a concise entry to `CHANGELOG.md`.
5. Describe the problem, resulting behavior, tests, skipped checks and remaining limits in the pull request. A green CI run does not close the final scientific acceptance checklist.

Use short subjects such as `fix: preserve ligand identity` or `docs: clarify Linux setup`. Details and commands belong in the commit body or pull request. Do not rewrite published history during routine maintenance.

The existing project and upstream licenses apply. Submitting code does not change the commercial-use boundary or authorize redistribution of third-party models, weights or data.
