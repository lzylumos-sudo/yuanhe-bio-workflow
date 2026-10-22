# Release process

Edit root VERSION, run the metadata generator, update changelog/compatibility notes, and run scientific-baseline, registry, protocol, package, API and frontend checks. Include exact Git tag/commit, tests collected/executed/skipped, unchanged versus revised scientific protocols, dependency locks and known limitations.

Do not upgrade a protocol for a UI-only release. Do not reuse a protocol ID after a scientific change. Do not label unvalidated hardware/model/backend combinations stable. Existing generated metadata is checked in for standalone distribution and must be regenerated, not edited manually.

Heavy GPU validation is manual/target-specific and separate from pull-request CI. The 0.11.0-dev.1 governance branch is not a public Release. No tag, merge, deployment, database migration or historical-data conversion is automatic.
