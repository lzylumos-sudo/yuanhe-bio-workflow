# Local upgrade status — 2026-10-19

Baseline: `c1a87a54d7c8447aac97820a7629b472ef87f28e`. Branch: `feature/screening-evidence-v2-20261019`.

## Delivered scope

- Root branding and all six module README navigation, with historical instructions preserved.
- Code/content map and tracked-file inventory, preserving existing paths and technical document history.
- Screening-policy corrections: QA-first representative, pose-bound strain/redocking evidence, missing-objective Pareto handling, strain warning/failure distinction and run-specific fingerprint cache clearing.
- Default strategy retained; optional Pareto/diversity configuration and versioned decisions added. Selection identities include evidence but exclude previous selection annotations.
- Source additions for weight sensitivity, evidence CSV/decision JSON and offline comparison. These are not evidence that a new scientific run or comparison was performed.
- Product value, customer segmentation, commercialization hypothesis, roadmap and public/private boundary. GitHub metadata remains a local proposal.

## Verification boundary

The user requested no execution during the upgrade. A local non-integration test command had already been started before that instruction; when cancellation was attempted it had completed with exit code 0. There were subsequent edits, so that result does not validate the final revision. No tests were started after the no-run instruction.

Final review is static: Python syntax parsing, YAML structure inspection, Markdown target existence, whitespace/diff inspection and file-copy integrity. No application launch, model training, generation, docking, screening run, MD, benchmark execution or report-rendering acceptance is claimed for this revision. Full regression, cross-platform and scientific acceptance remain open.

The existing test file for decision-policy regressions is retained as source; its existence is not evidence of final acceptance. Existing logs and candidate outputs are preserved rather than regenerated. The prior README's 9NU4 counts refer to historical evidence, not this revision.

## Known limitations requiring future work

- Weight sensitivity is neither model calibration nor predictive uncertainty.
- Raw contact counts can reward larger molecules; pose-quality and percentile conventions still need scientific review. They were not retuned for favorable outputs.
- A global pose-QA fallback still exists in legacy qualification paths; run provenance and stale-state cases need broader review before reuse across runs.
- Formal qualification and comparison agreement are policy consistency, not biological accuracy. Independent labels and held-out targets are required.
- Legacy evidence without pose IDs can become not-evaluated under the new policy. New evaluation must use a new run and retain old reports.
- Enterprise access control, customer isolation, calibrated ADMET, validated active learning and multi-target performance are not delivered.
- This is not a complete repository-history secret scan, IP ownership clearance or legal opinion.

No licenses, third-party notices, original research data or model weights were edited. No public push or pull request is part of this local revision.
