# Screening decision policy v2

Status: local development branch; static review after the user's no-run instruction. No new GPU or scientific calculations. Historical results remain unchanged.

## Compatibility and semantics

- Existing CLI, data directories and qualification thresholds are retained. New reports add fields and files; consumers should ignore unknown fields.
- Representative ordering is QA pass, flag, unknown, fail, followed by existing affinity/CNN/id tie breaks. This can legitimately change candidate identity and results in a future fresh run. It is not guaranteed to preserve historical selections.
- Strain and redocking qualification now require matching representative-pose provenance. Legacy records lacking identity become not-evaluated when reevaluated; existing reports are not rewritten.
- A measured blocking strain failure excludes every tier. A warning is retained as a warning. Missing strain does not become zero.
- Pareto objectives use raw affinity (minimize), region contacts (maximize) and strain (minimize). Their units remain separate. Ineligible or incomplete/nonfinite objective rows cannot enter the reported front. This is a computational tradeoff display, not a biological optimum. Contact count can reward size; it is not a validated mechanistic or potency score.
- Default `evidence_diversity` retains the prior allocation approach. Optional `pareto_diversity` applies evidence tier, then raw-objective Pareto layers, then diversity within a layer. Formal qualification still runs first. Incomplete objectives are ordered after complete layers within their tier.
- Weighted scores are within-pool percentiles, not cross-run values. A globally missing component is omitted from normalization; a missing component for one row contributes nothing. Coverage is relative to available weighted components, not completeness of all scientific requirements.

## Diagnostics and exports

`decision.json` records strategy, objective completeness, front membership and weight sensitivity. `candidate_evidence.csv` exposes pose identity, raw objectives, coverage and rank ranges. Sensitivity perturbs one weight by ±20% and removes one component at a time, holding the candidate pool and tiers fixed. It does not rerun docking, estimate model uncertainty or measure variability from seeds, receptor conformations or targets.

The offline `zvss.benchmark` module compares archived scoring/selection JSON files against docking-only, evidence ranking, formal-gate and random-list baselines. Without independent binary activity labels it cannot calculate or claim biological enrichment. With labels, comparison is retrospective and still needs leakage controls and target-independent validation. Selected list sizes and unfilled slots are explicit. Input files must belong to the same run and policy; hashes identify files but do not independently prove their scientific validity. Output must be a new directory.

## Deferred validation

Do not use new policy outputs to approve experimental purchases until the authorized acceptance work is completed. Planned checks: full module/workbench regressions; stale-pose provenance; deterministic repeat selection; default/optional strategy and missing-objective cases; report rendering; score-only/ablation comparison on independently labeled targets. These checks were not run after the user requested no execution. No model, dependency license, scientific threshold or historical result is altered here.
