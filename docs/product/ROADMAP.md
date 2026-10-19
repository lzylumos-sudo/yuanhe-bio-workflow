# Evidence-gated product roadmap

Planning assumptions: a small development team, periodic computational chemistry review and access to experimental collaborators. These are conditional targets, not committed delivery dates; resources and budgets have not been confirmed. No calculations or training are scheduled by this document.

| Window | Deliverable | Exit evidence | Commercial action |
|---|---|---|---|
| 0–3 months | Stable research workflow, candidate evidence package, complete rights inventory | Explicit unresolved/failed states; actual MurJ acceptance and reproducibility, measured operator time | Interview target users and agree pilot scope |
| 3–6 months | Reusable target profiles and reliable deployment | Several independent targets, identical-budget score-only comparisons, time/support cost | Limited pilots with experimental partners; test willingness to pay |
| 6–12 months | Multi-objective iterative selection prototype | Held-out targets, calibration and ablations; synthetic-accessibility and diversity accounting | Repeatable deployment/support pilot; decide service vs software emphasis |
| 12–18 months | Productized module or service with defensible scope | Repeated customer use, cost and margin records, rights clearance, external biological evidence where claimed | Commercial authorization and IP review; scale only the validated use case |

For a single part-time developer, scientific validation and product support can exceed these windows. Defer multi-tenancy, broad target coverage and large-model training until the first pilot demonstrates value. Failure to improve outcomes or reduce cost should trigger narrower scope, not broader marketing.

## Closed-loop design, not implemented capability

Reviewed structures and chemical constraints → molecular proposals → validity/deduplication/alerts → multi-objective evaluation → uncertainty and diversity-aware batch choice → limited high-cost calculations or experiment → versioned observations → next iteration.

Keep predicted and measured labels separate. Objectives include potency-related evidence, pose plausibility, properties/ADMET where models are licensed and validated, synthesis feasibility, diversity and cost. No single docking score drives all selection. SA estimates are not synthesis routes; ADMET predictions are not safety evidence. Unknown properties remain unknown. Fit uncertainty models only with data; score-weight sensitivity is not model uncertainty.

Baselines: random selection, docking-only, current workflow, and the same loop without the new acquisition rule. Freeze held-out targets/scaffolds/time splits, metrics and compute budgets before tuning. Report valid unique molecules, diversity, label coverage, calibration, cost per useful observation and prospective outcomes. Include independent pose checks and counter-screens to detect score exploitation. Archive failed proposals.

## Validation ladder

1. Engineering correctness and reproducible data handoff.
2. Computational benchmarks: redocking, seed repeats, enrichment only with independent labels, resource statistics and target generalization.
3. In-vitro assays with appropriate controls and synthesis/identity verification.
4. In-vivo evidence only when actually performed under a defined study.

MD stability is one diagnostic, not binding or efficacy proof. MurJ cannot serve as both the sole tuning example and generalization evidence. Do not relabel computational candidates as validated drugs.
