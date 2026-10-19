# Technical differentiation and priorities

Planning assessment, 2026-10-19. Existing capability is grounded in [repository evidence](repository_audit.md); all proposed features below remain unimplemented until separately delivered and tested.

## What the workflow actually improves

| Manual task | Existing assistance | Human work still required | Evidence needed for productivity claims |
|---|---|---|---|
| Copy structures and select chains/pocket | Preparation and staged configuration | Construct, protonation, biological assembly and binding-site interpretation | Paired setup time and error rate |
| Invoke two generators and track output counts | Pinned adapters, rejection/unique counts and core UI | Core/scaffold judgment and model suitability | Same budgets, seeds and chemistry acceptance rules |
| Convert and dock molecules | Batch GNINA calls with molecule/pose identities | Box/protocol and control selection | Failures per job, compute time and manual rescue minutes |
| Merge scores and choose MD inputs | Pose QC, qualification, diversity and identity-preserving handoff | Scientific interpretation of incomplete/conflicting evidence | Blinded review correctness and false handoff rate |
| Build membrane and recover interrupted MD | Reviewed setup, checkpoints, logs and analysis | Orientation, parameterization, convergence and biological relevance | Reproducible restart and review time, not just throughput |
| Assemble evidence | Source/input hashes and export implementation | Final acceptance and experiment selection | Complete export roundtrip and independent reproducibility |

These mechanisms plausibly remove transcription and shell coordination, but no controlled user study establishes hours saved or better scientific decisions. Proposed study: at least three independent researchers perform matched manual versus assisted tasks in counterbalanced order on held-out cases; freeze inputs and compute budget, measure active human minutes, failed/repeated jobs, unrecoverable provenance errors and review correctness. Report paired differences with uncertainty and learning effects; do not advertise a multiplier before measurement.

## Market context, not performance ranking

Schrodinger LiveDesign already integrates modeling, data and collaboration; Insilico Chemistry42 advertises generation, retrosynthesis, ADMET and simulation; XtalPi describes AI/physics/robotics workflows. Their product pages establish stated scope, not comparative accuracy. Yuanhe has no controlled benchmark against these systems, no demonstrated independent algorithm advantage and no equivalent enterprise or experimental infrastructure. The defensible initial segment is a small structure-focused research team with peptide-binding evidence and experimental access, seeking auditable candidate hypotheses.

Sources checked 2026-10-19: [LiveDesign](https://www.schrodinger.com/platform/products/livedesign/), [Chemistry42](https://insilico.com/chemistry42), [XtalPi](https://www.xtalpi.com/en/). Do not infer competitor pricing or internal implementation from these pages.

## Three focused directions

Ratings are planning judgments, not measured product scores. Durations assume two full-time developers (one scientific/ML, one engineering), part-time medicinal chemistry and experimental collaboration; they overlap the roadmap and are not additive commitments.

| Direction | Scientific novelty potential | Engineering difficulty | Commercial value hypothesis | Potential barrier | Validation cost | Prototype period | Reuse |
|---|---|---|---|---|---|---|---|
| Reproducible evidence and benchmark runner | Low as a method; high enabling value | Medium | Reduces review/rework for small teams | Curated failure cases, provenance and deployment reliability | Low-medium compute; researcher study | 4-8 weeks | High: stage records, hashes, controls, exports |
| Structure/peptide-constraint evaluation and generation interface | Medium-high only if new method wins ablations | High | Better candidates for a narrow target class | Validated constraint representation and owned data; current contact count is not a barrier | Medium-high; multi-target evaluation and assays | 3-6 months | Medium-high: contacts, geometry, cores, adapters |
| Cost-aware multi-objective active learning | Medium-high if calibrated and empirically efficient | High | Better use of scarce compute/assay budget | Proprietary longitudinal data and demonstrated sample efficiency | High; repeated independent benchmarks and wet-lab loop | 4-9 months | Medium: screening, identities, sampling backends |

The first direction is the immediate product foundation, not a claim of algorithmic novelty. Validate structure-aware objectives before using them as rewards. Start AL with a frozen library and inexpensive surrogate; a new foundation model is neither necessary nor justified now.

## Next three concrete functions

1. **Acceptance and benchmark runner:** one manifest captures code/environment/input hashes, seeds, runtime, expected controls, failures and export verification. Reproduce MurJ without changing thresholds, retain failed cases, and report held-out 1DRF controls. Gate: fresh installation plus identical candidate identities/qualification from archived deterministic inputs; GPU floating-point tolerance specified separately. No retrospective threshold tuning.
2. **Multi-objective candidate evidence table:** typed values/units, hard validity and pose gates, diversity, route evidence, ADMET applicability and unknown states; return a Pareto set with human review. Gate: invalid/missing evidence never becomes a pass; fixed-pool ablation versus current ranking across targets. No claim that SA score proves synthesis.
3. **Budgeted active-learning sandbox:** fixed public/licensed pool, train/validation/test separation, simple calibrated surrogate ensemble, diverse query batches and a durable iteration ledger. Gate: replicated improvement over random and greedy baselines at matched query and compute budgets; no test-label access during acquisition.

Do not prioritize a large model, broad commercial API or additional UI polish before these gates. Research/Enterprise separation is described in [boundary proposal](public_private_boundary.md); [closed-loop design](closed_loop_design.md) defines the evaluation safeguards.
