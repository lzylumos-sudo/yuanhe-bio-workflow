# Incremental AI closed-loop optimization design

Design only. No model training, new algorithm implementation, proprietary objective formula, private dataset or unpublished compound is included. All numerical budgets below are proposed pilot settings to freeze before execution, not achieved results or guaranteed adequate sample sizes.

## State and evidence flow

Reviewed target -> generator proposal -> chemical/pose validity gates -> multi-objective evidence -> calibrated uncertainty and applicability -> diverse acquisition batch -> permitted computational oracle or assay -> immutable observations -> next iteration. Human review approves objectives, protocol changes and experimental orders. Keep computational labels and measured assay labels in separate tasks; a docking-trained surrogate predicts a docking result, not biological activity.

Reuse current adapters, source snapshots, candidate/pose identities and stage manifests. Add an iteration record with parent model, dataset version/hash, training split, objective versions, acquisition decision, budget spent, failures and accepted observations. Never overwrite an old prediction with a later measurement. Unavailable properties retain null plus reason, not zero; use versioned schemas and units.

## Objectives and guardrails

| Dimension | Initial evidence | Limitation / decision rule |
|---|---|---|
| Validity and identity | RDKit sanitization, valence, charge/stereo review, unique canonical identity | Hard gate before scoring; record every attempted sample denominator |
| Novelty and diversity | Exact identity plus scaffold/fingerprint distance to a specified reference set | Novel to a dataset is not patent novelty; avoid rewarding trivial out-of-domain structures |
| Structure compatibility | Pose clash/strain checks, coordinate-consistent contact geometry, redocking controls | Contacts are hypotheses; require ablations and independent evaluation before rewarding them |
| Synthesis feasibility | Existing forward-reconstruction/route contracts, stock snapshot and chemist review | SA surrogate is only triage; route validity is not a laboratory yield guarantee |
| ADMET | Licensed endpoint-specific predictors with provenance and applicability flags | Train/calibrate per endpoint; do not silently impute unavailable/toxicology labels |
| Binding prediction | Docking/consensus measurements with consistent units and protocol | Not the sole objective; score gains need held-out oracle and assay confirmation |
| Cost and uncertainty | Wall time, GPU/CPU use, ensemble dispersion, calibration and domain distance | Ensemble variance is not automatically calibrated epistemic uncertainty |

Prefer explicit hard constraints plus a Pareto front over a hidden scalar. Normalize any display scores using development data only; predeclare tradeoffs and report raw endpoints. Keep missingness and uncertainty visible so an unmeasured candidate cannot dominate measured candidates by omission. Expensive MD applies to a small diverse subset after cheaper QC; trajectory stability is not a binding-affinity label.

## Phased implementation

**P0: offline evidence, weeks 1-8.** Freeze molecule standardization, scaffold/target splits and benchmark manifests. Add multi-objective reporting without changing generation. Use 3-5 targets with available controls as a pilot and at least five prespecified random seeds where feasible. MurJ is a case study with missing activity labels, not an enrichment benchmark.

**P1: fixed-pool active learning, months 3-6.** Proposed bounded experiment: 1,000 licensed/public pool molecules, 100 initial diverse labels, batches of 25, at most 10 iterations, five seeds, fixed oracle settings and hard compute/time cap. Re-estimate sample size and cost using a small development-only pilot, then freeze the protocol. Use a simple fingerprint-based surrogate or small ensemble before neural training. Compare random, diversity-only, greedy mean prediction and uncertainty-aware selection at equal oracle calls and total compute. Reserve an exploration fraction chosen on development data, with scaffold diversity in every batch.

Train on observed development labels; use a distinct validation set for model/hyperparameter selection and uncertainty calibration. Evaluate regression uncertainty with coverage at nominal intervals, interval width, residual plots and OOD slices; for classification use calibration curves/Brier scores plus discrimination. Report calibration breakdown when chemistry shifts; abstain or request measurements rather than producing confident extrapolations. Assay noise and model uncertainty must be distinguished.

**P2: structure-constrained proposals, months 6-12.** Evaluate explicit geometric constraints and scaffold restrictions with existing generators first. Compare unconstrained generator, constraints-only, AL-only and combined variants. Hold sampling and scoring budgets equal; retain rejected molecules and negative outcomes. Introduce reinforcement learning only after reward leakage and reward-hacking tests pass; constrain divergence/diversity and preserve an independent evaluation oracle. No large-model training is authorized by this plan.

**P3: measured feedback, months 9-18 conditional on partner.** Import assay records with concentration, units, protocol, replicate, controls, purity and censored values. Track assay drift and use target-specific tasks. Human approval selects synthesis/assays; measurements update models only after QC. A computational loop may precede this but must not be described as experimentally autonomous discovery.

## Prevent self-rewarding false improvements

- Separate acquisition scoring from final evaluation. An improved surrogate or docking reward is not success unless a frozen independent endpoint improves.
- Remove duplicate/near-duplicate leakage by scaffold clustering and target-family separation; add a temporal split when dates exist. Audit overlap with upstream training data where known; flag unknown pretraining overlap.
- Freeze test sets and reference hypervolume points before optimization. Never choose final objectives, thresholds or cherry-picked seeds after viewing test outcomes.
- Audit molecular size/lipophilicity/charge trends, strained poses, reactive groups and repeated scaffolds; run adversarial score-extreme controls. Reject invalid chemistry even when the reward is high.
- Version oracle updates. A changed engine or protocol creates a new label source; do not mix incompatible targets/units as one training task.
- Stop on budget exhaustion, no prespecified improvement, diversity collapse, invalid-rate increase or calibration failure. Thresholds are chosen in the pilot and preregistered; this document does not invent universal scientific cutoffs.

## Benchmark and decision rule

Generation: attempted/valid/unique counts, novelty relative to declared corpus, scaffold diversity and distribution metrics. Use [MOSES](https://github.com/molecularsets/moses) and [GuacaMol](https://github.com/BenevolentAI/guacamol) as reference protocols after dataset-license review; their scores do not prove target activity.

Optimization: feasible-set rate, Pareto hypervolume with fixed units/reference, independent top-k outcomes, query count to a prespecified target and area under the learning curve. Report per-target/per-seed outcomes and bootstrap uncertainty with target-level grouping; do not treat correlated poses as independent samples. Include failures and full compute cost. Report computational-oracle efficiency separately from assay hit rate. Promote a method only if a held-out improvement is reproducible and does not trade away validity, diversity or synthesis feasibility beyond prespecified limits. If intervals include no benefit, report inconclusive rather than superior.
