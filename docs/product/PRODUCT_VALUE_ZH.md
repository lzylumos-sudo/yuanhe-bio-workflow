# Yuanhe Bio: product value, innovation and commercial incubation

2026-10-19 assessment, translated for the English-language repository. The historical filename is retained for existing links. Positioning is based on local code and historical records. Commercial conclusions are hypotheses to validate, not market research results or return-on-investment promises.

## 1. The problem: moving from many computed molecules to defensible next experiments

The central problem is turning fragmented computational results into credible, reviewable candidate decisions under limited staffing, compute and experimental budgets. Tools produce scores, but researchers still handle structure preparation, molecular identity, file conversion, poor conformers, repeated runs, membrane context, deduplication, interpretation and handoff. A high score, computational feasibility, formal screening qualification and experimental activity are distinct.

| Specific problem | Existing Yuanhe approach | Evidence needed to establish value |
|---|---|---|
| Repeated manual work connecting tools | Shared workbench, stage interfaces, generation-to-screening and MD handoff | Manual time, handoff count and recovery time for the same task |
| High docking scores hide pose/evidence problems | QA, interaction regions, review, strain, membrane context and qualification gates | Problematic-candidate rate under independent criteria; internal gates cannot circularly prove greater accuracy |
| Similar molecules dominate rankings | Chemical diversity and candidate-budget allocation | Scaffold coverage, subsequent activity and cost under equal budgets |
| Results become hard to interpret or reproduce after personnel changes | Molecule/pose identity, configuration, logs, versions and exclusion reasons | Independent reproduction success and review time |

**Severity depends on the customer.** A student docking a few molecules occasionally may have little willingness to pay. Small teams repeatedly studying targets, lacking computational specialists and handing candidates to experimental colleagues face a more concentrated need. Established pharmaceutical companies have platforms and experts with high switching costs, making them less accessible first customers. The repository contains no customer interviews, paid pilots or measured time savings; claims such as "80% savings" or "several-fold hit-rate improvement" are unsupported.

Adoption requires Yuanhe to complete a real customer task with less manual effort, clearer evidence and acceptable total cost. Customers can assemble tools themselves or use established platforms. Preference must be earned within a defined use case, not asserted as necessity.

## 2. Strengths and innovation: a candidate-decision layer

**The project implements candidate-decision rules and research-evidence management in addition to connecting tools.** `03_virtual_screening/src/zvss/select/qualify.py` defines formal qualification; `select/pipeline.py` allocates candidates and records exclusions; `scoring/aggregate.py` combines pose evidence; `report/render.py` produces reviewable results. Workbench bridges connect actual stages while retaining human review.

| Level | Defensible positioning | Limits |
|---|---|---|
| Engineering | Multi-stage interfaces, traceable handoff, failure states and review workflows | Integration is not a proprietary generative model |
| Workflow/product innovation | Explicit candidate decisions from multiple evidence types, distinguishing formal/preliminary/exploratory candidates | An original implementation does not establish world-first novelty |
| Use-case differentiation | Research involving membrane context and peptide-binding information, such as MurJ | One MurJ case does not establish generality across membrane proteins or replacement of peptide function |
| Methodological potential | Structural constraints, multi-objective evaluation and budget-aware active learning | Comparative evidence of algorithmic uniqueness or superior performance is not yet sufficient |

Virtual screening can contribute to innovation, but the sequence "docking → screening → MD" alone is not strong novelty. A new decision strategy needs a clear rationale and gains on equal-budget independent tests, ablations and subsequent experiments. Pareto optimization, weighted scores and diversity selection are existing methods and cannot become original merely by renaming them.

The strongest current value is making candidate generation, rejection and further investment traceable for researchers. Existing records also indicate that shortlisted molecules recovered only a small fraction of peptide contacts; they cannot be presented as established small-molecule replacements of peptide function.

## 3. Initial customers and reasons to buy

| Priority | Customers and users | Buyer and rationale | Requirements |
|---|---|---|---|
| First | Academic groups/small biotechnology teams with structures and experimental facilities but no computational specialist | PI or project lead; reviewable screening deliverables and training | Real target, reasonable budget and defined follow-up experiments |
| Second | Small CROs or computational service teams | Business lead; less repetitive work, consistent deliverables and reusable staff workflows | Demonstrated time savings/stability, licensing and confidentiality |
| Third | Small biotech companies with ongoing projects | R&D lead; local deployment, standardized methods and maintenance | Stronger product support and data management |
| Later | Large pharmaceutical/platform teams | Platform lead; target-specific solutions or integrable decision modules | High security, validation, procurement and integration requirements |

Students and research enthusiasts can provide feedback and adoption but should not be the main revenue assumption. Academic institutional identity does not replace assessment of the actual use and applicable license terms.

## 4. Commercial sequence: reviewable deliverables before stable software

**Proposed first product: a candidate-prioritization assessment and evidence package for one defined target.** Agree on input structures, library size, review requirements, compute budget, delivery schedule and exclusions. Deliver a reviewable candidate list, per-candidate evidence, rejection reasons, methods, raw-result index and suggested next validation. Allow zero qualified candidates; do not promise active drugs or a minimum hit count. This is a proposed business model, not an already marketed capability commitment.

Suggested sequence: focused service pilots → deployment/maintenance for repeat customers → standardized commercial licensing → validated methods/APIs. A full multi-tenant SaaS is not the recommended initial investment. Service revenue can reveal actual demand, but extensive founder effort per project creates a capacity ceiling; record labor and progressively standardize.

Pricing requires evidence. The minimum quote must account for compute, failed reruns, human review, deployment, maintenance, storage, third-party licenses and contingency. Interviews and paid pilots should test willingness to pay. Per-target/agreed-scale project fees plus maintenance are easier to explain; actual terms still require customer and licensing agreement.

## 5. Commercial value and investment assessment

**Suitable for early incubation with a defined validation plan; insufficient evidence of a defensible, high-growth software company.** Inspectable code, an actual workflow and concrete deliverables make it more assessable than a concept alone. Missing evidence includes multi-target advantages, wet-lab results, paid retention, measurable gross margins and complete commercial rights to models/data. Current records do not support valuation or investment-return forecasts.

For users, a trial is reasonable if it can reduce organization costs, reliably preserve evidence and prevent handoff errors on their tasks; procurement should depend on pilot results. For investors, useful next evidence is paying customers, controlled delivery costs, repeat business and methodological advantage, rather than repository stars or tool count. Without that evidence, staged small-scale validation and research collaboration are more appropriate than claiming an established general AIDD moat.

Proposed validation goals, not existing achievements: interview 10–15 target teams, complete three clearly scoped pilots, seek at least one paid pilot and observe continued use. If most value requires extensive customization, begin as a specialist service. Compare the same target, inputs, budget and delivery standard; record failures and negative results, not only the best case.

## 6. Three priority value propositions

1. **Candidate-decision evidence packages.** The strongest existing foundation. Clarify report fields and policy boundaries, then improve per-candidate review, human-decision records and purchase/synthesis states. Evaluate independent review time and omissions.
2. **Validated target-specific research solutions.** Start with MurJ, specify membrane context, conformation and peptide-contact constraints, and add external targets and controls. Offer bounded, reusable methods rather than claiming applicability to every protein.
3. **Budget-aware multi-objective iterative optimization.** Later integrate synthesis, properties, activity feedback and uncertainty, comparing gains under equal budgets. Model calibration, experimental feedback and licensed data accumulation may support deeper differentiation. Currently planned, not trained.

## Market references and evidence sources

- [Schrödinger Virtual Screening Web Service](https://www.schrodinger.com/platform/products/virtual-screening-web-service/) illustrates an existing automated virtual-screening product form; it does not establish Yuanhe customers or profitability.
- [Optibrium StarDrop](https://optibrium.com/product/stardrop/) provides a reference for multi-parameter optimization; "multi-objective" alone is not a unique value proposition.
- Read [PolyForm Noncommercial 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0) alongside LICENSE_SCOPE.txt. Existing licenses, models, data and third-party terms require separate review. This is not a formal legal opinion.
- Local evidence: [maturity](../09_MATURITY_AND_APPLICATION_SCOPE.md), [historical acceptance](../10_INTEGRATED_ACCEPTANCE.md), [open work](../11_RELEASE_MAINTENANCE.md). This assessment did not rerun historical scientific results.
