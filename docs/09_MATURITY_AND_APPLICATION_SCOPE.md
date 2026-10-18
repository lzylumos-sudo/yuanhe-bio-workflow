# Application scope and maturity assessment

Assessment date: 2026-10-18. This is a research workflow under validation, not a commercially validated drug-discovery product. A successful calculation, software test, or short MD trajectory does not demonstrate functional replacement of a peptide, biological activity, or experimental hit-rate improvement.

## Intended initial application

The proposed initial users are teams with a target structure, peptide-binding evidence and experimental capability who need a reproducible workflow for generating and prioritizing small-molecule hypotheses. Candidate evidence, exclusions, provenance and reproducible inputs are deliverables. Functional replacement, oral availability and clinical equivalence are not established outputs.

The MurJ example uses a phage lysis-protein complex, not a marketed peptide drug. It is an antibacterial discovery hypothesis. No experimentally validated MurJ small molecule has been produced by this project.

## Evidence and unresolved gaps

- The previous RTX 5090 UI acceptance demonstrates one engineering path through membrane preparation and short MD. See `08_MEMBRANE_MD_ACCEPTANCE.md` for the exact conditions.
- A fresh bounded generation run, `20261018_091310_dea50e56`, completed all 16 attempts with one accepted molecular structure and 15 rejected attempts. This is not a hit-rate estimate or a new full-pipeline acceptance.
- Reconciliation found deployment-only REINVENT token/configuration fixes and screening identity handoff fixes absent from the repository. They are now included in source. Generation adapter checks passed (14 tests) on the compute environment.
- A missing `RegionConfig.usable` property prevented evaluation of configured contact regions. The correction and receptor/title handoff checks passed 19 targeted tests. The broader non-integration screening suite passed 39 tests with three skipped tests. Real-data integration acceptance remains separate.
- Historical zero formal candidates include missing evidence and an implementation defect. They must not be interpreted as proof that all candidates failed a complete scientific evaluation. Historical outputs remain preserved; corrected results require a new recorded evaluation.
- The current de novo path uses a peptide-derived receptor pocket; it does not yet establish full transmission of peptide pharmacophore/functional constraints into candidate selection. Selection of generated scaffold atoms is not evidence of peptide-function retention.
- Membrane-orientation qualification, additional targets, controlled comparisons, independent MD replicas and clean-environment reproducibility remain unfinished maturity workstreams. These corrections do not close them.

## Commercial and scientific differentiation

Peptide-inspired generation already exists in Peptide2Mol. Integrated discovery environments already exist in LiveDesign and Orion. Integration can provide engineering value, but tool count and a graphical interface alone do not establish novelty or competitive superiority.

The testable differentiation is: preserve traceable peptide constraints, reduce manual intervention and failures, and improve decision quality within a fixed computational budget. Compare the complete workflow with conventional screening, the underlying model alone, and a version without the added peptide constraints. Prespecify rules and include negative outcomes. Computational comparisons need subsequent functional assays before supporting biological replacement claims.

There is no customer interview, paid pilot, repeat-purchase or unit-economics evidence in this repository. Commercial potential is a hypothesis. Project cost must include compute, failed/repeated runs, human review, deployment/support, storage and third-party rights. The original-code and third-party/model licensing scopes require separate review; no licensing terms are changed by this assessment.

## Sources

- [Peptide2Mol original preprint](https://arxiv.org/abs/2511.04984)
- [9NU4 experimental structure and primary literature](https://www.rcsb.org/structure/9NU4)
- [LiveDesign product description](https://www.schrodinger.com/platform/products/livedesign/)
- [Orion product description](https://www.eyesopen.com/orion)
- [GLP-1R signaling study](https://www.nature.com/articles/s41467-023-41893-4)

Vendor descriptions establish their stated product scope, not independent comparative performance. None of these sources validates this project's molecules or business model.
