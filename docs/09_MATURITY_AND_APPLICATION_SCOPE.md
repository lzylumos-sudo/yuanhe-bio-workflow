# Application scope and maturity assessment

Assessment date: 2026-10-18. This is a research workflow under validation, not a commercially validated drug-discovery product. A successful calculation, software test, or short MD trajectory does not demonstrate functional replacement of a peptide, biological activity, or experimental hit-rate improvement.

## Intended initial application

The proposed initial users are teams with a target structure, peptide-binding evidence and experimental capability who need a reproducible workflow for generating and prioritizing small-molecule hypotheses. Candidate evidence, exclusions, provenance and reproducible inputs are deliverables. Functional replacement, oral availability and clinical equivalence are not established outputs.

The MurJ example uses a phage lysis-protein complex, not a marketed peptide drug. It is an antibacterial discovery hypothesis. No experimentally validated MurJ small molecule has been produced by this project.

## Current evidence and remaining scientific questions

The current engineering and computational evidence is in [integrated acceptance](10_INTEGRATED_ACCEPTANCE.md). It supersedes the earlier progress statements about unfinished membrane-orientation qualification, second-target docking controls, clean installation and repeat infrastructure. The historical records remain available.

The fresh 9NU4 run completed 16/16 generation attempts, eight unique REINVENT candidates and nine successful docking jobs. Corrected archived-pose assessment yielded two formal screening candidates. The formal-to-membrane handoff preserves pose identity and qualification provenance. Independent modules, a clean installation, second-target crystal-pose recovery and prespecified MD-repeat quality checks are implemented and accompanied by execution records.

Peptide-contact evidence is now explicit and reusable. In this run, every shortlisted molecule recovered only one of 22 peptide-contacting residues; contact-guided sorting did not change the ordering. This measured result does not support a peptide-function replacement claim or improved hit rate. Adding code around a metric cannot turn that negative result into biological evidence.

A complete computational path is valuable for generating and reviewing hypotheses. Validated functional replacement still requires a relevant assay, appropriate controls, independent replication and comparison with existing approaches. The repository contains no customer interview, paid pilot or unit-economics validation. These are external scientific and commercial evidence requirements, not software tests that can legitimately be marked passed.

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
