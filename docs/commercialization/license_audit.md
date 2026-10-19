# Preliminary license and intellectual-property audit

2026-10-19; baseline 841652e. Technical due diligence, not formal legal advice. No license or attribution is changed. Review scope: reachable source history and published primary references; not employment agreements, signed assignments, full runtime SBOMs or patent opinions.

## Current grant and commercialization

The root LICENSE is PolyForm Noncommercial 1.0.0; its digest matches deployment/license_provenance.json. LICENSE_SCOPE.txt limits coverage to contributions whose holders can license them and preserves third-party and earlier rights. Its effective-version wording remains 0.3.0-rc1 while the app is 0.4.0-rc1: record this provenance ambiguity for rights-holder review; do not silently amend a legal notice.

The standard license permits noncommercial purposes and specified organizations; the exact text controls eligibility. Source visibility alone does not grant unrestricted commercial use. Rights holders may separately authorize what they own, but cannot relicense third-party code, weights or data merely by selling a Yuanhe license. Dual licensing is a future option only after chain-of-title and dependency compatibility checks. [PolyForm primary text](https://polyformproject.org/licenses/noncommercial/1.0.0)

## Historical grants and file scope

`git show adfdfcb:03_virtual_screening/pyproject.toml` declares ZVSS 0.9.0 project code MIT. This declaration exists before 9d71299 introduced the root noncommercial license and explicit preservation of earlier rights. f83cbed includes a subsequent tracked change in this history; inspect content versus file-mode changes separately. The first reachable commit has no root LICENSE: absence does not imply public-domain or unrestricted permission for other subtrees.

The historical ZVSS subtree has 96 files, enumerated in [historical scope](historical_mit_scope.txt). This is a review universe, not a finding that every data file or third-party item was MIT-licensed. Compare original project code, copied upstream material and later protectable changes file by file before pricing exclusive rights. Root LICENSE_SCOPE.txt expressly preserves the earlier grant. MIT grants include commercial reuse subject to notice conditions; changing later terms does not simply withdraw rights to previously licensed copies. [MIT text](https://opensource.org/license/mit)

The history contains one author identity across 15 commits. Commit authorship is not proof of exclusive ownership, absence of contributors, or employer/university rights. Obtain contributor origin statements, university/lab policies, funding/contract terms and assignments before granting exclusivity. Earlier local distributions and external repositories may predate this Git history.

## Dependency rights matrix

Code rights, model rights, training-data rights, hosted use and redistribution are separate checks. Risk H means commercial delivery blocked pending evidence/review; M means obligations need package-specific review; L means relatively straightforward at source level, not an unconditional clearance.

| Component and evidence | Code | Weights / data | Commercial deployment | Redistribution | Risk |
|---|---|---|---|---|---|
| Peptide2Mol; pinned 2c5af962 in deployment/model_dependencies.json | Upstream advertises MIT; local notices record license hash | PMT_major.ckpt hash recorded; release page supplies checkpoints but explicit asset-specific grant and training-data chain not established here | Code permission does not clear checkpoint/data | Preserve MIT notices; obtain asset terms before bundling | H |
| Pocket2Mol nested dependency | Local notices record MIT and hash; vendor tree excluded | Any restored model or dataset needs separate inventory | Exact restored version/dependencies unresolved | Preserve notices and inspect nested dependencies | H |
| REINVENT4 / LibInvent; pinned 3cf206ea | Apache-2.0 per upstream and original publication | libinvent.prior SHA recorded; exact Zenodo record 20701824 could not be retrieved during this audit; do not substitute another DOI | Code generally permits commercial use under terms; optional licensed scorers and model rights remain separate | Include license/NOTICE as applicable; record modifications and exact prior license | H for model delivery; M for code |
| GNINA | Apache plus GPL notices; upstream explains OpenBabel requires GPL for its current combined build | CNN model/training inputs require per-artifact provenance | GPL is not a noncommercial prohibition; review exact executable and linked components | Do not claim Apache-only for a build retaining OpenBabel; satisfy applicable source/license obligations | H for proprietary bundled runtime |
| GROMACS | LGPL-2.1 per upstream; verify installed build/version | Force fields and parameter sets are distinct assets, not cleared by engine license | Commercial use can be permitted under LGPL terms | Preserve notices, modifications/source and applicable linking/relinking obligations | M |
| RDKit | BSD-style upstream license with additional bundled notices | SA-score data or other contributed assets need separate attribution | Permissive code subject to notices | Audit exact wheel and included third-party files | L/M |
| OpenMM / OpenMMForceFields / OpenFF / ParmEd / AmberTools | Mixed packages: inspect exact installed versions; OpenMM commonly has permissive core, AmberTools is not one blanket proprietary Amber grant | Force fields, charge models and parameter data separate | Not cleared from package names | Exported parameter files and container redistribution require package-level review | H pending solved runtime inventory |
| MDAnalysis / ProLIF / PyMOL / PoseBusters / Meeko / OpenBabel | Actual environment manifests/imports and bundled Meeko notice establish use; GPL/LGPL/permissive and distribution variants must be resolved individually | Reference databases and viewer assets separate | Importing a package versus spawning a separate tool changes integration analysis, not automatically rights | No blanket noncommercial relicense of a combined copyleft derivative | H pending integration review |
| PyTorch / CUDA / numerical and web stack | Version-specific notices; CUDA has separate vendor terms | Downloaded models not covered by framework license | Confirm deployment hardware/runtime terms | A container is a distribution artifact needing its own SBOM | M/H |
| PDB 9NU4 / 1DRF, PPM3 outputs, papers | Coordinates/provenance under data/structures; papers not bundled | Structural data, server-use terms, annotations and publication figures have distinct rights | Public access is not exclusive ownership of structures or compounds | Verify data source terms and citations; do not add papers or private outputs | M |

Primary references checked on the audit date: [Peptide2Mol source](https://github.com/BLUE-Flowing/Peptide2Mol), [checkpoint release](https://github.com/BLUE-Flowing/Peptide2Mol/releases/tag/v1.0), [REINVENT4 source](https://github.com/MolecularAI/REINVENT4), [original REINVENT4 paper](https://doi.org/10.1186/s13321-024-00812-5), [GNINA source/licensing explanation](https://github.com/gnina/gnina), [GROMACS license overview](https://www.gromacs.org/about.html), [RDKit license](https://github.com/rdkit/rdkit/blob/master/license.txt). Current upstream pages are not substitutes for pinned-version artifact licenses; two pinned license-page fetches failed. Existing notices are evidence leads, not independent verification of excluded vendor files.

## Risk register and required decisions

| Priority | Risk | Required evidence / owner | Gate |
|---|---|---|---|
| High | Historical MIT permissions and uncertain title | Maintainer plus IP counsel: historical file map, assignments, institution/funding agreements | No exclusive license for uncleared assets |
| High | Weight/training-data terms unresolved | Model integrator: exact asset license, source, hash, acquisition date and dataset obligations; seek upstream clarification if absent | No commercial model bundle or unqualified commercial inference promise |
| High | Combined GPL/LGPL and proprietary distribution | Release engineer plus counsel: solved SBOM, process/import/link map, source-offer and notice plan | No proprietary runtime shipment until cleared |
| High | Prior public disclosure / invention ownership | Inventor and counsel: dated disclosure ledger covering Git, presentations, demos and papers | Review before publishing new invention details |
| Medium | Runtime version drift and omitted transitive assets | Release engineer: lockfiles, binary hashes, force-field and viewer inventory | Per-release review, not one permanent approval |
| Medium | Public security scan incomplete | Maintainer/security reviewer: authenticated secret alerts, release/deployment review, broader scanner | No clean-bill-of-health claim |
| Medium | No signed commercial terms or legal entity verified | Founder and counsel: actual contracting party, support scope, customer-data ownership | Do not invent an incorporated company or commercial availability |

Public code cannot be made secret again simply by moving it. Disclosure may affect patent novelty; national rules and exceptions differ, so no universal grace period is assumed. Public copyright licensing is separate from patentability and freedom to operate. The PolyForm text also includes a patent grant within its terms. Record disclosures and obtain jurisdiction-specific review before publication. [WIPO patent FAQ](https://www.wipo.int/en/web/patents/faq_patents)

## Proposed rights process, not a license change

Maintain one asset register per release: file/artifact, origin/version/hash, holder, code license, model license, data license, usage purpose, distribution form, required notices, evidence and reviewer. Start with current models and the actual deployed environment. Retain source-available research licensing for eligible contributions; commercial exceptions apply only to cleared owned contributions. Separate executables are a useful architectural boundary, not a legal safe harbor. Leave LICENSE, LICENSE_SCOPE and all third-party notices untouched until explicit approval.
