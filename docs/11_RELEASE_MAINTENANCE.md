# Release acceptance and maintenance

Baseline agreed on 2026-10-18: five acceptance categories below define the current closeout scope. They are verification categories, not five proven defects. A category is complete only when its required evidence is available. Unknown, blocked and untested results must remain distinct from failures and passes.

## Supported product boundary

The proposed deliverable is a reproducible computational research workbench and reusable docking, screening, membrane MD and peptide-contact modules. The demonstrated compute environment is the recorded Linux / RTX 5090 D configuration, with browser control from the local computer. Other platforms or target classes are not supported merely because the UI can accept their inputs.

Release acceptance does not establish peptide-function replacement, experimental activity, equilibrium convergence, regulatory readiness, customer demand or commercial success. A software delivery decision must state its tested environments, input classes and exclusions. The existing project and third-party licenses remain unchanged; passing tests grants no additional usage rights.

## Fixed acceptance checklist

| ID | Category | Evidence required to close | Current evidence and open work |
|---|---|---|---|
| A1 | Reproducible installation and deployment | Install the candidate release in a clean supported environment using recorded dependencies and pinned model sources; start the UI; verify compute discovery; reproduce a bounded reference case with stated comparison tolerances. Record versions, commands, hashes and results. | Clean chemistry installation and independent module entrypoints passed. Automated CI installs the recorded dependencies. Full clean deployment and reference-case acceptance of the final release remain open. |
| A2 | Representative cases and failure recovery | Freeze a case matrix before execution: the complete reference case, a second full supported target case, module-level controls, rejected inputs, zero-qualified-candidate results, process interruption and checkpoint recovery. Record expected outcomes and observed results. Do not count docking-only controls as a second full pipeline case. | MurJ reached MD; independent 1DRF docking controls and regression cases passed. Full second-target acceptance and final-version recovery checks remain open. The current server-access loss has not yet been accepted as a successful recovery. |
| A3 | Input, output and module contracts | Document supported input classes, units, identities, coordinate frames, qualification gates and versioned artifacts. Demonstrate that exported modules consume valid upstream outputs and reject incompatible or altered inputs. Declare unsupported cases explicitly. | Module contracts, formal-candidate handoff, integrity checks and four independent packages exist. Review final documentation against the tested support matrix and verify final exported artifacts. |
| A4 | Data protection, dependency rights and operations | Review access boundaries, credential handling, storage/export/backup behavior, dependency and model usage permissions for the intended delivery, installation and recovery instructions, known limitations and maintenance responsibilities. Attach findings and unresolved decisions. | Loopback UI / SSH use, source provenance, license notices and operational documents exist. Their existence alone is not a security assessment or commercial-use clearance. Complete the scoped review before a commercial delivery decision; do not silently change licenses. |
| A5 | Final real-UI acceptance | Deploy the candidate release, run the bounded supported workflow through the real browser UI, inspect outputs, complete prespecified repeat-quality analysis, export raw records and verify all included file hashes. Record the execution commit separately from any later documentation commit. | Two 5 ns repeats were complete at the last remote observation. Third-repeat status, aggregate verdict, final export and activation of the prepared release remain unverified after loss of remote access. |

Case inputs and measurable tolerances must be frozen before the corresponding test. Do not relax scientific criteria after seeing a result. A correctly rejected candidate or zero passing molecules may be a successful software test; it is not a positive scientific result.

## Evidence already available

Commit titles were normalized on 2026-10-19 without changing any commit file tree. Historical execution hashes below remain the actual original evidence identifiers; [the mapping](history/commit_renames_2026-10-19.csv) relates them to the renamed history. New maintenance changes are separate commits.

- Source revision `b9cf13214c90f8986c2ca578aa34ae9df730b6f0` passed [main-branch CI](https://github.com/lzylumos-sudo/yuanhe-bio-workflow/actions/runs/39093935326).
- Original JUnit artifacts for the same executable code at `056c448e27194ea7d794824a2d2860a43c717214` record 63 generation passes, 55 screening passes and 80 workbench passes, with one historical-output test skipped: 198 passed, one skipped, no failures. [Run and artifacts](https://github.com/lzylumos-sudo/yuanhe-bio-workflow/actions/runs/39093619141).
- [Integrated execution record](10_INTEGRATED_ACCEPTANCE.md) contains actual observations and remaining execution limits. Historical test totals must not be added to current totals.
- Full raw outputs and model weights are not in GitHub. An executed-source archive is not a trajectory archive. Until the final export is copied and verified, do not claim that remote raw data have been backed up locally.

## Maintenance procedure

1. Inspect the actual branch, working tree, deployment and run state before editing. Preserve existing inputs, raw results, failed attempts and checkpoints. Never duplicate an unknown remote job merely because its SSH connection was lost.
2. Record a defect only with specific evidence: affected revision and category, expected and observed behavior, reproduction steps, impact and logs or artifact hashes. Use the repository defect template. Speculative concerns are review questions, not confirmed failures.
3. Fix the evidenced cause and run checks appropriate to the change. Record the exact test revision and results. Scientific recalculation must preserve original records and use an explicit new run or audited revision.
4. Update the corresponding acceptance evidence and known limitation. Push reviewed changes, verify CI, synchronize the external-drive checkout by fast-forward, and confirm matching commits and a clean working tree. Deployment is a separate action and must be recorded separately.
5. Close a category only against its stated evidence. Any newly discovered defect must link to a reproduction and retest; any proposed expansion of product scope must be identified as a scope change, not silently added as another unexplained release gate.

## Release decision

The software can be considered for commercial delivery within its explicitly tested and authorized support scope only after A1-A5 are closed, no known blocking defects remain, and the release record identifies evidence, limitations and operational responsibilities. This is not a guarantee of zero undiscovered defects. The current status is **release candidate; final acceptance open**.

This document establishes the maintenance process. It does not create scheduled monitoring, a support subscription or an unattended deployment task.
