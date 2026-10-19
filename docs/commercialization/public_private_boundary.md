# Public research edition and future enterprise boundary

Local planning only. No license change, file removal, private repository creation, code migration, history rewrite or publication is authorized by this document.

| Area | Recommended boundary | Reason / dependency |
|---|---|---|
| Documentation, input/output schemas, public examples | Public documentation; existing code remains source-available | Helps reproducibility; preserve third-party attribution |
| Existing workflow adapters, qualification and reports | Maintain current source-available scope | Do not promise exclusive ownership of upstream or prior MIT code |
| New target-specific scoring/optimization policies | Review before any publication | May embody future commercial know-how; local implementation is not publication |
| Proprietary models, active-learning strategy, private training data | Future private implementation proposal | Requires provenance and contributor/IP review; no implementation migrated now |
| Customer structures, candidates, results and credentials | Private customer-controlled storage | Never place in public examples or demo sites |
| Enterprise scheduling, API, access control, audit and isolation | Planned commercial layer | These capabilities are not certified as implemented |

Reuse interfaces rather than duplicate scientific stages: shared input/pose/evidence contracts → replaceable backend adapter → execution environment → immutable results and decision review. Keep UI separate from computational executors. New enterprise services should consume stable contracts and be independently deployable; current refactoring is not required. Linux GPU execution, Windows/WSL launch and macOS review are distinct deployment paths; each needs its own actual acceptance.

A repository license does not automatically license weights or training data. Prior MIT grants cannot be assumed withdrawn. Contributor, school/laboratory and partner rights remain to be confirmed by rights holders. Root LICENSE_SCOPE.txt already records these issues; preserve it verbatim. A full historical, dependency/model/data and distribution review remains required before commercial bundling. Source availability is not technical protection of secrets. Deleting a public file does not make its history confidential again.
