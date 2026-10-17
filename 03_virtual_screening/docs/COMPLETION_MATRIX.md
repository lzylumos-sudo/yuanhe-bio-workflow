# Historical completion matrix

DONE means implemented; VERIFIED means evidence existed in the recorded historical environment; BLOCKED means an external prerequisite was missing. These labels are not acceptance of the current release or another host.

| Area | Historical result | Evidence boundary |
|---|---|---|
| MurJ receptor audit | Verified | murj2: 4,824 deposited atoms; chain A 4,598 = 3,785 selected + 813 fusion atoms; selected range 3-500, 498 residues; 104 fusion residues removed; original complex retained |
| WT restoration | Not evaluated | K5F/S12I/M13A retained; no rotamer reconstruction |
| Formal qualification | Verified software gate | gq5: zero qualified primary, two preliminary, one exploratory; missing local-region/redocking/membrane-source evidence retained |
| PoseBusters scratch repair | Verified historical execution | Native file loading failed in a non-ASCII fallback directory; scratch capability entered the cache key. Real pose passed 22 checks; translated/clashing, stretched-bond and overlapping-atom controls failed appropriate checks |
| GNINA on that WSL host | Blocked | RTX 5070 detected; CUDA libraries remained incomplete, specifically libcusparse.so.12 |
| Redocking adapter | Implemented, real execution pending there | Planned, attempted, succeeded and supporting counts separated; RMSD uses unchanged receptor frame |
| MurJ trial | Blocked at that time | Reviewed local box and executable GNINA missing; no target poses or candidates fabricated |
| Handoff equivalence | Verified offline | handoff1: 2/2 bundles, eight checks, coordinate agreement <0.001 A; old t6 bundle retained missing protocol/units evidence |
| Membrane evidence | Proxy only | No sourced orientation or actual membrane construction in the screening module |

| Implemented area | Historical evidence |
|---|---|
| Discovery and pose import | t4: 313 files, 60/60 poses, no import failures; metadata and duplicate pose sources excluded |
| Score binding and units | Value/identity-based binding; CNNaffinity remains pK, not kcal/mol |
| Identity | CMP/STA/JOB/POS and deterministic selection identifiers |
| QA and contacts | 60 poses; later PoseBusters run checked 22 criteria per pose |
| Mapping | 498 SIFTS mappings, mutations at 5/12/13, 13 missing native sites |
| Ranking and diversity | t4 tiers B=18, C=2; missing components explicitly listed |
| Properties/strain | All 20 compounds assessed; failed/nonphysical energy remains not_evaluated |
| Budgets | <=200 redocking, <=50 deep, <=12 handoff, <=2 poses per parent |
| Candidates and exclusions | Reason codes, raw rankings and qualification reports |
| Reports | Markdown/JSON plus candidate/ranking/exclusion CSV |
| Handoff | Three historical bundles; checksums and actual blockers |
| CLI and read-only web | doctor/run/serve/runs/report and local HTTP smoke checks |
| Cache | Input/configuration/code keys; selective reruns preserve other stage summaries |
| Historical test suite | 29 passed: 26 unit and three integration |

The t4/t5 generic-target runs had identical candidate/ranking bytes and selection version `SEL6830370C0C42E694`. Their earlier primary labels predate stricter qualification. They establish neither MurJ screening nor membrane MD nor successful target redocking. Current workbench construction code is a later addition and still requires target-machine acceptance.
