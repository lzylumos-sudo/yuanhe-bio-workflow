# Historical screening test record

Historical environment: Windows, Python 3.13.5, separate virtual environment and a non-ASCII project path. The then-current suite recorded 29 passes: 26 unit and three integration tests using a two-compound subset of actual generic-target GNINA output. This is not a rerun of the current English release.

Coverage included deterministic identities and noise-tolerant coordinate digests, configuration limits, canonical chemistry, cLogP explicitly distinguished from unavailable logD, symmetry-aware RMSD without fitting away translation, score units and missing values, actual GNINA log parsing, SIFTS mappings, deterministic repeated selection and explicit blocked handoff contracts.

The 9NU4 mapping test recorded 498 mapped residues from a 614-residue construct sequence, mutations at author positions 5/12/13 and missing native positions 1-2 and 501-511. Coordinate-bearing residue counts and full construct sequence counts are different quantities.

| Run | Historical observation |
|---|---|
| t4 | 20 parents x 3 modes, 60 imported poses from 313 files; QC 56 pass/4 flag; tier B=18/C=2; two primary and one exploratory under the earlier rules; three handoff bundles |
| t5 | Same selection version and identical candidate/ranking bytes as t4 |
| murj1 | Structural mapping only; 4,598-atom pre-crop receptor, membrane proxy and blocked MD; no MurJ pose screening |

Local HTTP smoke tests returned 200 for health, run listing, home, run/report pages and run API; unknown runs returned 404. Browser automation was not performed in that record.

The original report had not exercised successful PoseBusters execution because of native path loading. The later ASCII-scratch repair and actual positive/negative checks are recorded in `acceptance/POSEBUSTERS_RUN.md` and supersede that specific limitation. Multiseed GNINA, membrane MD and target pose screening were still unverified in that historical environment. Use the root validation record for current engineering checks.
