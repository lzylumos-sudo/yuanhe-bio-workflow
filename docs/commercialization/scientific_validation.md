# Scientific validation roadmap

Proposed protocol and acceptance gates, not new experimental results. Current record: [integrated acceptance](../10_INTEGRATED_ACCEPTANCE.md). Freeze protocols before new runs; obtain compute/experimental budgets separately.

## Evidence levels

E1 computational prediction: generated molecule, docked pose or trajectory. E2 computational benchmark: independently reproducible comparison against prespecified controls/baselines. E3 in-vitro evidence: authenticated compound and assay QC with replicates. E4 in-vivo evidence: appropriate approved study and pharmacology/safety context. Software tests are engineering evidence outside this ladder. No E1/E2 result is a validated drug, and E3 alone does not establish clinical efficacy.

| Work package | Protocol and metrics | Gate / missing evidence |
|---|---|---|
| MurJ / 9NU4 reproduction | Preserve PDB/hash, chain A/B selection, construct, protonation, pocket/box, model/weight/environment hashes, seeds and rejected attempts | Recover existing server state read-only before any continuation; verify last complete checkpoint and source revision. If unavailable, archive interruption and start a distinct approved run |
| Generation quality | All attempted, valid, unique, novel and scaffold counts; chemistry/charge/stereo and route checks | Compare baseline generator and reviewed-core expansion at equal attempt budgets; current 1/16 acceptance is a small case, not a general rate estimate |
| GNINA pose/ranking | Reproduce 1DRF control with same seeds and symmetry-aware receptor-frame heavy-atom RMSD; top-ranked and best-pose reported separately | Preserve recorded 2/3 versus 3/3 at 2 A. Add targets with known ligands; do not retune box or threshold to pass |
| Screening | Deterministic identity/pose handoff, rejection reasons, seed consistency, held-out actives/inactives when licensed | Report EF1%, PR-AUC/ROC-AUC and uncertainty only where labels support them; no MurJ enrichment claim without labels |
| Membrane MD | Reviewed orientation, topology, force field, charge, lipid/water/ions, equilibration, production and checkpoint continuity | Finish/inspect three velocity repeats and final export before aggregate quality verdict. Shared starting coordinates are not independent conformations |
| MD interpretation | Energy/temperature/density, constraints, ligand RMSD relative to original pose and production reference, contact persistence, membrane properties and block variability | Short stability alone is not affinity or convergence. Longer/adaptive sampling requires a scientific question and budget, not an arbitrary duration |
| Seed reproducibility | At least five pilot generation/docking seeds and the prespecified three MD repeats; retain every outcome | Report distributions and failed runs; same seed is not guaranteed bitwise across hardware |
| Resource accounting | Per-stage wall/CPU/GPU hours, peak memory/VRAM, retries, storage and active human time | Compare at equal oracle calls AND total cost; include setup and failure costs |
| Baselines/ablations | Manual toolchain, underlying generator, library screen, no-contact/no-AL variants | Freeze inputs/budgets; independent test evaluation, no self-scoring-only success |
| Target generalization | Pilot 3-5 targets across membrane and soluble classes, with independent known-pose/activity controls where available | Family/scaffold split and pretraining overlap audit; target-level reporting, not only pooled mean. 1DRF alone is not full workflow generalization |
| Synthesis | Versioned building blocks, routes, forward identity/stereo checks, chemist review and prospective synthesis QC | Separate predicted accessibility, reviewed route and actual obtained/pure material |
| Wet-lab interface | Compound/lot identity, purity, assay version, units, concentration series, controls, replicates, confidence/censoring and raw-data provenance | For MurJ, plan mechanism-relevant assays with qualified partner; antibacterial activity alone does not prove on-target MurJ inhibition. In-vivo work requires separate study design/approval |

## Reproducibility packet

Each benchmark archive should contain a signed-off protocol, machine-readable manifest, immutable inputs with allowed redistribution, full commands, versions/hashes, driver/hardware, seeds, budget, raw outputs, failed attempts, analysis script revision, unit-labelled tables and a clean-room reproduction procedure. License-restricted artifacts are referenced by lawful acquisition route/hash rather than redistributed. Customer or unpublished compounds are excluded from the public packet.

## Decision gates

G1 engineering: clean install, no lost identity/provenance, interruption/recovery and complete export. G2 scientific baseline: controls reproduced and failures explained. G3 differentiation: held-out target/budget comparison with uncertainty and negative results. G4 prospective assay: independently QC-checked measurements. G5 product use: external users reproduce results within a measured support budget. Keep separate verdicts; an open gate cannot be closed by passing another one.
