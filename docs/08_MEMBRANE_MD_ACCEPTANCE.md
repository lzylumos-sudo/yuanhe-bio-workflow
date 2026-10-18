# Membrane MD acceptance: Ubuntu / RTX 5090 D, 2026-10-18

Release 0.3.0-rc2 closes the previously unexecuted screening-to-MD path for one real MurJ case. It is a bounded engineering acceptance, not general target validation, a converged affinity calculation, or experimental evidence. The existing upstream run was reused after hash verification, not regenerated. Its Peptide2Mol output remains explicitly partial. See the [earlier upstream record](06_RTX5090_ACCEPTANCE.md) and [independent module contracts](07_MODULE_CONTRACTS.md).

## Case and chain of custody

- Run: `20261018_054111_724ddc30`; handoff: `H_9422d94a82ca`; system: `MD_a2937bd0c6b3`.
- Molecule: `CMPab5e7a712d0bfbe3fca5`; screening compound: `CMPCF93BB6EAD737125`; pose: `POSDE0BE9245FD96FF8`.
- Isomeric SMILES: `CCc1nc(CC[C@@H]2CN(Cc3ccccc3)C[C@@H](O)[C@H]3O[C@H]23)cn1C`.
- Selected pose SHA-256: `29f7444628beda299c652a8fa6ce101814f0533c7a1e021f398100de766c0ab8`.
- Prepared receptor SHA-256: `0cba0d5a841112e467460955b0c62ebac8e871c07b2fdcc3eaa7ebe4d5f2864e`.
- 9NU4 MurJ chain A residues 3-500, 498 residues; original construct mutations retained. The ligand is neutral in its screened state; comprehensive protonation/tautomer enumeration was not performed.
- Formal screening candidates: **0**. MD-assessed exploratory molecules: **1**. Scientifically validated molecules: **0**. The UI required explicit exploratory acknowledgement; qualification remains `preliminary`.

The SDF stereochemistry is checked against the screened identity. Protein and ligand receive the same rigid transformation. PPM3 orientation was fitted to the exact receptor, with 3,785 matched heavy atoms and 0.00049822 angstrom fit RMSD. Center: [166.2042503493889, 164.11886069466075, 138.85715195277874] angstrom; membrane normal: [0.3434574946087956, 0.14688643067548426, -0.9276105464474116]. Source inputs, outputs and alignment evidence are retained in the delivery.

## Actual system and protocol

126,496 atoms: 7,799 protein atoms, 55 ligand atoms, 462 POPE lipids (PDB name POP), 20,255 waters, 55 sodium and 72 chloride ions. Amber14/Lipid17, TIP3P and GAFF 2.11 / AM1-BCC; 310 K and nominal 0.15 mol/L salt. Final net charge is approximately 7.55e-11 e. Pure POPE is an idealized model, not a complete bacterial membrane composition.

Energy minimization converged after 1,321 steps (maximum force 951.64008 kJ/mol/nm; threshold 1,000). This was followed by 0.1 ns NVT, 0.5 ns NPT with 1,000 kJ/mol/nm^2 position restraints, 0.5 ns NPT with 100 kJ/mol/nm^2 restraints, a separate 20 ps checkpoint continuity test, and **1 ns unrestrained production** at 2 fs/step. Equilibration uses semi-isotropic C-rescale pressure coupling; production uses semi-isotropic Parrinello-Rahman. Temperature coupling is V-rescale. Hydrogen bonds and rigid waters use constraints/SETTLE.

A single GROMACS rank with 12 OpenMP threads offloaded PP, PME, bonded forces and coordinate updates to the RTX 5090 D. Production speed was 396.133 ns/day, 218.109 seconds wall time. The entire GROMACS stage including its analysis took approximately 522 seconds; this excludes system construction, debugging, packaging and transfer. All five final phase logs contain zero LINCS warnings and zero fatal errors.

Software: Ubuntu 22.04; GROMACS 2026.4 CUDA/sm120; OpenMM 8.5.1; OpenMM ForceFields 0.15.1; OpenFF Toolkit 0.18.0; ParmEd 4.3; MDAnalysis 2.9; RDKit 2026.3.1. The builder uses CUDA double precision with a CUDA 13.1 runtime compatible with driver 595.71.05.

## Observations, not efficacy claims

101 trajectory frames cover 0-1 ns at 10 ps intervals. Protein-backbone least-squares alignment is applied to the entire solute; the ligand is not fitted independently.

| Observation | Result |
|---|---:|
| Final backbone RMSD | 1.296 angstrom |
| Maximum backbone RMSD | 1.640 angstrom |
| Final ligand RMSD | 3.325 angstrom |
| Maximum ligand RMSD | 3.651 angstrom |
| Written-sample mean temperature | 310.055 K |
| Written-sample mean pressure | 2.559 bar |
| Written-sample pressure SD | 114.658 bar |
| Mean density | 1030.909 kg/m^3 |
| Phosphorus-layer thickness proxy, first / last | 41.290 / 41.280 angstrom |

Thermodynamic CSV means use 1,001 written samples at 1 ps spacing. GROMACS's accumulated energy statistics use 5,001 points and report mean pressure -0.02717 bar; these are different sampling summaries, not interchangeable estimates. Neither establishes long-time convergence. The ligand moves substantially relative to the protein; the result does not demonstrate a stable binder. Distance-contact analysis uses all protein atoms (including hydrogen) versus ligand heavy atoms, periodic distances and a 4.5 angstrom cutoff; it is not a hydrogen-bond analysis. CA RMSF, raw thermodynamics, first/last structures and PNG/PDF figures accompany the raw trajectory.

## Executed acceptance and fixes

Real browser clicks selected the candidate, reviewed membrane orientation, launched parameterization, previewed the built bilayer, reviewed the system, started stage six, inspected first/last trajectory frames and plots, and generated the full run ZIP. Export contained 920 files and zero omissions. Terminal operations deployed code, diagnosed failures and verified artifacts; they did not substitute for the scientific launch clicks.

A separate checkpoint probe executed 5,000 steps, then resumed the same checkpoint with `-cpi -append` to 10,000 steps. Production itself was uninterrupted. During production, the local SSH tunnel was deliberately disconnected; the browser showed connection failure, the same remote parent process remained running, and a new tunnel restored the same run. This validates client disconnection for this case, not every possible server outage.

Five unsuccessful system-construction attempts remain in the run. Corrections cover rigid-water geometry during cross-engine export, stable membrane insertion relaxation, compatible CUDA runtime, lipid residue naming, bounded lipid charge rounding with corresponding 1-4 exception products, and GROMACS 2026 trajectory topology compatibility. See module contracts for exact policies. Input/output hashes, engine binary hash, phase TPR hashes, source snapshots and continuation revisions prevent silently mixing changed inputs with checkpoints.

The actual membrane environment passed **57 targeted regression tests** in 1.27 seconds: MD contracts, studio integration, workflow regression and screening revision. These include real OpenMM/ParmEd conversion checks and mocked orchestration tests; they are separate from the real UI MD evidence. A broader legacy collection was interrupted after 19 passed and 1 skipped during an unrelated expensive smoke test, so full-suite completion is not claimed. Docking, screening and MD source archives are independently exportable; packaging/import checks do not constitute new scientific runs.

## Remaining scope

The single-case integration is complete, but universal correctness and publication readiness are not established. Windows/WSL/RTX 5070, a second host, other target classes, mixed-lipid models and long independent replicas remain untested. Formal screening still needs its configured target-region and additional qualification evidence. Peptide2Mol needs a complete successful generation run rather than partial-output recovery. Target-specific scientific work requires justified protonation/construct/membrane choices, appropriate controls, replica/convergence analysis and experimental or otherwise adequate validation. Increasing candidate count alone does not repair missing qualification evidence.
