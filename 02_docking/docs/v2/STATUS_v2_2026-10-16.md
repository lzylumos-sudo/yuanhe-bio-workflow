# Historical implementation status: 2026-10-16

This record concerns the earlier standalone v2 prototype. Missing dependencies, absent target validation and untested features were not counted as passing.

Implemented and tested: relocatable paths, atomic writes/checksums/tamper detection, state transitions, hash-bound approval ledger and queue gates, receptor inventory, ligand CSV provenance and isolated failures, box/protocol hashes and argument-array commands, ranking/aggregation and explicit unavailable enrichment metrics.

Real 1DRF preparation repaired GLN12, added 1,511 hydrogens and recorded 0.27 A backbone drift. Water/heteroatom removal required explicit approval; missing loops were only reported. No pKa prediction or validated retained-metal/cofactor route was provided. Ligand state generation used coarse rdkit_rules_v0 acid/amine rules, not MolScrub; 100 parents produced 132 states with no recorded Gasteiger charge failures. These rules do not replace a validated protonation workflow.

Vina 1.2.7 arm64 with Meeko 0.8.0 executed in independent job directories with raw output/hashes. One interrupted run reused 11 verified jobs and quarantined one incomplete job before rerunning the remainder. GNINA/CNN/GPU was not executed on that arm64 host; generated command arrays alone were not acceptance.

FOL redocking used three conformers, two states and three seeds: 18 jobs. Top-1 RMSD ranged 1.88-3.10 A, with only 1/18 <=2 A. Any-pose success <=2 A occurred in 6/18, minimum 0.93 A. This is a Vina single-target result, not GNINA or enrichment validation.

The prototype did not implement its own web UI, SSH setup, WSL installer, PoseBusters integration, full report/export route, covalent/flexible/metal docking or MD. Later main-module features do not retroactively change this record. Historical baseline inventories and generated run artifacts are not bundled in GitHub.

The development environment used Python 3.11 and RDKit 2026.03.4. exFAT lacks the symlink/POSIX behavior required by many environments; keep high-frequency runtime data and environments on the compute host filesystem.
