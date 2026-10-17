# MD installation evidence

Python: 3.11.17
Installed prefix: `/home/lzy/.cache/nexus-structural-biology-v1/environments/md`

Actual CPU test: PASS; OpenMM 1,000 harmonic-oscillator steps, 2 ps, finite energy. This is an installation test, not a protein-ligand simulation or scientific validation.
GPU: CUDA plugin detected; explicit GPU MD testing has not been performed. See the preparation-platform caveat below for the first molecular smoke run.

| Package | Installed version |
|---|---|
| openmm | 8.6.1 |
| openmmforcefields | 0.16.0 |
| openff-toolkit | 0.18.0 |
| ambertools | 26.0 |
| rdkit | 2026.3.1 |
| openbabel | 3.2.1 |
| pdbfixer | 1.12.0 |
| ParmEd | 4.3.1 |
| MDAnalysis | 2.10.0 |
| mdtraj | 1.11.1 |
| numpy | 2.4.6 |
| pandas | 2.3.3 |
| matplotlib | 3.11.2 |
| biopython | 1.88 |
| PyYAML | 6.0.3 |
| conda-pack | 0.0.0 |
| parmed | 4.3.1 |

Install completed and native CPU test passed before module automation development.
Molecular preparation, simulation, analysis and 10 ns acceptance are recorded separately; this report does not imply those checks passed.

## New molecular workflow test

`results/public_1fkj_smoke_v1`: PASS; 23,394 atoms, CPU four threads, fresh preparation and GAFF/AM1-BCC parameterization, minimization, NVT, NPT, 2,000 production steps = 4 ps, four trajectory frames. Final checkpoint hash/XML time and numerical observables passed. All seven requested analysis outputs were created.

Actual AmberTools MMPBSA.py ran on four production frames with igb=5/mbondi2, 0.15 M salt, entropy excluded. This is a software validation, not a binding affinity result.

Original OpenMM versus exported Amber isolated-solute energy difference: 9.0206922e-06 kJ/mol (tolerance 0.01), PASS. All committed build-file hashes remained unchanged after this check.

Approximate CPU throughput: 9.09 ns/day; 10 ns estimated 26.4 hours from this very short benchmark. Longer-run throughput may differ.

The 10 ns test is NOT yet completed by this report. No MD_TEST_PASS.md has been generated for the 4 ps test. Production MD and Amber-export energy checks explicitly used CPU. The v1 PDBFixer preparation stage did not explicitly pin or record a platform, so its GPU use is unknown; no claim of an entirely CPU-only v1 preparation is made. This was corrected before the v2 10 ns run: PDBFixer now explicitly uses CPU and records it in protein_preparation.json.

## Corrected CPU preparation and input contracts

The corrected PDBFixer preparation explicitly records CPU. Actual GNINA output passed protein-H rebuilding, precise ligand mapping, and heavy-coordinate preservation; a single embedded public 1FKJ complex.pdb passed without any chemical sidecar. Evidence: `installation_evidence/input_contract_tests.json`. These are preparation acceptance tests, not additional production trajectories.

Independent tests reject missing/duplicate/corrupt embedded metadata and corrupt committed trajectory/checkpoint records. The first HADDOCK interface sample contains phosphorylated NEP and was correctly rejected as unsupported chemistry; its residue is not silently changed. A standard-residue PPI sample is required for automatic MD compatibility.

## Actual HADDOCK single-PDB to MD integration

`results/integration_haddock_single_file`: PASS; only the newly docked 1BRS complex.pdb was copied, without SDF or manifest sidecars. Explicit CPU preparation, 36,028-atom system, CPU two-thread minimization/NVT/NPT and 2,000 production steps = 4 ps completed. Four trajectory frames, checkpoint/time audit and all applicable analysis passed. Ligand RMSD is an explicitly labeled N/A figure; MMGBSA.csv and optional metadata consistently record NOT_APPLICABLE.

A PDB four-column residue-number parser warning affects only 2,835 solvent HOH atoms; protein identities are unaffected. MD uses the original indexed System, and analysis restores exact atom-indexed bonds and residue labels from topology.json. This warning does not establish a protein identity mismatch.

## Actual GNINA single-PDB full MD integration

`results/integration_gnina_single_file_full`: PASS. A single actual GNINA 3PTB/BEN output PDB, including its checked embedded chemical definition, completed CPU preparation, a 29,981-atom system, CPU two-thread minimization, 1 ps NVT, 1 ps NPT and 2,000 production steps = 4 ps. Four trajectory frames and checkpoint/time audit passed; real AmberTools MM/GBSA computed four frames without entropy. The source, copied example and staged PDB have identical SHA-256 hashes.

The original GNINA receptor fixture omitted one calcium ion and 62 crystallographic waters when selecting protein ATOM records. Its compatible export later removed 1,579 receptor H/D records, retaining all protein and ligand heavy-atom positions. This MD stage removed zero existing protein hydrogens, regenerated hydrogens at pH 7.0, and did not remove or infer calcium. This is an execution test of the defined calcium-free trypsin construct, not a claim about native calcium-containing trypsin stability. See `installation_evidence/gnina_md_input_scope.json`.

## Visual analysis and independent 10 ns acceptance QA

Eight actual PNGs from the public 1FKJ short run and actual HADDOCK 1BRS integration were individually displayed and reviewed: RMSD, RMSF, ligand RMSD and contact maps. Units, chain/residue labels and occupancy colorbars are readable and unclipped; the PPI ligand figure explicitly states NOT APPLICABLE. Distance labels use ASCII A for Angstrom, and RMSF uses C-alpha index. Hashes and observations: `installation_evidence/visual_qa.json`.

The independent `scripts/validate_10ns.py` passed actual short-trajectory hash, DCD timing and CPU binary checkpoint reload checks, then correctly rejected 2,000 steps / 4 ps as insufficient for 10 ns. Exit status was 2 and no PASS or prior validation marker changed. Evidence: `installation_evidence/10ns_validator_rejection_test.json`. The script only publishes 10 ns PASS after completion; writing/testing the validator is not completion of the ongoing long simulation. No frozen MD core source changed for these additions.

## Actual Rosetta single-PDB preparation integration

`results/integration_rosetta_single_file_v2`: actual Rosetta-generated public 5T35 VHL-MZ1-BRD4 model, copied as one PDB with no sidecars, passed independent CPU preparation. It contains 451 protein residues across four chains and 69 ligand heavy atoms. The original 3,641 protein hydrogens were rebuilt at pH 7.0; all existing protein and ligand heavy-atom coordinates were retained exactly during preparation. Exact SDF/PDB heavy-atom mapping and stereochemistry checks passed.

The first attempt was correctly rejected because the upstream atom_mapping included hydrogens. The Rosetta exporter corrected that metadata to list only heavy atoms while preserving original SDF indices, all coordinates, and the complete SDF chemical definition. Both attempts are retained. Preparation success is not proof of subsequent parameterization or MD completion; each later stage must supply its own evidence.

## Actual Rosetta single-PDB full 4 ps integration

`results/integration_rosetta_single_file_v2`: PASS. Independent GAFF 2.11 / AM1-BCC in the MD runtime produced a 215,834-atom solvated system with 129 ligand atoms and net ligand Amber charge -0.00000000 e. The existing saved CPU two-thread configuration completed minimization, 1 ps NVT, 1 ps NPT, and 2,000 production steps = 4 ps, followed by real AmberTools MM/GBSA and all applicable analysis. Independent verification reloaded the CPU binary checkpoint and checked each DCD/CSV hash, frame/time continuity and analysis artifact hash.

Evidence: `installation_evidence/integration_observer_integration_rosetta_single_file_v2/completion_audit.json`. Only the upstream complex.pdb was copied. This is a short interoperability test; the module 10 ns validation status and MD_TEST_PASS.md are not modified by this observer.
