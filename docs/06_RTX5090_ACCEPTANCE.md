# RTX 5090 D acceptance, 2026-10-18

Status: partial real computational acceptance. Membrane MD remains unexecuted.

The authorized small run used one 9NU4 MurJ construct, one generated core and one GPU. Peptide2Mol started 16 attempts: 12 completed (one valid, eleven failed decoding) and four were interrupted at the attempt ceiling. The valid molecule was recovered by validating receptor, pocket, model and SDF hashes. It was not regenerated or replaced by historical output. LibInvent sampled eight records, producing five unique molecules. These five and the original seed yielded six GNINA jobs and 54 poses. Screening displayed three preliminary candidates and zero formally qualified candidates. No ligand parameterization, bilayer, equilibration, production trajectory or MD analysis was produced. The requested 20 ns was a maximum authorization, not a completed or finalized production length.

## Scientific evidence and limits

- MurJ chain A residues 3-500: 498 residues, 3,785 original heavy atoms. The peptide reference has 29 residues. Construct mutations were retained.
- Receptor preparation at pH 7.4 retained every original heavy-atom coordinate, adding hydrogens and terminal OXT. Ligand protonation/tautomer enumeration and synthetic feasibility were not evaluated.
- GNINA 1.3.3, exhaustiveness 16, nine modes, seed 17, 24 CPU threads, device 0. Real per-job GPU process observations and original scores were retained. The docking box in angstroms has center [159.115,153.366,136.322], size [24.324,21.754,22.565].
- Ranking representative poses have MMFF94 strain values 37.977, 56.7549, 40.8381, 48.5638, 282.3218 and 54.8898 kcal/mol. The 25 kcal/mol warning and 50 kcal/mol failure thresholds were unchanged. Three were flagged and three failed; a flag does not prove chemical invalidity.
- Strain fixes preserve heavy coordinates, relax regenerated hydrogens and compare with a converged local minimum of the same molecule. This heuristic is not binding free energy. All 54 poses were independently checked; the lowest diagnostic strain was about 31.97 kcal/mol.
- Formal qualification lacks configured target-region, multi-seed redocking and integrated membrane evidence. Multi-seed rescreening exceeds this run's bounded scope. Official PPM3 was separately run on the actual receptor (33.9 +/- 0.2 angstrom hydrophobic thickness); this is orientation evidence, not a lipid bilayer or MD validation.
- The generated oxygen-containing three-membered ring has unvalidated reactivity and synthesis. A conservative MD hold was recorded rather than treating preliminary display candidates as approved MD inputs.

The first screening attempt unexpectedly ignored disabled redocking, completing five extra redock tasks before the sixth was stopped through the UI. This is a recorded scope deviation, not part of planned acceptance. The fixed gate now checks disabled stages and zero parent budgets before probing or invoking GNINA.

## Repairs

1. Keep the virtual-environment interpreter path without resolving its symlink into the base interpreter; stream model stdout/stderr.
2. Bound vendor sampling by started attempts, including decoding failures; preserve an attempt audit.
3. Explicitly recover verified partial generation with input/model/output hashes and honest UI labeling.
4. Condition LibInvent v1 on its supported achiral vocabulary, then rejoin decorations to the original stereochemical scaffold with identity and order checks. Conditional NLL applies to the achiral model conditioning graph, not a new stereochemical model.
5. Merge and deduplicate generation products before bounded GPU docking.
6. Honor disabled redocking; preserve PDB coordinates while validating author-residue mapping against official mmCIF metadata.
7. Use the same representative pose for strain and ranking; record convergence and avoid labeling descriptive metadata as nonconverged.
8. Add audited screening-only restart: preserve old source snapshots, verify upstream artifacts, reject unrelated source changes, archive the previous screen, and record the allowed revision.
9. Fix named HTTP route captures, Linux input filename case, and expose SMILES, qualification and pose IDs in the UI.

Vendor sources and weights are not committed. After restoring the exact locally supplied vendor revisions, run `python deployment/apply_vendor_patches.py`. The manifest records before/after source hashes, refuses unknown revisions and permits idempotent reapplication. Restore instructions and license boundaries still apply.

## Verification

- Remote workbench/generation suite: 113 passed, 1 skipped.
- Final remote screening/restart suite: 36 passed, 6 skipped. Skips require independent historical reference/GNINA caches; they are not presented as passes.
- Front-end syntax: 15/15 ES modules parsed.
- Vendor patch preimage/postimage and idempotence checks passed. Bounded empty-decoding mock and original-stereo rejoining mock passed; actual eight LibInvent products were separately inspected.
- GROMACS 2026.4 CUDA sm_120 build: 105/106 CTest groups passed initially; the complex group failed automatic thread decomposition. With official regression arguments `-nt 2 -ntomp 2`, the complex group passed. Build tests do not validate this MurJ system.
- Actual browser: input preparation, model/core controls, docking, screening, cancellation, selected structure, refresh restore and export/download event exercised. Double-clicking rescreen created one task. A running screen survived closing the SSH tunnel; after reconnecting, its completed results were visible.
- Membrane production, trajectory display/analysis and checkpoint continuity were not exercised. Proposed GPU/12-thread/checkpoint settings in the membrane module remain unvalidated for the actual system.

## Reproducibility boundaries

The run ID is `20261018_054111_724ddc30`; two precursor failure/interruption runs are preserved. Each run has its original source snapshot. Screening repair revisions are stored separately and old screening attempts are archived. The final repository uses the portable English directory layout; original computational runs used the preserved deployment layout. Do not attribute all upstream results to the final application commit.

Private local delivery includes full raw runs, CSVs, failed attempts, PPM evidence, logs, environments, hashes and UI captures. `selected_smiles.csv` is header-only; preliminary display molecules are exported separately. Large models, trajectories, local connection configuration and credentials are excluded from Git.

Official references: [GROMACS membrane tutorial](https://tutorials.gromacs.org/docs/membrane-protein.html), [GROMACS installation](https://manual.gromacs.org/2026.4/install-guide/index.html), [GNINA releases](https://github.com/gnina/gnina/releases), [PPM3](https://oprlm.org/ppm_server3_cgopm).
