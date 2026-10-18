# Integrated research acceptance, 2026-10-18

Run `20261018_091310_dea50e56` uses the original 9NU4 input, receptor chain A and peptide chain B, on the rented RTX 5090 D server. Research actions were initiated through the real browser UI. Shell access was used for implementation, installation, tests, artifact inspection and deployment. The execution source archive and versioned operation records distinguish the actual run from the subsequently packaged release.

## Generation and screening

| Step | Observed result |
|---|---|
| Prepared complex | 3,785 receptor heavy atoms; 258 pocket atoms |
| Peptide2Mol | 16 requested and completed attempts; 1 accepted, 15 rejected |
| REINVENT4 | 16 samples; 8 unique molecules; 0 adapter failures |
| GNINA | 9 molecules, 81 poses; 0 failed docking jobs |
| Deep screening | 27 redocking jobs, 3 seeds per molecule |
| Display shortlist | 7 compounds; 2 others excluded |
| Corrected formal qualification | 2 formal and 5 preliminary candidates |

The formal compounds are `CMP609F157C45053346` and `CMP82C801E536A70C1F`. Their existing strain results remain warnings rather than clean passes without alerts: the configured warning threshold is 25 kcal/mol and failure threshold is 50 kcal/mol. No threshold was relaxed after seeing the results.

The initial screening record reported zero formal candidates. An implementation error passed an invalid atom-map shape to RDKit RMSD comparison, causing all archived seed comparisons to be unavailable. The corrected comparison explicitly imports `rdMolAlign`, checks chemical identity and computes symmetry-aware heavy-atom RMSD in the fixed receptor coordinate frame without fitting the ligand independently. UI operation `1792289700350724415_8d646707` reevaluated the archived poses against each exact selected pose. Its SHA-256 is `27336b989a6657a82b7122b06414a19a61b0c80c5673b891ff78fbbe60abe144`. The original stage records remain unchanged. Qualification provenance accompanies the MD input.

The released redocking plan uses the same representative pose used for strain and selection. Cache accounting now includes physically verified outputs. A configured target-region predicate no longer silently becomes unavailable because of a missing property. Unknown or nonfinite evidence does not become a pass.

## Membrane orientation and handoff

An archived PPM3 result for this exact MurJ construct is transferred through matched receptor heavy atoms, file hashes and a proper rigid rotation. There are 3,785 matched atoms, coverage 0.999736 and fit RMSD 0.000498 angstrom. This is a sourced computational orientation for the specified construct, not an experimental membrane measurement. Other constructs do not inherit it automatically.

The selected formal candidate is `CMP609F157C45053346`, pose `POS625B279EF1A5579B`, upstream molecule `CMPe9b9c2fbd82571bb3892`. The UI-built parent is `MD_29eb12b2b0c8`, with 125,887 atoms, a reviewed idealized POPE bilayer, TIP3P water, 0.15 M salt, 310 K and GAFF-2.11 / AM1-BCC ligand parameters. Protein and ligand undergo the same rigid transform; molecular identity and copied-input hashes are checked before parameterization. The exported independent MD package validated the actual formal-candidate handoff in a clean environment.

Three UI-created velocity replicas use seeds 41, 97 and 17: `MD_258e6adedc11`, `MD_512d6c2c2c58`, `MD_a79cd96d2ad1`. Each has minimization, 0.1 ns NVT, 1 ns staged NPT, a separate 20 ps checkpoint continuity test and 5 ns production. Initial coordinates are shared; these are independent velocity repeats, not independent starting conformations. The prespecified `short-repeat-quality/1` protocol is frozen in each replica. The final per-replica results and quality verdict are recorded in the release evidence manifest.

The first replica completed 5 ns / 501 stored trajectory frames, with mean temperature 310.004 K, final backbone RMSD 1.535 angstrom and final ligand RMSD 2.464 angstrom. Its checkpoint test passed and no constraint warning was found. RMSD references the first production frame, so it does not by itself measure displacement during equilibration from the original docking pose.

## Modules, controls and ablation

The docking, screening, membrane MD and peptide-contact evidence archives each include file-level checksums and licensing. A new Python 3.12 virtual environment was installed from official registry wheels without system site packages; `pip check` passed. Docking and screening distributions installed in that environment; all four entrypoints imported successfully. Actual stage-six calculations run the exported MD core in isolated Python with no UI import path.

The independent docking module was invoked through the UI on 1DRF / crystal folate with seeds 17, 41 and 97. At a fixed 2 angstrom recovery threshold, best-pose recovery was 3/3 and top-ranked recovery was 2/3. Top-ranked RMSDs were 2.794, 1.919 and 1.919 angstrom; best RMSDs were 0.811, 0.815 and 0.865 angstrom. This is a second-target docking control, not a second full peptide-to-MD biological validation. The original operation's `clean_environment:false` field was a hardcoded metadata error; its actual request identifies the new interpreter and separate installation records establish isolation. The released detector measures the interpreter instead of hardcoding this field.

The peptide-contact module found 22 reference-contact residues at 4.5 angstrom. All seven shortlisted molecules shared one of these residues (recall 1/22). Geometric controls returned self Jaccard 1 and translated-peptide Jaccard 0. Adding contact-recall-first ordering did not change the affinity ordering. This is a measured negative ablation result: it provides no evidence that the current candidate pool preserves peptide function or improves prioritization. No active/inactive labels were available, so enrichment and activity-prediction accuracy are not reported.

## Reproducibility and software checks

- Full workbench suite in the clean chemistry environment: 76 passed, 4 skipped. Three skipped MD dependency tests were separately exercised in the actual membrane environment, where all six MD contract tests passed. The remaining skip concerns undistributed historical generated results.
- Screening suite: 55 passed, no skips after replacing an old machine-specific fixture with portable archived output. The offline test reproduces selection twice and verifies that excluded/zero-formal results do not silently create an MD handoff.
- Repeat quality, input tampering, peptide coordinate-frame checks and generic input binding: 12 targeted tests passed.
- The actual browser was disconnected from its local SSH tunnel and reconnected while production advanced from 2,244 to 2,478 ps. The remote computation continued. This is separate from each GROMACS checkpoint test.
- The real trajectory viewer was operated from 0 to 5 ns, and browser error logs were empty after recovery. Full exports preserve raw data; a compact review export is separate.
- Confirmed regenerable caches and unused website dependency installations were removed: 1,947 entries / 735,994,564 bytes. Lockfiles, unique environment images, weights, inputs, raw outputs and failures were preserved.

## Acceptance boundaries

This release supplies a bounded, traceable research workflow and modular components for the demonstrated Linux GPU setting. It is not evidence of a validated peptide-drug replacement, experimental activity, equilibrium convergence, broad target generalization, a paid customer outcome or regulatory readiness. The tested software features and measured negative results remain useful deliverables. The MurJ / phage-protein example is an antibacterial discovery hypothesis, not replacement of an established peptide medicine. See `09_MATURITY_AND_APPLICATION_SCOPE.md` for intended application and evidence boundaries.
