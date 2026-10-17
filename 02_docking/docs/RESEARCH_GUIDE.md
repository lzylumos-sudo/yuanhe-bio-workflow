# Research screening guide

## Receptor review

Choose author chains and model for PDB/mmCIF or externally prepared PDBQT. Review nonprotein components, alternate conformations, breaks, metals and nonstandard residues. Explicitly choose water/ion/cofactor retention. `keep_hetero` uses author numbering such as A:201; `altloc` selects a consistent conformation. Missing side chains require `repair_missing_sidechains=true`; missing loops are not rebuilt. Supported Meeko residue templates may be supplied, but added hydrogens do not establish exact pH-dependent microstates.

Define the box manually, from an actual reference SDF, or from reviewed pocket residues. Sizes are full side lengths in A. Reference coordinates are not re-embedded. Inspect the receptor/box before confirming. External PDBQT requires a corresponding prepared PDB and retains external/unknown protonation provenance.

## Protocol and chemical states

The standalone research_gnina preset uses a rigid receptor, Vina empirical score, CNN rescore, crossdock_default2018_1_3, exhaustiveness 16, nine poses and seed 20261012. It selects a previously tested CPU/GPU CNN path and does not silently substitute Vina or disable CNN. The unified workbench uses its own frozen run seed; record the actual protocol.

States default to preserve. Explicit enumeration example:

```json
{"ligand_states":"enumerate","ph_low":6.4,"ph_high":8.4,"max_states":8,"fragment_policy":"reject"}
```

Pinned Molscrub 0.2.2 rules enumerate acid/base and tautomer states. Canonical isomeric SMILES deduplicate and stably sort states; truncation is not population ranking. Multi-fragment inputs are rejected unless largest-fragment policy is explicit. Unspecified stereochemistry remains a warning.

research_refine is repeated docking: up to eight ETKDGv3 inputs with 0.5 A pruning, MMFF/UFF relaxation and up to three retained conformers; seeds 11/23/47, exhaustiveness 32. Conformer energies are not binding energies.

Each state/seed selects a representative across input conformers: affinity for Vina, CNN pose score for GNINA. Seed representatives are summarized by median; the parent uses its best state median. Keep every raw state/seed score. UI sorting does not overwrite the original ranking.

## Second stage and export

score_only scores existing same-frame poses; redock generates inputs and searches again; refinement minimizes and can change coordinates while retaining the original. Every new run saves receptor/box snapshots, original inputs, resolved settings, executable hashes, per-ligand commands, logs and outcomes. Cache reuse requires complete matching checksums.

Manual downstream export requires reviewer identity and pose approval. SDF is the chemical master; PDB is a same-frame preview. This module does not thereby parameterize or simulate MD.

Covalent docking, unvalidated metal/cofactor handling, flexible receptors, long-loop rebuilding and reliable MOL2 handling are outside the supported scope. No target activity labels or enrichment validation are supplied. A single FOL redocking example does not establish general screening performance.
