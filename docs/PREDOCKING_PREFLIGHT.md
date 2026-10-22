# Pre-docking toxicity and synthesizability screening

The current policy is synthesis-only. See the [complete configuration and model guide](../02_docking/docs/PREDOCKING_V2.md).

Execution: Peptide2Mol -> REINVENT4 or imported candidates -> standardization and exact deduplication -> structural alerts -> real ADMET-AI -> SA and SCScore -> bounded AiZynthFinder search -> independent gates and diversity ranking -> selected ligand preparation -> existing Vina/GNINA -> screening -> MD.

## Configuration and UI

New research tasks use `00_workbench/configs/research_synthesis.json`, engine 3.0.0 and protocol `murj-discovery@1.2.0`. Configure endpoints, thresholds, SA/SC limits, mode and route budget in the existing Pre-Docking Screening panel. Advanced JSON sets model paths, template/building-block assets, versions/licenses, search depth/time and route reviews. Imported CSV/SDF candidates enter the same gates. No procurement configuration or service remains.

The internal `availability` stage/directory and serialized key are retained for task navigation and historical result readers; their new contents assess synthesis routes only. Existing task artifacts and old immutable protocol snapshots are retained. Schema-1 proofs cannot qualify new docking; explicitly re-screen a new task using the current defaults. Existing results are not relabeled or deleted.

## Route qualification

A complete real search is required, with an exact target, at least one reaction, all terminal nodes in the versioned planning pool and successful replay of every SMARTS reaction. A pool is a computational reference set. It does not certify actual experimental synthesis. Stereo or template selectivity ambiguities require a route-hash-bound chemistry review. Missing models, incomplete trees, invalid templates and exhausted budgets remain review/error states.

All three modes retain the mandatory route gate. SA and SCScore thresholds are separate difficulty filters. Neither score can replace a route. Source searches and asset hashes, every rejected molecule, all endpoint predictions and identity provenance persist in JSON, CSV, SQLite and SDF exports. Docking rechecks raw search hashes, policy, selected IDs and exact structure identity.

## Commands and recovery

`python -m yuanhe_preflight.cli --help` describes candidate imports and explicit configuration. Run with `PYTHONPATH=02_docking/src` in an installed chemistry environment. `deployment/install_predocking_assets.py --help` describes pinned official model installation. The workbench uses the same stages and existing task cancellation/retry/recovery. Evidence revisions preserve prior attempts.

The final wet-lab handoff requires `synthesis_reviewed` and `synthesis_evidence` (project-local path and SHA256), plus existing identity, structure, MD and manual-review checks. This is technical review of a proposed synthesis, not procurement confirmation or experimental proof.
