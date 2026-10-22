# Dual RTX PRO 5000 Blackwell research configuration

The original upgrade delivered code changes and static acceptance only. It did not execute unit or integration tests, model inference, Vina, GROMACS, checkpoint trials, GUI acceptance or performance tests. That delivery does not establish dual-GPU validation or scientific validity. Subsequent engineering test results must be recorded separately from target-host scientific validation.

## Shared configuration and defaults

`00_workbench/configs/research_synthesis.json` is the current authoritative definition; `00_workbench/research/config.py` validates UI/API/CLI settings. Historical dual-Blackwell and smoke snapshots are retained unchanged; new tasks use the current schema. The `smoke` profile still retains one 500 ns MD production trajectory and does not automatically start MD.

| Field | Default |
|---|---|
| Raw Peptide2Mol attempts | 5000 |
| Reviewed core limit | 200 |
| Total REINVENT4 attempts | 200000 |
| Unique generated parent limit | 50000 |
| External parent limit | 50000 |
| Initial screening parent limit | 100000 |
| Refined screening parent limit | 2000 |
| Detailed review limit | 50 |
| MD candidate limit | 8 |
| Production trajectories per candidate | 1 × 500 ns |
| Concurrent MD GPUs | 2 |
| Per-MD-job stage timeout | 60 days (execution guard, not an ETA) |
| Wet-lab target count | 4–8 |
| Vina workers / threads / total budget | 24 / 2 / 48 |
| Initial receptor conformers / states / seeds / exhaustiveness / modes | 1 / 1 / 1 / 8 / 5 |
| Refined receptor conformers / states / seeds / exhaustiveness / modes | 2 / 2 / 2 / 32 / 10 |

All counts are limits or targets. If only one receptor or chemical state is available for refinement, use one and retain the shortfall; never duplicate structures to fill quotas. Additional states require review, file hashes and independent pre-docking evidence through `ligand_states`.

## Actual call chain

Workbench creation → `workflow.create` → `research.config.resolve` → frozen configuration, hashes, environment and asset manifests.

Peptide2Mol → human core review → global REINVENT4 quota → parent standardization in `research.qualification.prepare_risk` → shared `yuanhe_preflight` risk screening → ADMET-AI toxicity -> SA/SCScore -> local AiZynthFinder route validation → two-stage official Vina in `research.docking.run` → independent Y-MAVS per receptor → structural review → existing membrane building → `research.md_queue` → existing standalone GROMACS/analysis entry points → human wet-lab review.

Current rules include curated SMARTS and full RDKit PAINS/Brenk; AI endpoints and their limits are in the pre-docking guide. Missing models, templates, planning building blocks or complete routes remain pending. SA/SCScore alone never qualifies a route. Vina does not start without qualified molecules.

The MD `build` wrapper binds the actual docking receptor. `build_prepared`, force fields, charges, membrane models, MDP, run, checkpoint_probe and analyze retain the integration baseline. The existing independent 20 ps checkpoint self-check is a backend recovery check, not a new MD screening stage.

## Reuse, recovery and evidence

Configuration cannot change in place. `/api/research/run/<id>/derive` accepts only changes to Top N, MD count, wet-lab count and selection mode. Generation and pre-docking evidence are validated before reuse. Vina reuse keys include input SDF, chemical state, receptor, search box, seed, parameters, software and adapter hashes. All cached outputs are revalidated; only cache misses invoke Vina.

Completed MD is reused through original run files and hashes without copying or rewriting raw trajectories. Incomplete MD resumes within the same run from its original checkpoint. Completed generation batches can be reused by hash. If the number of raw attempts at interruption is unknown, preserve the budget reservation and block resampling pending human audit; unlimited retries must not consume new budget.

GPUs use exclusive UUID locks; each MD subprocess sees one GPU. Priority does not kill other calculations. Insufficient resources reduce concurrency or block execution. ETA remains null without measurements. The research panel displays actual production time, Performance and checkpoint state from MD logs.

## Target-machine acceptance commands (not executed in the original upgrade)

Run from the repository root in the configured chemistry environment:

```sh
python -m pytest 00_workbench/tests/test_research_dual_blackwell.py 02_docking/tests/preflight -q
python 00_workbench/run_research.py preview --config 00_workbench/configs/research_synthesis.json
python deployment/preflight_doctor.py
python 00_workbench/run_webui.py --help
```

Run creation accepts the complete workbench configuration: `python 00_workbench/run_research.py create --config <complete-config.json>`. Preview does not create a run; creation does not execute computation. Starting a stage requires explicit `start --run-id <ID> --stage <stage>`. First supply actual execution paths, reviewed cores, receptor, search box and synthesis evidence.

## Outstanding work and target-host validation

- The original delivery had no runtime acceptance results or measured rates. Blackwell CUDA, models, GROMACS/RESP and resource detection require target-machine validation.
- Initial screening combines Vina scores with diversity round-robin selection over actual scaffolds; refinement retains Y-MAVS. Initial screening does not yet rank all candidates by conformational strain/key-contact multi-objective criteria.
- REINVENT4 currently divides the total quota evenly across cores and executes them sequentially. Yield-adaptive allocation and dual-GPU core parallelism are not implemented.
- CPU workers obey quotas; online tuning based on measured throughput, cross-node clustering and automatic preemption are not implemented.
- Incomplete generation batches do not automatically resume sampling. Preserve budget and raw output for audit. Failed Vina caches remain for inspection and explicit handling.
- Some stages retain candidate metadata lists. Pose files are saved per job, but not all historical JSON/CSV reads use database pagination.
- Delivery ZIPs include hash manifests for structures, reviews and original trajectories. Large trajectories remain in their original directories and must accompany delivery according to the manifest.
- Chemical states require explicit import and review; additional protonation/tautomeric states are not automatically predicted.
- Tests supplied in the original upgrade were not executed then. Results from another historical task do not certify this integrated version.

This research workflow upgrade is integrated into the original project but still needs target-environment validation and the outstanding capabilities above. It does not establish completion of all P1/P2 work.
