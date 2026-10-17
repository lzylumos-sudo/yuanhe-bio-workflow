# Architecture

CLI calls `stages.run_all(ctx)`; each stage writes SQLite records and `runs/<id>/stages/<stage>/summary.json`, with error details when applicable. The standalone FastAPI/Jinja view reads existing records and reports. Report regeneration merges all stage summaries so selective reruns retain other stages' evidence.

## Layers

- `ingest/`: discover/classify inputs, parse SDF/PDB poses, normalize score semantics and bind scores to poses, then generate identity records.
- `receptor/`: component policy and geometry checks; PDB/UniProt mapping through SIFTS, `_struct_ref_seq` offsets, then alignment fallback.
- Computation: `qa/pose_qa.py`, `mech/interactions.py`, `motion/interference.py`, `redock/gnina.py`, `strain/strain.py`, `membrane/environment.py`, `scoring/aggregate.py`, `select/pipeline.py`.
- Delivery: `md/package.py`, `md/adapters.py`, `report/render.py`.

| Identity | Inputs |
|---|---|
| CMP compound | Canonical isomeric SMILES hash |
| STA chemical state | SMILES, formal charge and protonation description |
| JOB docking job | Receptor version, box and protocol digest |
| POS pose | Job, record index, coordinate digest and rank |

Coordinate digests round at 1e-3 A to suppress formatting noise. Identity and score units must not be inferred from record order.

Execution (`pending/running/completed/failed/cached/skipped`), assessment (`pass/flag/fail/needs_input/not_evaluated/unsupported`) and selection (`eligible/excluded/exploratory/selected/not_selected`) are separate state families. An unperformed check never becomes pass.

## Screening evidence

1. Chemical identity and pose quality: PoseBusters plus explicit box, clash and bond checks; unavailable optional checks are unsupported.
2. Local mechanism: mapped contacts; missing reference ligand means contact Jaccard is not evaluated.
3. Conformational interference and deterministic pose clustering.
4. Multiseed GNINA redocking in the unchanged receptor frame; missing engine blocks evidence.
5. Score aggregation: normalize only over available components and list omitted weights.
6. Molecular properties and membrane geometry proxy; no claim of actual membrane orientation.
7. Strain and diversity: failed parameterization or nonphysical energies are not evaluated, never zero strain.

Evidence tiers precede within-tier percentile-weighted ranking. Requirements come from budget configuration. Compound IDs break ties. Pareto objectives use available components. Multiple poses are not multiple hits: select one representative per parent, using affinity, CNN and pose identity tie-breaking; other poses inform consistency.

Budgets cap redocking at 200 parents, deep assessment at 50, MD handoff at 12 and poses per parent at two. Greedy max-min diversity also enforces scaffold, similarity and membrane-warning limits. Every exclusion has a reason code.

Sorted paths/IDs, resolved YAML and cache keys `(stage,input_hash,config_hash,code_version)` support deterministic selection. Capability fingerprints and host/CNN settings participate in invalidation. This does not promise GPU bitwise reproducibility. Workbench live counts track unique molecules, QC progress and formal/preliminary selection separately.
