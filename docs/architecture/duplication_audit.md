# Duplication audit

| Implementations | Class | Finding and action |
|---|---|---|
| Product version literals | A | Conflicting UI/package versions. Establish one platform authority and generated distribution metadata. Preserve historical version evidence. |
| workbench stage config merges | A | Scheduler and worker merge reviews separately. New runs freeze a single effective stage config and verify it before execution. |
| sbvs / sbvs2 | B/C | sbvs retains CLI, presets, preparation, runner and exports; sbvs2 has job/protocol/runtime components. Official workbench uses sbvs2 Vina. Add a public docking facade; retain both CLIs and scientific preparation variants. |
| Vina / GNINA | B/C | Different scoring semantics. New workbench already enforces Vina. Historical GNINA fields stay separate and absent values remain absent. |
| ZVSS / Y-MAVS | B | Y-MAVS is zvss.membrane, a subsystem of ZVSS, not a parallel complete platform. |
| research / scientific / automation | B | Budget/qualification/queue vs scientific backend vs standalone workflow. No safe wholesale merge. |
| OpenMM soluble / GROMACS membrane | D | 300 K / 1 atm versus 310 K / 1 bar semiisotropic; distinct construction, thermostat, barostat, cutoffs and analysis. Preserve protocols. |
| OpenMM builder / ParmEd / GROMACS | D | Consecutive transformation and execution steps, not redundant engines. |
| Preparation and chemical identity | B | Reviewed chemical-state preparation, generation normalization and parent identity have different scientific purposes. Preserve atom maps and stereochemistry. |
| Static RMSD / trajectory RMSD | B/D | Symmetry-aware pose comparison and protein-aligned trajectory RMSD have different references and atom selections. Never collapse by metric name. |
| ZLIS static contacts / MD contacts | B/D | Different geometry, atom definitions and time aggregation. Report definitions separately. |
| Scheduler implementations | B | Different jobs, databases and process leases. Keep compatibility; unified registry does not replace execution locks. |
| Configuration loaders | A/D | Shared governance and provenance can unify; engine-specific schema validation and units must remain. |
| archive/docking_prototype, historical pepwedge | E/C | Retain bytes and history. Lack of a static import is insufficient proof no external user invokes them. No files removed or moved. |

No deletion or physical directory migration is authorized by the audit. The safe initial changes are governance, shared entrypoints, effective configuration snapshots and truthful runtime capability reporting.
