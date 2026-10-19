# Code and content organization

This pass preserves all existing paths, public interfaces, scripts, historical documents, licenses and research outputs. Organization is through entry points, module boundaries and a tracked-file inventory rather than moving live scientific code. A path inventory is not a security audit or a certification that every function is correct.

| Layer | Maintained entry / source | Responsibility |
|---|---|---|
| Workbench | `00_workbench/run_webui.py`, `webui/workflow.py`, `pipeline_bridge.py` | Review, scheduling and stage handoff |
| Generation | `01_molecule_generation/murj_pipeline.py`, `src/` | Preparation, model adapters, core and chemical checks |
| Docking | `02_docking/` module guide and protocols | Receptor/ligand handling and GNINA execution |
| Screening | `03_virtual_screening/src/zvss/` | Evidence, qualification, diversity, selection and reports |
| Interaction analysis | `04_interaction_analysis/src/zlis/` | Interaction tables and structure-based visual outputs |
| Dynamics | Workbench `scientific/`; separate `05_molecular_dynamics/` | GROMACS membrane workflow versus independent OpenMM path |
| Deployment | `deployment/`, platform launchers | Dependency restoration and environment-specific installation |
| Evidence | `docs/04_VALIDATION.md`, `docs/10_INTEGRATED_ACCEPTANCE.md`, `docs/11_RELEASE_MAINTENANCE.md` | Historical validation and open acceptance requirements |
| Product | This directory and `docs/commercialization/` | Value, audience, roadmap and publication boundaries |

`TRACKED_FILE_INVENTORY.csv` lists tracked files from the baseline checkout plus this revision's added source/docs. It excludes untracked external-drive weights, environments, candidates and results. `docs/FILE_MAP.csv` is a historical migration/hash record and is intentionally not rewritten as a current inventory.

Changes to scientific selection are confined to the screening module and documented in its policy guide. Other modules receive navigation updates only; no global algorithm refactor, directory migration or cross-platform compatibility claim is made. Static parsing and link checks do not replace runtime or platform acceptance.
