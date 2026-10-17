# Architecture and code guide

The project's contributions are workflow organization, file contracts, scientific review gates, failure recovery, structure and chemistry checks, provenance and the research interface. Peptide2Mol, REINVENT4, GNINA, RDKit, OpenMM and GROMACS are upstream tools, not original models or engines developed here.

Read in this order:

1. Root README: scope, workflow and limitations.
2. `00_workbench/webui/workflow.py`: states, upstream hashes, frozen settings, process ownership and retry archives.
3. `00_workbench/pipeline_bridge.py`: real subprocess adapters and output checks.
4. `00_workbench/webui/studio.py` and `studio_worker.py`: UI operations, diagnostics, receptor review and export.
5. `00_workbench/scientific/`: membrane construction, staged GROMACS and trajectory analysis.
6. `01_molecule_generation/src/`: structure preparation, core review, chemical identity and model provenance.
7. `02_docking/src/sbvs2/docking/` and `03_virtual_screening/src/zvss/`: docking records, pose/score binding and screening.
8. `04_interaction_analysis/src/zlis/` and `05_molecular_dynamics/modules/`: optional standalone tools.
9. Module tests and validation records: the precise scope of observed behavior.

Engineering tests verify software behavior. Model performance, candidate quality and biological conclusions require independent computational and experimental evidence. Mark real execution as accepted only when matching inputs, commands, environment, raw outputs and logs exist for the same version.
