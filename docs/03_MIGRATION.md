# Naming and migration

The repository uses numbered English module directories. Package names, command-line arguments, scientific artifact identifiers and upstream model names retain their interfaces. Runtime references, launchers, structure selectors and model configuration are updated together with paths.

| Area | Repository path |
|---|---|
| Unified interface | `00_workbench` |
| Generation v2 | `01_molecule_generation` |
| Main GNINA backend | `02_docking` |
| Virtual screening | `03_virtual_screening` |
| ZLIS | `04_interaction_analysis` |
| Standalone OpenMM | `05_molecular_dynamics` |
| Earlier docking prototype | `archive/docking_prototype` |
| Structural inputs | `data/structures`, lowercase filenames |
| Launchers | `start_linux.sh`, `start_windows.cmd`, `start_windows.ps1` |

`FILE_MAP.csv` retains original file SHA-256 values and repository destinations. Source path labels are normalized to English; they are descriptive labels, not byte-exact historic filenames. Edited files naturally have different hashes. Old development prompts, temporary diagnostics, full-text papers, duplicate v1 generation copies and historical run outputs are excluded.

The historical `pepwedge` package supports legacy analysis pages; it is not the production model backend. This release is prepared in a separate checkout. Existing external-disk run records are not rewritten to suggest they can resume under changed source hashes. Resume old runs with their original code, or create new runs using this release.
