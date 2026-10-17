# Validation record

The earlier repository organization was checked on macOS / Apple Silicon on 2026-10-17: 302 Python files parsed; generation 57 passed; docking v2 29 passed, 2 skipped; workflow 15 passed; screening 24 passed, 3 skipped, 3 deselected; ZLIS 40 passed; standalone OpenMM contracts 76 passed. Total: 241 passed, 5 skipped, 3 deselected. These are historical checks of that revision, not current-release results.

The UI/membrane integration candidate was subsequently checked on 2026-10-18: workflow/UI 43 passed, legacy geometry 9 passed, generation 57 passed, screening 28 passed. Total: 137 passed; integration execution was excluded. Browser checks covered core selection, receptor/box preview, membrane preview and trajectory controls using clearly synthetic fixtures. Those fixtures are not target experiments.

Current English-release results are recorded separately in `../deployment/english_validation.json`; do not add overlapping historical totals together.

Not performed for this release: actual Peptide2Mol/REINVENT4 target inference, GNINA target batch docking, GROMACS membrane parameterization or dynamics, standalone OpenMM production, Windows launcher acceptance, or migration to a second computer. Adapter test doubles do not establish real inference. Dependency availability does not establish scientific validity.

The test interpreter uses an isolated local environment with RDKit, PyYAML and gemmi. Run suites separately because multiple modules use the same package names. Historical reports may reference artifacts excluded from GitHub; absence of those files is not evidence of a new successful run.
