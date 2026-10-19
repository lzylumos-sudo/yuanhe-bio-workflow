# Platform and validation matrix

The table separates recorded execution from portable engineering checks. A supported input format is not evidence that every target or platform is scientifically validated.

| Capability | Environment | Evidence / status |
|---|---|---|
| Workbench UI preview | Python 3.11+; local browser; POSIX host | Standard-library server. Preview does not install engines or run calculations. Windows uses WSL. |
| Recorded generation, GNINA and screening | Linux x86-64 / RTX 5090 D | [Integrated record](10_INTEGRATED_ACCEPTANCE.md); exact MurJ input and recorded conditions only. |
| Membrane/GROMACS workflow | Recorded Linux environment | Earlier bounded run completed; final three-repeat acceptance and export remain open. |
| RTX 5070 / Windows WSL2 | Planned target | Installation recipe exists; target execution is unverified. |
| Generation, docking, screening and MD contracts | CI Linux / local macOS checks | Unit/contract and public-fixture checks; no replacement for engine or GPU validation. |
| ZLIS interaction analysis | Recorded macOS environment | Optional module with separate environment. Windows remains unverified. Its UI environment uses Python 3.11; do not infer support for later Python versions from workbench preview. |
| Standalone OpenMM | Ubuntu / WSL target; macOS contracts | Historical module evidence is separate from current integrated GROMACS acceptance. |

## Run portable checks

From the repository root, using the chemistry test environment:

```bash
python deployment/check_english.py
python deployment/check_repository.py
(cd 01_molecule_generation && PYTHONPATH=. python -m pytest tests -q)
(cd 00_workbench && PYTHONPATH=. python -m pytest tests -q)
(cd 02_docking && PYTHONPATH=src python -m pytest tests -q)
(cd 03_virtual_screening && PYTHONPATH=src python -m pytest tests -m 'not integration' -q)
(cd 04_interaction_analysis && PYTHONPATH=src python -m pytest tests -q)
(cd 05_molecular_dynamics && python -m pytest tests -q)
python -m pytest deployment/tests -q
node --experimental-vm-modules 00_workbench/webui/check_js.mjs
```

Run modules separately: several historical modules use common Python package names. Install the recorded chemistry dependencies, OpenMM/ParmEd and the optional ZLIS test dependencies before running the combined set. See the CI workflow for the precise installation commands. A skipped dependency-dependent test must not be counted as passed. Workbench HTTP regressions require permission to bind a loopback socket.

## Explicit target acceptance tools

These commands are intentionally outside automatic test collection:

- `02_docking/scripts/validate_stop_guard.py`: dedicated legacy `/opt/sbvs` runtime and `SBVS_PROJECT_ROOT`; verifies that stop is refused while a computation lock is held.
- `05_molecular_dynamics/scripts/validate_native.py`: real OpenMM CPU installation smoke, not production MD.
- `05_molecular_dynamics/scripts/validate_amber_export.py RUN_ID`: existing complete run under the module's `results/`, OpenFF/Amber parameterization environment and exported topology required.
- `NEXUS_NATIVE_TESTS=1`: enables the separately documented Ubuntu preparation test, without changing its scientific acceptance rules.

Do not launch these against an unknown active deployment or substitute a simulation result for a missing validation record.
