# Yuanhe Bio Workflow

A traceable research workbench connecting protein-peptide complexes to Peptide2Mol seeds, reviewed molecular cores, REINVENT4 expansion, GNINA docking, virtual screening and membrane molecular dynamics.

**Release 0.3.0-rc2: single-case Ubuntu / RTX 5090 D integration tested.** The real UI now connects screened poses to membrane construction, equilibration, 1 ns production MD, checkpoint testing, trajectory analysis and full export. Upstream results were reused with hash verification; Peptide2Mol remains a verified partial-output run. Windows/WSL2 with RTX 5070 remains untested.

**Results:** six molecules / 54 docking poses upstream, three preliminary display candidates, **zero formal candidates**; one exploratory molecule completed MD, with **zero scientifically validated molecules**. See [upstream acceptance](docs/06_RTX5090_ACCEPTANCE.md), [actual membrane MD acceptance and limitations](docs/08_MEMBRANE_MD_ACCEPTANCE.md), and [reusable docking, screening and MD modules](docs/07_MODULE_CONTRACTS.md).

```mermaid
flowchart LR
 A[Complex and construct review] --> B[Peptide2Mol seeds]
 B --> C[Interactive core review]
 C --> D[REINVENT4 expansion]
 D --> E[GNINA batch docking]
 E --> F[Virtual screening]
 F --> G[Candidate and membrane review]
 G --> H[Build membrane system]
 H --> I[GROMACS and checkpoint recovery]
 I --> J[Trajectory analysis and export]
```

| Directory | Purpose |
|---|---|
| `00_workbench` | Unified UI, guided review, scheduling, membrane builder, GROMACS and analysis |
| `01_molecule_generation` | Preparation, model adapters, core and chemistry checks |
| `02_docking` | GNINA protocols, docking records and standalone tools |
| `03_virtual_screening` | Pose QC, evidence, diversity, qualification and ranking |
| `04_interaction_analysis` | Optional ZLIS interaction visualization |
| `05_molecular_dynamics` | Optional standalone OpenMM workflow |
| `data/structures` | Original structural inputs and provenance |
| `deployment` | Installation, runtime configuration and validation tools |
| `docs` | Architecture, setup and validation records |
| `archive` | Historical docking prototype |

## Start

After initial dependency installation, double-click `start_windows.cmd` and select the configured WSL distribution. On Linux use `bash start_linux.sh`. For a UI-only preview with Python 3.11+, run:

```bash
python3 00_workbench/run_webui.py --no-browser
```

Open the displayed local address, normally http://127.0.0.1:8777. The UI starts without launching calculations. Configure environments once, then prepare, review, execute, stop, resume, inspect trajectories and export records through the interface.

- [Setup](docs/01_SETUP.md)
- [UI workflow](USER_GUIDE.txt)
- [RTX 5070 installation and MurJ validation](deployment/validation_5070/INSTALL_AND_VALIDATE.txt)
- [Architecture](docs/02_PROJECT_GUIDE.md)
- [Validation scope](docs/04_VALIDATION.md)
- [Third-party sources](THIRD_PARTY_NOTICES.md)

## Scientific boundaries

The current input is MurJ / 9NU4. Review chains, construct mutations, residue ranges, protonation and the docking box. Peptide2Mol uses the peptide-defined receptor pocket; explicit peptide-pharmacophore retention is not claimed. REINVENT4 expands a reviewed molecular core, not a protein backbone.

Counts distinguish requested samples, valid unique molecules, docking poses, preliminary candidates and formally qualified candidates. Unperformed checks never become passes. The screening membrane orientation is still a geometric proxy; the separate membrane builder requires a sourced and reviewed orientation. Idealized single-lipid membranes are not full bacterial membrane models.

The built-in MD path includes parameterization, minimization, staged equilibration, production, checkpoint guards and analysis. It has not been executed on the target Windows machine in this release. Short trajectories and docking scores do not establish affinity, activity or convergence.

Weights, installed environments, VHDX images and historical run outputs are excluded from GitHub. Restore pinned models separately. An external disk does not remove host driver, WSL, environment or target-machine acceptance requirements.

## License

Original project contributions in this release use the unmodified [PolyForm Noncommercial License 1.0.0](LICENSE), subject to [scope and existing rights](LICENSE_SCOPE.txt). This is source-available software, not an unrestricted commercial-use open-source license. Commercial uses outside the license's permitted purposes require separate authorization. Third-party materials and prior grants retain their own terms; see [notices](THIRD_PARTY_NOTICES.txt).
