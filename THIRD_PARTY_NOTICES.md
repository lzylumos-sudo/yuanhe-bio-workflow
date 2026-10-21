# Third-party sources and attribution

This repository contains workflow orchestration, adapters and analysis code. Third-party model source trees, weights and installed environments are restored separately.

| Component | Source | Integration |
|---|---|---|
| Peptide2Mol | https://github.com/BLUE-Flowing/Peptide2Mol | Separate model process; upstream code uses MIT |
| REINVENT4 / LibInvent | https://github.com/MolecularAI/REINVENT4 | Separate sampling process; Apache-2.0 code |
| GNINA | https://github.com/gnina/gnina | External docking executable |
| RDKit | https://www.rdkit.org | Chemical identity, conformers and structure processing |
| OpenMM | https://openmm.org | Membrane construction and optional standalone dynamics |
| GROMACS | https://www.gromacs.org | Minimization, equilibration and production |
| ProLIF / PyMOL / PoseBusters | Module environment manifests | Interaction analysis, visualization and quality checks |

Pinned commits and weight hashes are in `deployment/model_dependencies.json`. Preserve upstream licenses and authorship when restoring dependencies. Source links identify provenance; they are not a complete license audit. See `LICENSE_SCOPE.txt` and `THIRD_PARTY_NOTICES.txt` for exclusions from the project license.

Full-text papers are not redistributed. Original PDB coordinates and provenance are in `data/structures/manifest.json`. These structures were not determined by this team and do not establish experimental activity of generated compounds.

Current production docking uses official AutoDock Vina 1.2.7. New preparation/MD/analysis dependencies include Meeko, RDKit, Psi4, cdsgroup RESP, AmberTools, OpenMM, ParmEd, MDAnalysis and GROMACS. Their component-specific license and distribution obligations remain; see `docs/16_RESEARCH_DEPENDENCY_LICENSES.md` at the repository root. Historical GNINA notices remain applicable to historical code and artifacts.

## Pre-docking risk and acquisition filters (2026-10-21)

The curated risk configuration uses individually selected RDKit Brenk patterns,
an explicitly adapted NIH epoxide alert, and PAINS catechol_A(92). RDKit attribution
and its BSD license are retained in `docs/preflight_sources/rdkit_license.txt`.
The nitrosamine SMARTS is project-defined and informed by FDA's structural-risk
discussion; it is not an FDA classifier or a measured toxicity result.

AiZynthFinder 4.4.1 is an optional MIT-licensed local dependency. filelock is an
existing dependency; current upstream is MIT (check each installed version).
Software notices and source hashes are in `docs/preflight_sources/`. No model
weights, reaction-template datasets, raw-material inventories or vendor catalogs
are redistributed here. Their licenses remain separate and must be recorded for
each installation/import; see `deployment/preflight_licenses.json`.
