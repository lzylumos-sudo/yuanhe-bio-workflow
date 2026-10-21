# Third-party sources and licenses

The complete distribution is not uniformly Apache-2.0. GNINA includes its documented GNU/Apache portions and dependencies; preserve and inspect licenses/ and sources.json. NVIDIA user-space libraries retain independent terms; drivers do not migrate with this project.

Package versions/builds/sources/licenses are recorded in env/package_provenance.json, pip_sources.json and vendor.lock.json when the runtime is restored. Ubuntu Base originates from https://cdimage.ubuntu.com/ubuntu-base/releases/24.04/release/. RCSB PDB/CCD provenance is recorded in data/sources.json. Localized 3Dmol.js 2.5.5 assets retain their upstream license and source hashes. Review each component's actual terms before redistribution.

Current production docking uses official AutoDock Vina 1.2.7. New preparation/MD/analysis dependencies include Meeko, RDKit, Psi4, cdsgroup RESP, AmberTools, OpenMM, ParmEd, MDAnalysis and GROMACS. Their component-specific license and distribution obligations remain; see `docs/16_RESEARCH_DEPENDENCY_LICENSES.md` at the repository root. Historical GNINA notices remain applicable to historical code and artifacts.
