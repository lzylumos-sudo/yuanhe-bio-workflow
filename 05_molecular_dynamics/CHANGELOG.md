# 2026-10-16: v2 input and task management

- Unified input, decisions, preparation, construction, OpenMM, analysis and acceptance through run_md.sh and scripts/run_workflow.py; legacy entry points retained.
- Input resolution honors explicit paths, limits discovery and blocks unattended ambiguity.
- Fixed-column atom inventory preserves stable identity, components, links, missing data and protonation variants, with bidirectional normalization maps.
- Ligand recovery validates templates, bounds RDKit inference time, preserves original charge sites/hydrogens/stereochemistry/coordinates and requires exact candidate approval.
- Backend preserves histidine evidence, requires terminal decisions and validates parameter/configuration/checkpoint contracts.
- Analysis mapping tracks original atoms/residues; rebuilt hydrogens are explicitly new.
- State management adds fingerprints, locks, decisions, stage logs, integrity and resume guards.
- Setup checks imports, tools and a minimal CPU/checkpoint exercise rather than Python-file existence alone. The explicit lock header was repaired while preserving 338 package URLs.
- Ubuntu validation script checks a short actual molecular run, interruption/resume and completed-run idempotence; native MD was not executed in this development round.

The R2 input regression preserves 3,094 atoms and a ligand with 92 atoms, 98 edges and +2 charge sites, with candidate confirmation required. HXT remains undecided rather than being silently renamed; four ring-hydrogen-supported HID residues are retained. See V2_VALIDATION.json.

Preparation, AM1-BCC parameterization, construction, full simulation, actual recovery and final mapping still require target Ubuntu acceptance. Old v1 passes do not become v2 passes, and Mac platform blocking is not a target-runtime scientific failure.
