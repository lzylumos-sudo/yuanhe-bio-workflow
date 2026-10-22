# Membrane Protein MD Protocol

Authority: `00_workbench/scientific/membrane.py`. OpenMM constructs the membrane system; HF/RESP and GAFF2 parameterize the ligand; ParmEd converts the topology; GROMACS runs production. These are different roles and must not be merged into the soluble route.

Defaults: 310 K, 1 bar semiisotropic pressure, 2 fs, one 500 ns trajectory (250000000 steps). The default builder lipid is **POPE**, not POPC. Alternative selected lipids remain explicit protocol deviations. Salt 0.15 M, padding 1.2 nm. Source force-field XMLs are `amber14-all.xml` and `amber14/tip3p.xml`, ligand `gaff-2.11`. The exact resolved lipid XML/version is unknown until installed files are inspected; AMBER naming alone does not verify Lipid17. New built systems record parameter file hashes and installed library versions.

The unchanged MDP uses steepest descent for minimization, V-rescale temperature control, 100 ps NVT heating from 50 K, two restrained NPT stages, C-rescale equilibration and Parrinello-Rahman production pressure control, semiisotropic coupling, PME, LINCS hydrogen-bond constraints and 1.0 nm cutoffs. The exact recipe is stored in the registry and compared as text for all phases/seeds/temperatures in regression tests.

Topology checks preserve ligand RESP charges, total neutrality, atom order and 1-4 exception scaling. Completion checks require the planned final step, final coordinates and production trajectory; partial output or engine exit zero alone is insufficient. The new protocol metadata adds an exact duration/timestep/step-count consistency gate. The original checkpoint continuity probe, TPR/input hashes and `-cpi -append` policy remain.

Analysis: protein-backbone least-squares alignment; first production frame reference; original PBC/whole-molecule treatment; ligand not independently fitted; 4.5 Angstrom contact occupancy is not hydrogen-bond occupancy. Original windows, PBSA membrane review and raw/derived output definitions remain unchanged. System-specific reviewed orientation, pose/state IDs, charge evidence and complete parameters are authoritative for each run.

Validation in this refactor is static/contract and exact MDP regression only. No real membrane building, topology conversion, checkpoint or 500 ns run was performed here. Existing method notes and dependency citations are in `docs/15_VINA_RESP_RESEARCH.md`; no new scientific claim is made.
