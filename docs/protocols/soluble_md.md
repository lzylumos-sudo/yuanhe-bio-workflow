# Soluble Protein MD Protocol

Authority: `05_molecular_dynamics/config/default.yaml`; full content is frozen in the registry. Backend: standalone OpenMM, distinct from membrane GROMACS.

Defaults: 300 K; 1 atm; 2 fs; one 500 ns production trajectory (250000000 steps). Protein `amber14-all.xml`, water `amber14/tip3p.xml`, ligand `gaff-2.11`; HF/6-31G* two-stage RESP with charge preservation. Do not reduce this to an unspecified AMBER label. Actual installed package and parameter file versions must be read from the run environment/system metadata.

The unchanged YAML also specifies pH 7, 1.2 nm solvent padding, 0.15 M ions, 0.9 nm nonbonded cutoff, 0.8 nm switch distance, 1/ps friction, constrained integration, 100 ps heating from 50 K, 1 ns NVT and 1 ns NPT. Production uses the original OpenMM implementation; this registration does not modify integrators, thermostat, barostat, 1-4 scaling or charge treatment.

Analysis: protein-backbone alignment; formal 400-500 ns window; 10 ns blocks; contacts at 4.5 Angstrom; hydrogen-bond 3.0 Angstrom / 150 degrees; SASA settings, PBC handling, stride and PBSA settings stay in the original analysis code/config. Static ZLIS measurements are not interchangeable with trajectory occupancies.

Inputs: original complex and selected components, chemical-state/atom mapping, coordinates, topology, parameters and explicit decisions. Outputs: systems, trajectory/checkpoints, analysis, acceptance and versioned run metadata. Checkpoint/state/atom order and charge checks remain unchanged. Coordinate-only restart is not represented as continuous checkpoint recovery.

Scope: regression/contract tests on this Mac; no new production trajectory, RESP calculation or target CUDA acceptance. See existing `05_molecular_dynamics/docs` and `docs/15_VINA_RESP_RESEARCH.md` for methods and scientific limitations. Parameter choices are engineering protocol definitions, not a claim of a unique literature standard.
