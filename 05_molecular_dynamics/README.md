# Nexus MD v2

Use `bash run_md.sh` for lossless input organization, chemical recovery, explicit scientific decisions and task management with the OpenMM backend. Target: x86-64 Ubuntu / WSL2 Ubuntu. Mac engineering tests do not validate Linux dynamics.

Read [v2 guide](V2_USER_GUIDE.md), [changelog](CHANGELOG.md) and V2_VALIDATION.json. No ChatGPT login is required to run this standalone module.

```bash
bash setup.sh
bash run_md.sh --input "/path/to/my_complex R2.pdb" --inspect-only
source runtime.sh
python scripts/validate_v2_ubuntu.py
```

New tasks use runs/ and default to one 100 ps production segment. The legacy run.sh below retains its older long-run defaults; do not confuse the entry points. Historical v1 acceptance does not certify v2.

---

# Universal MD Simulation — v1 archived instructions

Independent OpenMM simulation from a supplied `complex.pdb`. No docking project or Shared_Environment is imported. Default: 300 K, 1 atm, 2 fs, three independent 100 ns production replicas; 100 ps NVT and 100 ps NPT per replica. This installation defaults to CPU, four threads.

## Install and run

```bash
bash setup.sh
bash run.sh examples/1fkj/complex.pdb --config config/smoke.yaml --run-id public_1fkj_smoke
bash run.sh examples/1fkj/complex.pdb --config config/demo_10ns.yaml --run-id public_1fkj_10ns
# Default production protocol, 3 x 100 ns:
bash run.sh inputs/complex.pdb
# Resume a preserved run with exactly its stored settings and input:
bash run.sh --run-id public_1fkj_10ns --resume
```

Use Linux/WSL. Paths are derived from `PROJECT_ROOT`; the runtime defaults to `${XDG_CACHE_HOME:-$HOME/.cache}/nexus-structural-biology-v1/environments/md`. `NEXUS_ENV_PREFIX` overrides it. `setup.sh` restores the bundled conda-pack archive when available, then runs conda-unpack and a real CPU test. Otherwise it rebuilds from the explicit lock. Copy this entire module and rerun setup; copying an extracted conda prefix without relocation is unsupported. Heavy environments remain on the Linux filesystem. No WSL installation or GPU driver replacement is performed.

## Input chemistry

The command takes only `complex.pdb`. Nexus exporters embed the verified ligand SDF and manifest in hashed, bounded `REMARK 950` records; copying that single PDB preserves the chemistry. The MD module reads these records locally and validates the chemical graph and atom mapping. `examples/1fkj_single_file/complex.pdb` demonstrates this contract without sidecars. Ordinary external PDB files need a precise adjacent `ligand.sdf`, which is read automatically. If embedded and adjacent definitions both exist but disagree, the run is rejected. Bare PDB coordinates lack reliable bond orders, protonation and stereochemistry. A single nonstandard residue is selected automatically. An optional adjacent `manifest.json` specifies exact `ligand` resname/chain/resid/template_sdf and PDB-atom-name to 1-based SDF `atom_mapping`; `analysis.groups` assigns chains to biological roles. All preparation remains local.

Examples include public 1FKJ (FKBP12–CCD FK5), 1BRS (protein pair), and 1FAP (public ternary complex). Their fixture provenance and chemistry mappings are copied, while old trajectories and PASS files are excluded. For 1FAP, role names do not imply that this example is an E3 degradation system.

The current Amber export and analysis contract uses `protein.remove_water: true`; retaining input crystal waters is not supported by this version. Unresolved alternate conformers, unknown cofactors/metals, cross-residue covalent ligands, modified residues and ambiguous chemical mappings produce explicit errors. Missing residues require an explicit modeling decision. Protein hydrogens are rebuilt by default (`protein.rebuild_hydrogens: true`) to remove malformed or inconsistent docking-converter H records while preserving all input heavy atoms; the removed count is recorded. Set this false only for explicitly prepared, valid protonation inputs. PDBFixer is explicitly pinned to CPU. Standard protein protonation uses pH 7; it is not a site-specific pKa prediction. Private inputs are never sent to external services.

## Outputs and protocol

`results/<run-id>/inputs`, `md`, `analysis/replica_XX`, `config.yaml`, `manifest.json` and `acceptance.json` retain inputs, actual settings, topology, independent seeds, staged NVT/NPT, production DCD/CSV segments, SHA-256, binary checkpoints and portable XML states. Resume preserves velocities and never silently starts a new trajectory. Binary checkpoints may require the same OpenMM/hardware; XML recovery requires a separate explicit migration audit. Uncommitted segments can be regenerated from the last committed checkpoint.

Protein Amber ff14SB, TIP3P water/ions, GAFF 2.11 small molecules with AM1-BCC, PME, 0.9 nm cutoff, switching at 0.8 nm, 1.2 nm solvent padding, neutralization and 0.15 M NaCl. Minimization precedes velocity initialization. Heavy atoms are restrained during NVT/NPT and released in production; production time is reset to zero. Change parameters in a new module-relative YAML configuration, not in an active run.

Analysis outputs: `RMSD.png`, `RMSF.png`, `Ligand_RMSD.png`, `Hbond.csv`, `Contact_map.png`, `SASA.csv`, `MMGBSA.csv`, raw numeric tables and method metadata. Protein Cα least-squares alignment uses the first analyzed production frame. Ligand RMSD uses heavy atoms after protein alignment, without ligand self-fit or symmetry correction. Hydrogen bonds use donor–acceptor distance ≤3 Å and D–H–A angle ≥150°. Contacts use 4.5 Å. SASA uses Shrake–Rupley, 1.4 Å probe and 100 sphere points; pairwise buried area is not additive in ternary systems. PPI ligand RMSD is explicitly not applicable.

MM/GBSA runs the actual AmberTools `MMPBSA.py` against complete unconstrained Amber parameter exports, using igb=5, mbondi2 radii and 0.15 M salt. It reports a sampled single-trajectory end-state energy estimate in kcal/mol, **without entropy**. Default partition is all protein chains versus the small molecule; for a ternary complex this is not a cooperativity calculation. No small molecule yields an explicit not-applicable result. `mmgbsa.max_frames` controls evenly sampled frames; the exact times and backend output are retained. It must never be interpreted as an experimental binding free energy.

## Acceptance

`INSTALLATION_REPORT.md` records real installation and the native CPU test. `SHORT_RUN_TEST.md` is written for a completed short molecular test. `MD_TEST_PASS.md` is created only after a new 10 ns run passes actual production-step/time, final-checkpoint, finite-observable, trajectory-frame and analysis-file audits. Neither a short run nor 10 ns proves convergence, a correct pose, affinity, recruitment or degradation.

Production MD is explicitly CPU-limited in this build, and PDBFixer preparation is now explicitly pinned to CPU. The first v1 smoke test did not record its automatic preparation platform; that historical uncertainty is preserved in INSTALLATION_REPORT.md. CUDA availability is recorded separately from GPU execution. Hardware selection must be configured explicitly, including resource coordination with other local simulations.

Sources: [OpenMM documentation](https://docs.openmm.org/latest/userguide/), [OpenMMForceFields](https://github.com/openmm/openmmforcefields), [ParmEd OpenMM export](https://parmed.github.io/ParmEd/html/openmm.html), [Amber MM/PB(GB)SA tutorial](https://ambermd.org/tutorials/advanced/tutorial3/).

## Explicit installation and export checks

Use `python scripts/validate_native.py` in the installed OpenMM environment for the CPU smoke check. To compare an existing run's original system with its exported Amber topology, use `python scripts/validate_amber_export.py RUN_ID` in the complete OpenFF/Amber environment. These tools require real dependencies or run artifacts and are not executed during ordinary pytest collection.

## Independent completed 10 ns audit

After normal analysis, audit the public CPU example independently:

```bash
source runtime.sh
"$NEXUS_ENV_PREFIX/bin/python" scripts/validate_10ns.py --run-id public_1fkj_10ns_v2
```

The auditor only reads a completed run and reloads its binary checkpoint in a separate CPU Context without advancing simulation. It checks at least 5,000,000 steps at 2 fs, input/configuration/system/segment/analysis hashes, finite coordinates, global CSV step/time and actual DCD frame counts/times. Segment-local DCD times are offset using independently verified preceding steps before comparison to the global CSV. It does not automatically use a GPU.

Only full success writes installation_evidence/10ns_final_audit.json, validation_status.json and MD_TEST_PASS.md. Incomplete or 4 ps runs exit nonzero without overwriting prior successful evidence; the historical rejection record is installation_evidence/10ns_validator_rejection_test.json. Numerical completion does not establish convergence or affinity.


## Yuanhe Bio integration

See the [project code map](../docs/product/CODE_AND_CONTENT_MAP.md) for current module boundaries and the [product value assessment](../docs/product/PRODUCT_VALUE_ZH.md) for intended users. Existing technical instructions and historical acceptance records above retain their original scope; this documentation refresh does not revalidate execution.
