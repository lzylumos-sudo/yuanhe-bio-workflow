# Nexus MD v2 user guide

## Install and validate

On x86-64 Ubuntu or WSL2 Ubuntu, enter this module's actual mount location and run `bash setup.sh`. Environments stay on the Linux filesystem; source/configuration/results can remain on the external disk. GitHub excludes the historical environment archive. If a separately supplied archive exists, setup verifies it; otherwise it rebuilds from pinned URLs. The older explicit list lacks per-package hashes, so preserve the actual solved environment evidence.

Use `bash setup.sh --repair` to archive a broken environment before rebuilding. Do not move an extracted Conda prefix directly. NEXUS_ENV_PREFIX selects a Linux location. Setup verifies imports, tools and a small CPU/checkpoint operation.

After installation, `source runtime.sh` and run `python scripts/validate_v2_ubuntu.py`. This uses public 1FKJ with 1 ps NVT, 1 ps NPT and 10 ps production, actual interruption, binary checkpoint reload/resume, analysis and completed-run idempotence. Default wall-time limit is 1,200 seconds; --max-seconds changes the wall limit, not the simulation duration. An interruption not achieved is NOT_RUN, not passed. This acceptance has not been executed by the Mac-only engineering checks.

## Commands

```bash
bash run_md.sh
bash run_md.sh --input "/path/to/my_complex (R2).pdb"
bash run_md.sh --input "/path/to/Complex.pdb" --config config/v2_short.yaml
bash run_md.sh --input "/path/to/Complex.pdb" --inspect-only
bash run_md.sh --input "/path/to/Complex.pdb" --prepare-only
bash run_md.sh --input "/path/to/Complex.pdb" --smoke-test
bash run_md.sh --input "/path/to/Complex.pdb" --config config/v2_short.yaml --non-interactive
bash run_md.sh --resume "runs/TASK_ID" --decisions "/path/to/decisions.yaml"
```

The default v2 run is one 100 ps production, 100 ps NVT and 100 ps NPT on two CPU threads. Smoke mode uses one replica, two threads, 1+1 ps equilibration and 10 ps production. Set actual duration/resources in YAML. Legacy run.sh retains different defaults.

Explicit input takes precedence. Discovery is restricted to the selected input directory, never old run outputs. Each task has a unique directory; existing IDs are not overwritten. Unattended missing information writes pending.json and exits instead of waiting for a keyboard.

## Scientific decisions

Supply a complete complex PDB and intended settings. Embedded or adjacent SDF/manifest chemistry is preferred and checked for elements, bonds, atom mapping, charges and stereochemistry. With PDB alone, local inference considers CONECT, explicit hydrogens, charges and coordinates; a candidate is not accepted without review.

Confirm the actual construct and internal gaps; sequence numbering alone does not establish completeness. Review ligand chemistry and approve its exact hash or provide precise SMILES. Choose models, alternate conformations, retained components, special termini and whether input waters may be replaced. Protein pH is an assumption, not site-specific pKa prediction.

```yaml
construct_confirmed: true
termini:
  r0188: standard_charged
ligands:
  LIGAND_ID:
    accept_candidate_sha256: "REVIEWED_FULL_HASH"
```

This illustrates syntax only; IDs and hashes must come from the actual task. Check the generated decision schema before use. Preserve existing protein heavy-atom coordinates and HID/HIE/HIP or ring-hydrogen evidence. Simple side-chain completion is recorded. HXT is not silently renamed OXT; standard-terminus selection explicitly records removal/addition. Large backbone gaps and abnormal links are not invented.

## Supported scope

The enabled preset is AMBER ff14SB + GAFF 2.11/AM1-BCC + amber14 TIP3P with matching Na/Cl parameters. ff19SB/OPC, CHARMM36m/CGenFF, membranes and nonaqueous solvents are not validated in this standalone module. The unified workbench membrane path is separate.

Ordinary noncovalent protein/single-ligand and protein-only inputs are recognized; multichain roles can be supplied for analysis. Multiple ligands are preserved during inventory but blocked by the single-ligand simulation adapter. Metals, covalent ligands, cofactors, complex glycosylation and missing backbone require another supported route rather than silent simplification.

## Outputs and recovery

runs/TASK_ID stores original inputs/hashes, inventories, normalization attempts, chemical candidates, approvals, configuration, decisions, logs, topology, trajectories/checkpoints, mappings and acceptance. READY means input checks passed; PREPARED means system construction; COMPLETE requires actual MD, analysis and trajectory/checkpoint audits. Other states explicitly distinguish waiting, dependency blocks, failure and interruption.

Exit codes: 0 completes the requested mode; 2 failure; 3 decision required; 4 dependency/platform block; 5 active lock; 130 user interruption. Resume validates inputs, settings, decisions, code, dependencies and artifacts. Changed protonation/force field/structure requires a new task. Binary checkpoint recovery preserves velocities/state; incompatible checkpoints never silently become fresh-velocity restarts.

After a crash inspect .running/owner.json and confirm the corresponding host process has stopped before removing only that task lock. Preserve checkpoints and progress. The 2 GiB free-space guard is not a long-run storage estimate.

Analysis uses protein C-alpha alignment, ligand heavy-atom RMSD after protein fitting, C-alpha RMSF, contacts/hydrogen bonds/distances/SASA and original-identity maps. It handles periodicity and molecular integrity. C-alpha RMSD is not all-backbone RMSD. MM/GBSA excludes entropy; all protein chains versus one ligand in a ternary complex is not cooperativity. Short tests do not establish equilibrium, affinity or degradation activity.

The R2 regression fixture tests software contracts; simulated approval is not scientific review. If a deployment backup exists, its rollback tool checks current hashes before restoring and refuses to overwrite later changes. Stop tasks before rollback; preserve runs and trajectories.
