# Standalone OpenMM molecular dynamics

Current research default: **one 500 ns production trajectory per selected complex**,
GAFF2 2.11 + HF/6-31G* ESP / two-stage RESP, ff14SB, TIP3P and compatible ions.
Minimization, heating and equilibration are additional. GPU use is required in the
research configuration and the actual device/precision is recorded. No charge
fallback to AM1-BCC and no automatic CPU fallback are allowed.

```bash
bash run_md.sh --input /path/to/handoff/complex.pdb --config config/default.yaml \
  --run-id complex_A_single500 --non-interactive
bash run_md.sh --resume runs/complex_A_single500 --non-interactive
```

Use trusted embedded/adjacent ligand SDF chemistry and the selected pose mapping.
A QM optimized copy never replaces the protein-bound pose. Missing scientific
input decisions block the run and produce a pending decision file.

Core output: protein-aligned ligand heavy-atom RMSD mean and sample SD (ddof=1)
over **400–500 ns**, with raw series, block statistics, full RMSD/RMSF, distances,
contacts and hydrogen bonds. Short trajectories cannot pass this duration gate.
MM/PBSA uses actual AmberTools PB calculations on the same trajectory and records
energy components, uncertainty and omitted entropy. MM/GBSA remains an explicit
optional historical method; it is not renamed PBSA.

`config/default.yaml` separates 500 ns production, 1 ns transaction segments,
10 ps output and 100 ps inflight checkpoints. Resume restores a committed
checkpoint without assigning velocities or resetting production time. An
interrupted segment is replayed from its preceding committed checkpoint.
`config/smoke.yaml` is engineering-only. One trajectory provides no independent-
repeat uncertainty; no low-SD or elapsed-time rule establishes binding/convergence.

The workbench has a separate explicit-membrane GROMACS route; ordinary aqueous
OpenMM is not a MurJ membrane simulation. Review orientation, lipid model and
membrane-specific PB settings before using that route.

See [full methods, migration and commands](../docs/15_VINA_RESP_RESEARCH.md),
[license obligations](../docs/16_RESEARCH_DEPENDENCY_LICENSES.md), and
[archived instructions](docs/ARCHIVED_README_PRE_RESP.md). This revision was not
run or tested after the owner's code-only instruction; no acceptance results
from this task are retained. Historical validation files do not certify it.
