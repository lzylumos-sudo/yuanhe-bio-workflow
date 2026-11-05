# Public same-compound UI acceptance

This is a minimum engineering test of the imported-candidate main workflow. It
does not establish scientific convergence, binding affinity, or experimental
activity of a newly generated molecule. Acetylcholine is a public known M2
agonist, rather than a generated hit. Peptide2Mol and REINVENT are deliberately
not invoked by this imported-candidate path; their independent checks remain
separate evidence.

## Input and scope

- Public receptor: [deposited 4MQS PDB](https://files.rcsb.org/download/4MQS.pdb).
  SHA-256: `ec46471ca53e1131d70f0977e528e94658fd48a22a96f2ade5f578820b763267`.
- Public molecule: [acetylcholine CSV](public_minimal/m2-acetylcholine.csv).
  Identity and pharmacology: [Guide to Pharmacology ligand 294](https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId=294&tab=analogues).
- Receptor chain A, deposited residue range 19–456. Chain B is the crystallization
  nanobody, not an orthosteric peptide reference. Select the generic target profile.
- Review the unresolved SER215–PRO377 segment boundary and the artificial termini.
  Do not connect those fragments by a peptide bond or silently rebuild the missing loop.
- Retain the current default screening rules. Quaternary nitrogen and triple-bond
  alerts are warnings; other gates and candidate qualification remain unchanged.
- Use explicit exploratory acknowledgements for a preliminary candidate and the
  geometric membrane proxy. These acknowledgements do not turn the candidate into
  a formally qualified scientific result.

## Original deployment, followed only by UI execution

Follow [DEPLOY.md](../DEPLOY.md) in a new source directory and runtime prefix.
Record the exact Git commit and installer receipt. Do not add packages or alter
installed code to repair a failed module. Fix GitHub source, pull and deploy it,
then repeat the affected UI workflow. Use a separate fresh project and environment
for the final reproduction. Verified public download caches may be reused; task
outputs, prepared receptors, charges, trajectories and analysis must not be copied
from the first run into the reproduction.

All business actions below use the actual browser: file upload, form editing,
review buttons, stage execution, trajectory controls and downloads. No direct
business API request or command-line module execution is part of this acceptance.

## Minimum settings

In a new research run, upload the PDB and CSV and verify both appear in the form.
Set `target_profile=generic`, `receptor_chain=A`, `peptide_chain=B`, and
`receptor_range=[19,456]`. Keep the uploaded paths supplied by the UI.

Set all candidate/sample limits to 1, including `peptide_attempts`, `core_limit`,
`reinvent_attempts`, `generated_unique_limit`, `external_limit`, `initial_limit`,
`refine_limit`, `review_limit`, `md_limit`, `wetlab_min` and `wetlab_max`.
Set `vina_workers`, `md_concurrency` and `generation_workers` to 1. Initial docking
uses one receptor, one state and one seed; refinement uses one receptor, one state
and three seeds. Keep the default docking exhaustiveness and screening thresholds.
Set `research.md_production_ns=0.01` before freezing the run.

Run preparation, risk and pre-docking screening. Prepare the receptor in the UI
at pH 7.4, allowing missing side-chain completion and reviewing the recorded
termini. Review and apply the following orthosteric docking box:

| Setting | Value |
| --- | --- |
| Center X/Y/Z, angstrom | -6.359, -9.777, -2.459 |
| Side lengths X/Y/Z, angstrom | 20, 20, 20 |
| Pocket residues | 103,104,107,108,111,155,194,195,400,403,404,426,429,430 |

Run docking and screening. Continue the same compound into membrane construction;
do not replace it with another molecule to make the downstream stages pass.

| Membrane/MD setting | Value |
| --- | --- |
| Lipid | POPC, idealized single-lipid bilayer |
| Salt / temperature / padding | 0.15 M / 310 K / 1.2 nm |
| NPT total equilibration | 0.1 ns |
| Production | 0.01 ns, one trajectory |
| Trajectory output interval | 1 ps |
| Velocity seed | 17 |
| HF/RESP threads | 12 |
| Extra receptor hydrogen addition | Disabled for the already prepared receptor |
| Alignment | protein and backbone |
| PB start/end/stride | 0 / 0.01 ns / 1 saved frame |
| PB slab | 30 angstrom, reviewed geometric proxy; dielectric 2, poretype 1 |
| PB solute/solvent dielectric | 1 / 80 |

Use the current run's screening report for its membrane center and normal.
Preview the orientation, record that it is a receptor-PCA geometric proxy, and
review the construct, pose and PB assumptions. Explicitly acknowledge exploratory
simulation. This is not an experimentally sourced membrane placement.

After construction, preview the actual system and review its compound mapping and
parameters. Choose **Use for GROMACS stage**, then run the GROMACS stage. It performs
minimization, NVT, two NPT stages, the independent checkpoint-continuity probe,
production and analysis. Inspect the trajectory slider through the last frame,
RMSD/RMSF and thermodynamic plots, and the real AmberTools PB result. Download
`systems.json` and the reproduction archive through the UI.

## Acceptance evidence

Both runs passed on 2026-11-05. The final scientific source was
`7e0c66cfc8a847e0d678d07267039498d538516f`.

| Evidence | First acceptance | Independent reproduction |
| --- | --- | --- |
| Project | YH-20261105-a38b083b2b6c4dddb61c2b310c7eb7ac | YH-20261105-f15ce681f5544a488fa45f97f605e076 |
| System | MD_cccf49d60e7c | MD_9608e057161b |
| Compound | CMP84D51F4020B502D4 | CMP84D51F4020B502D4 |
| Production / saved frames | 0.01 ns / 11 | 0.01 ns / 11 |
| Finite trajectory / constraint warnings | yes / none | yes / none |
| Checkpoint continuity | passed | passed |
| Real PB analysis frames | 11 | 11 |
| UI last-frame slider and archive download | passed | passed |

The first run began at 05c3826 and completed after audited UI MD-only continuations
for source repairs. The second ran from new-project creation through completion
entirely at 7e0c66c in a separate original-flow installation. It independently
prepared the receptor, fitted HF/RESP charges, built the membrane and computed
the trajectory and PB results. Public download caches were reused, task outputs
were not. Reproduction means the workflow completes again, not bitwise identical
molecular coordinates.

[Machine-readable evidence and artifact SHA-256 values](public_minimal/mainflow-ui-acceptance.json)
record both completed UI downloads. The reproduction ZIPs retain inputs, settings,
reviews, logs and analyses; raw-binary inclusion was left unchecked. The full local
run evidence is retained separately. Later administrator-entrypoint and UI wording
changes do not change this scientific source.

The formal 400–500 ns core report must still say `insufficient_trajectory` for this
short test. It is not a software failure and must never be relabeled as a converged
scientific result. The target administrator must separately verify their machine's
deployment identity and execute the same small UI test before handing over its URL.
