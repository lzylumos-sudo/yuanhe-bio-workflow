# One pipeline, reusable modules

The workbench is the orchestrator. Scientific cores retain explicit input/output contracts, provenance and failure states. A module invocation and a UI stage must call the same engine; neither path may silently change candidate qualification.

| Stage | Independent component | Handoff |
|---|---|---|
| Generation | `01_molecule_generation/v2` adapters | Stable molecule IDs, isomeric SMILES, raw SDF, generation history |
| Docking | Installable `02_docking` / `sbvs2` | Receptor/box/protocol, molecule/state IDs, poses, scores, engine logs |
| Screening | Installable `03_virtual_screening` / `zvss` | Candidate compound/state/pose IDs, selected SDF, QC and qualification evidence |
| Peptide evidence | `scientific.peptide_contacts` | Receptor-frame check, peptide and candidate contacts, geometry controls, same-pool ranking ablation |
| Membrane MD | `00_workbench/scientific` / `scientific.md_cli` | Versioned, checksum-verified MD bundle, built system, simulation phases, analysis |

Docking and screening already have their own Python distributions. Export them with `python deployment/package_component.py --module docking --output /new/path/docking.zip` or `--module screening` into a separate output file; the archives include source, environment guidance and checksum manifests. The membrane core is exportable with `python deployment/package_md.py --output /new/path/yuanhe-md.zip`. This archive contains the same implementation called by the UI, dependencies and license notices. It has no UI dependency at import time. GPU engines and model weights remain external dependencies; packaging is not a claim that every host is validated.

## Screening to MD: `yuanhe.md-input/1`

The membrane wizard exports `md_handoffs/H_*/handoff.json` before construction. Each bundle contains copied `inputs/receptor.pdb` and `inputs/ligand.sdf`, SHA-256 digests, candidate identities, isomeric SMILES, unchanged qualification, explicit exploratory acknowledgement where needed, orientation and physical settings, execution configuration, originating run and stage hashes. The SDF chemistry must match the screening identity before parameterization. The pose coordinates are retained and both protein and ligand undergo the same rigid membrane transform.

External workflows may call `scientific.md_cli.export_bundle()` with the same explicit fields; no opaque global run directory is required. Absolute executable paths in `execution` are host-specific configuration and must be reviewed on another host. Preserve the original bundle and record any host migration in a new copy.

From the exported module's root, using the membrane environment:

```bash
python -m scientific.md_cli validate --bundle /path/to/H_bundle
python -m scientific.md_cli build --bundle /path/to/H_bundle
python -m scientific.md_cli run --bundle /path/to/H_bundle --system-id MD_ID --review-system
python -m scientific.md_cli analyze --bundle /path/to/H_bundle --system-id MD_ID
```

These commands are the independent-module interface, not substitutes for the recorded browser acceptance. The UI uses the same `build_prepared`, `run` and `analyze` functions. The bundle has an exclusive compute lock. Input mutation, changed simulation protocol, changed engine binary, failed minimization or incomplete production block continuation or analysis.

## State and scientific gates

Construction produces a system ID, force-field and charge records, atom groups, topology/coordinate hashes and software versions. ParmEd receives unconstrained bond and angle parameters; GROMACS applies hydrogen-bond constraints and water SETTLE. This preserves the water geometry during conversion rather than deleting the constraints from production. Nominally neutral lipids may contain tiny tabulation residuals: only residuals up to 1e-5 e per lipid are corrected on the phosphate atom, with every correction recorded. Coulomb exception products are updated simultaneously to preserve their original 1-4 scaling. Larger residuals are rejected. The complete exported system must then be neutral within 1e-3 e.

The GPU build environment pins CUDA 13.1 for the validated driver; other hosts need compatible driver/runtime versions. The membrane insertion relaxation is seeded, uses 0.5 fs steps with the same elapsed relaxation time as the library default, and allows 300 initial minimization iterations. These settings are recorded in system metadata.

System review precedes minimization, 0.1 ns NVT, two NPT restraint stages, a separate 20 ps checkpoint continuity test, and production. Production length is configurable. The historical acceptance used 1 ns; the current repeat protocol uses three 5 ns production trajectories with seeds 17, 41 and 97, sharing a prepared structure but drawing independent velocities. The independent checkpoint test does not splice or interrupt the production trajectory.

Analysis uses the actual production trajectory: protein-backbone aligned RMSD, ligand RMSD under the same fit, CA RMSF, periodic distance contacts, a membrane phosphorus-separation proxy, thermodynamic time series, first/last structures, CSV, PNG/PDF and UI trajectory frames. Contact occupancy is not a hydrogen-bond definition. Short MD does not establish activity, affinity, convergence or publishability.

`preliminary` remains `preliminary` after successful MD. A zero-formal-candidate screen may enter an explicitly acknowledged exploratory branch, but it does not gain a formal pass by completing that branch. New targets require appropriate structures, protonation, box, membrane model and acceptance criteria. The current membrane builder is a reviewed single-lipid model, not a universal biological-system builder.

## Reproduction and failures

Keep upstream immutable stage records, continuation source revisions, all unsuccessful attempts, operation requests/logs, handoff bundle, system inputs, phase checkpoints and raw trajectory/energy files. A successful stage is separate from scientific qualification. Reopening a completed run verifies artifacts rather than regenerating them. Export the full run for trajectory-level reproduction and the small review archive for inspection.

## Reference implementations

- [OpenMM membrane modeling](https://docs.openmm.org/latest/userguide/application/03_model_building_editing.html)
- [OpenMM ForceFields / GAFF templates](https://github.com/openmm/openmmforcefields)
- [ParmEd OpenMM topology conversion](https://parmed.github.io/ParmEd/html/openmmobj/parmed.openmm.load_topology.html)
- [ParmEd rigid-water export limitation](https://github.com/ParmEd/ParmEd/issues/1030)
- [GROMACS membrane-protein tutorial](https://tutorials.gromacs.org/docs/membrane-protein.html)

The insertion policy and tiny charge-rounding correction are project-specific implementation choices documented above, not claims that these sources prescribe them.

## Independent peptide evidence

Export with `python deployment/package_evidence.py --output /new/path/peptide-evidence.zip`. Its pure core accepts explicit receptor, original-frame receptor, peptide, candidate SDF paths and selected docking scores. Run `python -m scientific.peptide_contacts --manifest input.json --output /new/result`. The input manifest resolves paths relative to itself and preserves each candidate's compound and pose IDs. The screening orchestrator calls the same core after collecting selected poses; the UI also exposes a separately recorded comparison.

Contact recall is the fraction of peptide-contacting receptor residues also contacted by the candidate within 4.5 angstrom. The self-contact and translated-peptide controls should return Jaccard 1 and 0. The ablation compares affinity ranking with contact-recall-first ranking within the same candidate pool. This is a transparent distance metric, not a learned pharmacophore model or evidence of activity.

## Versioned qualification and repeat acceptance

Archived screening results are never silently rewritten. An explicit UI reassessment records original candidate/selection/redocking hashes, the exact selected pose hash and corrected fixed-frame symmetry-aware RMSD. Only a completed, internally consistent reassessment supplies effective qualification to the UI and subsequent MD handoff. Changes to the evidence cause refusal.

The independent stage-six process exports the scientific core and runs it in isolated Python without UI imports. It retains its request, source hashes, system phases and analysis. Repeat comparison verifies the recorded analysis and raw inputs, reports each prespecified criterion and preserves failed criteria. Short-MD-quality acceptance and convergence/activity validation are separate output fields. The protocol is frozen in each replica before execution and cannot be changed by the summary operation.
