# Independent module UI

Open /modules.html from **Independent modules** in the main workbench.
No research run or completed upstream stage is needed. Save an input snapshot,
then start the module. A saved job is not a completed calculation.

## Inputs and execution

| Entry | Independent input | Execution and evidence |
| --- | --- | --- |
| Complex preparation | Protein/peptide PDB, chains, residue range | Existing preparation; receptor and pocket artifacts |
| Peptide2Mol | Prepared pocket PDB, attempt count, seed | Pinned real model; raw output and actual attempt report |
| REINVENT4 | Reviewed scaffold SMILES with exactly one labelled [*:1], count, seed | Pinned real model; batch size one supported |
| Substructure risk | SDF, CSV or tab-separated SMILES table | Existing strict rules, exclusions and molecule records |
| Pre-docking availability | The same library formats | Both gates, including toxicity and complete-route evidence |
| Candidate docking | One ligand SDF, receptor PDB, reviewed box | Both gates recomputed before Vina; rejected ligands are not docked |
| Public docking control | No upload | Hash-pinned public 1DRF/FOL reference, three seeds and RMSD |
| Post-docking screening | Receptor PDB and scored pose SDF in the same frame | Generic protocol; missing evidence remains missing |
| Interaction analysis | Complex PDB, ligand selector, optional template SDF | Offline chemistry, tables, diagrams and scene files |
| MD system construction | Receptor PDB, ligand SDF, reviewed membrane orientation | Build, inspect, review, run/resume, then analyze |
| Prepared MD | Self-contained production TPR, device, threads, review | Run recorded duration, preserve checkpoints, audit completion |

SMILES tables require id and smiles headers. A .smi file is tab-separated,
for example an id/smiles header and a row containing PUBLIC_ETHANOL/CCO.
Each file is limited to 64 MiB; total decoded input is limited to 96 MiB.
Uploaded filenames are never used as server paths. Executables and environments
come from administrator settings, not uploaded files or browser input.

Post-docking screening binds its generic target to the exact uploaded
single-chain receptor. It does not borrow the demo receptor or infer missing
membrane, redocking, synthesis or biological evidence. Molecules with missing
requirements cannot become formal candidates.

The membrane builder retains the existing force fields, charge calculation,
equilibration and analysis. Independent duration settings do not modify the
formal research protocol. Short runs are engineering experiments, not a
replacement for the research trajectory requirement. A prepared TPR must
record a compressed trajectory and contain a positive finite dynamics run.
Use the same or a compatible GROMACS version that created the TPR.

## Records, cancellation and restart

Jobs are separate from research runs under 00_workbench/module_jobs/ and are
excluded from Git and release source archives. An administrator may set
YUANHE_MODULE_JOBS to an absolute persistent directory before starting the
service. Preserve that directory when redeploying.

Each job records file hashes, frozen parameters, administrator settings and a
source snapshot. Each action has its own log and result. The UI lists
downloadable artifacts, including failures. Changed source or input hashes
block execution: create a new job instead of modifying existing records.

Use **Stop job** to stop its worker and engine descendants. MD retains its
checkpoint. After a service restart the UI distinguishes a live worker from
interrupted records. Prepared MD retries use the same input, engine and
settings; partial output without a checkpoint blocks resuming.
Completed results cannot be silently overwritten.

Prepared MD success requires:

- Finished mdrun and the exact planned final step/time.
- A final checkpoint with the planned step.
- A readable compressed trajectory, final coordinates and energy file.
- No detected fatal or constraint warning.
- A completion record containing hashes of the audited files.

These are execution checks, not evidence of binding or convergence.
Membrane jobs retain their separate system review and trajectory analysis.

## Reference controls and strict candidate gates

The 1DRF control has its own yuanhe.reference-control/1 evidence. Only the
bundled public receptor/ligand hashes are accepted, even if another molecule
is renamed to match the reference filename. Its outputs explicitly set
candidate_qualified to false and do not create qualification sidecars.
The ordinary Vina candidate adapter and strict evidence validation are unchanged.

Independent jobs never write research-stage completion markers or count as
qualified research candidates. Calculation completion and scientific candidate
acceptance remain distinct states.

## REINVENT4 batch-size-one fix

The pinned LibInvent encoder used unrestricted tensor squeezing, removing
the batch dimension at size one. Its sampling loop also squeezed that
dimension from logits and target labels. The vendor patches preserve these
dimensions without increasing the requested sample budget.

deployment/vendor_patches/manifest.json records exact before/after hashes.
deployment/apply_vendor_patches.py checks all files before applying patches.
Readiness checks every corresponding installed Python module, including
encoder and sampler; an old installed copy cannot pass just because its
source checkout was updated.

For a fresh install, follow [DEPLOY.md](../DEPLOY.md). Installation applies
patches before installing the model environment. For an existing installation,
stop scientific jobs, preserve runtime settings and job directories, deploy a
fresh checkout and rebuild REINVENT4 from the patched pinned source using
the documented fresh prefix. Point runtime settings to that environment and
rerun diagnostics. Do not edit patch hashes or reuse an old installed package.

## Reproducible acceptance

- python -m pytest 00_workbench/tests/test_independent_modules.py -q:
  independent inputs, strict gates, tampering, MD review and incomplete output
  rejection. Command doubles in this suite are not real simulations.
- node 00_workbench/tests/run_module_ui.cjs: actual browser and real preparation,
  risk screening and interaction calculations. Set MODULE_TEST_GMX to enable
  the real prepared-MD browser test.
- deployment/vendor_patches/test_libinvent_batch.py: official prior and real
  CLI for one and two samples, plus encoder shape checks.
- deployment/validate_reference_control.py: real fixed public Vina control
  and verification that reference evidence cannot qualify a candidate.

CI retains logs and test artifacts. Local or CI success does not mean the
remote server has been redeployed. Verify its deployed revision and runtime
diagnostics before testing actual project inputs.
