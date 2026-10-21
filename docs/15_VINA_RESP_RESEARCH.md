# Vina / HF-RESP research migration

## Delivery scope

This is a source, configuration and documentation change. At the owner's request,
no further calculations, tests, builds or deployments are performed for this
revision. Task-generated acceptance results and test logs were deleted. No runtime
PASS, 500 ns trajectory, convergence or activity claim accompanies this revision.
Run validation on the intended compute computer before starting research jobs.
Historical project results remain historical and are not reclassified.

**Current MD default: one 500 ns production trajectory per selected complex.**
Minimization, heating and equilibration are additional. There is no requirement to
create three repeats. One trajectory cannot provide between-repeat uncertainty.

## User workflow

The administrator deploys this GitHub source and supplies the workbench URL. Users
operate the existing UI manually or with Codex; no separate compute-computer entry
is required. See [the UI guide](18_UI_RESEARCH_WORKFLOW.md). Commands below document
existing backends for administrators and maintainers, not required user steps.

## Docking and score migration

New workbench jobs, SBVS/SBVS2 jobs and screening redocking use official
[AutoDock Vina 1.2.7](https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7).
The new adapter checks the executable's version banner. It does not invoke GNINA,
a CNN model, a GPU docking mode or an automatic GNINA fallback.

Use `deployment/runtime.example.json`, `02_docking/configs/research_vina.yaml`
and `03_virtual_screening/configs/host.yaml` as new-task templates. Replace the
executable paths and review the receptor and search box. Center and full box
size are in angstroms, never nanometers. Defaults: scoring `vina`, exhaustiveness
32, nine modes, energy range 3 kcal/mol. Configure seeds explicitly; docking
multiseed sampling is separate from the single MD trajectory policy.

Workbench Vina scheduling enforces `docking_workers * docking_cpu <=
docking_total_threads`. SBVS2's process lease serializes jobs in its CPU slot;
SBVS uses one worker. Do not run separate independent queues with budgets that
jointly exceed the host allocation. Increase exhaustiveness and CPU budgets in
new frozen protocols, never mutate an active run.

`vina_affinity_kcal_mol` is a raw docking energy estimate, lower is better.
`cnn_pose_score` / `cnn_affinity_pK` are null and `cnn_status` is
`not_applicable`. No zero-filled CNN values or Vina-as-CNN aliases are created.
Vina aggregate ranking omits CNN weights and renormalizes available evidence.
An affinity cutoff is a protocol-specific filter, not a universal activity
threshold. Calibrate it for the target with suitable controls. Raw GNINA fields,
engine tags and historical license notices remain readable. Mixed-engine
ranking is refused: import and compare them in separate runs. Y-MAVS retains its
membrane evidence and qualification gates; Vina alone does not satisfy them.

Caches include engine, version, protocol, box, seeds and input hashes. Screening
redocking also verifies cached SDF hashes. Start a new run after migration;
never rename an old GNINA directory into a Vina run or reuse its score scale.
Old GNINA executable installers/launchers now refuse new execution. Historical
GNINA configuration examples are archival, not current launch presets.

Meeko 0.8.0 and RDKit prepare the ligand from a trusted chemical SDF. Mapped
PDBQT remarks carry chemical graph and atom indices into SDF reconstruction.
The adapter rejects missing mappings or changed stereochemical identity.
Docking charges are preparation artifacts, not MD charges. Keep compound_id,
state_id, pose_id, the selected pose SDF and handoff hashes together. Do not
substitute a different conformer after selection.

## Administrator reference: dependency configuration

The pinned core environment is `deployment/research-environment.yml`. The
standalone MD and membrane environment files include Psi4 and RESP as well.
Generation models retain their separate environments and licenses.

```bash
conda env create -f deployment/research-environment.yml
conda activate yuanhe-research
conda list --explicit > /path/to/evidence/environment-platform.lock
python -m pip freeze > /path/to/evidence/pip-freeze.txt
```

The Vina conda package must expose the actual `vina` executable; Python bindings
alone are not sufficient. The alternative official Linux x86-64 download is:

```bash
SBVS_PROJECT_ROOT="$PWD/02_docking" python 02_docking/scripts/install_vendor.py vina
```

The installer records the downloaded SHA-256 and whether a publisher digest was
available. A locally computed digest is not described as publisher verification.
Use a GPU-enabled GROMACS build for membrane production. The research Dockerfile
is a build recipe, not a tested image. NVIDIA drivers/container runtime and the
host's CUDA-compatible package solve still need target-host verification.

## HF / RESP / GAFF2

The production ligand path uses Psi4 1.11, cdsgroup RESP 1.0.0 and AmberTools
24.8. The supported ordinary organic protocol uses a separate gas-phase
HF/6-31G* optimized conformer, HF/6-31G* ESP, four van der Waals shells
(1.4, 1.6, 1.8, 2.0; density 1 point/A2) and two-stage hyperbolic RESP
(a=0.0005 then 0.001; b=0.1). Stage two refits alkyl CH2/CH3 groups with
equivalent hydrogens and freezes other atoms. If all atoms are fixed, the
second-stage status states that explicitly. This is a single-conformer charge
model; it is not an automatic multi-conformer RESP study.

The supplied protonation, formal charge, multiplicity, input atom order,
optimized conformer, grids, ESP values, constraints, fitted charges and package
versions are recorded. H/C/N/O/F/P/S/Cl closed-shell connected ligands are the
supported input class. Other elements, radicals, undefined stereochemistry,
metals, covalent ligands and ambiguous states require a separate protocol.
RESP failure blocks parameterization. AM1-BCC is never a fallback.

The final force field is ff14SB + GAFF2 2.11 + TIP3P/compatible NaCl; the
membrane builder adds OpenMM Amber14 Lipid17 parameters. GAFF1 is an explicit
option of the isolated parameterization CLI, not a separately validated full
protein/membrane preset. RESP charges are passed to the template generator;
charge-vector agreement is checked again on the final MD System. The independent
AmberTools path imports charge values using antechamber `-c rc`, then checks
MOL2, Amber and GROMACS topology charges and the original bound coordinates.
Topology export must retain Lorentz-Berthelot rules and Amber 1-4 scaling.

```bash
python 05_molecular_dynamics/modules/parameterization/hf_resp.py \
  --ligand /path/to/selected-pose.sdf --output /path/to/new-parameters --gaff gaff2
```

The bound SDF remains the coordinate source for complex assembly. The QM
optimized SDF must never replace it. Cached RESP evidence is hash checked.

## Single-trajectory production and recovery

OpenMM: 2 fs, minimization, a 100 ps heating ramp from 50 K within 1 ns NVT,
1 ns NPT, then one 500 ns production trajectory. The production clock starts
at zero exactly once. Output defaults to 10 ps, committed segments to 1 ns and
inflight checkpoint copies to 100 ps. A segment is committed only with closed
trajectory/log files and a checkpoint/hash manifest. Recovery replays an
interrupted segment from the preceding committed checkpoint; it never appends
unverified partial frames. The inflight copy is diagnostic, not an automatic
partial-DCD resume point. Lost work is bounded by the segment length.

```bash
cd 05_molecular_dynamics
bash run_md.sh --input /path/to/handoff/complex.pdb \
  --config config/default.yaml --run-id complex_A_single500 --non-interactive
bash run_md.sh --resume runs/complex_A_single500 --non-interactive
```

Missing scientific input decisions produce a pending file. Supply those decisions
on resume; do not bypass chemical mapping or preparation errors. Preserve the
whole run directory. Existing historical settings/checkpoints are immutable;
this release does not silently convert old AM1-BCC or short MD jobs into formal
RESP/500 ns jobs. GPU selection is explicit; research mode refuses CPU fallback.
OpenMM records the actual platform, device and precision. CUDA mixed precision is
the research default; changing precision/hardware can invalidate binary checkpoints.
`config/smoke.yaml` is a separate short engineering configuration, never formal MD.

For independent complexes, edit `deployment/gpu-queue.example.json` with absolute
paths, fixed run IDs and explicit resume commands, then on the compute host:

```bash
python deployment/gpu_queue.py /path/to/queue.json --gpus 0,1 \
  --output /path/to/queue-records
```

The queue grants one lease per GPU and sets process-local device index 0. Use the
same physical GPU identifiers and lock directory in all cooperating queues.
Re-invocation uses the recorded job's resume command. Each selected complex still
has exactly one trajectory; queue parallelism does not create repeats.

The workbench GROMACS handoff remains `yuanhe.md-input/1`. The bundle preserves
selected molecule/compound/pose identity and chemistry. Its execution commands are:

```bash
PYTHONPATH=00_workbench python -m scientific.md_cli validate --bundle /path/to/bundle
PYTHONPATH=00_workbench python -m scientific.md_cli build --bundle /path/to/bundle
PYTHONPATH=00_workbench python -m scientific.md_cli run --bundle /path/to/bundle \
  --system-id MD_ID_FROM_BUILD --review-system
# Identical run command resumes with -cpi and -append; do not regenerate velocities.
PYTHONPATH=00_workbench python -m scientific.md_cli analyze --bundle /path/to/bundle \
  --system-id MD_ID_FROM_BUILD
```

Keep `00_workbench` and `05_molecular_dynamics` sibling directories in source
exports: the membrane engine shares RESP and research-analysis code. See
`deployment/membrane-research.example.json` for parameters. Review membrane
orientation, construct, thickness and the stable TM alignment selection before
building. The supplied bilayer builder supports idealized single-lipid POPE,
POPC, DPPC or DOPC, not an arbitrary lipid mixture or a full bacterial membrane.
Record that model limitation; a mixed bacterial membrane requires an externally
prepared and separately validated compatible system. GROMACS uses NVT heating,
restrained NPT relaxation and semiisotropic pressure coupling; production uses
Parrinello-Rahman. Its production TPR encodes the full 500 ns. Periodic checkpoints
and append recovery do not redefine total production length.

## RMSD, stability and PB analysis

The core metric is ligand heavy-atom RMSD after the prespecified protein backbone
or TM core fit, in A, over **400–500 ns**. Bonds reconstruct PBC fragments, protein
chains and the intact ligand are assembled in the same image, then the protein
selection is aligned to the recorded reference. There is no ligand self-fit.
Outputs include the mean, sample SD (`ddof=1`), 10 ns blocks, full protein/ligand
RMSD, RMSF, contacts, distances and hydrogen-bond occupancies, with CSVs and plotting
code. A trajectory shorter than 500 ns returns `insufficient_trajectory`; it does
not borrow its last 100 ps/ns as a formal window. Low SD and elapsed 500 ns alone
never establish stability, convergence, binding strength or activity. The result
may require extension. With one trajectory, between-repeat comparison is N/A.

MM/PBSA actually invokes AmberTools MMPBSA.py with an `&pb` block and requires an
EPB result. Default sampling is 400–500 ns, every tenth saved frame (100 ps at the
10 ps output default), indi=1, exdi=80, salt=0.15 M, mbondi2 radii, no entropy.
It uses the same trajectory/topology and records selected frame indices, times,
parameters, raw energy terms, hashes, mean, frame sample SD and a correlation-
adjusted SEM. The SEM assumes stationarity and sufficient sampling; it is not
proof of convergence. Alternative post-equilibration windows may be configured
for sensitivity checks; full production length remains required in formal mode.

Membrane PB requires an explicit reviewed slab/pore model and source. The code
uses memopt=1, emem, mthick, mctrdz=0, poretype and bcopt=10; preserves the z normal
and centers from explicit bilayer phosphate coordinates. It refuses ordinary
soluble-protein PB settings for the membrane path. Slab thickness, dielectric,
pore treatment and geometry must be justified for the actual MurJ system.
Membrane PB runtime and parameter suitability remain for target-host validation.
End-state PB values support comparisons under the same protocol; they are not
experimental binding free energies or entropy-complete affinities.

## Optional maintainer tools (not executed for this delivery)

`deployment/validate_vina_reference.py` runs a public 1DRF/FOL reference redocking
with symmetry-aware RMSD in the unchanged receptor frame. `deployment/validate_research_gpu.py`
accepts a real HF-RESP ethanol parameter directory and checks GPU continuation,
trajectory analysis and real PB on a small solvated model. These are explicitly
invoked tools; they do not run during installation/import. Run them, module
regressions and UI/container checks on the intended computer. No results from
this task are included. A complete MurJ 500 ns run and membrane PB remain unrun.

License obligations remain component specific. See
[dependency and distribution notes](16_RESEARCH_DEPENDENCY_LICENSES.md).
