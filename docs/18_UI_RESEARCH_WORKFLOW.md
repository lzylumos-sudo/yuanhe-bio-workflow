# Using the research workbench through a UI link

An administrator deploys the workbench from GitHub and provides its URL. Users open that URL and operate the page manually, or authorize Codex to operate the same page. Computation and file storage run on the deployment machine; users do not need a separate environment or launch entry point. The public static demo only illustrates the workflow. Real tasks use the administrator-provided research workbench.

## Daily workflow

| Page step | Actions and required review |
|---|---|
| New run | Select the structure, chains and ranges. Create a new Vina / RESP run after migration and retain historical runs. |
| Tool 1 and complex preparation | Prepare the receptor and review protonation, construct and search box. Centers and side lengths are in angstroms. |
| Generation and Tool 2 | Review generated molecules and retained core atoms. |
| Docking stage settings | Official Vina 1.2.7; start at exhaustiveness 32; configure the CPU budget and seed array; save, then execute. |
| Tool 3 and screening | Inspect Vina affinity, qualification and poses; verify compound_id / state_id / pose_id. |
| Tool 4 | Select the final pose; review membrane orientation, lipids, GAFF2 / HF-RESP, GPU, one 500 ns trajectory, output and checkpoint intervals. |
| Tool 5 | Preview the built system and confirm before starting. Resume the same system with Start / Resume from checkpoint. |
| Tools 5 and 6 | Trigger analysis; inspect ligand RMSD over 400–500 ns, sample standard deviation, full curves, contacts, hydrogen bonds and actual PB results. |
| Export | Export reproduction records; include raw trajectories when needed. |

Controls and additional trajectories in Tool 7 are optional, not prerequisites for completing one 500 ns workflow. Docking seeds and MD trajectory counts are independent settings.

## Parameters and interpretation

- **Docking:** Vina affinity is in kcal/mol; lower is better. CNN is not applicable. Historical GNINA and new results are separated by engine. Concurrent tasks multiplied by threads per task must not exceed the total CPU budget.
- **Parameterization:** Optimize the isolated ligand at HF/6-31G*, calculate ESP, fit two-stage RESP charges, and use GAFF2. Preserve the bound pose and RESP charges in the final topology. Unsupported chemistry or failures block continuation.
- **Simulation:** One 500 ns production trajectory per selected complex, excluding equilibration. Save trajectory frames every 10 ps and GROMACS checkpoints every 2 minutes by default. Record the initial velocity seed, GPU and execution information; resuming preserves velocities and the time axis.
- **Core statistics:** Align the specified stable backbone/transmembrane core, then report the mean and sample standard deviation of ligand heavy-atom RMSD over 400–500 ns. Between-replica statistics are not applicable to one trajectory. Short trajectories must not substitute their last segment for the formal window.
- **PB:** Use the same trajectory, default window 400–500 ns and every tenth saved frame. PB window, sampling, dielectric and membrane parameters are adjustable before building; salt follows the selected system. MurJ requires review of membrane thickness, dielectric, pore model and geometric evidence. Do not automatically apply a soluble-protein model. Entropy is not calculated.
- **Interpretation:** Completion and low standard deviation do not prove strong binding, activity or convergence. Results may remain unstable and require extension. The built-in ideal single-lipid membrane is not a realistic mixed bacterial membrane.

## Recovery and responsibilities

Closing the browser does not change backend run identity. After stopping, wait for the checkpoint to be written, then reopen the original run and system to resume. If inputs change, checkpoints mismatch or dependencies are missing, retain logs for the administrator. Do not overwrite old runs or regenerate initial velocities to concatenate trajectories.

The administrator manages GitHub deployment, dependencies and models, GPU allocation, paths and access permissions. Workbench environment settings are administrator configuration. Existing standalone OpenMM and GROMACS commands remain backend maintenance interfaces; users do not need a terminal.

The original migration delivered source files only, without running the modified calculations, tests or builds. It did not establish runtime acceptance. See [method and configuration migration](15_VINA_RESP_RESEARCH.md) and [component licenses](16_RESEARCH_DEPENDENCY_LICENSES.md).
