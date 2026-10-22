# Dependency graph

```mermaid
flowchart TD
 UI[Research Workbench Web UI] --> API[webui.server existing API and versioned governance API]
 CLI[run_research.py / yuanhe] --> WF[webui.workflow]
 API --> WF
 WF --> CONFIG[research.config / protocol and configuration resolver]
 WF --> RUN[pipeline_bridge subprocess]
 RUN --> GEN[Peptide2Mol and REINVENT4 adapters]
 RUN --> QUAL[yuanhe_preflight / research.qualification]
 RUN --> DOCK[research.docking -> Vina public adapter]
 RUN --> SCREEN[ZVSS -> Y-MAVS membrane evidence]
 RUN --> QUEUE[research.md_queue]
 QUEUE --> BUILD[OpenMM membrane builder / HF-RESP / GAFF2]
 BUILD --> CONVERT[ParmEd topology conversion]
 CONVERT --> GMX[GROMACS membrane production]
 GMX --> ANALYSIS[MDAnalysis / MMPBSA]
 ANALYSIS --> RESULTS[Stage records / artifacts / reviewed handoff]
 SOLCLI[Standalone MD CLI] --> SOL[OpenMM soluble protocol]
 SOL --> SOLAN[Standalone trajectory analysis]
 ZLIS[ZLIS CLI / Web] --> STATIC[Static interactions / PyMOL]
 HIST[sbvs and sbvs2 legacy CLI] --> VINA[Official Vina adapter]
 PUBLIC[Public React demo] --> FIXTURE[Precomputed static fixtures; no GPU API]
```

`research` implements budget, identity, qualification, orchestration and queues; `scientific` implements membrane construction/execution/analysis and portable handoffs; `modules/automation` is the standalone OpenMM workflow. They are not interchangeable layers. Existing `webui.jobs`, docking schedulers and GPU queue have different persistence and lease boundaries. Preserve them until a cross-process scheduling migration is validated.

Full AST imports and subprocess/environment records are in `tracked_file_audit.json`. Dynamic paths include runtime.json, configured Python interpreters, Vina/GROMACS executables, reviewed receptor paths, model directories, result bundles and source-revision exceptions. File checksums do not prove dynamic third-party behavior.
