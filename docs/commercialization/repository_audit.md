# Local repository review — bounded static assessment

Date: 2026-10-19. Baseline `c1a87a5`, local external-drive branch `codex/rtx5090-workflow-20261018`; changes prepared separately. README release is `0.4.0-rc1`. Remote GitHub visibility and its current HEAD were not rechecked in this pass; local state must not be equated with the current public repository.

| Module | Classification for this assessment | Evidence and qualification |
|---|---|---|
| Workbench | B: implemented, final acceptance open | `00_workbench/webui/workflow.py`, `pipeline_bridge.py`, `scientific/`; historical execution in docs/10 |
| Molecular generation | B | `01_molecule_generation/src/`, model pins in deployment/model_dependencies.json; historical attempts do not establish general performance |
| GNINA docking | B | `02_docking/`, workbench bridge and historical acceptance; score is not affinity proof |
| Screening existing policy | B | `03_virtual_screening/src/zvss/`; formal qualification, exclusions, diversity and handoff exist |
| Screening v2 additions | B, static final review only | Changed source and tests; final revision not runtime-accepted |
| Interaction analysis | B | `04_interaction_analysis/src/zlis/`, historical module guide; not reexecuted |
| MD paths | B | Workbench GROMACS and standalone OpenMM are distinct; final aggregate acceptance remains open |
| Enterprise platform | D: planned | Multi-user isolation, commercial API and enterprise scheduling not established by this audit |
| Validated learning loop | D: planned | Roadmap only; no large-model training performed |
| Archive/historical prototypes | E for current acceptance purposes | Retained as historical material, not the current product entry point |

A means implemented and currently demonstrated by a relevant real test; no entire module is promoted to A on this pass. C means partial implementation, D planned, E legacy/problematic. B does not mean scientific validation. Per-feature evidence and final acceptance boundaries remain in docs/09–11.

All tracked Python source is included in static syntax review. The current content map indexes tracked files without moving any. Local untracked demo/video and business-plan documents remain untouched; model/environment/result directories are not copied into this change set. No credentials or private candidate content are reproduced in these reports. A full historical secrets scan and rights review remains open; absence of a reported secret is not proof of absence.

Hardcoded paths and deployment dependencies must be reviewed in their host context. Existing model restoration, platform launchers and environment locks are retained. macOS syntax checks cannot certify Windows/WSL or Linux/GPU execution. Public-demo isolation and remote final MD state are not revalidated.

See [implementation status](../product/IMPLEMENTATION_STATUS.md), [code map](../product/CODE_AND_CONTENT_MAP.md) and [license review scope](license_audit.md).
