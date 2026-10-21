# Code-only delivery, 2026-10-21

This branch integrates the previously local Y-MAVS work and preserves the local
research UI changes before applying the Vina / HF-RESP migration. The latest owner
instructions override the initial runtime-acceptance request:

- One MD trajectory per selected complex, with 500 ns formal production.
- Modify source, configuration and documentation only; testing will be done on
  another computer.
- Delete this task's generated acceptance results and test logs.
- Synchronize the changes to GitHub through a feature branch.
- The administrator deploys from GitHub and supplies a UI URL; users operate that
  existing page manually or with Codex, with no separate computer launch entry.

No new PASS or performance claim is made. No 500 ns production was completed for
this delivery. No retained task-generated acceptance dataset is published. The
commit uses `[skip ci]` so push/PR automation does not perform the tests the owner
asked to defer. Existing project historical results and earlier release records
are retained and explicitly distinguished from this migration.

The user-facing run/resume workflow is in
[18_UI_RESEARCH_WORKFLOW.md](18_UI_RESEARCH_WORKFLOW.md). Technical migration and
existing backend references are in
[15_VINA_RESP_RESEARCH.md](15_VINA_RESP_RESEARCH.md). The dependency/license review
is in [16_RESEARCH_DEPENDENCY_LICENSES.md](16_RESEARCH_DEPENDENCY_LICENSES.md).

Pending on the other computer: dependency solve and GPU build compatibility,
Vina preparation/export/handoff and reference redocking, HF/ESP/RESP charges in
both final topologies, short GPU/continuation/analysis checks, real PB including
the chosen membrane model, UI/container/module regressions, and formal production.
The built-in membrane is an idealized single-lipid model; a bacterial lipid mixture
needs an externally prepared compatible system and its own validation. Unknown
protonation/chemistry, membrane orientation and PB slab parameters remain explicit
scientific input decisions. This branch is a draft pending that review.
