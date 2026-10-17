# Workbench web interface

The current main page is static/index.html with studio.js and studio-tools.js. It connects run management, receptor/box review, core selection, screening budgets, membrane construction, GROMACS recovery, analysis and export. See the root USER_GUIDE.txt for daily operation.

Start `python run_webui.py --no-browser` from 00_workbench, or use the root platform launcher. The server binds to loopback by default and advances to an available port when necessary. Do not expose this unauthenticated local service publicly. Host/origin checks, request limits, path containment and explicit compute ownership protect local operations; they are not a multi-user authentication system.

The server uses the Python standard library; frontend modules use plain JavaScript and Canvas 2D with no build step. Scientific tasks and RDKit structure drawing require their configured environments. Starting a dependency-light server is not proof that calculations can run.

## Current workflow implementation

- workflow.py: frozen configuration, source/input hashes, stage states, shared compute lock, process identity, cancellation and attempt archives.
- studio.py / studio_worker.py: bounded preparation/review operations, isolated executors, diagnostics and archive generation.
- pipeline_bridge.py: real model, docking, screening and MD stage adapters.
- scientific/: membrane parameterization, staged GROMACS, checkpoint checks and analysis.
- server.py: HTTP routing, static files and error handling.

Run-scoped logs and unique-molecule counts update during execution. A completed process is not automatically a completed or scientifically qualified stage. Failed operations preserve their diagnostics.

## Historical analysis pages

console.html and legacy.html retain earlier WEDGE-Design views: overview, scaffolds, hotspots, pharmacophores, pocket, 3D scene, evidence, jobs and documents. They read historical result files rather than fabricate them. Empty output directories are expected in a fresh GitHub checkout. The legacy enumeration launch is archived; it is not a replacement for the actual model pipeline.

| View | Source artifacts |
|---|---|
| Overview | summary.json, figures and availability metadata |
| Scaffolds | ranked CSV, structure filters and pose files |
| Hotspots | interface_hotspots.csv and residue evidence |
| Pharmacophores | pharmacophore.json and axes.json |
| Pocket | pocket.json and cached coordinate-derived voxels |
| Structure | centralized PDB inputs and optional ligand poses |
| Evidence | report.md, prior checks, census, interfaces and anchors |
| Documents | repository Markdown files |

Canvas controls: drag rotates, right/Shift-drag pans, wheel zooms, double-click resets, hover identifies atoms/features and click selects details. Shared receptor C-alpha alignment places peptide structures and historical poses in a common frame; switching references transforms the whole scene. Canvas depth sorting approximates occlusion and does not provide full cartoon/surface rendering.

## Legacy read endpoints

GET /api/health, /api/meta, /api/overview, /api/hotspots, /api/pharmacophore, /api/axes, /api/pocket, /api/scaffolds, /api/filters and /api/poses return structured data. /api/pocket/points may return 202 while computing. /api/scene accepts ref, ligand and radius. /api/evidence/report, priors, census and interface expose the corresponding evidence. /api/docs lists documents and /api/docs/<name> reads one. /api/jobs lists jobs; /api/jobs/<id>?since=N retrieves incremental logs. Allowed file downloads are confined to explicit project artifact directories. Missing static assets return 404 instead of HTML.

## Verification

```bash
python -m pytest tests/test_workflow_regressions.py tests/test_studio_integration.py
node webui/test_pages.mjs
node --experimental-vm-modules webui/check_js.mjs
```

Run commands from the workbench directory with appropriate interpreters. The Node tests use a local DOM/Canvas shim and labeled snapshots; they check rendering and routing logic, not pixels or real scientific execution. The historical test_api.py expects historical outputs and can optionally launch a real legacy smoke job; it is not a clean-checkout production acceptance command.

Troubleshooting: inspect operation logs for blocked prerequisites; refresh after output changes; wait for initial cavity voxelization; use the actual displayed port; create a new run after source/input changes. Historical result caches are separate from immutable current-run records. Stop and confirm process exit before moving the disk.
