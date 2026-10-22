# Pre-docking screening: synthesis-only policy v3

New research runs use `murj-discovery@1.2.0` and
`00_workbench/configs/research_synthesis.json`. Historical records and frozen v1/v2 protocols remain readable and unchanged.
New docking requires schema-2 synthesis evidence; historical purchase proofs
cannot be reused as synthesis qualification. Re-screen in a new task. No Vina scoring, MD settings, or historical protocol
values are replaced.

## Execution and identity

Peptide2Mol -> REINVENT4 (or imported CSV/SDF) -> exact stereochemical identity,
parent/salt mapping and duplicate detection -> curated SMARTS plus full RDKit
PAINS/Brenk -> ADMET-AI -> SA and SCScore -> bounded
AiZynthFinder -> independent threshold gates and lexicographic/scaffold ranking
-> selected SDF/SMILES -> ligand preparation -> existing Vina/GNINA adapter
-> existing virtual screening -> existing MD. The production workbench uses
Vina 1.2.7; GNINA retains its existing evidence gate.

`pipeline_bridge.stage_risk/stage_availability/stage_docking` and the research
qualification adapter use this implementation. Vina verifies the immutable
selection manifest, screening policy, model evidence and synthesis-route proof.
A molecule omitted by the final cap cannot reuse an otherwise valid proof.
The compound/molecule ID and state ID survive the handoff; pose IDs remain
owned by the existing docking contract. Salt parents are recorded separately;
the route must match the exact target structure.

Full PAINS and Brenk matches default to REVIEW, not proven toxicity. Curated
severe reactive alerts BLOCK. Invalid structures, duplicates and rejected rows
remain in `records.json`, CSV, SQLite and readable SDF exports (invalid SMILES
cannot be encoded as molecules in SDF). Structured endpoint results distinguish
SUCCESS, ERROR, UNSUPPORTED and NOT_EVALUATED. Missing results require review or
abort. There is no score-only bypass of the complete synthesis-route gate.

## Installation

Use Python 3.12 and independent environments:

```sh
python3 -m venv .venv-admet
.venv-admet/bin/python -m pip install -r deployment/predocking-model-requirements.txt
python3 -m venv .venv-retro
.venv-retro/bin/python -m pip install -r deployment/predocking-retro-requirements.txt
.venv-admet/bin/python deployment/install_predocking_assets.py \
  --output /absolute/path/predocking-assets --with-retro \
  --tox-python /absolute/path/.venv-admet/bin/python \
  --retro-python /absolute/path/.venv-retro/bin/python
```

The installer checks official download SHA256 values and writes an asset
receipt and `predocking-config.json`. Apply this configuration in the new-run
pre-docking JSON field. Set resource limits for the deployment. Do not commit local paths or checkpoints.
Model paths are mandatory readiness requirements, not silently downloaded at
first inference. Local installation is not remote GPU deployment.

ADMET-AI 2.0.1/Chemprop 2.3.1 execute bundled weights for AMES, hERG, DILI
(probabilities), and LD50_Zhu (log10 inverse mol/kg). Higher LD50_Zhu means lower
predicted lethal molar dose; it is not a probability. Applicability-domain
assessment is UNSUPPORTED. CPU is the default; CUDA honors scheduler-visible
devices, but NVIDIA hardware acceptance is separate. Per-molecule caches include
weights, implementation/package versions and configuration; batch failures are
retried per molecule. Failed predictions are not cached as success.

SA uses RDKit's Ertl implementation and fragment data. SCScore uses the actual
Coley 1024-bit boolean network, pinned to upstream commit
`37090a6aa8220408f572de20224c33c221312d16`. Both increase with difficulty, but SA
1..10 and SC 1..5 are never added. Defaults are SA <= 6, SC <= 4. Independent
review thresholds for AMES/hERG/DILI are 0.8 and LD50_Zhu is 3.0. These are
configurable research triage thresholds, not clinical safety guarantees.

AiZynthFinder 4.4.1 uses the pinned USPTO expansion policy, reaction templates
and a public ZINC planning building-block pool. Isolate NumPy 1.26.4 and
PyArrow 19.0.1 from the model worker environment. Search only retained priority
candidates, default maximum 10. Results export tree nodes, precursors, scores,
step counts, search statistics and raw routes. The complete tree must have
`route_metadata.is_solved=true`; every terminal node must be in that versioned
planning pool. The raw search, model, templates and pool are hash-bound and
rechecked at docking. RDKit replays every reaction template against the exact
product and precursors. Stereochemistry and ambiguous template outcomes require
an explicit route-hash-bound human review. Missing steps or unresolved terminal
nodes cannot be approved away. A predicted route is not experimental synthesis.
Route labels: ROUTE_FOUND, NO_ROUTE_FOUND, NOT_EVALUATED, ERROR.

## Synthesis policy and modes

Procurement adapters, vendor catalogs, supplier queries, inventory labels and
purchase-based bypasses have been removed. No candidate structure is transmitted
to a supplier service. Planning-pool membership is a computational stopping
criterion, without a claim about current material inventory.

All modes apply independent risk, toxicity, SA/SC and complete-route gates.
Discovery alternates scaffolds across synthesis tiers; balanced alternates
within tiers; delivery prioritizes validated routes, then toxicity tiers and
separate SA/SC objectives. Deterministic IDs break ties. The cap is applied before
qualification proofs are issued. Candidates outside the route budget remain
NOT_EVALUATED/review; low difficulty scores alone do not prove a route.

## Workbench and persistence

New research run -> Pre-Docking Screening controls configure mode, endpoint
selection, thresholds, rule actions, CPU/CUDA, batch size, SC weights, expensive
search budget. Advanced pre-docking JSON configures per-rule overrides, route
models and search depth/time. The production
stage between generation and docking shows counts, pass rate, timing, rejection
reasons, individual predictions/synthesis and complete evidence. CSV/JSON/SDF
links refer to actual backend artifacts.

`GET /api/project/run/<id>/predocking/summary` returns actual persisted summary,
policy and export links. Existing rows and pre-docking detail APIs return the
full molecule record. Run cancellation kills the managed stage process group;
resume uses existing recovery/integrity checks and per-molecule model caches.
Explicit evidence revision preserves previous artifacts. No historical database
migration is required; each stage adds a local SQLite result store.

## Acceptance scope (2026-10-23 local macOS CPU)

1000 unique valid RDKit NCI molecules: 43.779 s total, 22.842 molecules/s;
342 reached real AI/SA/SC inference, zero model failures. Structural rules
blocked 15 and kept 643 for review. AI inference: 15.276 s; synthesis scoring:
0.717 s; bounded retrosynthesis: 20.052 s. With a search budget of two, both
priority molecules had complete validated routes and qualified for docking;
other molecules remained blocked or explicitly pending. Parent maximum RSS:
310,280,192 bytes; child maximum RSS: 4,956,389,376 bytes (separate maxima).
These are execution measurements, not held-out toxicology accuracy results.
Both selected molecules subsequently passed real Vina 1.2.7 against the prepared
1DRF receptor. The other 998 records were not docked. Screening and docking use
the same chemistry runtime: mismatched RDKit versions invalidate rule evidence
and require re-screening, rather than silently reusing changed rules.

Run `deployment/validate_predocking.py` with installed SCScore weights and
`--retro-config` / `--retro-python` to reproduce the real batch. Then use
`deployment/validate_synthesis_docking.py --screening <run>/availability
--vina <binary> --receptor <prepared.pdb> --box <box.json> --out <new-folder>`
for real selected-only Vina acceptance. Neither command injects route evidence
or model predictions. Real-model CI uses the 1000-molecule CPU path without
large retro assets; missing routes remain pending rather than passing.

The unchanged 127-file protected scientific baseline includes docking and MD
core settings. NVIDIA GPU and complete MD execution require target-host
acceptance; this local CPU test does not establish GPU or experimental results.

## Licenses and sources

- ADMET-AI: MIT, https://github.com/swansonk14/admet_ai ; bundled model artifacts,
  with upstream TDC dataset terms applicable separately.
- Chemprop: MIT, https://github.com/chemprop/chemprop .
- RDKit/SA implementation: BSD-3-Clause, https://github.com/rdkit/rdkit .
- SCScore: MIT, https://github.com/connorcoley/scscore ; retained MIT notice in
  `licenses/scscore-MIT.txt`. Proprietary Reaxys training data are not copied.
- AiZynthFinder: MIT, https://github.com/MolecularAI/aizynthfinder . USPTO weights
  and templates are CC-BY-4.0 from Zenodo 7797465 and 7341155. Public stock
  redistribution/commercial rights still depend on upstream supplier terms;
  the installer records this limitation and does not claim certification.
