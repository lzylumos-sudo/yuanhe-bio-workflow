# Pre-docking screening v2

New research runs use `murj-discovery@1.1.0` and
`00_workbench/configs/research_predocking_v2.json`. Archived v1 runs retain their
frozen policy and evidence. No Vina scoring, MD settings, or historical protocol
values are replaced.

## Execution and identity

Peptide2Mol -> REINVENT4 (or imported CSV/SDF) -> exact stereochemical identity,
parent/salt mapping and duplicate detection -> curated SMARTS plus full RDKit
PAINS/Brenk -> ADMET-AI -> SA and SCScore -> catalog matching -> bounded
AiZynthFinder -> independent threshold gates and lexicographic/scaffold ranking
-> selected SDF/SMILES -> ligand preparation -> existing Vina/GNINA adapter
-> existing virtual screening -> existing MD. The production workbench uses
Vina 1.2.7; GNINA retains its existing evidence gate.

`pipeline_bridge.stage_risk/stage_availability/stage_docking` and the research
qualification adapter use this implementation. Vina verifies the immutable
selection manifest, screening policy, model evidence and acquisition proof.
A molecule omitted by the final cap cannot reuse an otherwise valid proof.
The compound/molecule ID and state ID survive the handoff; pose IDs remain
owned by the existing docking contract. Salt parents are recorded separately;
neutral parent similarity does not establish exact supplier equivalence.

Full PAINS and Brenk matches default to REVIEW, not proven toxicity. Curated
severe reactive alerts BLOCK. Invalid structures, duplicates and rejected rows
remain in `records.json`, CSV, SQLite and readable SDF exports (invalid SMILES
cannot be encoded as molecules in SDF). Structured endpoint results distinguish
SUCCESS, ERROR, UNSUPPORTED and NOT_EVALUATED. Missing results require review or
abort. There is no score-only bypass of the existing acquisition evidence gate.

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
pre-docking JSON field. Set resource limits and authorized supplier catalogs for
the deployment. Do not commit credentials, catalogs, local paths or checkpoints.
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
and public ZINC stock. Isolate NumPy 1.26.4 and PyArrow 19.0.1 from the model
worker environment. Search only retained priority candidates, default maximum
10. Results export tree nodes, precursors, scores, step counts, search statistics
and raw routes. A ZINC stock node does not prove present supplier availability;
computed routes still require the existing route review and evidence contract.
Route labels: ROUTE_FOUND, NO_ROUTE_FOUND, NOT_EVALUATED, ERROR.

## Acquisition and modes

Authorized local catalogs retain the existing exact stereochemistry/isotope/
salt match and dated, hash-verified stock or quotation evidence. Only those
records can become PURCHASABLE_VERIFIED. Catalog-only hits are
PURCHASABLE_UNCONFIRMED. A predicted route becomes SYNTHESIZABLE_PREDICTED;
unknown and unmatched are distinct states.

The optional PubChem PUG REST adapter resolves exact isomeric SMILES to CID and
supplier SIDs, retaining source, SKU, linked CID, HTTPS catalog link, timestamp
and raw response hashes. It retries transient errors and caches responses.
Supplier-name filtering is explicit; it does not scrape vendor inventory or
claim free Enamine/Mcule APIs. Public lookup transmits candidate structures and
is disabled by default. Enable only for structures authorized for disclosure.
No price, inventory or delivery time is inferred. Live PubChem hits alone do
not satisfy procurement verification.

All modes apply the same required evidence and independent risk thresholds.
Discovery alternates scaffolds across acquisition tiers; balanced alternates
scaffolds within acquisition tiers; delivery prioritizes acquisition evidence,
then toxicity threshold tiers and separate SA/SC objectives. Deterministic IDs
break ties. The cap is applied before qualification proofs are issued.

## Workbench and persistence

New research run -> Pre-Docking Screening controls configure mode, endpoint
selection, thresholds, rule actions, CPU/CUDA, batch size, SC weights, expensive
search budget and public supplier lookup. Advanced pre-docking JSON configures
catalogs, per-rule overrides, route models and search depth/time. The production
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

1000 unique valid RDKit NCI molecules: 30.610 s total, 32.669 molecules/s;
342 reached real AI/SA/SC inference, zero model failures; structural stage
excluded 15 and retained 643 for review. AI 5.515 s, synthesis 0.701 s,
acquisition/search 17.903 s. Two bounded real searches returned routes. Final
qualified count was zero because there was no current supplier/route delivery
evidence. The remaining records were retained. Parent maximum RSS was
312,705,024 bytes and maximum child RSS was 5,849,382,912 bytes (separate maxima,
not a simultaneous total). No GPU, docking or MD ran in that benchmark.

Independent aspirin AiZynthFinder search: 1.848 s, 20 iterations, 31 nodes,
21 solved routes in the search tree, 7 exported routes, top score 0.997629,
one-step best route. SCScore matched the upstream current NumPy implementation
on 197 structures with maximum absolute difference 8.89e-16. This is numerical
conformance, not a new experimental synthesis benchmark. Live PubChem aspirin
lookup returned Sigma-Aldrich catalog records and preserved raw responses;
current inventory/quote verification remains unavailable.

The real-model CI runs the 1000-molecule CPU workflow with SC weights, without
supplier calls or costly retro assets. It is execution validation, not a
held-out toxicity accuracy study. No test fixture vendor or injected probability
is used by that acceptance workflow. Unit tests clearly label injected fixtures.
Full clinical/toxicology validation, target NVIDIA GPU deployment, and a real
procurement-qualified positive docking/MD campaign remain outside the completed
acceptance evidence. Do not report the entire requested P0 acceptance COMPLETE.

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
- PubChem API: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest ; linked supplier
  content retains the relevant provider's terms. No supplier keys are stored.

Additional local handoff acceptance executed actual ADMET-AI and SCScore,
selected ethanol from four mixed valid/invalid/duplicate/reactive input records,
generated its conformer and ran official Vina 1.2.7 against the repository's
prepared 1DRF receptor. Three Vina poses were exported. The acquisition proof
in this isolated test was a clearly marked simulated unit fixture; this proves
the executable handoff and exclusion behavior, not a real procurement-qualified
campaign or binding activity. The 1000-molecule run above used no such fixture.
