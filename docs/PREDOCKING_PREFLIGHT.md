# Pre-docking risk and acquisition screening

Execution order: REINVENT4 → structure standardization and exact deduplication → `risk` → `availability` → AutoDock Vina → existing screening/MD. Peptide2Mol seeds enter before risk and cannot be appended at the docking entry point. Both stages execute actual checks; UI labels or pass booleans are not qualification evidence. Without supplier catalogs, models and stock, the defaults produce no acquisition-pass labels.

## Configuration and entry points

- Rules: `02_docking/src/yuanhe_preflight/configs/risk_rules.v1.json`.
- Configuration template: `02_docking/src/yuanhe_preflight/configs/predocking.example.json`.
- Workflow configuration: `predocking` in run `config.json`. Stage settings live in `reviews/risk.json` and `reviews/availability.json`, then freeze as stage `review.json` and `config.snapshot.json` at execution. Merge effective configuration per stage before passing it to the actual module.
- Dependencies: `deployment/preflight-environment.yml`, separate from MD. The AiZynthFinder adapter targets 4.4.1. Actual assets are neither bundled nor automatically downloaded.
- License inventory: `deployment/preflight_licenses.json`; upstream licenses and interface-source hashes are in `docs/preflight_sources/`.

Run the UI from the repository root:

```sh
python 00_workbench/run_webui.py --no-browser
```

Edit `predocking` JSON under Substructure risk screening / Purchase and synthesis availability → Stage settings. Absolute paths are recommended for `rules_path`, `cache_dir`, catalogs, model/template/stock assets and the retrosynthesis Python executable. Configure rule actions, human exceptions, catalog validity, search time, maximum depth, iterations, workers and batch size. These settings do not include MD parameters.

For external supplier libraries, set risk-stage `imported_candidates` to SDF or a headered CSV/TSV/SMILES file. Receptor preparation and review remain necessary before docking. External libraries may skip generation but cannot skip risk/availability.

The standalone CLI performs only pre-docking screening, not generation, docking or MD:

```sh
PYTHONPATH=02_docking/src python -m yuanhe_preflight.cli \
  --input /absolute/library.sdf \
  --config /absolute/predocking.json \
  --out /absolute/new_preflight_run
python deployment/preflight_doctor.py --config /absolute/predocking.json
```

CSV requires `smiles` and may include provenance fields such as `id`, `core` and `seed_id`. `.smi/.smiles/.tsv` use headered tab-separated records. A plain SMILES list without supply fields is candidate input, not purchase evidence. Unspecified stereocenters are not interpreted as a specific enantiomer; SMILES/SDF identity disagreement after conformer generation requires review.

## Risk rules

Defaults are a small curated set with positive and negative examples, not a comprehensive toxicity predictor; `ALL` is not enabled.

| Category | Default rule/action |
|---|---|
| Toxicity risk | Project-defined N-nitrosamine pattern, exclude. FDA sources support the risk context, not authorship of this SMARTS. |
| Nonspecific reactivity | Brenk acyl halides and isocyanates: exclude. Peroxides and project-adapted NIH epoxide alert: review. |
| Assay interference | PAINS catechol_A(92): retain with a warning; not a toxicity conclusion. |
| Other medicinal-chemistry alerts | Brenk nitro: warn. |

Each rule has an ID, SMARTS, category, severity, action, source, version, positive example and negative example. Startup checks SMARTS and FilterCatalog and runs each example. Invalid rules fail the task rather than being ignored. Match indices refer to zero-based heavy-atom indices in standardized SMILES; mappings to original SDF atom indices are retained. SVG output shows indices and highlights.

The shared preflight standardizer sanitizes structures, removes atom-map numbers and produces canonical isomeric SMILES. It does not automatically desalt, neutralize, merge tautomers or remove isotopes. The research orchestration profile additionally normalizes parent identities before this shared stage, as described in RESEARCH_DUAL_BLACKWELL.md. Duplicates retain original rows and all provenance. Exclusion removes the whole molecule from progression, while retaining original inputs and excluded SDFs. No hit is reported only as "No match to current rules".

`risk_policies` looks like `{"BRENK_NITRO":"review"}`; actions are `exclude/review/warn`. Human exceptions in `risk_exceptions` record `rule_id`, canonical `structure_smiles`, `rules_sha256`, `reviewer`, `reason` and `reviewed_at`; exceptions cannot erase original alerts. To add review after risk completes, use Revise pre-docking evidence: archive the old stage, edit configuration and rerun the pre-docking stages without regenerating molecules. Runs with docking or MD results cannot have their evidence retroactively edited in place.

## Actual purchase evidence

`catalogs` (finished products) and `stock_catalogs` (retrosynthesis terminal materials) are arrays of descriptors, for example:

```json
{"path":"/absolute/vendor.csv","version":"vendor-2026-10-21","source":"supplier-authorized-export","license":"contract/reference","trusted":true}
```

SDF properties, CSV and headered TSV/SMILES are supported. Each valid supply record requires:

- `smiles`, `vendor`, `sku`, `supply_type` (`in_stock` or `made_to_order`).
- `catalog_date`, `checked_at` (ISO dates), with optional `valid_until`.
- `evidence_kind`: `stock_confirmation`, `dated_quote` or `orderable_confirmation`.
- `evidence_file`: a catalog-relative or absolute file path; `evidence_sha256`: SHA256 of the actual quote/stock confirmation file.

Catalog and check dates obey `evidence_max_age_days` (default 30). The system checks file existence, content hash, structure and dates; it does not replace human verification of supplier-document authenticity. A product-page URL, CAS, PubChem entry, similar molecule or old catalog alone does not qualify. A shared scaffold, opposite enantiomer, different isotope or different salt form is not an exact match. Optional `parent_smiles` records salt/parent mapping without granting purchase qualification for the target salt.

Unconfirmed and expired supply statuses are separate; in-stock and made-to-order evidence have different labels. Reconfirm before actual procurement.

## Local retrosynthesis and route completeness

Purchasing takes priority. AiZynthFinder subprocess search starts only without valid purchase evidence and actually loads `model_path`, `templates_path` and `stock_path`. Each asset requires a version and independent license field. Only local template-based models and local inchiset stock are allowed; arbitrary remote model configurations are rejected. `external_interfaces` must currently be empty: this version has no third-party interface client and does not upload structures.

Stock assets support actual search; `stock_catalogs` independently verify exact supply evidence for every terminal material. AiZynthFinder `in_stock=true` alone is insufficient.

Route checks cover exact root identity, complete molecule/reaction trees, absence of cycles, per-step template provenance, exact precursor sets from RDKit replay of retro reactions, and supply records for all terminal materials. Retain the full tree, step count, maximum depth, templates, metadata, software versions, raw search and asset SHA256. Potential stereochemistry issues require conservative review; backend-declared unresolved selectivity or missing steps do not automatically pass. SA/SCScore do not determine qualification.

`route_reviews` may include `route_sha256`, `target_smiles`, `reviewer`, `reason`, `reviewed_at` and `identity_stereo_selectivity_confirmed:true`. Human review cannot invent missing steps or material evidence. `search_artifact` can link an existing search result so a purchasable target can also display existing synthesis-route labels without forcing a new search.

Pass labels are "Purchasable | stock evidence", "Purchasable | made-to-order evidence", "Synthesizable | computed route, chemistry confirmation pending" and "Synthesizable | route reviewed by a person". They may coexist. Computed routes may enter docking but do not replace final acquisition review before wet-lab work or commissioned synthesis.

`model_missing`, `templates_missing`, `stock_missing`, `dependency_missing`, `search_timeout`, `no_route_found` and `route_needs_review` are recorded separately. No route found does not establish synthetic infeasibility. The process timeout includes search time plus model-startup grace. Concurrency is bounded by workers and batches by batch_size; each result is written independently. Cross-run cache keys include exact structure, rule hash, model/template/stock hashes and versions, catalog versions and hashes, search parameters and worker-code hash. File locks and atomic writes allow reuse only of complete search results after interruption; timeouts/errors are not cached as success. Supply validity is rechecked on route reuse, so expired evidence does not grant permanent passage.

## Docking gates, identities and export

`vina_official.run_vina_job` first calls `require_ligand`. The formal workflow, scheduler, new jobs, cache recovery, remote batch CLI and redocking use the same validation. Historical PDBQT-only Python entry points without evidence refuse execution. Legacy `sbvs` CLI/API requires `qualification_manifest` pointing to the new `availability/qualified.json`; every enumerated exact chemical state requires its own evidence.

Validation does not just read `decision=true`: it rechecks structural identity, original SDF hash, screening-record hash, current rules and unresolved blockers, catalog/stock versions, valid supply or complete routes, raw search and model assets. Changed structures, tampered evidence or expired evidence are rejected. Historical results do not acquire new qualification labels.

Each qualified record includes compound_id, state_id, source core, generation batch, original generated ID and evidence hashes. Existing zvss compound/state ID formulas and pose_id generation remain in use. Evidence flows into docking jobs, pose sidecars, the screening database, candidates.csv/json, `predocking_handoff.json` and final `predocking_evidence.json`/reports. This pre-docking integration does not change MD source or numerical parameters.

Main stage artifacts:

- `risk/raw_input.json`, `records.json`, `rules.snapshot.json`, `inputs/*.sdf`, `passed/excluded/review.sdf`, highlighted molecule SVGs.
- `availability/records.json`, `qualified.json`, `identity_evidence.json`, `evidence/*.json`, `ligands/*.sdf` and `.qualification.json` sidecars.
- `source_evidence/` retains supplier catalogs, supply files and raw search snapshots; stages also have CSV and count summaries.

Export reproduction records includes these stage files and hash manifests. Lightweight export retains the existing 20 MB cutoff, listing omissions in EXPORT_MANIFEST. Choose full export for complete catalogs/raw evidence. Model weights are excluded, but hashes and versions remain. Moving hosts requires reconfiguring local asset paths and rerunning preflight; editing a pass boolean cannot reuse old qualification. Continue with the actual qualified count when short; zero qualified molecules prevents Vina execution. Additional generation is not automatic.

## Acceptance scope

Tests are in `02_docking/tests/preflight/test_preflight.py`. Suppliers, quotes and routes are explicitly simulated unit fixtures, while RDKit matching, template replay, file hashes and gate decisions execute. Actual model searches, real stock/quotes and production-environment acceptance require the corresponding assets. No formal screening batch or MD was started by the original integration.
