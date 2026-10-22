# Configuration resolution

The integrated workbench uses `webui.workflow.validate_config` for API, UI and CLI. It invokes existing research, pre-docking and docking validators, plus registry/protocol integrity checks. The original default JSON is unchanged.

Research defaults -> selected research profile -> supplied task values -> allowed stage reviews. `yuanhe validate/run` additionally accepts `--target-config`, `--config`, `--override` in that order. Nested supplied layers merge with source/override records; existing stage reviews retain their original shallow semantics. Unknown top-level keys and malformed research values fail validation. Existing unit/range and chemistry checks remain in native validators.

`resolved_config.yaml` records Default, User Override, Resolved Config and protocol conformance. It uses valid YAML 1.2 JSON syntax. Existing `resolved_config.json` remains available. At stage dispatch, final review values are frozen in both formats. The actual worker checks the hash, agreement of formats and agreement with frozen run/review records, then executes that exact configuration. Each stage records its executing Python environment separately from the UI host.

The original scientific budget remains Peptide2Mol 5,000 attempts; 200 cores; REINVENT4 200,000 attempts; generated/external caps 50,000 each; initial 100,000 parents; refine 2,000; review 50; MD up to 8 candidates; one 500 ns trajectory each; concurrency 2. Initial Vina: 1 receptor x 1 state x 1 seed, exhaustiveness 8, num_modes 5. Refine: maxima 2 x 2 x 2, exhaustiveness 32, num_modes 10. These are budgets, not guarantees of unique valid molecules or completed jobs.

Scientific changes are not silently labelled the original protocol: registry comparison returns `custom@unknown`, the base protocol and differences. Operational resource settings may differ without inventing a new scientific recipe. Reviewed structure, pose, model hashes and system-specific MD parameters are recorded separately.

Units remain native and explicit: docking Angstrom / kcal/mol; MD ns, fs, K; soluble pressure atm; membrane pressure bar. `integration_steps(500,2)` is exactly 250000000 and rejects non-integral steps and non-finite values. No unit stripping, default force-field substitution or automatic duration reduction occurs.

Standalone scientific loaders remain authoritative for their engine schemas. This is not yet a universal schema replacement for every historical CLI; see the compatibility matrix.
