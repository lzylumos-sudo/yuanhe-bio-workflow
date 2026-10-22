# Data contracts and artifacts

New workbench runs write `run_manifest.json`, `artifact_registry.json`, `resolved_config.yaml`, validated per-stage `stage_input.json` / `stage_output.json`, execution environment and structured `events.jsonl`. Native soluble MD tasks in a full checkout also create protocol/software manifests. `yuanhe.contracts` enforces required identity, full Git SHA, hash and separated status fields before manifest writes. JSON schema declarations are under `configs/schemas`.

Run identity, software identity, scientific protocol identity and model/data identity are separate. Each artifact has an ID, type, run/stage, source, content hash, schema version, validation status and relative storage path. Recorded file hashes mean integrity provenance, not biological or numerical validation. Existing per-stage manifests retain detailed schemas and file hashes.

Parent compound -> chemical state -> conformer -> docking pose -> receptor conformation -> simulation system -> trajectory remain separate identities in original adapters, research cache inputs and system records. No canonicalization rule or atom mapping was changed. Vina job cache identity retains receptor/state hashes, seed, box, exhaustiveness, num_modes, actual executable hash/version and adapter hash. Missing GNINA/CNN fields are not replaced with zero.

Execution: created, validated, queued, running, completed, failed, cancelled, blocked. Scientific review: not_reviewed, review_required, accepted, rejected. A completed pipeline requests review; it does not self-approve. Existing native states are preserved and mapped only for new governance views. Historical manifests are never relabelled.

No historical migration is performed. A future export/migration must preserve original files and include old/new paths and hashes, migration version/time. Automatic reinterpretation of old results under the new protocols is forbidden.
