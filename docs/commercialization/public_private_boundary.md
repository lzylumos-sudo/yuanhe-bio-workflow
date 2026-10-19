# Research and Enterprise boundary proposal

Planning only, 2026-10-19. No repository creation, migration, deletion, visibility change or license change is authorized by this document. Public history and earlier grants remain effective according to their terms.

## Module decisions

| Asset | Recommendation | Reason and boundary |
|---|---|---|
| README, setup, scientific limitations, public benchmark protocols | Public documentation | Enables evaluation without disclosing unpublished research |
| Existing 00_workbench orchestration and reviewed local UI | Remain source-available under existing applicable terms | Reuse current behavior; do not duplicate into an enterprise fork |
| 01 generation adapters, identity/core validation | Remain source-available, with upstream exclusions | Interface work is distinct from a new proprietary generator |
| 02 docking wrappers and public controls | Remain source-available / retain third-party terms | Transparent execution and benchmark inputs; no blanket claim over engines |
| 03 screening | Further rights review; retain public files | Earlier MIT scope must be mapped; existing rules cannot regain secrecy |
| 04 interaction analysis and 05 standalone MD | Remain source-available subject to component terms | Reusable scientific tooling; optional engines retain separate environments |
| deployment packaging and shared data contracts | Public interfaces/source-available implementation | Publish compatibility contracts; review each binary/runtime separately |
| Public structural inputs and approved regression fixtures | Public with provenance | No customer structures or unpublished candidates mixed into fixtures |
| New original generator/scoring/AL acquisition implementation | Develop privately after explicit repository approval | Ownership, publication/patent decision and evidence before disclosure; no implementation included here |
| Enterprise identity, tenant isolation, billing and commercial API service | Future private implementation | Public API schema may be documented independently after review |
| Proprietary weights/training sets and prospective optimization results | Private by default after authorization | Separate access controls and provenance; no automatic upload from runs |
| Customer molecules, targets, assays, credentials | Private, customer-governed | Contracted use, retention/deletion and no cross-customer training without consent |
| Already disclosed results, candidate IDs and source history | Review, preserve current state | Cannot promise restoration of trade-secret protection |

## Minimal architecture

Research UI/CLI -> existing pipeline adapters -> versioned execution request -> local or Linux GPU executor -> immutable artifacts and evidence manifest. A future Enterprise control plane can call the same executor contract, adding authenticated projects, quotas, audit records and deployment policy. A future private model service implements the same generator interface without being imported into the public distribution. This is a proposed seam, not implemented tenancy or access control.

Reuse `webui/workflow.py` stage records, source snapshots and review gates; `scientific/artifacts.py` hashing/path checks; `scientific/md_cli.py` independent entrypoint; `deployment/package_component.py` modular exports; current model adapters. Keep current CLI paths and file formats until contract fixtures prove equivalence. Avoid a new microservice for every scientific function.

## Proposed contracts

An execution request carries schema_version, run_id, stage_id, target/input hashes, molecule/pose identity, coordinate-frame identity, protocol and environment digests, model/version/weight hash, seed, budget, expected outputs and explicit review approvals. Outputs carry artifact hashes, units, status, rejection reasons, timings and source revision. A missing result is unknown/failed, never a numeric zero or a pass. Keep candidate, pose, run, training-set and model identifiers separate.

Executor operations: validate, plan, execute, status, cancel, resume, collect. Generator operations: capabilities, propose and provenance; scoring backends report value, unit, uncertainty, applicability and evidence class. API version negotiation must fail on unsupported contracts. These operation names describe a proposed interface; no new service exists yet.

Retain a local process runner first. Introduce a durable job ledger only when queued work is needed; ensure idempotent submission, leases/heartbeats, bounded retry and cancellation propagation. A retry may repeat computation but cannot overwrite scientific evidence. Resume validates system/input/engine compatibility and last complete checkpoint. Add a queue broker only after measured concurrency demand. Separate tenant authorization from the job schema and from object-store paths; test cross-tenant denial, signed artifact access and audit integrity before any external API.

## Deployment paths and provenance

| Host | Near-term path | Acceptance gate |
|---|---|---|
| Linux workstation / GPU server | Existing environments and local runner; remote UI through authenticated tunnel | Repeat known case, interruption/recovery, exported hashes |
| Windows | Browser plus reviewed WSL2 Linux executor; optional remote GPU | Actual WSL/driver/target GPU run; no native parity claim |
| macOS | Browser, preview and supported CPU checks; remote Linux inference | Client-to-worker contract and artifact roundtrip, not CUDA equivalence |
| Enterprise private deployment | Later control plane with versioned shared executor packages | Rights-cleared runtime, tenant/auth tests, backup/restore and customer acceptance |

One release line owns shared scientific code and contracts; private extensions depend on pinned released interfaces. Do not maintain copied pipelines. Package runtime dependencies separately according to their licenses; process separation does not resolve every copyleft question. No customer telemetry or automatic training-data contribution is implied.

## Incremental implementation gates

1. Freeze current acceptance fixtures and document existing formats. No extraction until clean-run and negative cases are reproducible.
2. Propose adapters around existing entrypoints; show equivalent outputs/hashes and no lost review gates on public fixtures. Obtain approval before refactoring.
3. Add a single remote-worker transport and durable job records only after a local compatibility suite passes.
4. Build private optimization and enterprise layers only after owner approval, rights clearance and a concrete customer need.

Future algorithm formulas, training recipes, customer structures and undisclosed candidates require a separate disclosure decision. This public plan contains generic interfaces and priorities, not new proprietary implementation details.
