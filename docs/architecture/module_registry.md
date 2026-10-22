# Module and capability registries

`configs/module_registry.yaml` is parsed by `yuanhe.registry`. JSON syntax is used as the YAML 1.2 subset so the stdlib-only HTTP service needs no YAML dependency. It registers module identities, implementation paths, public symbols, backend dependencies, native environment selectors, capabilities, lifecycle status, protocol compatibility and old aliases.

Consumers: workbench inventory, Web UI System & Versions, `/api/v1/modules`, `yuanhe modules`, `yuanhe doctor`, configuration validation and CI. AST checks verify implementations and exported symbols without importing scientific packages in the UI process.

- `implemented`: source mapping exists.
- `installed`: declared interpreter/module, executable and model-file prerequisites were found on the probe host.
- `verified`: actual validation evidence for this environment and backend. Currently false unless a future evidence integration explicitly establishes it; old reports are not automatically promoted.
- `runnable`: known runtime readiness, not scientific success. Complex GPU backends conservatively require target-environment validation. Input and human review gates still apply.

`active`, `experimental`, `deprecated`, `archived` describe support lifecycle. They are independent of installation, task execution and scientific validation. New requests use the public `docking` facade; `sbvs` and `sbvs2` are retained. `ymavs` has parent `zvss`. No historical implementation was deleted or silently routed into a different algorithm.

A new backend must provide an existing implementation, precise dependency checks, input/output contracts, supported protocol IDs and tests. Do not introduce names ending in `2`, `new`, `final` or `latest`; existing identifiers remain compatibility aliases.
