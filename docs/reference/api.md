# API compatibility

New read-only governance endpoints: `GET /api/v1/versions`, `GET /api/v1/modules`, `GET /api/v1/protocols`. `POST /api/v1/validate` validates a prospective workbench configuration without creating a run. Existing `/api/project`, `/api/project/runs`, task/stage/review and studio endpoints remain intact.

Workbenches return registry-backed module status and platform identities. Task details return stored `manifest` and `resolved_configuration`; historical runs lacking these fields return null, not the current platform version. No endpoint or database was deleted or migrated.

Errors retain the legacy `error` string and add `error_code`, `message`, `run_id`, `stage`, `details`, `recoverable`, `suggested_action`. The HTTP layer separates input/client failures from runtime/server failures; engine-specific original logs retain their detailed chemical/parameterization diagnosis. Further normalization of every historical exception class is not implied.

Version 1 governance endpoints are additive. Breaking response changes require a new API namespace; existing aliases need a documented compatibility period. Host, Origin and body-size protection remains enabled and is tested through a real localhost server.

The public demo consumes precomputed public assets only. No API route connects it to the private compute service.

`GET /api/project` retains the legacy five-entry `modules` shape (`name`, `path`, `present`). The six registry modules are additive under `registered_modules`; old clients are not required to adopt a new payload shape. Stage contracts record the actual frozen configuration hash and input/output artifact identities.
