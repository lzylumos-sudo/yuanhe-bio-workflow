# Compatibility matrix

| Existing surface | After refactor | Evidence / limitation |
|---|---|---|
| Root start scripts / run_webui.py | unchanged launch | real localhost page inspected |
| Existing `/api/*` | retained | existing Host/Origin/body guard test passed |
| `run_research.py` | same workflow implementation | existing regression suite |
| New `yuanhe` CLI | repository-native facade | version/registry/validate/create/status checks |
| `sbvs` CLI | retained own workflow and schema | existing docking suite and independent package |
| `sbvs2` CLI/imports | retained | new `docking` facade re-exports original functions |
| Peptide2Mol / REINVENT4 | independent adapters behind lazy registry | original adapters byte-identical |
| ZVSS / Y-MAVS | same implementation, parent-child registration | screening suite; no algorithm changes |
| ZLIS | unchanged | 40 existing tests; original source byte-identical |
| Standalone OpenMM | same scientific core; additive new-run metadata | native contracts; no GPU acceptance |
| Membrane portable package | governance included alongside original core | package checksum and help checks |
| Historical research records | read unchanged; absent metadata stays unknown | no migration/backfill |
| Old checkpoints / running jobs | no deployment into their environment | no target-host resume attempted |
| Public demo | same static architecture, generated version footer | build and static-source tests; not redeployed |

Deprecated aliases have no scheduled removal in this release. Independent docking/screening native manifests, schedulers and configuration schemas are intentionally preserved; the integrated workbench supplies the unified parent-run contract. Converting every standalone store into that contract is additional migration work, not achieved by registering a name.
