# Official AutoDock Vina docking module

Current production engine: Vina 1.2.7, CPU, exhaustiveness 32, explicit seeds/box/thread budget. CNN fields are N/A. Historical GNINA results remain archived. See [migration and commands](../docs/15_VINA_RESP_RESEARCH.md). This source revision has no retained runtime acceptance results.

This module contains the main docking implementation, preparation/review utilities, protocol records and standalone classroom/research interface. The unified project starts from the root `start_windows.cmd`; standalone module entry points and old WSL runtime scripts have their own environment requirements.

Read `docs/USER_GUIDE.md`, `docs/RESEARCH_GUIDE.md`, `docs/TROUBLESHOOTING.md` and `docs/MIGRATION.md`. Model binaries, installed environments, VHDX images and generated runs are excluded from GitHub. Restore and validate them separately.

Historical acceptance records in `PROGRESS.md` and `docs/ACCEPTANCE_REPORT.md` describe specific earlier hosts and inputs. They do not establish current MurJ, RTX 5070 or second-machine acceptance. Scores, poses, failures and provenance must come from actual execution.

Portable tests: `PYTHONPATH=src python -m pytest tests -q` from this module. The legacy runtime stop-guard check is an explicit tool, `python scripts/validate_stop_guard.py`; it requires the dedicated `/opt/sbvs` runtime and `SBVS_PROJECT_ROOT`.

## Yuanhe Bio integration

See the [project code map](../docs/product/CODE_AND_CONTENT_MAP.md) for current module boundaries and the [product value assessment](../docs/product/PRODUCT_VALUE_ZH.md) for intended users. Existing technical instructions and historical acceptance records above retain their original scope; this documentation refresh does not revalidate execution.
