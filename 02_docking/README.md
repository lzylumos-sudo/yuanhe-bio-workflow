# GNINA/SBVS docking module

This module contains the main docking implementation, preparation/review utilities, protocol records and standalone classroom/research interface. The unified project starts from the root `start_windows.cmd`; standalone module entry points and old WSL runtime scripts have their own environment requirements.

Read `docs/USER_GUIDE.md`, `docs/RESEARCH_GUIDE.md`, `docs/TROUBLESHOOTING.md` and `docs/MIGRATION.md`. Model binaries, installed environments, VHDX images and generated runs are excluded from GitHub. Restore and validate them separately.

Historical acceptance records in `PROGRESS.md` and `docs/ACCEPTANCE_REPORT.md` describe specific earlier hosts and inputs. They do not establish current MurJ, RTX 5070 or second-machine acceptance. Scores, poses, failures and provenance must come from actual execution.

Portable tests: `PYTHONPATH=src python -m pytest tests -q` from this module. The legacy runtime stop-guard check is an explicit tool, `python scripts/validate_stop_guard.py`; it requires the dedicated `/opt/sbvs` runtime and `SBVS_PROJECT_ROOT`.
