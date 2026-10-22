# Testing and evidence

Existing module tests use pytest, including unittest classes. Do not use unittest-only discovery as a substitute for collection. CI retains all existing module suites and adds registry/version/protocol/default gates. Automatic CI does not launch 500 ns production.

```bash
python deployment/check_english.py
python deployment/check_repository.py
python deployment/check_governance.py
python deployment/check_scientific_baseline.py
python deployment/sync_platform_version.py --check
# Run from each numbered module; src path for DOCK, SCREEN and INTERACT:
PYTHONPATH=. python -m pytest tests -q -ra --junitxml=report.xml
# SCREEN excludes explicitly marked external integrations:
PYTHONPATH=src python -m pytest tests -m 'not integration' -q -ra
node --experimental-vm-modules 00_workbench/webui/check_js.mjs
```

Unit/contract tests include simulated dependency/error inputs. They do not constitute real inference/docking/MD validation. The lightweight end-to-end acceptance creates a real run and executes structure preparation only. Browser checks examine the actual service. GPU/native checks remain separately opted in and must identify hardware, environment, inputs and actual duration. Every skipped test has an explicit reason in `-ra` output/JUnit.

Historical `webui/test_api.py` expects undistributed generated artifacts; it is not a fresh-checkout acceptance suite. `webui/test_pages.mjs` uses recorded payload fixtures and a DOM shim, not a live browser. Both boundaries are retained in the validation report rather than hidden behind aggregate passes.
