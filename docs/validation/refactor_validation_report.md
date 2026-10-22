# Architecture refactor validation report

Date: 2026-10-22. Baseline: `e4906f2c76c711522654516b2ad15111707b1484` on main, no tags, clean cloned checkout. Changes are on `refactor/architecture-governance`; no main merge or Release. The baseline file hashes and 391-source AST audit were captured before code edits.

## Environment and scope

Existing macOS ARM64 Python 3.11.17 test environment reused. Installed versions observed: RDKit 2026.9.1, OpenMM 8.6.1, ParmEd 4.3.1, NumPy 2.4.6, MDAnalysis 2.10.0, Meeko 0.8.0 and PyYAML 6.0.2. No scientific dependencies were upgraded. Frontend used Node 24.19.0 and the existing frozen pnpm lock. No production GPU environment, model weights or GROMACS executable was available for full scientific acceptance.

## Actual Python test execution

`$PYTHON` below denotes that existing Python 3.11.17 interpreter. Each command ran in its listed directory, separately to avoid conflicting module names. Report paths below are normalized; original JUnit hashes and durations are in [test_results.json](test_results.json).

| Suite / directory | Command | JUnit records | Passed | Failed/errors | Skipped |
|---|---|---:|---:|---:|---:|
| generation / `01_molecule_generation` | `PYTHONPATH=. $PYTHON -m pytest tests -q -ra --junitxml=report.xml` | 63 | 63 | 0/0 | 0 |
| workbench / `00_workbench` | `PYTHONPATH=. $PYTHON -m pytest tests -q -ra --junitxml=report.xml` | 128 | 127 | 0/0 | 1 |
| screening / `03_virtual_screening` | `PYTHONPATH=src $PYTHON -m pytest tests -q -ra --junitxml=report.xml -m 'not integration'` | 71 | 71 | 0/0 | 0 |
| docking / `02_docking` | `PYTHONPATH=src $PYTHON -m pytest tests -q -ra --junitxml=report.xml` | 81 | 79 | 0/0 | 2 |
| interaction / `04_interaction_analysis` | `PYTHONPATH=src $PYTHON -m pytest tests -q -ra --junitxml=report.xml` | 40 | 40 | 0/0 | 0 |
| md / `05_molecular_dynamics` | `PYTHONPATH=. $PYTHON -m pytest tests -q -ra --junitxml=report.xml` | 80 | 79 | 0/0 | 1 |
| package / `deployment` | `PYTHONPATH=. $PYTHON -m pytest tests -q -ra --junitxml=report.xml` | 9 | 9 | 0/0 | 0 |

Total: **468 passed, 0 failed, 0 errors, 4 skipped JUnit records**. There are 471 collected test items: 468 passed plus 3 item skips. The remaining skipped record is a module-collection skip, not an executed test. No repeat run is added to this total. The final workbench run includes all 27 new governance tests.

Skipped coverage:

- Docking `tests/v2/test_prepare_dock.py`: missing `pdbfixer`; all seven source test functions are excluded before collection, including the real Vina path. This is a meaningful untested area.
- Docking molecule preparation: `molscrub` not installed.
- Native MD Ubuntu acceptance: opt-in environment unavailable.
- Workbench legacy generated-results fixture: historical generated data is not distributed.

Collection audit also identified standalone `*_acceptance.py` scripts, `webui/test_api.py`, browser scripts and external integrations that pytest does not run by filename or environment. They are not counted as passed. CI retains module pytest suites; unittest-only discovery is not used. Mocked command/failure tests remain contract evidence only.

Baseline run: generation 63 pass; screening 71; docking 79 plus 2 skipped records; interaction 40; MD 79 plus 1 skip; packaging 9. Baseline workbench had 99 pass, 1 skip and one localhost-binding sandbox failure. The final workbench suite ran with localhost access and passed 127 tests plus 1 skip. During development, new target-file validation exposed incomplete temporary fixtures (9 failures and 17 setup errors); fixtures now copy the real target YAML and the full suite passed afterward. No scientific assertion or validation threshold was weakened.

## Other executed checks

- `python deployment/check_english.py`: pass, no content-policy errors.
- `python deployment/check_repository.py`: 878 files checked, 0 errors; repository, links, packaging and version integrity passed.
- `python deployment/check_governance.py`: 6 modules, 13 backends, all source-symbol/protocol/schema references valid.
- `python deployment/check_scientific_baseline.py`: **127 original scientific/configuration/data/history files byte-identical**.
- `python deployment/sync_platform_version.py --check`: no stale generated version metadata.
- `git diff --check`: pass.
- `node --experimental-vm-modules 00_workbench/webui/check_js.mjs`: 19 JavaScript modules passed.
- `node 00_workbench/webui/test_pages.mjs`: 15 fixture/DOM checks passed; this is not a live browser/scientific test.
- In `06_websites/01_public_demo/yuanhe-public`: `node node_modules/typescript/bin/tsc`, `node node_modules/vite/bin/vite.js build`, `node --test tests/*.test.mjs`: typecheck/build passed; 3 source-contract tests passed. Vite warns about the existing large NGL chunk. Frozen dependency install preserved the lockfile. The pnpm build wrapper initially attempted a mismatched local-store reinstall; the exact TypeScript/Vite package-script executables succeeded directly.
- `python -m pip wheel --no-deps --no-build-isolation --no-cache-dir` for docking and ZVSS: both wheels built as `0.11.0.dev1` (PEP 440 normalization). Independent export/package acceptance: 9/9 passed.

## Real lightweight execution and UI

A new real workbench run `YH-20261022-fae062126f964a9589a9ebdf0c8b451b` executed only `prepare` on the repository's 9NU4 input. It produced 3,785 receptor atoms and 258 pocket atoms. Both stage contracts validated, resolved JSON/YAML matched, artifact hashes were recorded, and all seven downstream stages remained `created`. The run correctly records the baseline Git SHA plus `git_dirty=true` because it was executed before the refactor commits; it has not been relabelled.

The actual localhost HTTP service was opened in the in-app browser. Product/subtitle, dark workbench, registry-backed availability, System & Versions and the actual run's stored version/configuration were inspected. Missing dependencies appeared as Not Installed; the interface did not present local scientific backends as validated. Old `/api/project.modules` retains its five-entry payload and new `registered_modules` adds the six functional identities. No public service was deployed and no GPU endpoint was exposed.

## Scientific invariance

Original Peptide2Mol/REINVENT4 adapters, chemistry identity/preparation, official Vina implementation, ZVSS/Y-MAVS, ZLIS and native OpenMM/RESP/trajectory implementations remain unchanged in protected paths. Original research budgets and all native default files are unchanged. Facades re-export the same implementations. Added changes concern provenance, validation, dispatch consistency and version identity.

The original membrane MDP function is text-identical and regression-tested across its phases and temperature inputs. Formal membrane production remains **one 500 ns trajectory per candidate**, 2 fs, **250000000 steps**, 310 K, 1 bar semiisotropic. Soluble defaults remain independently 300 K and 1 atm. No force-field, charge, lipid, integration or analysis default was changed. POPE is the actual default membrane lipid; no unverified Lipid17 version is invented.

## Not validated / remaining work

No real new Peptide2Mol/REINVENT4 inference, real Vina docking, RESP computation, OpenMM-to-GROMACS topology conversion, native checkpoint continuity test, complete trajectory analysis or 500 ns production was performed. Local OpenMM import and contract tests do not establish CUDA, target-host, topology or scientific correctness. A target-host evidence run requires the existing authorized environments and actual inputs, and heavy compute remains explicitly excluded from this task.

The integrated workbench and native soluble task provenance are implemented. Independent historical tool schemas/schedulers remain preserved; universal per-field configuration adapters for every standalone tool, export-only historical-format migration, exhaustive backend-specific error classes and per-capability target-host validation are not complete. This report does not claim the entire broader platform acceptance checklist is satisfied. See [migration](../migration/architecture_refactor.md) for exact implemented scope and rollback.
