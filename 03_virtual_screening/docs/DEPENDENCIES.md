# Dependencies

The following versions were recorded in the historical Windows environment, not solved or certified for the current WSL release:

| Package | Recorded version | Role |
|---|---|---|
| RDKit | 2026.9.1 | Chemistry and conformers |
| NumPy | 2.5.4 | Geometry/statistics |
| SciPy | 1.18.1 | Numerical helpers |
| pandas | 3.0.6 | Tables |
| gemmi | 0.7.5 | Structure parsing |
| MDAnalysis | 2.10.0 | Structure/trajectory support |
| ProLIF | 2.2.3 | Optional interaction analysis |
| PoseBusters | 0.6.5 | Optional pose checks |
| FastAPI | 0.143.3 | Read-only web API |
| Uvicorn | 0.54.0 | ASGI server |
| Pydantic | 2.14.1 | Validation |
| PyYAML | 6.0.3 | Configuration |
| Jinja2 | 3.1.6 | Templates |
| Matplotlib | 3.11.2 | Optional plots |
| requests | 2.34.2 | Reference retrieval support |
| psutil | 7.2.2 | Host probes |
| pytest | 9.1.1 | Historical test runner |

Use `pyproject.toml` and the target installation evidence for current dependencies. Earlier Windows PEP 517 temporary-directory failures required `--no-build-isolation`; do not disable isolation without a matching diagnosed cause.

`zvss.compat` finds a writable ASCII scratch directory and mirrors RDKit data there when required by native path handling. If unavailable, relevant optional checks remain unsupported.

Dependencies retain their own licenses: RDKit BSD-3-Clause, gemmi MPL-2.0, MDAnalysis GPL-2.0-or-later, ProLIF Apache-2.0, PoseBusters/FastAPI/PyYAML MIT and Uvicorn BSD-3-Clause, subject to actual installed versions. Importing a GPL component can create obligations; optional or local execution is not a blanket exemption. Review combined-distribution obligations and the root license scope. GNINA is an external executable restored separately.
