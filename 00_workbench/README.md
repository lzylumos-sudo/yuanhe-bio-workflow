# Unified research workbench

Start from the repository root with `python3 00_workbench/run_webui.py --no-browser`, or use the platform launcher.

`webui/workflow.py` manages the pipeline; `pipeline_bridge.py` invokes real stages. `webui/studio.py` supplies guided preparation and review tools; `scientific/` implements membrane construction and GROMACS analysis. `webui/static/index.html` is the main page.

Historical pages and `pepwedge` preserve earlier analysis interfaces. Their old output datasets are not bundled and are not model-generated candidates. Follow the root README, `USER_GUIDE.txt` and `docs/01_SETUP.md`. `README.historical.md` records the earlier prototype, not acceptance of the current release.
