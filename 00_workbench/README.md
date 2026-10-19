# Unified research workbench

Start from the repository root with `python3 00_workbench/run_webui.py --no-browser`, or use the platform launcher.

`webui/workflow.py` manages the pipeline; `pipeline_bridge.py` invokes real stages. `webui/studio.py` supplies guided preparation and review tools; `scientific/` implements membrane construction and GROMACS analysis. `webui/static/index.html` is the main page.

Historical pages and `pepwedge` preserve earlier analysis interfaces. Their old output datasets are not bundled and are not model-generated candidates. Follow the root README, `USER_GUIDE.txt` and `docs/01_SETUP.md`. `README.historical.md` records the earlier prototype, not acceptance of the current release.


## Yuanhe Bio integration

See the [project code map](../docs/product/CODE_AND_CONTENT_MAP.md) for current module boundaries and the [product value assessment](../docs/product/PRODUCT_VALUE_ZH.md) for intended users. Existing technical instructions and historical acceptance records above retain their original scope; this documentation refresh does not revalidate execution.
