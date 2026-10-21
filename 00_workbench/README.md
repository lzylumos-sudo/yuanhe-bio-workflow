Current research migration: official Vina 1.2.7 CPU docking; GAFF2 + HF/RESP; one 500 ns trajectory per selected complex; core analysis and actual PB at 400–500 ns. See docs/15_VINA_RESP_RESEARCH.md at the repository root. Previous short-run defaults below are historical. This revision is code-only; validate on the intended computer.

# Unified research workbench

Users open the workbench UI URL supplied by the administrator and operate the existing interface, manually or with Codex. See [the UI workflow](../docs/18_UI_RESEARCH_WORKFLOW.md). Administrators retain the existing deployment entrypoint; no additional client launcher is needed.

`webui/workflow.py` manages the pipeline; `pipeline_bridge.py` invokes real stages. `webui/studio.py` supplies guided preparation and review tools; `scientific/` implements membrane construction and GROMACS analysis. `webui/static/index.html` is the main page.

Historical pages and `pepwedge` preserve earlier analysis interfaces. Their old output datasets are not bundled and are not model-generated candidates. Follow the root README, `USER_GUIDE.txt` and `docs/01_SETUP.md`. `README.historical.md` records the earlier prototype, not acceptance of the current release.


## Yuanhe Bio integration

See the [project code map](../docs/product/CODE_AND_CONTENT_MAP.md) for current module boundaries and the [product value assessment](../docs/product/PRODUCT_VALUE_ZH.md) for intended users. Existing technical instructions and historical acceptance records above retain their original scope; this documentation refresh does not revalidate execution.
