# ZVSS: Zhiyu Virtual Screening Studio

ZVSS consumes existing GNINA batch outputs and produces traceable candidates, exclusions and downstream task bundles. It connects to standalone OpenMM and ZLIS through explicit file contracts. The unified workbench additionally provides a separate membrane/GROMACS path.

Install the package into the chemistry environment with `python -m pip install -e ./03_virtual_screening` from the repository root. Run from this module directory:

```bash
zvss --root . doctor --with-tools
zvss --root . run --input /path/to/gnina/results --target configs/targets/ecmurj.yaml --budget configs/budgets/quick.yaml
zvss --root . serve --port 8770
zvss --root . report --run-id YOUR_RUN_ID
zvss --root . runs
```

`--root` is a global option and precedes the subcommand. The standalone web view is read-only; daily integrated execution is in the main workbench.

`src/zvss/` contains ingestion, receptor audit, quality control, mechanism, motion, redocking, scoring, strain, selection, membrane proxy, handoff and reporting. Target definitions are in `configs/targets/`; budgets in `configs/budgets/`. Each run stores resolved configuration, logs, stage summaries and reports. Runtime databases and generated outputs are not committed.

Historical generic-target runs t4/t5 imported 60 poses from 20 compounds and produced identical selection/ranking files. They do not validate MurJ docking. The later qualification audit distinguishes qualified primary, preliminary and exploratory candidates; historical primary labels must not be interpreted as passing current strict qualification.

Missing dependencies produce explicit blocked/unsupported states. Use an ASCII scratch directory for RDKit/PoseBusters where necessary. For the exact 9NU4 construct, screening can transfer the archived PPM3 orientation after hash, mapping and rigid-fit checks. Other inputs retain proxy or missing-evidence status. This screening orientation check does not build a membrane or establish experimental orientation. MurJ pose-level conclusions require real target poses. See the architecture, interface and historical acceptance records for details.
