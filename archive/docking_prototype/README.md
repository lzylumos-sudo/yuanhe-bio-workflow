# Archived docking prototype v2

This earlier implementation is retained for reference and is not the current workbench backend. Historical baseline records counted 2,259 files, excluding runtime images; those generated artifacts are not bundled here.

```bash
export SBVS_PYTHON="$HOME/.sbvs2_env/bin/python"
./run.sh doctor
./run.sh new-project my_project
```

source/sbvs2 contains paths, hashing, states, approvals, receptor inventory, ligand preparation/review, docking protocols/engines, ranking, project and CLI modules. source/tests contains unit tests. Projects separate originals, preparation, run records, logs, work products, results and analysis; runs are never overwritten. Configuration values are starting hypotheses requiring target calibration.

scripts/redock_fol.py creates a 1DRF/FOL test project; supplying project/run IDs resumes only after hash checks. See docs/STATUS.md for historical Vina results and explicitly untested GNINA capability.
