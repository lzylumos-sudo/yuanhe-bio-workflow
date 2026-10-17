# Standalone docking user guide

For the unified workflow use the repository-root launcher and guide. The standalone classroom interface uses this module's `START_WINDOWS.cmd` and normally port 5902; it requires its separately restored runtime.

Enter CCO or multiline `SMILES Name`, then start screening. Missing names receive internal IDs; duplicate names do not overwrite files. Inspect actual status, success/failure counts and parent/pose rankings. Toggle box and residue labels, save at least two distinct actual pose images, then export CSV, manifest and HTML report. The 20/100-molecule demonstration libraries have no activity labels and are not DHFR inhibitor validation sets.

```bash
./run.sh doctor
./run.sh demo --preset class_quick
./run.sh demo --preset class_full
./run.sh serve --host 127.0.0.1 --port 5902
./run.sh prepare-receptor --input data/raw/1DRF.pdb --chains A --out data/prepared/my_review
./run.sh prepare-ligands --input data/demo/ligands_20.smi --out data/prepared/my_ligands
./run.sh make-box --reference data/reference/FOL_crystal.sdf --padding 5 --out data/prepared/my_box.json
./run.sh screen --config class_quick --ligands data/demo/cco.smi
./run.sh report --run runs/YOUR_RUN_ID
./run.sh resume --run runs/YOUR_RUN_ID
./run.sh rescore --run runs/YOUR_RUN_ID --engine gnina --operation score_only --top-k 10 --device gpu
./run.sh export-md --run runs/YOUR_RUN_ID --selection data/manual_selection.json
```

Replace placeholders with actual records. PowerShell `SBVS.ps1` forwards the same arguments. Research preparation first returns REVIEW_REQUIRED: inspect the plan, bind its input hash and confirmed policy, then supply `--policy`. The classroom policy applies only to its explicit teaching case.

Raw outputs and CSV are authoritative. Keep empirical affinity, CNN pose score and CNN affinity/pK separate. Aggregate seed representatives by their median while retaining an actual output pose; no synthetic median-coordinate pose is created.
