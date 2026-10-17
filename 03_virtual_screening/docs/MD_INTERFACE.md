# Downstream handoff contract

Each `runs/<run_id>/md_handoff/<role>_<order>_<compound_id>/` contains `complex.pdb`, `ligand.sdf`, `manifest.json`, `task.json`, `README.md`, `adapter04.json`, `adapter05.json` and `checksums.sha256`.

The ligand is HETATM/LIG, chain L, residue 1. REMARK 900-903 records run, compound, state and pose identity. The manifest maps PDB atom names to one-based SDF atom indices and records receptor identity, coordinates, chemical identity and checksums. Downstream tools must not guess atom order.

```json
{"format":"zvss-md-manifest-1","compound_id":"CMP...","state_id":"STA...","pose_id":"POS...","ligand":{"resname":"LIG","chain":"L","resid":"1","atoms":9,"smiles":"...","atom_mapping":{"C1":1}},"receptor":{"path":"receptor.pdb","version":"REC...","atoms":4598},"checksums":{"complex.pdb":"sha256:...","ligand.sdf":"sha256:..."}}
```

The atom count above illustrates the schema; it is not a fixed target for all preparations.

`adapter04` targets standalone `05_molecular_dynamics` using `bash run_md.sh --input <complex.pdb>` with adjacent SDF/manifest chemistry. It verifies file equivalence and records command templates. The standalone runner is Linux-only and lacks a membrane protocol, so relevant bundles retain `blocked_platform_linux_required` and `blocked_missing_membrane_capability`. `verified_offline: true` means contract checks only, not MD execution.

`adapter05` targets `04_interaction_analysis` with `bash run.sh analyze <pdb> --ligand L:LIG1 --template ligand.sdf`. It records expected interaction tables, diagrams and scene artifacts, but retains `executed: false` until execution actually occurs.

Formal handoff gates must remain enforced. The newer unified workbench membrane wizard is a separate construction/execution path in `00_workbench/scientific/`, with explicit exploratory approval for preliminary candidates. Its availability does not retroactively mark historical standalone bundles executed or scientifically qualified.
