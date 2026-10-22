# Module inventory

Baseline: `e4906f2c76c711522654516b2ad15111707b1484`, main; clean checkout, no tags. Refactor branch: `refactor/architecture-governance`. All 824 tracked files hashed; all 391 Python sources AST-parsed before implementation. The machine-readable audit records file hashes, top-level definitions, imports, subprocess and environment accesses. AST edges are an inventory, not proof of runtime reachability.

| ID | Path | Responsibility | Callers / entry | Input | Output | Parameters | Status |
|---|---|---|---|---|---|---|---|
| GEN | 01_molecule_generation | Pocket-conditioned Peptide2Mol; REINVENT4 LibInvent prior scaffold decoration | pipeline_bridge.stage_peptide2mol / stage_reinvent4 | PDB pocket, checkpoint, scaffold, seed, attempts | SDF, candidates, failures, model manifests | configs/project.json; workbench research config | active; model execution requires external environments |
| DOCK | 02_docking | Official Vina preparation, docking and pose export; historical GNINA readers | research.docking.run; sbvs.runner; sbvs2.docking.vina_official | chemical states, reviewed receptor/box, seed, protocol | mapped poses, engine-specific scores, job hashes | research config; sbvs YAML; sbvs2 Protocol | active Vina; deprecated GNINA execution |
| SCREEN | 03_virtual_screening | ZVSS evidence screening with Y-MAVS membrane sub-analysis | pipeline_bridge.stage_screening; research.screening.run; zvss.cli | docking poses, target profile, evidence, budget | decisions, ranking, handoff and resolved YAML | configs/stages.yaml, targets/*.yaml, budgets/*.yaml | active; membrane proxy is not activity evidence |
| INTERACT | 04_interaction_analysis | ZLIS static interaction analysis and visualization | zlis.cli; zlis.web; zlis.pipeline | PDB and ligand definitions | interaction tables, images, reports | configs; CLI arguments | active independent tool; no claim all runs invoke it |
| MD | 05_molecular_dynamics; 00_workbench/scientific | Standalone soluble OpenMM; OpenMM membrane construction -> ParmEd -> GROMACS production | modules.automation.workflow; scientific.membrane; research.md_queue | pose/state, topology, RESP parameters, protocol | checkpoint, trajectory, per-phase records and analysis | config/default.yaml; scientific.membrane.mdp/build | distinct active scientific routes; target execution unverified in this refactor |
| WORKBENCH | 00_workbench | Persistent task orchestration, review, provenance and resource leases | webui.server -> webui.workflow -> pipeline_bridge | project-local PDB and configuration | run directories, stage records, reviews, handoffs | webui.workflow.defaults + research.config | active; historical pepwedge remains separate |

## Dependencies and evidence

GEN: PyTorch, Peptide2Mol, REINVENT4, RDKit. DOCK: official Vina, Meeko, RDKit. SCREEN: RDKit, MDAnalysis, PoseBusters, optional executors. INTERACT: RDKit, PyMOL, optional PLIP. MD: OpenMM, GROMACS, ParmEd, OpenFF, AmberTools, Psi4, MDAnalysis, MMPBSA.py. No dependency upgrade is part of this change.

Every numbered module has tests. Baseline commands/results are recorded in the validation report. Stored `deployment/validation_*`, `deployment/public_minimal/VALIDATION.json` and module test reports are historical evidence only; they do not establish this branch or a complete 500 ns run. No model weights or production trajectories are present in this source checkout.

## Tracked implementation and configuration inventory

### 00_workbench

114 tracked files.

- `00_workbench/pepwedge/__init__.py`
- `00_workbench/research/__init__.py`
- `00_workbench/run_evidence.py`
- `00_workbench/run_export_poses.py`
- `00_workbench/run_filters.py`
- `00_workbench/run_pipeline.py`
- `00_workbench/run_research.py`
- `00_workbench/run_webui.py`
- `00_workbench/scientific/__init__.py`
- `00_workbench/tests/test_md_contract.py`
- `00_workbench/tests/test_peptide_contacts.py`
- `00_workbench/tests/test_research_dual_blackwell.py`
- `00_workbench/tests/test_screening_revision.py`
- `00_workbench/tests/test_smoke.py`
- `00_workbench/tests/test_studio_integration.py`
- `00_workbench/tests/test_target_profile.py`
- `00_workbench/tests/test_validation.py`
- `00_workbench/tests/test_workflow_regressions.py`
- `00_workbench/webui/__init__.py`

### 01_molecule_generation

25 tracked files.

- `01_molecule_generation/src/__init__.py`
- `01_molecule_generation/tests/run_all.py`
- `01_molecule_generation/tests/test_adapters.py`
- `01_molecule_generation/tests/test_audit_regressions.py`
- `01_molecule_generation/tests/test_chemistry.py`
- `01_molecule_generation/tests/test_cores.py`
- `01_molecule_generation/tests/test_sampling_accounting.py`
- `01_molecule_generation/tests/test_structure.py`

### 02_docking

234 tracked files.

- `02_docking/configs/class_full.yaml`
- `02_docking/configs/class_quick.yaml`
- `02_docking/configs/receptor_policy.example.yaml`
- `02_docking/configs/research_gnina.yaml`
- `02_docking/configs/research_refine.yaml`
- `02_docking/configs/research_vina.yaml`
- `02_docking/configs/two_stage.yaml`
- `02_docking/configs/v2_default.yaml`
- `02_docking/env/environment.yml`
- `02_docking/env/sbvs2.environment.yml`
- `02_docking/pyproject.toml`
- `02_docking/scripts/run_docking.py`
- `02_docking/src/sbvs/__init__.py`
- `02_docking/src/sbvs/engines/__init__.py`
- `02_docking/src/sbvs/preparation/__init__.py`
- `02_docking/src/sbvs/web/__init__.py`
- `02_docking/src/sbvs2/__init__.py`
- `02_docking/src/sbvs2/analysis/__init__.py`
- `02_docking/src/sbvs2/docking/__init__.py`
- `02_docking/src/sbvs2/docking/run_gnina.py`
- `02_docking/src/sbvs2/docking/run_vina.py`
- `02_docking/src/sbvs2/ligands/__init__.py`
- `02_docking/src/sbvs2/receptor/__init__.py`
- `02_docking/src/yuanhe_preflight/__init__.py`
- `02_docking/tests/e2e/browser_acceptance.py`
- `02_docking/tests/e2e/pages_acceptance.py`
- `02_docking/tests/integration/artifact_acceptance.py`
- `02_docking/tests/integration/corruption_acceptance.py`
- `02_docking/tests/integration/enumeration_acceptance.py`
- `02_docking/tests/integration/recovery_acceptance.py`
- `02_docking/tests/integration/research_acceptance.py`
- `02_docking/tests/preflight/test_preflight.py`
- `02_docking/tests/unit/test_science.py`
- `02_docking/tests/v2/conftest.py`
- `02_docking/tests/v2/fixtures/1DRF.pdb`
- `02_docking/tests/v2/fixtures/FOL_crystal.sdf`
- `02_docking/tests/v2/test_core.py`
- `02_docking/tests/v2/test_ligands.py`
- `02_docking/tests/v2/test_prepare_dock.py`
- `02_docking/tests/v2/test_receptor.py`
- `02_docking/tests/v2/test_scheduler.py`

### 03_virtual_screening

137 tracked files.

- `03_virtual_screening/configs/budgets/full.yaml`
- `03_virtual_screening/configs/budgets/handoff_validation.yaml`
- `03_virtual_screening/configs/budgets/quick.yaml`
- `03_virtual_screening/configs/budgets/redock_trial.yaml`
- `03_virtual_screening/configs/host.yaml`
- `03_virtual_screening/configs/stages.yaml`
- `03_virtual_screening/configs/stages_pareto.yaml`
- `03_virtual_screening/configs/targets/ecmurj.yaml`
- `03_virtual_screening/configs/targets/generic.yaml`
- `03_virtual_screening/configs/targets/membrane_protein_template.yaml`
- `03_virtual_screening/pyproject.toml`
- `03_virtual_screening/src/zvss/__init__.py`
- `03_virtual_screening/src/zvss/chem/__init__.py`
- `03_virtual_screening/src/zvss/ingest/__init__.py`
- `03_virtual_screening/src/zvss/jobs/__init__.py`
- `03_virtual_screening/src/zvss/md/__init__.py`
- `03_virtual_screening/src/zvss/mech/__init__.py`
- `03_virtual_screening/src/zvss/membrane/__init__.py`
- `03_virtual_screening/src/zvss/motion/__init__.py`
- `03_virtual_screening/src/zvss/qa/__init__.py`
- `03_virtual_screening/src/zvss/receptor/__init__.py`
- `03_virtual_screening/src/zvss/redock/__init__.py`
- `03_virtual_screening/src/zvss/report/__init__.py`
- `03_virtual_screening/src/zvss/scoring/__init__.py`
- `03_virtual_screening/src/zvss/select/__init__.py`
- `03_virtual_screening/src/zvss/strain/__init__.py`
- `03_virtual_screening/src/zvss/web/__init__.py`
- `03_virtual_screening/src/zvss/web/static/__init__.py`
- `03_virtual_screening/src/zvss/web/templates/__init__.py`
- `03_virtual_screening/tests/data/1drf/box.json`
- `03_virtual_screening/tests/data/1drf/poses.sdf`
- `03_virtual_screening/tests/data/1drf/provenance.json`
- `03_virtual_screening/tests/data/1drf/receptor_prepared.pdb`
- `03_virtual_screening/tests/data/9nu4_updated.cif`
- `03_virtual_screening/tests/data/FIXTURES.txt`
- `03_virtual_screening/tests/data/P0AF16.fasta`
- `03_virtual_screening/tests/data/gnina_1drf_seed17.log`
- `03_virtual_screening/tests/test_budget_strain_fix.py`
- `03_virtual_screening/tests/test_chem.py`
- `03_virtual_screening/tests/test_config.py`
- `03_virtual_screening/tests/test_decision_policy.py`
- `03_virtual_screening/tests/test_handoff_identity.py`
- `03_virtual_screening/tests/test_ids.py`
- `03_virtual_screening/tests/test_membrane_aware_ymavs.py`
- `03_virtual_screening/tests/test_membrane_orientation.py`
- `03_virtual_screening/tests/test_pipeline_e2e.py`
- `03_virtual_screening/tests/test_redock_coordinates.py`
- `03_virtual_screening/tests/test_score_semantics.py`
- `03_virtual_screening/tests/test_sifts_mapping.py`
- `03_virtual_screening/tests/test_ui_integration_regressions.py`

### 04_interaction_analysis

36 tracked files.

- `04_interaction_analysis/env/environment.yml`
- `04_interaction_analysis/src/zlis/__init__.py`
- `04_interaction_analysis/tests/conftest.py`
- `04_interaction_analysis/tests/data/1HSG.pdb`
- `04_interaction_analysis/tests/data/3PTB.pdb`
- `04_interaction_analysis/tests/helpers.py`
- `04_interaction_analysis/tests/test_diagram.py`
- `04_interaction_analysis/tests/test_pipeline.py`
- `04_interaction_analysis/tests/test_thresholds.py`

### 05_molecular_dynamics

87 tracked files.

- `05_molecular_dynamics/config/default.yaml`
- `05_molecular_dynamics/config/demo_10ns.yaml`
- `05_molecular_dynamics/config/demo_10ns_cpu10.yaml`
- `05_molecular_dynamics/config/integration.yaml`
- `05_molecular_dynamics/config/smoke.yaml`
- `05_molecular_dynamics/config/v2_short.yaml`
- `05_molecular_dynamics/environment.yml`
- `05_molecular_dynamics/modules/__init__.py`
- `05_molecular_dynamics/modules/analysis/__init__.py`
- `05_molecular_dynamics/modules/automation/__init__.py`
- `05_molecular_dynamics/modules/md/__init__.py`
- `05_molecular_dynamics/modules/parameterization/__init__.py`
- `05_molecular_dynamics/modules/utils/__init__.py`
- `05_molecular_dynamics/run_md.sh`
- `05_molecular_dynamics/scripts/run_md.py`
- `05_molecular_dynamics/scripts/run_workflow.py`
- `05_molecular_dynamics/tests/fixtures/complex_R2.pdb`
- `05_molecular_dynamics/tests/test_backend_native_v2.py`
- `05_molecular_dynamics/tests/test_backend_v2.py`
- `05_molecular_dynamics/tests/test_chemistry_v2.py`
- `05_molecular_dynamics/tests/test_embedded_chemistry.py`
- `05_molecular_dynamics/tests/test_integrity.py`
- `05_molecular_dynamics/tests/test_inventory_v2.py`
- `05_molecular_dynamics/tests/test_mapping_v2.py`
- `05_molecular_dynamics/tests/test_workflow_v2.py`

### 06_websites

33 tracked files.

- `06_websites/01_public_demo/yuanhe-public/pnpm-lock.yaml`
- `06_websites/01_public_demo/yuanhe-public/pnpm-workspace.yaml`
- `06_websites/01_public_demo/yuanhe-public/tests/data.test.mjs`
- `06_websites/01_public_demo/yuanhe-public/tests/site.spec.ts`

### deployment

58 tracked files.

- `deployment/Dockerfile.research`
- `deployment/membrane-environment.yml`
- `deployment/preflight-environment.yml`
- `deployment/public_minimal/run_public.py`
- `deployment/public_minimal/run_public.sh`
- `deployment/research-environment.yml`
- `deployment/tests/test_environment_exports.py`
- `deployment/tests/test_packages.py`

### archive

40 tracked files.

- `archive/docking_prototype/configs/default.yaml`
- `archive/docking_prototype/source/sbvs2/__init__.py`
- `archive/docking_prototype/source/sbvs2/analysis/__init__.py`
- `archive/docking_prototype/source/sbvs2/docking/__init__.py`
- `archive/docking_prototype/source/sbvs2/docking/run_vina.py`
- `archive/docking_prototype/source/sbvs2/ligands/__init__.py`
- `archive/docking_prototype/source/sbvs2/receptor/__init__.py`
- `archive/docking_prototype/source/tests/conftest.py`
- `archive/docking_prototype/source/tests/fixtures/1DRF.pdb`
- `archive/docking_prototype/source/tests/fixtures/FOL_crystal.sdf`
- `archive/docking_prototype/source/tests/test_core.py`
- `archive/docking_prototype/source/tests/test_ligands.py`
- `archive/docking_prototype/source/tests/test_prepare_dock.py`
- `archive/docking_prototype/source/tests/test_receptor.py`

### data

5 tracked files.


### .github

4 tracked files.

- `.github/ISSUE_TEMPLATE/config.yml`
- `.github/workflows/engineering.yml`
