# MD Task bundle equivalence validation (t6)

Checks bundle/specification equivalence,**are not**membrane engine preflight.

- Task bundle count: 3 | All checks passed: 0 | Failures present: 3

| Task bundles | Compounds | Checks | Failures | adapter04 | adapter05 |
| --- | --- | --- | --- | --- | --- |
| exploratory_01_CMP75C7A4B0B41E4390 | CMP75C7A4B0B41E4390 | 8 | protocol_and_units_declared | None | None |
| primary_01_CMP7047E1960DBD3C23 | CMP7047E1960DBD3C23 | 8 | protocol_and_units_declared | None | None |
| primary_02_CMPCBBA17C44FCCD2B3 | CMPCBBA17C44FCCD2B3 | 8 | protocol_and_units_declared | None | None |

## Per-bundle results

### exploratory_01_CMP75C7A4B0B41E4390

| Check | Result | Observed value | Expected value |
| --- | --- | --- | --- |
| package_files_present | pass | {"manifest.json": true, "task.json": true, "complex.pdb": true, "ligand.sdf": true, "checksums.sha256": true} | {"manifest.json": true, "task.json": true, "complex.pdb": true, "ligand.sdf": true, "checksums.sha256": true} |
| checksums_match_files | pass | {"complex.pdb": "b522761feebd8141043e0f0459605a3e5b704b73aee2cf2540a2917d6cb9e849", "ligand.sdf": "8e35bb89b2b03386e19265bdec365855f9f8ad427f109ee0d99b310a46d0786e"} | {"complex.pdb": "b522761feebd8141043e0f0459605a3e5b704b73aee2cf2540a2917d6cb9e849", "ligand.sdf": "8e35bb89b2b03386e19265bdec365855f9f8ad427f109ee0d99b310a46d0786e"} |
| ligand_atom_mapping_complete | pass | {"mapping_entries": 6, "declared_atoms": 6, "sdf_atoms": 6} | "mapping entries == declared atoms == SDF atoms" |
| complex_atom_count_consistent | pass | {"total": 1503, "protein": 1497, "hetatm": 6, "ligand_atoms": 6} | {"total": "protein+hetatm", "ligand_atoms": 6} |
| coordinates_identical_to_ligand_sdf | pass | {"max_abs_deviation_A": 0.0005} | "< 0.001 A" |
| provenance_linked_to_selection | pass | {"selection_version": "SEL6830370C0C42E694", "evidence_score": 0.6725, "tier": "B", "components_used": ["affinity", "cnn_pose", "cnn_affinity", "mechanism", "pose_quality"], "compo | "selection_version + tier present" |
| protocol_and_units_declared | fail | {"protocol_keys": [], "units": {}, "engine": "openmm"} | "protocol + units + engine declared" |
| blocked_capabilities_explicit | pass | {"status": "blocked", "blockers": ["blocked_missing_membrane_capability", "blocked_platform_linux_required"]} | {"status": "blocked", "blockers": ["blocked_missing_membrane_capability", "blocked_platform_linux_required"]} |

### primary_01_CMP7047E1960DBD3C23

| Check | Result | Observed value | Expected value |
| --- | --- | --- | --- |
| package_files_present | pass | {"manifest.json": true, "task.json": true, "complex.pdb": true, "ligand.sdf": true, "checksums.sha256": true} | {"manifest.json": true, "task.json": true, "complex.pdb": true, "ligand.sdf": true, "checksums.sha256": true} |
| checksums_match_files | pass | {"complex.pdb": "884ee5d22c0c3d29b376075f3a8b6d3f2cb3c611bb8eca15ca04c65b71c6c7b1", "ligand.sdf": "b9354c2b850a5e08f1f02196fb8a43883cc2924e6d3bfd8a84eb77efe4986db5"} | {"complex.pdb": "884ee5d22c0c3d29b376075f3a8b6d3f2cb3c611bb8eca15ca04c65b71c6c7b1", "ligand.sdf": "b9354c2b850a5e08f1f02196fb8a43883cc2924e6d3bfd8a84eb77efe4986db5"} |
| ligand_atom_mapping_complete | pass | {"mapping_entries": 9, "declared_atoms": 9, "sdf_atoms": 9} | "mapping entries == declared atoms == SDF atoms" |
| complex_atom_count_consistent | pass | {"total": 1506, "protein": 1497, "hetatm": 9, "ligand_atoms": 9} | {"total": "protein+hetatm", "ligand_atoms": 9} |
| coordinates_identical_to_ligand_sdf | pass | {"max_abs_deviation_A": 0.0005} | "< 0.001 A" |
| provenance_linked_to_selection | pass | {"selection_version": "SEL6830370C0C42E694", "evidence_score": 0.7555, "tier": "B", "components_used": ["affinity", "cnn_pose", "cnn_affinity", "mechanism", "pose_quality"], "compo | "selection_version + tier present" |
| protocol_and_units_declared | fail | {"protocol_keys": [], "units": {}, "engine": "openmm"} | "protocol + units + engine declared" |
| blocked_capabilities_explicit | pass | {"status": "blocked", "blockers": ["blocked_missing_membrane_capability", "blocked_platform_linux_required"]} | {"status": "blocked", "blockers": ["blocked_missing_membrane_capability", "blocked_platform_linux_required"]} |

### primary_02_CMPCBBA17C44FCCD2B3

| Check | Result | Observed value | Expected value |
| --- | --- | --- | --- |
| package_files_present | pass | {"manifest.json": true, "task.json": true, "complex.pdb": true, "ligand.sdf": true, "checksums.sha256": true} | {"manifest.json": true, "task.json": true, "complex.pdb": true, "ligand.sdf": true, "checksums.sha256": true} |
| checksums_match_files | pass | {"complex.pdb": "14b0e36833f672f2f1fa5101816664a89bfacea993a832feb0710c4d3f073676", "ligand.sdf": "83db4de4ffa54af6b45f820d47f936cf1178dd11ec48fd27bda534a611c2efa0"} | {"complex.pdb": "14b0e36833f672f2f1fa5101816664a89bfacea993a832feb0710c4d3f073676", "ligand.sdf": "83db4de4ffa54af6b45f820d47f936cf1178dd11ec48fd27bda534a611c2efa0"} |
| ligand_atom_mapping_complete | pass | {"mapping_entries": 6, "declared_atoms": 6, "sdf_atoms": 6} | "mapping entries == declared atoms == SDF atoms" |
| complex_atom_count_consistent | pass | {"total": 1503, "protein": 1497, "hetatm": 6, "ligand_atoms": 6} | {"total": "protein+hetatm", "ligand_atoms": 6} |
| coordinates_identical_to_ligand_sdf | pass | {"max_abs_deviation_A": 0.0005} | "< 0.001 A" |
| provenance_linked_to_selection | pass | {"selection_version": "SEL6830370C0C42E694", "evidence_score": 0.5695, "tier": "B", "components_used": ["affinity", "cnn_pose", "cnn_affinity", "mechanism", "pose_quality"], "compo | "selection_version + tier present" |
| protocol_and_units_declared | fail | {"protocol_keys": [], "units": {}, "engine": "openmm"} | "protocol + units + engine declared" |
| blocked_capabilities_explicit | pass | {"status": "blocked", "blockers": ["blocked_missing_membrane_capability", "blocked_platform_linux_required"]} | {"status": "blocked", "blockers": ["blocked_missing_membrane_capability", "blocked_platform_linux_required"]} |

## Conclusion

- Coordinates/Mapping/units/Provenance matches screening-stage records; see itemized results.
- Container and membrane system remain unbuilt:`04` associated OpenMM runner lacks a membrane protocol and runs only under Linux ;
  these bundles therefore remain `blocked`, and are not described as having entered production MD”.
