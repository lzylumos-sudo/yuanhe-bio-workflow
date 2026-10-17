# PoseBusters Actual execution record

- PoseBusters Version: `0.6.5`
- ASCII Runtime temporary directory: `C:\Users\user\AppData\Local\Temp\dsh-iKCsDM\zvss-scratch`
- Actual pose source: `F:\MemPilot\runs\public20_engineering\entries\P000005_S01_C01_seed20261015\input.sdf`
- Receptor: `F:\MemPilot\runs\public20_engineering\protocol_inputs\receptor_prepared.pdb`

| Case | Atom count | Checks | Passed | Failures |
| --- | --- | --- | --- | --- |
| real_pose | 12 | 22 | 22 | None |
| translated_8A_into_protein | 12 | 22 | 20 | minimum_distance_to_protein, volume_overlap_with_protein |
| stretched_bond_2.2A | 12 | 22 | 20 | bond_lengths, internal_energy |
| all_atoms_overlapped | 12 | 21 | 18 | bond_lengths, bond_angles, protein-ligand_maximum_distance |

## Decision

- Number of failed checks for actual pose: 0(should be 0)
- Failures for pose translated into protein: minimum_distance_to_protein, volume_overlap_with_protein
- Failures for stretched-bond pose: bond_lengths, internal_energy
- Failures for overlapping atoms: bond_lengths, bond_angles, protein-ligand_maximum_distance

Failures come directly from PoseBusters boolean check columns without postprocessing or clipping.
