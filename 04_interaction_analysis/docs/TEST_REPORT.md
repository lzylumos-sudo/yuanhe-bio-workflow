# ZLIS historical test and installation record

Historical host: macOS 26.5.1 arm64, Python 3.11.17; RDKit 2026.9.1, ProLIF 2.2.3, MDAnalysis 2.10.0, Biopython 1.88, OpenMM 8.6.1, pandas 3.0.6, openpyxl 3.1.5, NumPy 2.4.6, gemmi 0.7.5, CairoSVG 2.7.1, pytest 8.3.4 and open-source PyMOL 3.1.0. Optional PLIP 3.0.0 used a separate interpreter. This is not a target-Windows acceptance record.

The earlier run recorded 40 passing tests: 16 interaction-threshold tests, 17 pipeline tests and 7 diagram tests. Coverage included donor/acceptor roles, hydrogen-bond distance 3.49/3.55 A and angle 125/135 degrees, salt-bridge 4.49/4.6 A boundaries and charge negatives; aromatic, cation-pi, halogen, metal and water candidates; atom identity, unchanged coordinates, template errors, unknown protonation, insertion codes, model/chain separation, relocated paths, scene reopening and diagram bounds/layout round trips.

| System | Ligand | Recorded interactions |
|---|---|---|
| 3PTB trypsin | A:BEN1, benzamidinium template | 3 hydrogen bonds, 4 salt bridges, 4 water candidates, 2 hydrophobic, 5 vdW |
| 1HSG HIV protease | B:MK1902, RCSB default requiring confirmation | 4 hydrogen bonds, 13 water candidates, 7 hydrophobic, 9 vdW |

3PTB hydrogen bonds: ASP189 N2-OD2 2.87 A / 147.0 degrees; ASP189 N1-OD1 2.92 A / 151.3 degrees; GLY219 N1-O 2.82 A / 144.6 degrees. 1HSG assignments included A:ASP25, B:ARG8, B:GLY27 and B:ASP29 with distinct donor roles.

PLIP and ZLIS were not automatically reconciled. PLIP reported two ASP189 hydrogen bonds for 3PTB and omitted GLY219 and the salt bridge; for 1HSG it reported six bonds versus four, with a different ASP25 chain assignment. These disagreements require manual review and are not resolved by this record.

Remaining limitations: Windows launchers were not executed; earlier diagram checks were programmatic rather than visual; browser drag/save was not automated; hydrogen bonds and salt bridges may coexist on one atom pair; ligand protonation and HID histidines are assumptions; water bridges have no donor/acceptor assignment; aromatic and hydrophobic criteria may overlap. Alternate conformer selection defaults to highest occupancy, and a batch multi-file GUI is not implemented. exFAT can regenerate AppleDouble metadata; exclude it from releases.
