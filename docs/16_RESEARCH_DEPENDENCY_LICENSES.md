# Research dependency licenses and distribution

The project retains its existing PolyForm Noncommercial policy and exclusions in
LICENSE_SCOPE.txt. Replacing GNINA does not remove other dependencies' obligations
or relicense historical code. No binaries, environments, neural weights, private
candidates or task-generated simulation results accompany this source change.

| Component / pinned core version | License evidence | Distribution obligations to retain |
|---|---|---|
| AutoDock Vina 1.2.7 | [Official Apache-2.0 license](https://github.com/ccsb-scripps/AutoDock-Vina/blob/v1.2.7/LICENSE) | License/notice and modification notices; upstream authorship |
| Meeko 0.8.0 | [Official LGPL license](https://github.com/forlilab/Meeko/blob/develop/LICENSE); installed metadata LGPL-2.1-only | LGPL source/relinking and notice duties as applicable to distributed builds |
| RDKit 2026.03.1 | [BSD-3-Clause](https://github.com/rdkit/rdkit/blob/master/license.txt) | Copyright, terms and disclaimer; inspect bundled components |
| Psi4 1.11 | [License sources](https://github.com/psi4/psi4/tree/master); package metadata LGPL-3.0-only AND BSD-3-Clause AND MIT | Component notices and applicable source/relinking duties; bundled math/QM dependencies keep their terms |
| cdsgroup RESP 1.0.0 | [BSD-3-Clause](https://github.com/cdsgroup/resp/blob/master/LICENSE) | Copyright, terms and disclaimer |
| AmberTools 24.8: antechamber, parmchk2, tleap, PBSA, MMPBSA.py, cpptraj | [AmberTools licensing](https://ambermd.org/AmberTools.php); installed package GPL-3.0-or-later AND LGPL-3.0-or-later AND BSD-3-Clause AND MIT | Inspect each distributed component, include corresponding licenses/source obligations; not the proprietary Amber/pmemd license |
| OpenMM 8.6.1 | [Official license](https://github.com/openmm/openmm/blob/master/docs-source/licenses/Licenses.txt); package LGPL-3.0-or-later | Preserve MIT/LGPL component notices and applicable source/relinking terms |
| openmmforcefields 0.15.1 / OpenFF Toolkit 0.17.1 | [OpenMM forcefields](https://github.com/openmm/openmmforcefields), [OpenFF Toolkit](https://github.com/openforcefield/openff-toolkit) | Preserve software and separately sourced force-field/data notices; do not assume code license covers every parameter set |
| ParmEd 4.3.1 | [LGPL-2.1-or-later](https://github.com/ParmEd/ParmEd/blob/master/LICENSE) | Notice and applicable source/relinking obligations |
| MDAnalysis 2.10.0 | [Official license](https://github.com/MDAnalysis/mdanalysis/blob/develop/LICENSE); installed metadata LGPL-3.0-or-later | Apply the pinned version's terms; older MDAnalysis versions had different license metadata |
| GROMACS, selected host build | [LGPL-2.1-or-later policy](https://github.com/gromacs/gromacs/blob/main/COPYING) and distributed COPYING | Preserve license/source obligations; record exact host version/build and linked FFT/GPU libraries |
| CUDA and NVIDIA runtime | Vendor distribution terms | Required only by selected GPU builds; driver visibility is not permission to redistribute all libraries |

The environment recipe is not an exhaustive bill of materials. Export the solved
platform lock and retain `conda-meta` license metadata and each package's license
files for a binary/container release. Linking, redistribution and modification
conditions require evaluation for that exact release. Citation of scientific
methods is separate from software license compliance. Public demo assets retain
React/Vite/NGL and structural-data notices. Historical GNINA notices remain for
historical code/results; they are not a new production runtime requirement.
