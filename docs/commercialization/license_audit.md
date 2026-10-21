# Preliminary license and ownership review

Technical inventory, not a legal opinion or permission to distribute a commercial bundle. No license or attribution has been changed.

## Confirmed local facts

Root LICENSE is PolyForm Noncommercial 1.0.0. LICENSE_SCOPE.txt explicitly preserves prior MIT grants for ZVSS, third-party terms and pre-existing rights. The notices identify Peptide2Mol, REINVENT4, GNINA, RDKit, OpenMM, GROMACS, ProLIF, PyMOL and PoseBusters. deployment/model_dependencies.json records source commits and hashes for Peptide2Mol/REINVENT4 weights. Hashes establish identity, not permission. Local license-related history includes commits 1873f36 and c40cae4; this pass has not reconstructed a complete file-by-file historical grant map.

| Risk | Issue | Before commercial distribution |
|---|---|---|
| High | Weights and training/data rights are separate from code | Obtain terms for exact pinned artifacts, permitted uses and redistribution |
| High | Contributor, school/lab and partner ownership not established | Confirm employment/research agreements and contributor authorization |
| High | Original project's noncommercial license does not license third-party dependencies | Review full runtime distribution and any derivative/copy-left obligations |
| High | Public history and prior permissive grants | Map affected versions/files; do not promise exclusive ownership of previously granted code |
| Medium | Incomplete software bill of materials and deployment provenance | Inventory versions, license texts, notices, linking/process boundaries and distribution modes |
| Medium | Patent novelty and public disclosures | Specialist review before publishing potentially patentable methods; deletion cannot reverse disclosure |

## Component scope requiring separate confirmation

| Component | Local evidence | Code vs weights/data/commercial use |
|---|---|---|
| Peptide2Mol | Notices state MIT code; pinned source/weight hash | Weight and training-data conditions not cleared by the MIT statement |
| REINVENT4 / LibInvent | Notices state Apache-2.0 code; pinned prior | Exact prior terms and upstream data provenance require review |
| GNINA | External binary integration | Exact build, bundled neural models and redistribution obligations need version-specific review |
| RDKit | Chemistry dependency | Confirm installed version and bundled notices; no claim to upstream IP |
| GROMACS / OpenMM | Two distinct execution paths | Review actual builds, dependency obligations and distribution packaging |
| ProLIF / PyMOL / PoseBusters / AmberTools and other dependencies | Environment/module manifests | Must be included in full runtime review; their terms are not replaced by the root license |

Dual licensing may be possible for rights actually controlled by the project, subject to prior grants and third-party obligations. It is not established for the entire runtime or weights. The original rights holder's commercial authorization cannot grant rights it does not own. Do not assume a separate process alone removes all licensing obligations.

[PolyForm standard text](https://polyformproject.org/licenses/noncommercial/1.0.0) defines permitted purposes and exceptions; the standard license and unchanged local scope govern. Full upstream/version-specific verification is deferred, so this document is not a completed commercialization clearance.
