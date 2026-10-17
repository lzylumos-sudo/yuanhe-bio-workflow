# MurJ v2: historical implementation record

This implementation was separated from the old `pepwedge` results. Historical SDF files failed RDKit readback and must not be used as v2 chemical inputs. Current structural input is `../data/structures/9nu4.pdb`; provenance is recorded in its manifest.

| Dependency | Pinned identity |
|---|---|
| Peptide2Mol | `2c5af9628aeeab5ed2dff5dc5326ce2463316291` |
| PMT checkpoint SHA-256 | `7dc82c4f895efea83764ce25c2fae5787f9f25b7935590e5119bde7199540129` |
| REINVENT4 | `3cf206ea5b385ac955dc44bbae6063b55ee31e0a` |
| LibInvent prior | Zenodo version DOI `10.5281/zenodo.20701824`; hash in the dependency manifest |

`inputs/raw/` records input references; `configs/project.json` defines parameters and evidence boundaries. `structure.py` and `hypotheses.py` prepare structures; `adapters/` invokes actual models; chemistry, core and synthesis modules track identity and evidence; geometry, diversity and export modules apply gates. Runs preserve logs and staged artifacts, and only successful stages receive `stage_success.json`.

Models communicate through files and run in separate environments. The earlier Apple Silicon machine lacked NVIDIA/CUDA. An observed CPU sampling startup was not a completed inference validation. No substitute molecules may be emitted after model failure. Model commands preserve hashes, source identity and stdout/stderr; bad hashes, missing weights, invalid input or model errors fail.

The strict standalone export requires traceable real Peptide2Mol seeds and LibInvent events; valid single-component RDKit chemistry without dummy/query atoms; unique standardized identity; fully supported synthesis evidence; at least one conformer passing a local hypothesis and severe-clash gate; and complete foreign-key/SDF readback checks including charges, bond orders, stereochemistry and coordinates. Deliver a smaller qualified subset with reasons rather than duplicate molecules or count conformers as compounds.

Historical environment paths, CPU timing, incomplete sampling and then-unfixed Hydra quoting issues were development observations. They are not statements about current code. Use the current setup and validation records. Environments must remain on a filesystem supporting their symlinks and absolute prefixes, rather than inside an exFAT model directory.
