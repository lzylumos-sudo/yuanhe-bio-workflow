# Historical GNINA runtime investigation

This record concerns an earlier WSL2 Ubuntu-22.04 installation on an RTX 5070 host. It does not certify the current project runtime.

GNINA 1.3.3 `gnina.cuda12.8.static` was downloaded according to the upstream vendor lock: 2,056,131,000 bytes, recorded SHA-256 prefix/suffix `3340c1f4...f72172`. The old `/opt/sbvs/env` environment was absent on that machine. WSL driver libraries existed, but complete CUDA user-space dependencies did not.

| Failure | Historical action | Outcome |
|---|---|---|
| libcudnn.so.9 missing | NVIDIA cuDNN 9.8.0.87 CUDA12 archive; hash recorded | Resolved |
| libcudart.so.12 missing | Matching NVIDIA CUDA runtime wheel library | Resolved |
| libcublas.so.12 missing | Matching cuBLAS and cuBLAS-Lt libraries | Resolved |
| libcusparse.so.12 missing | Remaining CUDA runtime not installed | Blocked |

The filename does not mean a fully statically linked executable. A WSL driver interface alone does not supply the CUDA runtime. The planned repair was to install matching user-space dependencies, verify `gnina --version`, then execute the five-parent, three-seed redocking trial with real logs. Until actual execution succeeds, redocking remains not_evaluated/needs_input and cannot satisfy formal qualification. Current instructions scope library paths to GNINA rather than exporting them globally into unrelated model environments.
