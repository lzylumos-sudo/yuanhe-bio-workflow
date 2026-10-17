# Historical build and acceptance record

Date: 2026-10-13. Dedicated WSL distribution GNINA-SBVS-1e0ea3f4; project UUID 1e0ea3f4-3fef-4108-ae52-4af02d120ed8. This records an earlier host, not current-release acceptance. Referenced runtime/results were excluded from GitHub.

1. Ubuntu Base 24.04.5 rootfs hash matched; WSL2 import/start succeeded with an external exFAT-hosted VHDX. An unrelated Ubuntu-22.04 distribution was retained unchanged.
2. Vina 1.2.7, Python 3.11.17, RDKit 2026.9.1 and Meeko 0.8.0 were used. 1DRF FOL A:187 had 32 heavy atoms; CCD chemistry preserved crystal coordinates. GLN12 repair was explicit, with no loop rebuilding or residue deletion.
3. Three fresh CCO browser runs took 1.973/1.891/1.902 seconds and produced actual affinity -2.546, consistent across UI/CSV/raw output; no prefilled results.
4. The classroom batch completed 20/20 compounds and 180 poses; actual page rendering and two distinct pose images were checked. Demonstration libraries had no activity labels.
5. GNINA 1.3.3 hash matched. Initial missing CUDA libraries and invalid alias failures were retained. crossdock_default2018_1_3 executed on GPU and CPU --no_gpu; a trace recorded 768 successful /dev/dxg ioctl calls. Multiseed 11/23/47 and second-stage operations were exercised independently of Vina.
6. Custom receptor review, mmCIF planning, external-PDBQT unknown provenance, Molscrub 0.2.2 states and parent summaries were checked. Vina exported 180 chemically consistent poses. FOL same-frame symmetry-aware top-1 RMSD was 1.8458 A. Manual MD export gates were tested with labeled automated samples, not an actual MD calculation.
7. Cancellation, invalid entries, service interruption, completed-job reuse, corruption rejection/retry, offline operation, path forwarding, same-host registration and safe stop/restart were tested. A 13,558,968,320-byte runtime tar was exported and hashed, but restoration, copied-image registration and a second physical host were not validated.

Historical evidence names include reports/browser_acceptance.json, class_report_images.json, gnina_* records, gpu_device_trace.json, research_acceptance.json, artifact_acceptance.json, enumeration_acceptance.json, FOL_redocking.json and recovery/corruption/offline/path/runtime-export reports. Their presence in an old record is not a claim that those files ship here.

Unsupported scope included covalent docking, general metals/complex cofactors, flexible receptors, long-loop reconstruction and reliable MOL2 chemistry. No enrichment validation was performed. CPU CNN ran on a GPU-equipped host, not a separate GPU-less machine. An instructor endpoint returned HTTP 401; no authenticated case or remote-platform acceptance was obtained. Preserve raw results and inspect evidence before making new claims.
