# Troubleshooting

- Port 5902 occupied: verify project UUID before acting. Choose another explicit port; the launcher does not kill unknown processes.
- WSL E_ACCESSDENIED: resolve host permissions/policy; do not silently switch to an unrelated distribution.
- Missing libcudnn.so.9 or CUDA12 libraries: inspect installation logs and scoped runtime paths. Do not install a Linux NVIDIA driver inside WSL.
- Invalid CNN alias: use the actual executable help and pinned protocol. The recorded model is crossdock_default2018_1_3, not an assumed default1.3 alias. Preserve failed attempts.
- Meeko template failure: inspect missing atoms and nonstandard residues; do not delete residues to hide errors. The 1DRF classroom repair explicitly adds missing atoms without rebuilding loops.
- INVALID_SMILES / UNSUPPORTED_ELEMENT: failures.csv records rejected entries while valid compounds continue.
- REVIEW_REQUIRED: research preparation requires a hash-bound policy; this is distinct from the classroom shortcut.
- Protocol/environment changed: create a new run from saved inputs rather than reusing mismatched cache.
- GNINA device not verified: execute CNN validation on the current hardware. Detection alone is insufficient.
- not_tested / blocked: neither means passed.
- Historical instructor endpoint HTTP 401: authentication was required; no case was accessed and that backend was not used for computation.
