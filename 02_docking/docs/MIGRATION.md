# Migration and recovery

These standalone runtime scripts require the separately retained WSL image; the image is not distributed through GitHub.

1. Stop submissions and finish or explicitly cancel active work.
2. Run PREPARE_MIGRATION.ps1, optionally with -CancelActive. It verifies process identity/UUID, stops this project's processes, syncs and terminates only its distribution.
3. Confirm other disk users have released files before ejecting. Do not copy an active VHDX.
4. On a new Windows x86-64/WSL2 host, run REGISTER_ON_THIS_PC.ps1 against the retained image. Registration uses path and UUID and avoids overwriting conflicting distributions.
5. Rerun doctor, the CCO probe and actual GNINA CPU/GPU CNN checks. Hardware changes invalidate old device acceptance.

Project paths are resolved through the current directory and wslpath; Linux prefixes remain inside ext4. Windows, WSL registration and GPU drivers are host prerequisites. Update WSL explicitly if import-in-place is unsupported. Never use `wsl --unregister` as a supposedly safe detach operation; it may delete distribution data.

EXPORT_RUNTIME.ps1 safely stops this project, exports a tar snapshot and records SHA-256. Restore into a new location with `wsl --import UNIQUE_NAME NEW_DIRECTORY SNAPSHOT --version 2`, verify the project UUID and retain the original until recovery succeeds. Exporting a snapshot is not proof of restoration.

After interruption, matching orphan processes and run states are reconciled. Resume reuses only checksum-valid completed jobs; incomplete or corrupted jobs remain failed/interrupted and are retried explicitly. Cancellation is not completion.

Historical checks covered same-host paths, forwarding, existing-image registration, stop/restart and offline operation. A second physical computer, copied-image registration, Windows ARM and WSL1 were not accepted. Frozen absolute paths may require new runs after relocation; do not rewrite evidence to imply automatic portability.
