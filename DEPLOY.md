# Administrator download and deployment

Always open [the current administrator entrypoint](https://github.com/lzylumos-sudo/yuanhe-bio-workflow/blob/main/DEPLOY.md) when preparing a new deployment. Its [README download button](https://github.com/lzylumos-sudo/yuanhe-bio-workflow#yuanhe-bio-platform) and [machine-readable release record](deployment/current_release.json) are generated from the single root [VERSION](VERSION).

1. Download the named `yuanhe-bio-platform-<version>.zip` release asset and `SHA256SUMS.txt` from that current entrypoint. Keep the versioned filename. A previously saved ZIP or bookmark to an older release does not update itself.
2. Verify the downloaded checksum with `sha256sum -c SHA256SUMS.txt` (Linux) or `shasum -a 256 -c SHA256SUMS.txt` (macOS). Extract into a new directory. The named package includes `SOURCE_MANIFEST.json`; GitHub's automatic source ZIP does not include this delivery manifest.
3. Inside the extracted directory run `python3 deployment/verify_release.py`. It checks every packaged file against the source manifest and confirms the software version. Preserve its full source commit with the installation record. `python3 yuanhe version` reports the software version separately from Git identity; without a Git checkout, use the verified delivery manifest for the packaged source commit.
4. Compare the extracted `VERSION` with the **current** entrypoint before installation. Stop if it names a newer version. Historical packages remain available only for reproducibility and controlled recovery.
5. Follow the [complete workbench deployment guide](deployment/WORKBENCH_DEPLOYMENT.md): install environments, restore licensed models and evidence, configure paths, run diagnostics, and complete small real-machine acceptance before starting production research.

For a Git deployment, use the current release's immutable tag from `deployment/current_release.json`, then record `git rev-parse HEAD`. `main` is the integration branch; tags identify exact deployment builds. Never replace an active computation's checkout or frozen run records.

## Current software, historical evidence

The current build integrates version governance and deployment-readiness fixes. Its development label is intentional: engineering CI does not certify a complete GPU installation or a 500 ns trajectory on a new host. The formal research setting remains one 500 ns trajectory per candidate; the requested quarter-size batch selects two candidates.

Older product/package labels, archived scientific results and development branches are indexed in [version history](docs/releases/HISTORY.md). The repeated text beside files in GitHub's code browser is each file's last commit message, not a product-version list. Use the current download button for installation and the history index for old records.
