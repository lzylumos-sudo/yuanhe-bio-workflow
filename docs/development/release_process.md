# Release process

Edit root VERSION, run the metadata generator, update changelog/compatibility notes, and run scientific-baseline, registry, protocol, package, API and frontend checks. Include exact Git tag/commit, tests collected/executed/skipped, unchanged versus revised scientific protocols, dependency locks and known limitations.

Do not upgrade a protocol for a UI-only release. Do not reuse a protocol ID after a scientific change. Do not label unvalidated hardware/model/backend combinations stable. Existing generated metadata is checked in for standalone distribution and must be regenerated, not edited manually.

Heavy GPU validation is manual/target-specific and separate from pull-request CI. The original 0.11.0-dev.1 governance branch was not a public Release. Publish current administrator builds only after protected merge and verification of the merged main commit. Never move an existing tag or rewrite historical run identities.

## Administrator publication

1. Advance root `VERSION` for the build; run `deployment/sync_platform_version.py`. The generator updates package metadata, the README current-download block and `deployment/current_release.json` together. CI rejects stale generated copies.
2. Preserve former identities in `docs/releases/HISTORY.md`; keep original evidence documents unchanged. Merge through protected PR checks and verify the distinct post-merge main run.
3. Tag the verified merged commit as `v<VERSION>`. Run `python3 deployment/package_release.py --ref <commit> --output /outside/repository/yuanhe-bio-platform-<VERSION>.zip`; it packages only that committed tree, never local ignored environments or run data.
4. Extract the ZIP into a fresh directory and run its `deployment/verify_release.py`. Save the ZIP SHA-256 as `SHA256SUMS.txt`. Create a GitHub Release with the named ZIP and checksum assets. Mark development versions as prereleases; never label them stable to obtain GitHub's Latest badge.
5. Download the published assets again and verify their checksums and source manifest. Check the README download URL and release tag match root VERSION. The permanent administrator entrypoint is `/blob/main/DEPLOY.md`; do not send older fixed ZIP links as the current download.

GitHub file-list commit messages describe the last change to each path. They are not version labels and must not be rewritten just to make directory rows look uniform.
