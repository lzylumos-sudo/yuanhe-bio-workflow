# Runtime and configuration upgrade

Use one immutable administrator release and its matching source manifest, frontend and configuration templates. The current development build is 0.11.0-dev.5; historical tags and research tasks retain their original versions and configurations.

## Administrator settings

For a generic machine, copy `deployment/site.example.json` outside the source directory. For the previously inspected campus machine, copy `deployment/site.campus-upgrade.json`. It selects a new `/home/user/yuanhe-runtime-0.11.0-dev.5` prefix, reuses the existing public model bundle and CUDA GROMACS executable, binds `10.71.121.57`, and uses **8778** for validation while the existing service keeps **8777**.

The campus template contains observed paths, not a guarantee that they still exist. Preflight checks them and refuses an existing destination. If the destination already exists from a failed attempt, preserve its logs and choose another new prefix. Do not overwrite it or reuse its evidence as a successful install.

The original installer bounds uv downloads to four concurrent requests, with a 60-second connection timeout, 120-second read timeout and five HTTP retries. These defaults address transient connection failures observed while fetching REINVENT build dependencies alongside large CUDA wheels. They do not change package indexes, locked versions, hashes or TLS checks. Explicit administrator `UV_CONCURRENT_DOWNLOADS`, `UV_HTTP_CONNECT_TIMEOUT`, `UV_HTTP_TIMEOUT` and `UV_HTTP_RETRIES` settings take precedence and the effective values are logged. A persistent network failure still stops installation and prevents service handoff; preserve the evidence and retry in a new prefix after resolving connectivity. See the [official uv environment settings](https://docs.astral.sh/uv/reference/environment/).

`model_source` supplies only verified local weights and, when available, cached Git objects for the pinned upstream commits. Restoration reconstructs clean vendor code inside the new source directory, then the installer applies this release's patches there. It never copies the old vendor working tree or imports its installed environment. If the old vendor directories have no usable Git metadata, restoration fetches the same pinned public source commits. REINVENT4 is installed from the new directory into the new prefix, and generated runtime paths point to that new environment.

## Install without changing the running service

1. Download the named administrator ZIP and SHA256SUMS from the release linked by `deployment/current_release.json`, verify checksums, and extract into a new source directory. Run `python3 deployment/verify_release.py`. Alternatively, clone the exact new release tag into a new directory and record `git rev-parse HEAD`. Do not update a source directory used by running workers.
2. Copy and review the chosen site JSON outside that directory. Run `python3 deployment/setup_workbench.py --config /home/user/yuanhe-site-0.11.0-dev.5.json`. This is a read-only plan; fix any failed checks before proceeding. Source configuration hashes must match the recorded protocol; a cross-platform checkout must preserve the repository's line endings.
3. Repeat with `--execute`. The installer creates new environments and runtime settings from current defaults, validates dependencies and records its source identity. It does not switch the service or copy old frozen scientific settings into new runs.
4. Start the exact launch script printed in the new prefix's `evidence/deployment.json`. The campus template uses port 8778. Keep old job/checkpoint directories separate.
5. From the new source, run `python3 deployment/verify_service.py --url http://10.71.121.57:8778 --expected-commit FULL_SOURCE_COMMIT --out /home/user/yuanhe-service-0.11.0-dev.5.json`, replacing the commit with the verified package manifest's `git_commit` or exact Git HEAD. This checks source/version and seven served frontend file hashes; version text alone is insufficient.

## Compatibility and acceptance

Old unversioned runtime settings remain readable; unknown future schemas fail explicitly. Runtime normalization derives legacy REINVENT4 interpreter paths and uses the selected module environment for tools. New installations generate the current runtime schema. Keep old runtime files as evidence; do not copy old run configurations or change their source hashes to make a task resume under different code.

To inspect a legacy runtime migration without overwriting it, run `python3 deployment/migrate_runtime.py --input /path/to/old-runtime.json --output /path/to/new-runtime-review.json`. The destination must not exist. The export records the original hash and explicitly leaves execution validation unperformed. For a full reinstall, prefer the current site template and installer instead of copying this file over newly generated runtime settings.

Verify diagnostics and public small-input UI paths. Confirm that a second conflicting task reports a recoverable resource-busy state, explicitly **not queued**, that a successful refresh clears resolved errors, and that zero valid samples are displayed as zero. Model completion is not molecular validity, and molecular validity is not candidate qualification. Keep scientific thresholds and required manual reviews unchanged. Validate short MD and checkpoint recovery separately before starting long scientific runs.

After acceptance, stop new submissions to the old service, wait for or safely checkpoint its workers, and have the administrator switch the service command to the new source/runtime. To retain port 8777, use the new source's `start_linux.sh --host 10.71.121.57 --port 8777` under the service account after the old listener stops. Confirm the same exact commit and frontend hashes again. Keep the old environment and service definition for rollback; never let both versions write the same task directory.

The release and configuration do not claim that the target machine was redeployed, that a 500 ns trajectory completed, or that generated candidates have experimental activity.
