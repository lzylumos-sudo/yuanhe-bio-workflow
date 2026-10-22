# Version identities

| Identity | Authority | Meaning |
|---|---|---|
| Platform software | root `VERSION` | SemVer; current channel and download are recorded in `deployment/current_release.json` |
| Git revision | execution checkout | Full 40-character SHA, branch, dirty flag; unknown for archives without Git |
| Scientific protocol | `configs/protocol_registry.yaml` | Immutable scientific recipe independent of UI/software version |
| Model/data/tool | actual run evidence | Declared sources kept separate from measured hashes and executable versions |
| Run | new task directory | `YH-YYYYMMDD-<UUID4>`; UUID4 is deliberately used instead of adding a ULID dependency |

The previous runtime reported 0.4.0-rc1, Web package 1.0.0, docking distribution 0.1.0, screening distribution 0.10.0 and standalone automation frontend 2.0.0. None was a repository tag. The original governance branch used **0.11.0-dev.1**, chosen after inspecting the existing 0.10.0 distribution. The current integrated version is read from root `VERSION`; development builds do not advertise stable acceptance. Automation configuration/schema versions remain 2; algorithm and chemical-standardization versions do not change.

`python deployment/sync_platform_version.py` derives component Python metadata and the static site's package/version JSON from `VERSION`. Generated copies exist only to preserve independent packaging and deployment. `--check` rejects drift. For example, Setuptools may normalize the current SemVer development identifier to PEP 440 `0.11.0.dev2` in wheel metadata; it represents the same version. Live UI/API/CLI import the authority. Old release narratives, upstream engine versions and archived results are not rewritten.

Changing only UI code does not change protocol IDs. Scientific parameter changes produce a custom protocol identity until an explicitly reviewed revision is registered. A model filename or declared Git commit is not a verified installed model version. Missing provenance is `unknown`; missing Git dirty information is null, never falsely clean.
