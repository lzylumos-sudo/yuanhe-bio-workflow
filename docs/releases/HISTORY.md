# Version history

For new installations, use the [current administrator download](../../DEPLOY.md). This page distinguishes historical software identities from the current build. Old scientific results keep their original parameters, source hashes and acceptance status.

| Identity | Role | Deployment status |
|---|---|---|
| Version in root [VERSION](../../VERSION) | Current integrated administrator build; [release metadata](../../deployment/current_release.json) | Use the current download entrypoint; target-machine acceptance required |
| 0.11.0-dev.2 | Earlier package at `d0bdfd51b936e8b16422e4bf2dae898897ad8c0d`; later main commits retained this label | Historical tag; lacks subsequent synthesis-only deployment integration |
| 0.11.0-dev.1 | Original architecture-governance branch, source `a6bbe7a2399b0603894dcd6feb97f86efa61b512` | Historical development snapshot, superseded by the integrated build |
| 0.4.0-rc1 | Legacy workbench/UI identity, including deployment fix source `dc7c2251b32afb1864188db823e2f0a8629a4d64` | Historical main-branch snapshot; no longer the current package |
| 0.10.0 / 0.1.0 / 2.0.0 / 1.0.0 | Former screening / docking / standalone automation / public-site package identities | Historical component labels, not competing current platform releases |

The old labels above were not repository release tags when version governance was introduced. Do not infer a downloadable release or scientific qualification from a label alone. The integrated line uses one software authority; scientific protocols, model/data versions, schema versions, Git commits and Run IDs remain separate identities.

## Historical acceptance records

- [Integrated acceptance](../10_INTEGRATED_ACCEPTANCE.md): original target observations, including incomplete repeat/export acceptance.
- [Release acceptance and maintenance](../11_RELEASE_MAINTENANCE.md): A1-A5 evidence requirements; historical observations are not revalidated by this integration.
- [Original governance validation](../validation/refactor_validation_report.md): tests for the original development branch, not substituted for current CI.
- [Public minimal GPU acceptance](../../deployment/public_minimal/README.md): bounded CLI docking/OpenMM evidence; not full browser/model/membrane acceptance.
- [Changelog](../../CHANGELOG.md): changes by build and historical maintenance date.

Git history and original source hashes are preserved. Repeated file-list commit messages do not indicate duplicate releases, and changing those messages would change historical commit identities. New deployment builds have explicit versioned tags, release notes, source manifests and checksums.
