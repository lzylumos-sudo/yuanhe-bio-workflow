. (Join-Path $PSScriptRoot 'scripts\windows_common.ps1')
$distro=Find-ProjectDistro
if(!$distro){throw 'Project registration not found'}
& (Join-Path $PSScriptRoot 'PREPARE_MIGRATION.ps1')
$target=Join-Path $PSScriptRoot ('runtime/downloads/runtime_snapshot_'+(Get-Date -Format 'yyyyMMdd_HHmmss')+'.tar')
& wsl.exe --export $distro $target
if($LASTEXITCODE -ne 0){throw 'Project snapshot export failed'}
@{project_uuid=$script:Project.uuid;file=(Split-Path $target -Leaf);sha256=(Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash;created=(Get-Date -Format o)}|ConvertTo-Json|Set-Content -Encoding UTF8 ($target+'.json')
