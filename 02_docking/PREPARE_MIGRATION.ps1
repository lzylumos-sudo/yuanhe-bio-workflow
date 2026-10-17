param([switch]$CancelActive)
. (Join-Path $PSScriptRoot 'scripts\windows_common.ps1')
$distro=Find-ProjectDistro
if(!$distro){throw 'Project registration not found'}
$linux=Get-LinuxRoot $distro
$args=@('stop')
if($CancelActive){$args+='--cancel-active'}
& wsl.exe -d $distro -u sbvs --exec env ('SBVS_PROJECT_ROOT='+$linux) /opt/sbvs/env/bin/python ($linux+'/scripts/service_control.py') @args
if($LASTEXITCODE -ne 0){throw 'Active tasks remain. Wait for completion or rerun with explicit -CancelActive'}
& wsl.exe -d $distro -u sbvs --exec sync
if($LASTEXITCODE -ne 0){throw 'Sync failed'}
& wsl.exe --terminate $distro
if($LASTEXITCODE -ne 0){throw 'Project terminate failed'}
Write-Host 'Only this project has stopped. Other programs or MD tasks may still use this disk; check them before unplugging.'
