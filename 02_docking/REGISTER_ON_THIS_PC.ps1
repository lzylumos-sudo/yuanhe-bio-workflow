. (Join-Path $PSScriptRoot 'scripts\windows_common.ps1')
Start-Transcript -Path (Join-Path $PSScriptRoot 'logs/install/register.log') -Append
if($env:PROCESSOR_ARCHITECTURE -ne 'AMD64'){throw 'Windows x86-64 required'}
$vhd=Join-Path $PSScriptRoot 'runtime/wsl/ext4.vhdx'
if(!(Test-Path -LiteralPath $vhd)){throw 'Project ext4.vhdx missing'}
$distro=Find-ProjectDistro
if(!$distro){
 $help=(& wsl.exe --help | Out-String)
 if($help -notmatch 'import-in-place'){throw 'Update host WSL to support import-in-place; runtime does not need reinstall'}
 $names=@(Get-ChildItem 'HKCU:\Software\Microsoft\Windows\CurrentVersion\Lxss' | ForEach-Object {(Get-ItemProperty $_.PSPath).DistributionName})
 $distro='GNINA-SBVS-'+$script:Project.uuid.Substring(0,8)
 if($names -contains $distro){$distro+='-'+[guid]::NewGuid().ToString().Substring(0,6)}
 & wsl.exe --import-in-place $distro $vhd
 if($LASTEXITCODE -ne 0){throw 'Registration failed; existing distributions preserved'}
 $verified=Find-ProjectDistro
 if($verified -ne $distro){throw 'Registered VHDX identity verification failed'}
}
@{distro=$distro;project_uuid=$script:Project.uuid;vhdx=[IO.Path]::GetFullPath($vhd);registered=(Get-Date -Format o)}|ConvertTo-Json|Set-Content -Encoding UTF8 (Join-Path $PSScriptRoot 'runtime/local_machine.json')
Invoke-Sbvs $distro @('doctor')
$linux=Get-LinuxRoot $distro
& wsl.exe -d $distro -u sbvs --exec env ('SBVS_PROJECT_ROOT='+$linux) ('PYTHONPATH='+$linux+'/src') /opt/sbvs/env/bin/python ($linux+'/scripts/cco_smoke.py')
if($LASTEXITCODE -ne 0){throw 'Migration CCO smoke failed'}
Stop-Transcript
