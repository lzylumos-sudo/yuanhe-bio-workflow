param([int]$Port=5902,[switch]$NoBrowser)
. (Join-Path $PSScriptRoot 'scripts\windows_common.ps1')
$distro=Find-ProjectDistro
if(!$distro){throw 'Run REGISTER_ON_THIS_PC.ps1 first'}
if(Test-ProjectHealth $Port){if(!$NoBrowser){Start-Process ('http://localhost:'+$Port+'/')};exit 0}
$connections=@(Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue)
if($connections.Count){$connections|Select-Object LocalAddress,LocalPort,OwningProcess|Format-Table;throw 'Port occupied; no process was stopped'}
$linux=Get-LinuxRoot $distro
& wsl.exe -d $distro -u sbvs --exec env ('SBVS_PROJECT_ROOT='+$linux) /opt/sbvs/env/bin/python ($linux+'/scripts/service_control.py') start --port $Port
if($LASTEXITCODE -ne 0){throw 'Service start failed'}
for($i=0;$i -lt 60;$i++){if(Test-ProjectHealth $Port){if(!$NoBrowser){Start-Process ('http://localhost:'+$Port+'/')};exit 0};Start-Sleep -Seconds 1}
throw 'Health check timed out. See logs/web.log'
