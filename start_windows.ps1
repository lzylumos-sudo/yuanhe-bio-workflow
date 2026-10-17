param([string]$Distro = $env:YUANHE_WSL_DISTRO, [string]$Python = '')
$ErrorActionPreference = 'Stop'
if (-not (Get-Command wsl.exe -ErrorAction SilentlyContinue)) { throw 'WSL2 is not installed.' }
if (-not $Distro) {
  Write-Host 'Available WSL distributions:'
  wsl.exe --list --verbose
  $Distro = Read-Host 'Enter the existing distribution name (no distribution will be installed or changed)'
}
if ([string]::IsNullOrWhiteSpace($Distro)) { throw 'A WSL distribution name is required.' }
$LinuxRoot = (& wsl.exe -d $Distro -- wslpath -u $PSScriptRoot).Trim()
if ($LASTEXITCODE -ne 0 -or -not $LinuxRoot) { throw 'Cannot translate the project path in this distribution.' }
Write-Host 'Open http://localhost:8777 in your Windows browser after the server starts.'
Write-Host 'Ctrl+C stops the UI. Running scientific jobs must be stopped from the UI first.'
& wsl.exe -d $Distro -- bash -lc 'set -e; cd -- "$1" && export YUANHE_UI_PYTHON="$2"; exec bash start_linux.sh' yuanhe $LinuxRoot $Python
exit $LASTEXITCODE
