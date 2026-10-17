# ZLIS launcher (Windows PowerShell). UNTESTED on a real Windows machine in this release.
$ErrorActionPreference = "Stop"
$Here = Split-Path -Parent $MyInvocation.MyCommand.Path
$EnvName = if ($env:ZLIS_ENV) { $env:ZLIS_ENV } else { "zlis" }
$roots = @($env:CONDA_PREFIX, "$env:USERPROFILE\miniconda3", "$env:USERPROFILE\miniforge3", "$env:USERPROFILE\anaconda3",
           "$env:LOCALAPPDATA\miniconda3", "C:\ProgramData\miniconda3", "C:\ProgramData\miniforge3")
$py = $null
foreach ($r in $roots) {
  if ($r -and (Test-Path "$r\envs\$EnvName\python.exe")) { $py = "$r\envs\$EnvName\python.exe"; break }
}
if (-not $py) {
  $conda = Get-Command conda -ErrorAction SilentlyContinue
  if ($conda) {
    $base = (& conda info --base).Trim()
    if (Test-Path "$base\envs\$EnvName\python.exe") { $py = "$base\envs\$EnvName\python.exe" }
  }
}
if (-not $py) { Write-Error "ZLIS environment '$EnvName' not found. Run env\setup_windows.ps1"; exit 2 }
$env:PYTHONPATH = "$Here\src"
$env:PYTHONWARNINGS = "ignore"
$env:PYTHONUTF8 = "1"
& $py -m zlis @args
exit $LASTEXITCODE
