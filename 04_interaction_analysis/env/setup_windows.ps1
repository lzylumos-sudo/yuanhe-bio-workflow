# Creates the 'zlis' conda environment on Windows. UNTESTED on a real Windows machine in this release.
$Here = Split-Path -Parent $MyInvocation.MyCommand.Path
if (-not (Get-Command conda -ErrorAction SilentlyContinue)) { Write-Error "conda not found. Install Miniforge first."; exit 2 }
conda env create -f "$Here\environment.yml" --yes
if ($LASTEXITCODE -ne 0) { conda env update -n zlis -f "$Here\environment.yml" --prune }
Write-Host "Done. Start with Start_Windows.cmd"
