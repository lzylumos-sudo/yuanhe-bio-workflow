param([switch]$InstallGnina)
$ErrorActionPreference='Stop'
if($env:PROCESSOR_ARCHITECTURE -ne 'AMD64'){throw 'Windows x86-64 required'}
if(!(Test-Path (Join-Path $PSScriptRoot 'project.json'))){throw 'Project identity missing; do not install over unknown directories'}
$p=Get-Content (Join-Path $PSScriptRoot 'project.json') -Raw -Encoding UTF8|ConvertFrom-Json
Start-Transcript -Path (Join-Path $PSScriptRoot 'logs/install/setup_windows.log') -Append
& wsl.exe --version
if($LASTEXITCODE -ne 0){throw 'Working WSL2 host required'}
$drive=[IO.DriveInfo]::new([IO.Path]::GetPathRoot($PSScriptRoot))
if($drive.AvailableFreeSpace -lt 15GB){throw 'At least 15 GB free space required'}
$probe=Join-Path $PSScriptRoot 'runtime/write_probe.tmp';[IO.File]::WriteAllText($probe,'ok');Remove-Item -LiteralPath $probe
Get-NetTCPConnection -LocalPort 5902 -State Listen -ErrorAction SilentlyContinue | Select-Object LocalAddress,OwningProcess
$vhd=Join-Path $PSScriptRoot 'runtime/wsl/ext4.vhdx'
if(Test-Path $vhd){throw 'Existing runtime preserved. Use REGISTER_ON_THIS_PC.ps1; installation/upgrade is not automatic'}
$base='https://cdimage.ubuntu.com/ubuntu-base/releases/24.04/release/';$file='ubuntu-base-24.04.5-base-amd64.tar.gz';$dl=Join-Path $PSScriptRoot 'runtime/downloads'
& curl.exe -fL --retry 2 --max-time 120 ($base+'SHA256SUMS') -o (Join-Path $dl 'SHA256SUMS')
if($LASTEXITCODE -ne 0){throw 'Official checksum download failed'}
$lines=@(Get-Content (Join-Path $dl 'SHA256SUMS')|Where-Object{$_ -match ([regex]::Escape($file)+'$')})
if($lines.Count -ne 1){throw 'Pinned official rootfs missing'}
$hash=($lines[0] -split '\s+')[0];$tar=Join-Path $dl $file
if(!(Test-Path $tar)){
 & curl.exe -fL --retry 2 --connect-timeout 20 --max-time 600 -C - ($base+$file) -o ($tar+'.part')
 if($LASTEXITCODE -ne 0){throw 'Rootfs download failed'}
 if((Get-Item ($tar+'.part')).Length -lt 10000000 -or (Get-FileHash ($tar+'.part') -Algorithm SHA256).Hash.ToLower() -ne $hash){throw 'Rootfs verification failed'}
 Move-Item -LiteralPath ($tar+'.part') -Destination $tar
}
if((Get-FileHash $tar -Algorithm SHA256).Hash.ToLower() -ne $hash){throw 'Stored rootfs checksum mismatch'}
@{url=($base+$file);sha256=$hash;publisher_hash_verified=$true}|ConvertTo-Json|Set-Content -Encoding UTF8 (Join-Path $PSScriptRoot 'env/rootfs.lock.json')
$names=@(Get-ChildItem 'HKCU:\Software\Microsoft\Windows\CurrentVersion\Lxss'|ForEach-Object{(Get-ItemProperty $_.PSPath).DistributionName})
if($names -contains $p.distro){throw 'Distro name already registered; existing registration preserved'}
& wsl.exe --import $p.distro (Split-Path $vhd -Parent) $tar --version 2
if($LASTEXITCODE -ne 0){throw 'WSL2 import failed'}
$linux=(& wsl.exe -d $p.distro -u root --exec wslpath -u $PSScriptRoot.Replace('\','/')).Trim()
& wsl.exe -d $p.distro -u root --exec bash ($linux+'/scripts/bootstrap_runtime.sh') $linux $p.uuid
if($LASTEXITCODE -ne 0){throw 'Runtime bootstrap failed; inspect logs before resuming'}
& wsl.exe -d $p.distro -u sbvs --exec env MAMBA_ROOT_PREFIX=/opt/sbvs/mamba /opt/sbvs/vendor/micromamba install -y -p /opt/sbvs/env -c conda-forge pdbfixer=1.12
if($LASTEXITCODE -ne 0){throw 'PDBFixer installation failed'}
& wsl.exe -d $p.distro -u sbvs --exec env ('SBVS_PROJECT_ROOT='+$linux) /opt/sbvs/env/bin/python ($linux+'/scripts/install_vendor.py') vina
if($LASTEXITCODE -ne 0){throw 'Vina installation failed'}
& wsl.exe -d $p.distro -u sbvs --exec bash ($linux+'/scripts/refresh_locks.sh') $linux
& wsl.exe -d $p.distro -u sbvs --exec env ('SBVS_PROJECT_ROOT='+$linux) ('PYTHONPATH='+$linux+'/src') /opt/sbvs/env/bin/python ($linux+'/scripts/install_demo.py')
if($LASTEXITCODE -ne 0){throw 'Classroom setup/smoke failed; not ready'}
if($InstallGnina){
 & wsl.exe -d $p.distro -u sbvs --exec bash ($linux+'/scripts/gnina_dependencies.sh')
 if($LASTEXITCODE -ne 0){throw 'GNINA runtime library install failed'}
 & wsl.exe -d $p.distro -u sbvs --exec env ('SBVS_PROJECT_ROOT='+$linux) /opt/sbvs/env/bin/python ($linux+'/scripts/install_vendor.py') gnina
 if($LASTEXITCODE -ne 0){throw 'GNINA failed; classroom remains available'}
}
Stop-Transcript
