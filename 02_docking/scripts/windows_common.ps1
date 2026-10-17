$ErrorActionPreference='Stop'
$script:ProjectRoot=Split-Path -Parent $PSScriptRoot
$script:Project=Get-Content -LiteralPath (Join-Path $script:ProjectRoot 'project.json') -Raw -Encoding UTF8 | ConvertFrom-Json
function Find-ProjectDistro {
 $expected=[IO.Path]::GetFullPath((Join-Path $script:ProjectRoot 'runtime\wsl\ext4.vhdx'))
 $matches=@()
 foreach($key in Get-ChildItem 'HKCU:\Software\Microsoft\Windows\CurrentVersion\Lxss'){
  $entry=Get-ItemProperty $key.PSPath
  if(!$entry.BasePath){continue}
  $base=$entry.BasePath -replace '^\\\?\\',''
  $vhd=if($entry.VhdFileName){$entry.VhdFileName}else{'ext4.vhdx'}
  $actual=[IO.Path]::GetFullPath((Join-Path $base $vhd))
  if($actual -ieq $expected){$matches+= $entry}
 }
 if($matches.Count -gt 1){throw 'Ambiguous project VHDX registration'}
 if($matches.Count -eq 0){return $null}
 $entry=$matches[0]
 if($entry.Version -ne 2){throw 'Project requires WSL2'}
 $id=(& wsl.exe -d $entry.DistributionName -u root --exec cat /opt/sbvs/project_uuid)
 if($LASTEXITCODE -ne 0 -or $id.Trim() -ne $script:Project.uuid){throw 'VHDX project UUID mismatch'}
 return $entry.DistributionName
}
function Get-LinuxRoot($distro){
 $p=(& wsl.exe -d $distro -u sbvs --exec wslpath -a -u $script:ProjectRoot.Replace('\','/'))
 if($LASTEXITCODE -ne 0){throw 'wslpath failed'}
 return $p.Trim()
}
function Invoke-Sbvs($distro,$arguments){
 $linux=Get-LinuxRoot $distro
 & wsl.exe -d $distro -u sbvs --exec env ('SBVS_PROJECT_ROOT='+$linux) bash ($linux+'/run.sh') @arguments
 if($LASTEXITCODE -ne 0){throw ('SBVS exit '+$LASTEXITCODE)}
}
function Test-ProjectHealth($port=5902){
 try{$h=Invoke-RestMethod -Uri ('http://127.0.0.1:'+$port+'/api/health') -TimeoutSec 3;return ($h.project_uuid -eq $script:Project.uuid -and $h.ready)}catch{return $false}
}
