param([Parameter(ValueFromRemainingArguments=$true)][string[]]$CommandArgs)
. (Join-Path $PSScriptRoot 'scripts\windows_common.ps1')
$distro=Find-ProjectDistro
if(!$distro){throw 'Run REGISTER_ON_THIS_PC.ps1 first'}
Invoke-Sbvs $distro $CommandArgs
