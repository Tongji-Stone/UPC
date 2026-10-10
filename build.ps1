param([switch]$SkipCalculations)
$ErrorActionPreference = 'Stop'
$python = Join-Path $PSScriptRoot '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $python)) { $python = 'python' }
if ($SkipCalculations) { & $python (Join-Path $PSScriptRoot 'build_paper.py') }
else { & $python (Join-Path $PSScriptRoot 'reproduce.py') --compile }
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
& $python (Join-Path $PSScriptRoot 'verify_translation.py')
exit $LASTEXITCODE
