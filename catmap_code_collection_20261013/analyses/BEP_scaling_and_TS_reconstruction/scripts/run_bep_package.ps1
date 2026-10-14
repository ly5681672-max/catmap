$ErrorActionPreference = 'Stop'
$python = 'G:\anaconda\python.exe'
$root = Split-Path -Parent $PSScriptRoot
& $python (Join-Path $PSScriptRoot 'verify_catmap_bep.py')
if ($LASTEXITCODE -ne 0) { throw 'Official CatMAP BEP scaler verification failed.' }
& $python (Join-Path $PSScriptRoot 'analyze_bep_reconstruction.py')
if ($LASTEXITCODE -ne 0) { throw 'BEP analysis failed.' }
& $python (Join-Path $PSScriptRoot 'render_bep_figures.py')
if ($LASTEXITCODE -ne 0) { throw 'Figure rendering failed.' }
& $python (Join-Path $PSScriptRoot 'validate_bep_package.py')
if ($LASTEXITCODE -ne 0) { throw 'Package validation failed.' }
& $python (Join-Path $PSScriptRoot 'write_bep_reports.py')
if ($LASTEXITCODE -ne 0) { throw 'Final report/hash generation failed.' }
