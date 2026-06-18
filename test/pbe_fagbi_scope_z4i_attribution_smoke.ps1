param(
    [string]$Cxx = '',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4i_unit')
)

$ErrorActionPreference = 'Stop'
$RepoRoot = Split-Path -Parent $PSScriptRoot
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
if ([string]::IsNullOrWhiteSpace($Cxx)) {
    $clang = 'C:\msys64\ucrt64\bin\clang++.exe'
    $Cxx = if (Test-Path -LiteralPath $clang) { $clang } else { 'g++' }
}

$exe = Join-Path $WorkRoot 'titan_fmm_m2l_scope_z4i_attribution_smoke.exe'
$sources = @(
    'test\titan_fmm_m2l_scope_z4i_attribution_smoke.cpp',
    'src\tabi\tabi_fmm_m2l_prototype.cpp',
    'src\titan_executor.cpp',
    'src\tabi\tabi_near.cpp',
    'src\tabi\tabi_treecode.cpp',
    'src\tabi\tabi_treecode_internal.cpp',
    'src\tabi\tabi_treecode_geometry.cpp',
    'src\tabi\tabi_treecode_series.cpp',
    'src\tabi\tabi_treecode_policy.cpp',
    'src\tabi\tabi_treecode_leaf.cpp',
    'src\tabi\tabi_treecode_operator.cpp'
)
$buildArgs = @('-std=c++23','-O2','-mavx2','-mfma',
    '-DTITAN_RUNTIME_DISPATCH_HOTSPOTS=1','-Isrc') + $sources + @('-o',$exe)

Push-Location $RepoRoot
try {
    & $Cxx @buildArgs
    if ($LASTEXITCODE -ne 0) { throw 'failed to build Z4I attribution smoke' }
    & $exe
    if ($LASTEXITCODE -ne 0) { throw 'Z4I attribution smoke failed' }
} finally {
    Pop-Location
}
