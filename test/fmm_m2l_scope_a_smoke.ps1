param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_a'),
    [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
$RepoRoot = Split-Path -Parent $PSScriptRoot
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

if ([string]::IsNullOrWhiteSpace($Cxx)) {
    $Clang = 'C:\msys64\ucrt64\bin\clang++.exe'
    if (Test-Path $Clang) {
        $Cxx = $Clang
    } else {
        $Cxx = 'g++'
    }
}

$ExePath = Join-Path $WorkRoot 'titan_fmm_m2l_prototype_smoke.exe'
$Args = @(
    '-std=c++23',
    '-O2', '-mavx2', '-mfma',
    '-Isrc',
    'test\titan_fmm_m2l_prototype_smoke.cpp',
    'src\tabi\tabi_fmm_m2l_prototype.cpp',
    'src\tabi\tabi_fmm_m2l_fallback.cpp',
    'src\tabi\tabi_fmm_m2l_cartesian.cpp',
    'src\tabi\tabi_fmm_m2l_local.cpp',
    'src\tabi\tabi_near.cpp',
    'src\tabi\tabi_treecode.cpp',
    'src\tabi\tabi_treecode_internal.cpp',
    'src\tabi\tabi_treecode_geometry.cpp',
    'src\tabi\tabi_treecode_series.cpp',
    'src\tabi\tabi_treecode_policy.cpp',
    'src\tabi\tabi_treecode_leaf.cpp',
    'src\tabi\tabi_treecode_operator.cpp',
    'src\titan_executor.cpp',
    '-o',
    $ExePath
)

Push-Location $RepoRoot
try {
    & $Cxx @Args
    if ($LASTEXITCODE -ne 0) {
        throw "failed to build titan_fmm_m2l_prototype_smoke"
    }

    & $ExePath
    if ($LASTEXITCODE -ne 0) {
        throw "titan_fmm_m2l_prototype_smoke failed"
    }
} finally {
    Pop-Location
}
