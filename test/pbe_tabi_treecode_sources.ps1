function Get-PbeTabiTreecodeSources {
    param([Parameter(Mandatory = $true)][string]$RepoRoot)

    @(
        (Join-Path $RepoRoot 'src\tabi\tabi_treecode.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_treecode_internal.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_treecode_geometry.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_treecode_series.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_treecode_policy.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_treecode_leaf.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_treecode_operator.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_fmm_m2l_prototype.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_fmm_m2l_fallback.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_fmm_m2l_cartesian.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_fmm_m2l_local.cpp')
    )
}
