function Get-PbeTabiTreecodeSources {
    param([Parameter(Mandatory = $true)][string]$RepoRoot)

    @(
        (Join-Path $RepoRoot 'titan_pbe\tabi_treecode.cpp')
        (Join-Path $RepoRoot 'titan_pbe\tabi_treecode_internal.cpp')
        (Join-Path $RepoRoot 'titan_pbe\tabi_treecode_geometry.cpp')
        (Join-Path $RepoRoot 'titan_pbe\tabi_treecode_series.cpp')
        (Join-Path $RepoRoot 'titan_pbe\tabi_treecode_policy.cpp')
        (Join-Path $RepoRoot 'titan_pbe\tabi_treecode_leaf.cpp')
        (Join-Path $RepoRoot 'titan_pbe\tabi_treecode_operator.cpp')
    )
}
