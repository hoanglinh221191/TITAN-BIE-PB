[CmdletBinding()]
param(
    [string]$TitanExe = '.\titan.exe',
    [string]$RepoRoot = '',
    [string]$FixtureRoot = '',
    [string]$WorkRoot = '',
    [string]$ReportPath = '',
    [string[]]$Fixtures = @('barnase','1a63','1sh1'),
    [ValidateSet('sphere-patch','grid-ses','analytic-ses')][string]$PrimarySurface = 'sphere-patch',
    [ValidateRange(0.25,4.0)][double]$SesGridSpacing = 1.25,
    [ValidateRange(1,86400)][int]$CaseTimeoutSeconds = 1800,
    [ValidateRange(256,8192)][int]$RasMemoryBudgetMb = 1024,
    [switch]$Quick,
    [switch]$Force,
    [switch]$RerunFailed,
    [string[]]$RerunStages = @(),
    [switch]$RecoverVerifiedProfiles,
    [switch]$RecoveryOnly,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$Inv = [Globalization.CultureInfo]::InvariantCulture
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

function Split-Values([string[]]$Values) {
    return @($Values | ForEach-Object { ([string]$_).Split(',') } |
        ForEach-Object { $_.Trim().ToLowerInvariant() } |
        Where-Object { -not [string]::IsNullOrWhiteSpace($_) })
}
function D($Value,[double]$Fallback=[double]::NaN) {
    $parsed = 0.0
    if ($null -ne $Value -and [double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float,$Inv,[ref]$parsed)) { return $parsed }
    return $Fallback
}
function I($Value,[int64]$Fallback=0) {
    $parsed = 0L
    if ($null -ne $Value -and [int64]::TryParse([string]$Value,[ref]$parsed)) { return $parsed }
    return $Fallback
}
function P($Object,[string]$Name,$Default='') {
    if ($null -eq $Object) { return $Default }
    $property = $Object.PSObject.Properties[$Name]
    if ($null -eq $property) { return $Default }
    return $property.Value
}
function Xml([string]$Value) { return [Security.SecurityElement]::Escape($Value) }
function Get-SolverContractHash([string]$Root) {
    $relativePaths = @(
        'src/tabi/tabi_solver.cpp',
        'src/tabi/tabi_surface.cpp',
        'src/tabi/tabi_analytic_ses.cpp',
        'src/tabi/tabi_fmm_m2l_preconditioner.cpp',
        'src/tabi/tabi_fmm_m2l_prototype.cpp',
        'src/tabi/tabi_fmm_m2l_fallback.cpp',
        'src/tabi/tabi_preconditioner.cpp',
        'test/pbe_fmm_m2l_cli_helpers.ps1'
    )
    $parts = foreach ($relative in $relativePaths) {
        $path = Join-Path $Root ($relative -replace '/', [IO.Path]::DirectorySeparatorChar)
        if (-not (Test-Path -LiteralPath $path)) { throw "Missing solver-contract file: $path" }
        "$relative=$((Get-FileHash -Algorithm SHA256 -LiteralPath $path).Hash.ToLowerInvariant())"
    }
    $sha = [Security.Cryptography.SHA256]::Create()
    try {
        $bytes = [Text.Encoding]::UTF8.GetBytes(($parts -join "`n"))
        return ([BitConverter]::ToString($sha.ComputeHash($bytes)) -replace '-','').ToLowerInvariant()
    } finally { $sha.Dispose() }
}

if ([string]::IsNullOrWhiteSpace($RepoRoot)) { $RepoRoot = Split-Path -Parent $PSScriptRoot }
if ([string]::IsNullOrWhiteSpace($FixtureRoot)) { $FixtureRoot = Join-Path $PSScriptRoot 'pbe_inputs_public' }
if ([string]::IsNullOrWhiteSpace($WorkRoot)) { $WorkRoot = Join-Path $PSScriptRoot 'tmp\fagbi_graph_preconditioner_ablation' }
$RepoRoot = [IO.Path]::GetFullPath($RepoRoot)
$FixtureRoot = [IO.Path]::GetFullPath($FixtureRoot)
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
if (-not (Test-Path -LiteralPath $manifestPath)) { throw "Missing fixture manifest: $manifestPath" }
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$manifest = @(Import-Csv -LiteralPath $manifestPath)
$fixtureNames = Split-Values $Fixtures
if ($Quick) { $fixtureNames = @('1sh1') }
$rerunLabels = Split-Values $RerunStages
if ($RecoveryOnly -and -not $RecoverVerifiedProfiles) {
    throw '-RecoveryOnly requires -RecoverVerifiedProfiles.'
}

$matrixPath = Join-Path $WorkRoot 'fagbi_graph_preconditioner_matrix.csv'
$tracePath = Join-Path $WorkRoot 'fagbi_graph_preconditioner_residual_trace.csv'
$decisionPath = Join-Path $WorkRoot 'fagbi_graph_preconditioner_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fagbi_graph_preconditioner_summary.dat'
$plotPath = Join-Path $WorkRoot 'fagbi_graph_preconditioner_convergence.svg'
if ([string]::IsNullOrWhiteSpace($ReportPath)) {
    $reportName = if ($PrimarySurface -eq 'analytic-ses') {
        'TITAN_FAGBI_ANALYTIC_SES_GRAPH_PRECONDITIONER_REPORT.txt'
    } else {
        'TITAN_FAGBI_GRAPH_PRECONDITIONER_ABLATION_REPORT.txt'
    }
    $ReportPath = Join-Path $RepoRoot $reportName
}
$reportPath = [IO.Path]::GetFullPath($ReportPath)
$titanSha = (Get-FileHash -Algorithm SHA256 -LiteralPath $TitanExe).Hash.ToLowerInvariant()
$solverContractSha = Get-SolverContractHash $RepoRoot
$protocol = 'surface_points=8;theta=0.45;tol=1e-4;restart=128;maxiter=720;nt=4;salt=0.15;near=analytic_singular;graph_cap=64'
if ($PrimarySurface -in @('grid-ses','analytic-ses')) {
    $protocol += ";surface=$PrimarySurface;ses_spacing=$($SesGridSpacing.ToString('G15',$Inv))"
}

$stages = @(
    [pscustomobject]@{label='diagonal'; precond='diagonal'; expected='diagonal'; block=32; overlap=64; ilu=64},
    [pscustomobject]@{label='leaf_block_jacobi'; precond='fmm-block-jacobi'; expected='fmm_block_jacobi'; block=64; overlap=64; ilu=64},
    [pscustomobject]@{label='fmm_near_graph_ras'; precond='fmm-overlap-schwarz'; expected='fmm_overlap_schwarz'; block=32; overlap=64; ilu=64},
    [pscustomobject]@{label='distance_knn_graph_ras'; precond='distance-knn-overlap-schwarz'; expected='distance_knn_overlap_schwarz'; block=32; overlap=64; ilu=64},
    [pscustomobject]@{label='distance_knn_ilu0'; precond='distance-knn-ilu0'; expected='distance_knn_ilu0'; block=32; overlap=64; ilu=64},
    [pscustomobject]@{label='fmm_near_graph_ilu0'; precond='fmm-nearfield-ilu0'; expected='fmm_nearfield_ilu0'; block=32; overlap=64; ilu=64}
)

$existing = @{}
if (Test-Path -LiteralPath $matrixPath) {
    foreach ($row in @(Import-Csv -LiteralPath $matrixPath)) { $existing[[string]$row.row_id] = $row }
}
$rows = [Collections.Generic.List[object]]::new()
$plannedIds = @{}
foreach ($fixtureName in $fixtureNames) {
    foreach ($stage in $stages) {
        $surfaceToken = if ($PrimarySurface -in @('grid-ses','analytic-ses')) {
            "$(($PrimarySurface -replace '-','_'))_h$(([string]$SesGridSpacing).Replace('.','p'))_"
        } else { '' }
        $plannedIds["$fixtureName`_$surfaceToken$($stage.label)"] = $true
    }
}
foreach ($oldRow in $existing.Values) {
    if (-not $plannedIds.ContainsKey([string]$oldRow.row_id)) {
        $rows.Add($oldRow)
    }
}
foreach ($fixtureName in $fixtureNames) {
    $fixtureRows = @($manifest | Where-Object fixture -eq $fixtureName)
    if ($fixtureRows.Count -ne 1) { throw "Expected one public fixture row for '$fixtureName'" }
    $fixture = $fixtureRows[0]
    $pqr = Join-Path $FixtureRoot $fixture.pqr
    if (-not (Test-Path -LiteralPath $pqr)) { throw "Missing fixture PQR: $pqr" }
    $pqrSha = (Get-FileHash -Algorithm SHA256 -LiteralPath $pqr).Hash.ToLowerInvariant()
    foreach ($stage in $stages) {
        $isRas = [string]$stage.label -in @('fmm_near_graph_ras','distance_knn_graph_ras')
        $stageMemoryBudgetMb = if ($isRas) { $RasMemoryBudgetMb } else { 256 }
        $stageProtocol = if ($isRas) {
            "$protocol;ras_memory_budget_mb=$stageMemoryBudgetMb"
        } else { $protocol }
        $surfaceToken = if ($PrimarySurface -in @('grid-ses','analytic-ses')) {
            "$(($PrimarySurface -replace '-','_'))_h$(([string]$SesGridSpacing).Replace('.','p'))_"
        } else { '' }
        $rowId = "$fixtureName`_$surfaceToken$($stage.label)"
        $runDir = Join-Path $WorkRoot $rowId
        $old = $existing[$rowId]
        if ($RecoverVerifiedProfiles) {
            $recoveryProfilePath = Join-Path $runDir 'pbe_tabi_profile.csv'
            $recoveryArgPath = Join-Path $runDir 'arg.in'
            if ((Test-Path -LiteralPath $recoveryProfilePath) -and
                (Test-Path -LiteralPath $recoveryArgPath) -and
                (Get-Item -LiteralPath $recoveryProfilePath).LastWriteTimeUtc -ge
                    (Get-Item -LiteralPath $TitanExe).LastWriteTimeUtc) {
                $recoveryProfile = @(Import-Csv -LiteralPath $recoveryProfilePath)[-1]
                $recoveryStatus = [string](P $recoveryProfile 'status' 'failed')
                $recoveryIterations = I (P $recoveryProfile 'iterations')
                $recoveryEffective = [string](P $recoveryProfile 'preconditioner_effective' '')
                $normalizedPqr = $pqr.Replace('\','/')
                $requiredArgLines = @(
                    "-pbe-test-pqr=$normalizedPqr",
                    '-pbe=tabi',
                    '-pbe-max-iter=720',
                    '-pbe-rms=0.0001',
                    '-pbe-dphi=0.0001',
                    '-nt=4',
                    '-pbe-tabi-sum=treecode',
                    '-pbe-tabi-fmm-m2l=single-matvec-analytic',
                    '-pbe-tabi-surface-points=8',
                    "-pbe-tabi-primary-surface=$PrimarySurface",
                    "-pbe-tabi-ses-grid-spacing=$($SesGridSpacing.ToString('G15',$Inv))",
                    '-pbe-tabi-opening-theta=0.45',
                    '-pbe-tabi-near-integration=analytic_singular',
                    '-pbe-tabi-taylor-order=4',
                    "-pbe-tabi-precond=$($stage.precond)",
                    "-pbe-tabi-precond-block-panel-cap=$($stage.block)",
                    "-pbe-tabi-precond-overlap-panel-cap=$($stage.overlap)",
                    "-pbe-tabi-precond-ilu-neighbor-cap=$($stage.ilu)",
                    '-pbe-tabi-precond-ilu-drop-tol=0',
                    "-pbe-tabi-precond-memory-budget-mb=$stageMemoryBudgetMb",
                    '-pbe-tabi-gmres-restart=128',
                    '-pbe-tabi-matvec-audit=sampled',
                    '-pbe-tabi-matvec-audit-sample-count=256'
                )
                $argLines = @(Get-Content -LiteralPath $recoveryArgPath)
                $argMatches = @($requiredArgLines | Where-Object { $argLines -notcontains $_ }).Count -eq 0
                $profileComplete = $recoveryStatus -in @('ok','failed') -and
                    $recoveryIterations -gt 0 -and $recoveryIterations -le 720 -and
                    $recoveryEffective -eq $stage.expected
                if ($argMatches -and $profileComplete) {
                    $recoveryMemoryPath = Join-Path $runDir 'pbe_memory_profile.csv'
                    $recoveryMemory = if (Test-Path -LiteralPath $recoveryMemoryPath) {
                        @(Import-Csv -LiteralPath $recoveryMemoryPath)[-1]
                    } else { $null }
                    $residual = D (P $recoveryProfile 'true_residual_rel' (P $recoveryProfile 'final_residual'))
                    $audit = D (P $recoveryProfile 'tree_vs_direct_matvec_rel' (P $recoveryProfile 'matvec_audit_rel_rms'))
                    $finite = (Test-FmmCliFinite $residual) -and (Test-FmmCliFinite $audit)
                    $failure = if ($recoveryIterations -ge 720 -and $finite) { 'not_converged' } elseif ($recoveryStatus -ne 'ok') { 'runtime_failure' } elseif (-not $finite -or $audit -gt 1e-2) { 'accuracy_gap' } elseif ($residual -gt 1e-4) { 'not_converged' } else { 'pass' }
                    $old = [pscustomobject]@{
                        row_id=$rowId; fixture=$fixtureName; protein_name=$fixture.protein_name; pdb_id=$fixture.pdb_id
                        surface_protocol=$PrimarySurface; ses_grid_spacing_A=$(if($PrimarySurface -in @('grid-ses','analytic-ses')){$SesGridSpacing}else{'NA'})
                        preconditioner_label=$stage.label; preconditioner_requested=$stage.precond
                        preconditioner_effective=$recoveryEffective; graph_cap=64; panels=I (P $recoveryProfile 'panels')
                        status=$recoveryStatus; failure_class=$failure; iterations=$recoveryIterations
                        true_residual=$residual; final_residual=D (P $recoveryProfile 'final_residual')
                        audit_rms=$audit; wall_seconds=D (P $recoveryProfile 'total_seconds')
                        peak_rss_mb=(D (P $recoveryMemory 'peak_working_set_bytes') 0.0)/1MB
                        setup_seconds=D (P $recoveryProfile 'preconditioner_setup_seconds') 0.0
                        apply_seconds=D (P $recoveryProfile 'preconditioner_apply_seconds') 0.0
                        matvec_seconds_avg=D (P $recoveryProfile 'matvec_seconds_avg')
                        nnz_blocks=I (P $recoveryProfile 'preconditioner_sparse_nnz_blocks')
                        avg_degree=D (P $recoveryProfile 'preconditioner_sparse_avg_degree')
                        max_degree=I (P $recoveryProfile 'preconditioner_sparse_max_degree')
                        valid_blocks=I (P $recoveryProfile 'preconditioner_valid_blocks')
                        valid_panels=I (P $recoveryProfile 'preconditioner_valid_panels')
                        quality_status=[string](P $recoveryProfile 'preconditioner_quality_status')
                        residual_history=[string](P $recoveryProfile 'residual_history')
                        protocol=$stageProtocol; solver_contract_sha256=$solverContractSha
                        titan_sha256=$titanSha; pqr_sha256=$pqrSha; run_dir=$runDir
                    }
                    Write-Host "recover verified: $rowId"
                }
            }
        }
        # A profile without a matrix row has no trustworthy executable/contract
        # provenance.  Do not recover it from timestamps alone.
        # Raising the RAS budget does not change a RAS matrix that was already
        # constructed successfully under the old 256-MiB cap. Reuse only those
        # exact-effective rows; fallback rows must be rebuilt with the new cap.
        $oldEffectiveExact = $null -ne $old -and
            [string](P $old 'preconditioner_effective') -eq [string]$stage.expected
        $protocolCompatible = $null -ne $old -and (
            [string](P $old 'protocol') -eq $stageProtocol -or
            ($isRas -and $oldEffectiveExact -and [string](P $old 'protocol') -eq $protocol)
        )
        $same = $null -ne $old -and
            [string](P $old 'titan_sha256') -eq $titanSha -and
            [string](P $old 'solver_contract_sha256') -eq $solverContractSha -and
            [string](P $old 'pqr_sha256') -eq $pqrSha -and
            $protocolCompatible -and
            (-not $isRas -or $oldEffectiveExact)
        $rerunStage = $rerunLabels -contains [string]$stage.label
        $reuse = -not $Force -and -not $rerunStage -and $same -and
            (-not $RerunFailed -or [string](P $old 'failure_class') -eq 'pass')
        if ($reuse) {
            Write-Host "resume: $rowId"
            $rows.Add($old)
            continue
        }
        if ($RecoveryOnly) {
            Write-Host "recovery skip (not verified complete): $rowId"
            continue
        }
        Write-Host "run: $rowId (timeout ${CaseTimeoutSeconds}s)"
        if (Test-Path -LiteralPath $runDir) { Remove-Item -LiteralPath $runDir -Recurse -Force }
        New-Item -ItemType Directory -Force -Path $runDir | Out-Null
        $arg = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot -PqrPath $pqr `
            -Backend tabi -Mode treecode -FmmMode single-matvec-analytic `
            -SurfacePoints 8 -MaxIter 720 -Restart 128 -Tolerance 1e-4 -Threads 4 `
            -SaltM 0.15 -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
            -OpeningTheta 0.45 -AuditSampleCount 256 -AuditMode sampled `
            -Preconditioner $stage.precond -PreconditionerBlockPanelCap $stage.block `
            -PreconditionerOverlapPanelCap $stage.overlap -PreconditionerIluNeighborCap $stage.ilu `
            -PreconditionerIluDropTol 0.0 -PreconditionerMemoryBudgetMb $stageMemoryBudgetMb `
            -NearIntegration analytic_singular `
            -TaylorOrder 4 -PrimarySurface $PrimarySurface -SesGridSpacing $SesGridSpacing
        $bare = @(Get-Content -LiteralPath $arg | Where-Object { $_ -match '^-[A-Za-z0-9-]+$' })
        if ($bare.Count -gt 0) { throw "Bare CLI options in $arg" }
        $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir -ArgPath $arg `
            -CaseTimeoutSeconds $CaseTimeoutSeconds
        $profilePath = Join-Path $runDir 'pbe_tabi_profile.csv'
        $profile = if (Test-Path -LiteralPath $profilePath) { @(Import-Csv -LiteralPath $profilePath)[-1] } else { $null }
        $status = if ($run.timed_out) { 'timeout' } elseif ($run.exit_code -ne 0) { 'failed' } elseif ($null -eq $profile) { 'profile_missing' } else { [string](P $profile 'status' 'failed') }
        $effective = [string](P $profile 'preconditioner_effective' '')
        $residual = D (P $profile 'true_residual_rel' (P $profile 'final_residual'))
        $audit = D (P $profile 'tree_vs_direct_matvec_rel' (P $profile 'matvec_audit_rel_rms'))
        $finite = (Test-FmmCliFinite $residual) -and (Test-FmmCliFinite $audit)
        $effectiveOk = $effective -eq $stage.expected
        $iterations = I (P $profile 'iterations')
        $failure = if ($status -eq 'timeout') { 'timeout' } elseif ($null -ne $profile -and $iterations -ge 720 -and $finite) { 'not_converged' } elseif ($status -ne 'ok') { 'runtime_failure' } elseif (-not $effectiveOk) { 'preconditioner_fallback' } elseif (-not $finite -or $audit -gt 1e-2) { 'accuracy_gap' } elseif ($residual -gt 1e-4) { 'not_converged' } else { 'pass' }
        $rows.Add([pscustomobject]@{
            row_id=$rowId; fixture=$fixtureName; protein_name=$fixture.protein_name; pdb_id=$fixture.pdb_id
            surface_protocol=$PrimarySurface; ses_grid_spacing_A=$(if($PrimarySurface -in @('grid-ses','analytic-ses')){$SesGridSpacing}else{'NA'})
            preconditioner_label=$stage.label; preconditioner_requested=$stage.precond
            preconditioner_effective=$effective; graph_cap=64; panels=I (P $profile 'panels')
            status=$status; failure_class=$failure; iterations=$iterations
            true_residual=$residual; final_residual=D (P $profile 'final_residual')
            audit_rms=$audit; wall_seconds=$run.wall_seconds; peak_rss_mb=$run.peak_working_set_mb
            setup_seconds=D (P $profile 'preconditioner_setup_seconds') 0.0
            apply_seconds=D (P $profile 'preconditioner_apply_seconds') 0.0
            matvec_seconds_avg=D (P $profile 'matvec_seconds_avg')
            nnz_blocks=I (P $profile 'preconditioner_sparse_nnz_blocks')
            avg_degree=D (P $profile 'preconditioner_sparse_avg_degree')
            max_degree=I (P $profile 'preconditioner_sparse_max_degree')
            valid_blocks=I (P $profile 'preconditioner_valid_blocks')
            valid_panels=I (P $profile 'preconditioner_valid_panels')
            quality_status=[string](P $profile 'preconditioner_quality_status')
            residual_history=[string](P $profile 'residual_history')
            protocol=$stageProtocol; solver_contract_sha256=$solverContractSha
            titan_sha256=$titanSha; pqr_sha256=$pqrSha; run_dir=$runDir
        })
        # Persist every completed case immediately. Long cohort runs may outlive
        # an interactive shell; the next invocation can then resume only rows
        # whose executable, solver contract, PQR, and protocol hashes match.
        @($rows.ToArray() | Sort-Object fixture,preconditioner_label) |
            Export-Csv -LiteralPath $matrixPath -NoTypeInformation -Encoding UTF8
    }
}
$matrix = @($rows.ToArray() | Sort-Object fixture,preconditioner_label)
foreach ($row in $matrix) {
    if ($row.PSObject.Properties.Name -notcontains 'surface_protocol') {
        $row | Add-Member -NotePropertyName surface_protocol -NotePropertyValue 'sphere-patch'
    }
    if ($row.PSObject.Properties.Name -notcontains 'ses_grid_spacing_A') {
        $row | Add-Member -NotePropertyName ses_grid_spacing_A -NotePropertyValue 'NA'
    }
    if ($row.PSObject.Properties.Name -notcontains 'solver_contract_sha256') {
        $value = if ([string](P $row 'titan_sha256') -eq $titanSha) {
            $solverContractSha
        } else { '' }
        $row | Add-Member -NotePropertyName solver_contract_sha256 -NotePropertyValue $value
    }
    $profilePath = Join-Path ([string]$row.run_dir) 'pbe_tabi_profile.csv'
    if (Test-Path -LiteralPath $profilePath) {
        $profile = @(Import-Csv -LiteralPath $profilePath)[-1]
        $row.setup_seconds = D (P $profile 'preconditioner_setup_seconds') 0.0
        $row.apply_seconds = D (P $profile 'preconditioner_apply_seconds') 0.0
        $row.nnz_blocks = I (P $profile 'preconditioner_sparse_nnz_blocks')
        $row.avg_degree = D (P $profile 'preconditioner_sparse_avg_degree')
        $row.max_degree = I (P $profile 'preconditioner_sparse_max_degree')
    }
    if ((I $row.iterations) -ge 720 -and
        (Test-FmmCliFinite (D $row.true_residual)) -and
        (Test-FmmCliFinite (D $row.audit_rms))) {
        $row.failure_class = 'not_converged'
    }
}
$matrix | Export-Csv -LiteralPath $matrixPath -NoTypeInformation -Encoding UTF8

$trace = [Collections.Generic.List[object]]::new()
foreach ($row in $matrix) {
    $pointsByIteration = @{'0' = 1.0}
    $progressPath = Join-Path ([string]$row.run_dir) 'pbe_tabi_longrun_progress.csv'
    if (Test-Path -LiteralPath $progressPath) {
        foreach ($event in @(Import-Csv -LiteralPath $progressPath | Where-Object {
                $_.phase -in @('gmres_restart','gmres_true_probe',
                    'gmres_early_true_confirmation','solve_done') })) {
            $iteration = I $event.gmres_iter -1
            $residual = D $event.residual
            if ($iteration -ge 0 -and (Test-FmmCliFinite $residual) -and
                $residual -gt 0) {
                $pointsByIteration[[string]$iteration] = $residual
            }
        }
    }
    if ($pointsByIteration.Count -eq 1) {
        $fallback = @(Convert-FmmCliResidualHistory ([string]$row.residual_history))
        $lastIndex = [Math]::Max(1,$fallback.Count-1)
        foreach ($point in $fallback) {
            $iteration = [int][Math]::Round((I $row.iterations) *
                ([double]$point.iteration / $lastIndex))
            $pointsByIteration[[string]$iteration] = [double]$point.residual
        }
    }
    foreach ($iteration in @($pointsByIteration.Keys | ForEach-Object {[int]$_} |
            Sort-Object)) {
        $residual = [double]$pointsByIteration[[string]$iteration]
        $trace.Add([pscustomobject]@{fixture=$row.fixture; protein_name=$row.protein_name
            surface_protocol=$row.surface_protocol; ses_grid_spacing_A=$row.ses_grid_spacing_A
            preconditioner_label=$row.preconditioner_label; iteration=$iteration
            residual=$residual; log10_residual=[Math]::Log10($residual)})
    }
}
$trace.ToArray() | Export-Csv -LiteralPath $tracePath -NoTypeInformation -Encoding UTF8

$colors = @{diagonal='#555555'; leaf_block_jacobi='#D55E00'; fmm_near_graph_ras='#0072B2';
    distance_knn_graph_ras='#56B4E9'; distance_knn_ilu0='#CC79A7'; fmm_near_graph_ilu0='#009E73'}
$dashes = @{diagonal=''; leaf_block_jacobi=''; fmm_near_graph_ras='';
    distance_knn_graph_ras='13 8'; distance_knn_ilu0='3 7'; fmm_near_graph_ilu0=''}
$labels = @{diagonal='Diagonal'; leaf_block_jacobi='Leaf block Jacobi'; fmm_near_graph_ras='FMM-near RAS';
    distance_knn_graph_ras='Distance-kNN RAS'; distance_knn_ilu0='Distance-kNN ILU(0)'; fmm_near_graph_ilu0='FMM-near ILU(0)'}
$plotStages = @($stages | Sort-Object @{Expression={
    if ($_.label -eq 'fmm_near_graph_ilu0') { 4 }
    elseif ($_.label -eq 'distance_knn_ilu0') { 5 }
    else { [array]::IndexOf(@($stages.label),$_.label) }
}})
$surfaceOrder = @(@('sphere-patch','grid-ses','analytic-ses') | Where-Object {
    $surface = $_
    @($matrix | Where-Object surface_protocol -eq $surface).Count -gt 0
})
$width=1800; $height=210+455*$surfaceOrder.Count; $left=145; $right=40; $top=140; $bottom=85
$columnGap=45; $rowGap=86
$panelWidth=($width-$left-$right-2*$columnGap)/3
$panelHeight=($height-$top-$bottom-([Math]::Max(0,$surfaceOrder.Count-1))*$rowGap)/[Math]::Max(1,$surfaceOrder.Count)
$svg=[Collections.Generic.List[string]]::new()
$physicalHeight=55*$surfaceOrder.Count
$svg.Add("<svg xmlns=`"http://www.w3.org/2000/svg`" width=`"180mm`" height=`"${physicalHeight}mm`" viewBox=`"0 0 $width $height`">")
$svg.Add('<rect width="100%" height="100%" fill="white"/>')
for($li=0;$li -lt $stages.Count;$li++){
    $stage=$stages[$li];$label=$stage.label;$name=Xml $labels[$label];$color=$colors[$label]
    $legendRow=[Math]::Floor($li/3);$legendColumn=$li%3
    $legendX=125+$legendColumn*555;$legendY=30+$legendRow*44
    $dash=if([string]::IsNullOrWhiteSpace($dashes[$label])){''}else{" stroke-dasharray=`"$($dashes[$label])`""}
    $svg.Add("<line x1=`"$legendX`" y1=`"$legendY`" x2=`"$($legendX+50)`" y2=`"$legendY`" stroke=`"$color`" stroke-width=`"4.5`"$dash/>")
    if($label -in @('distance_knn_graph_ras','distance_knn_ilu0')){
        $svg.Add("<line x1=`"$($legendX+25)`" y1=`"$($legendY-7)`" x2=`"$($legendX+25)`" y2=`"$($legendY+7)`" stroke=`"$color`" stroke-width=`"3`"/>")
        $svg.Add("<line x1=`"$($legendX+18)`" y1=`"$legendY`" x2=`"$($legendX+32)`" y2=`"$legendY`" stroke=`"$color`" stroke-width=`"3`"/>")
    } elseif($label -in @('fmm_near_graph_ras','fmm_near_graph_ilu0')){
        $svg.Add("<circle cx=`"$($legendX+25)`" cy=`"$legendY`" r=`"5.5`" fill=`"white`" stroke=`"$color`" stroke-width=`"3`"/>")
    }
    $svg.Add("<text x=`"$($legendX+64)`" y=`"$($legendY+8)`" font-family=`"Arial`" font-size=`"24`">$name</text>")
}
$panelIndex=0
for($si=0;$si -lt $surfaceOrder.Count;$si++){
    $surface=$surfaceOrder[$si]
    $rowY=$top+$si*($panelHeight+$rowGap)
    $surfaceLabel=switch($surface){
        'grid-ses' { "Voxel-derived watertight SES (h = $($SesGridSpacing.ToString('G4',$Inv)) A)" }
        'analytic-ses' { "Analytic patch-ray SES (h = $($SesGridSpacing.ToString('G4',$Inv)) A)" }
        default { 'Sphere-patch surface' }
    }
    $svg.Add("<text x=`"82`" y=`"$($rowY+$panelHeight/2)`" text-anchor=`"middle`" font-family=`"Arial`" font-size=`"24`" font-weight=`"bold`" transform=`"rotate(-90 82,$($rowY+$panelHeight/2))`">$(Xml $surfaceLabel)</text>")
    for($fi=0;$fi -lt $fixtureNames.Count;$fi++){
        $fixture=$fixtureNames[$fi]
        $points=@($trace|Where-Object {$_.fixture -eq $fixture -and $_.surface_protocol -eq $surface})
        if($points.Count -eq 0){continue}
        $x0=$left+$fi*($panelWidth+$columnGap);$x1=$x0+$panelWidth;$y0=$rowY;$y1=$rowY+$panelHeight
        $svg.Add("<line x1=`"$x0`" y1=`"$y1`" x2=`"$x1`" y2=`"$y1`" stroke=`"black`" stroke-width=`"1.3`"/>")
        $svg.Add("<line x1=`"$x0`" y1=`"$y0`" x2=`"$x0`" y2=`"$y1`" stroke=`"black`" stroke-width=`"1.3`"/>")
        $maxIter=1000.0;$minIter=1.0;$xLogDen=[Math]::Log10($maxIter/$minIter);$minLog=-6.0;$maxLog=0.0
        foreach($tick in @(0,-2,-4,-6)){
            $y=$y0+($y1-$y0)*(($maxLog-$tick)/($maxLog-$minLog))
            if($tick -eq -4){
                $svg.Add("<line x1=`"$x0`" y1=`"$y`" x2=`"$x1`" y2=`"$y`" stroke=`"#777777`" stroke-width=`"1.8`" stroke-dasharray=`"5 6`"/>")
                if($fi -eq $fixtureNames.Count-1){
                    $svg.Add("<text x=`"$($x1-8)`" y=`"$($y-10)`" text-anchor=`"end`" font-family=`"Arial`" font-size=`"18`" fill=`"#555555`">convergence threshold</text>")
                }
            }else{
                $svg.Add("<line x1=`"$x0`" y1=`"$y`" x2=`"$x1`" y2=`"$y`" stroke=`"#dddddd`"/>")
            }
            $svg.Add("<text x=`"$($x0-12)`" y=`"$($y+8)`" text-anchor=`"end`" font-family=`"Arial`" font-size=`"22`">$tick</text>")
        }
        foreach($tick in @(1,10,100,1000)){
            $x=$x0+($x1-$x0)*([Math]::Log10($tick/$minIter)/$xLogDen)
            $svg.Add("<line x1=`"$x`" y1=`"$y1`" x2=`"$x`" y2=`"$($y1+6)`" stroke=`"black`"/>")
            $svg.Add("<text x=`"$x`" y=`"$($y1+31)`" text-anchor=`"middle`" font-family=`"Arial`" font-size=`"22`">$tick</text>")
        }
        foreach($stage in $plotStages){
            $label=$stage.label
            $line=@($points|Where-Object preconditioner_label -eq $label|Sort-Object {[int]$_.iteration})
            if($line.Count -eq 0){continue}
            $coords=foreach($p in $line){
                $displayIteration=[Math]::Max($minIter,[double]$p.iteration)
                $x=$x0+($x1-$x0)*([Math]::Log10($displayIteration/$minIter)/$xLogDen)
                $v=[Math]::Max($minLog,[Math]::Min($maxLog,[double]$p.log10_residual))
                $y=$y0+($y1-$y0)*(($maxLog-$v)/($maxLog-$minLog))
                ('{0:F1},{1:F1}' -f $x,$y)
            }
            $dash=if([string]::IsNullOrWhiteSpace($dashes[$label])){''}else{" stroke-dasharray=`"$($dashes[$label])`""}
            $strokeWidth=if($label -in @('distance_knn_graph_ras','distance_knn_ilu0')){3.2}elseif($label -in @('fmm_near_graph_ras','fmm_near_graph_ilu0')){5.2}else{4.2}
            $svg.Add("<polyline fill=`"none`" stroke=`"$($colors[$label])`" stroke-width=`"$strokeWidth`" stroke-linejoin=`"round`" stroke-linecap=`"round`"$dash points=`"$($coords -join ' ')`"/>")
            foreach($marker in $line){
                $markerIteration=[Math]::Max($minIter,[double]$marker.iteration)
                $markerX=$x0+($x1-$x0)*([Math]::Log10($markerIteration/$minIter)/$xLogDen)
                $markerV=[Math]::Max($minLog,[Math]::Min($maxLog,[double]$marker.log10_residual))
                $markerY=$y0+($y1-$y0)*(($maxLog-$markerV)/($maxLog-$minLog))
                if($label -in @('distance_knn_graph_ras','distance_knn_ilu0')){
                    $svg.Add("<line x1=`"$markerX`" y1=`"$($markerY-6)`" x2=`"$markerX`" y2=`"$($markerY+6)`" stroke=`"$($colors[$label])`" stroke-width=`"2.8`"/>")
                    $svg.Add("<line x1=`"$($markerX-6)`" y1=`"$markerY`" x2=`"$($markerX+6)`" y2=`"$markerY`" stroke=`"$($colors[$label])`" stroke-width=`"2.8`"/>")
                } elseif($label -in @('fmm_near_graph_ras','fmm_near_graph_ilu0')){
                    $svg.Add("<circle cx=`"$markerX`" cy=`"$markerY`" r=`"5.1`" fill=`"white`" stroke=`"$($colors[$label])`" stroke-width=`"2.8`"/>")
                }
            }
        }
        $protein=Xml $(switch($fixture){'barnase'{'Barnase'}'1a63'{'1A63'}'1sh1'{'1SH1'}default{([string]$fixture).ToUpperInvariant()}})
        $letter=[char](65+$panelIndex);$panelIndex++
        $svg.Add("<text x=`"$x0`" y=`"$($y0-16)`" font-family=`"Arial`" font-size=`"26`" font-weight=`"bold`">$letter</text>")
        $svg.Add("<text x=`"$($x0+38)`" y=`"$($y0-16)`" font-family=`"Arial`" font-size=`"25`">$protein</text>")
    }
}
$svg.Add("<text x=`"27`" y=`"$($top+($height-$top-$bottom)/2)`" text-anchor=`"middle`" font-family=`"Arial`" font-size=`"24`" transform=`"rotate(-90 27,$($top+($height-$top-$bottom)/2))`">log10 relative residual</text>")
$svg.Add("<text x=`"$($left+($width-$left-$right)/2)`" y=`"$($height-22)`" text-anchor=`"middle`" font-family=`"Arial`" font-size=`"24`">GMRES iteration (log scale)</text>")
$svg.Add('</svg>')
Set-Content -LiteralPath $plotPath -Value $svg -Encoding ASCII

$requiredRows=@($matrix|Where-Object {$fixtureNames -contains [string]$_.fixture -and $surfaceOrder -contains [string]$_.surface_protocol})
$iluRows=@($requiredRows|Where-Object {$_.preconditioner_label -in @('distance_knn_ilu0','fmm_near_graph_ilu0')})
$requiredIluCount=2*$fixtureNames.Count*$surfaceOrder.Count
$requiredIluConverged=$iluRows.Count -eq $requiredIluCount -and @($iluRows|Where-Object { [string]$_.preconditioner_effective -notin @('distance_knn_ilu0','fmm_nearfield_ilu0') -or -not (Test-FmmCliFinite (D $_.true_residual)) -or (D $_.true_residual) -gt 1e-4 }).Count -eq 0
$requiredIluAuditPass=$requiredIluConverged -and @($iluRows|Where-Object { -not (Test-FmmCliFinite (D $_.audit_rms)) -or (D $_.audit_rms) -gt 1e-2 }).Count -eq 0
$graphSourcesIdentical=$true
foreach($surface in $surfaceOrder){
    foreach($fixture in $fixtureNames){
        $k=@($requiredRows|Where-Object {$_.surface_protocol -eq $surface -and $_.fixture -eq $fixture -and $_.preconditioner_label -eq 'distance_knn_ilu0'})
        $f=@($requiredRows|Where-Object {$_.surface_protocol -eq $surface -and $_.fixture -eq $fixture -and $_.preconditioner_label -eq 'fmm_near_graph_ilu0'})
        if($k.Count -ne 1 -or $f.Count -ne 1 -or [string]$k[0].residual_history -ne [string]$f[0].residual_history){$graphSourcesIdentical=$false}
    }
}
$timeouts=@($matrix|Where-Object failure_class -eq 'timeout').Count
$decision=if($timeouts -gt 0){'graph_ablation_runtime_limited'}elseif(-not $requiredIluConverged){'graph_ablation_ilu_convergence_gap'}elseif(-not $requiredIluAuditPass){'graph_ablation_operator_accuracy_gap'}elseif($graphSourcesIdentical){'graph_ablation_graph_sources_indistinguishable'}else{'graph_ablation_ready'}
@([pscustomobject]@{decision=$decision;fixtures=($fixtureNames -join ';');rows=$matrix.Count;ilu_rows=$iluRows.Count;ilu_converged=$requiredIluConverged;ilu_audit_pass=$requiredIluAuditPass;graph_source_histories_identical=$graphSourcesIdentical;protocol=$protocol})|Export-Csv -LiteralPath $decisionPath -NoTypeInformation -Encoding UTF8
@("decision $decision","fixtures $($fixtureNames -join ',')","rows $($matrix.Count)","matrix $matrixPath","trace $tracePath","figure $plotPath")|Set-Content -LiteralPath $summaryPath -Encoding ASCII
$surfaceLines=[Collections.Generic.List[string]]::new()
foreach($surface in $surfaceOrder){
    $surfaceName=switch($surface){
        'grid-ses' { "Grid-SES h=$($SesGridSpacing.ToString('G4',$Inv)) A" }
        'analytic-ses' { "Analytic patch-ray SES h=$($SesGridSpacing.ToString('G4',$Inv)) A" }
        default { 'Sphere-patch' }
    }
    $surfaceLines.Add("[$surfaceName]")
    foreach($fixture in $fixtureNames){
        $fmmIlu=@($requiredRows|Where-Object {$_.surface_protocol -eq $surface -and $_.fixture -eq $fixture -and $_.preconditioner_label -eq 'fmm_near_graph_ilu0'})
        $leaf=@($requiredRows|Where-Object {$_.surface_protocol -eq $surface -and $_.fixture -eq $fixture -and $_.preconditioner_label -eq 'leaf_block_jacobi'})
        $diag=@($requiredRows|Where-Object {$_.surface_protocol -eq $surface -and $_.fixture -eq $fixture -and $_.preconditioner_label -eq 'diagonal'})
        if($fmmIlu.Count -eq 1){
            $surfaceLines.Add(('{0}: panels={1}; diagonal_iter={2}; leaf_block_iter={3}; ilu_iter={4}; ilu_residual={5:G6}; ilu_audit={6:G6}' -f
                $fixture,(I $fmmIlu[0].panels),(I $(if($diag.Count -eq 1){$diag[0].iterations}else{0})),
                (I $(if($leaf.Count -eq 1){$leaf[0].iterations}else{0})),(I $fmmIlu[0].iterations),
                (D $fmmIlu[0].true_residual),(D $fmmIlu[0].audit_rms)))
        }
    }
    $surfaceLines.Add('')
}
$report=@('TITAN GRAPH PRECONDITIONER ABLATION','====================================','',
    "Decision: $decision","Fixtures: $($fixtureNames -join ', ')",
    "Protocols represented: $($surfaceOrder -join ', ')",
    "ILU convergence gate: $requiredIluConverged",
    "ILU operator-audit gate: $requiredIluAuditPass",
    "kNN/FMM-near residual histories identical: $graphSourcesIdentical",'',
    'Measured surface comparison:')
$report += $surfaceLines.ToArray()
$report += @('Interpretation: each surface generator changes the discretized operator and preconditioner quality; no surface is treated as merely a larger copy of another. ILU convergence alone is not an operator-accuracy claim, so sampled direct-operator audit remains mandatory. Identical kNN/FMM-near ILU histories mean this matrix does not establish superiority of either graph source.','',
    "Matrix: $matrixPath","Residual trace: $tracePath","Figure: $plotPath")
$report|Set-Content -LiteralPath $reportPath -Encoding ASCII
Write-Host "Graph ablation decision: $decision";Write-Host "Matrix: $matrixPath";Write-Host "Figure: $plotPath"
if(-not $AllowFailedDecision -and $decision -ne 'graph_ablation_ready'){throw "Graph ablation decision: $decision"}
