param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4j'),
    [int]$CaseTimeoutSeconds = 2400,
    [switch]$Quick,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
$matrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4j_matrix.csv'
$channelCsv = Join-Path $WorkRoot 'fagbi_scope_z4j_channel_summary.csv'
$levelCsv = Join-Path $WorkRoot 'fagbi_scope_z4j_level_summary.csv'
$decisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4j_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fagbi_scope_z4j_summary.dat'
$plotPath = Join-Path $WorkRoot 'fagbi_scope_z4j_comparison.svg'
$rootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4J_CARTESIAN_M2P_REPORT.txt'

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $manifestPath)) { throw "Missing public fixture manifest: $manifestPath" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
$manifest = @(Import-Csv -LiteralPath $manifestPath)
$auditTargets = if ($Quick) { 128 } else { 256 }

function Get-Z4jFixture([string]$Name) {
    $rows = @($manifest | Where-Object fixture -eq $Name)
    if ($rows.Count -ne 1) { throw "Expected exactly one fixture row for $Name" }
    return $rows[0]
}

function Get-Z4jDouble([object]$Value, [double]$Default = [double]::NaN) {
    $parsed = 0.0
    if ($null -eq $Value -or -not [double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
        return $Default
    }
    return $parsed
}

function Get-Z4jErrorSq([object[]]$Rows) {
    $sum = 0.0
    foreach ($row in $Rows) {
        $rms = Get-Z4jDouble $row.error_rms 0.0
        $sum += $rms * $rms * [Math]::Max(0, [int]$row.samples)
    }
    return $sum
}

function Invoke-Z4jCase([pscustomobject]$Fixture, [double]$Theta) {
    $thetaLabel = $Theta.ToString('0.00', [Globalization.CultureInfo]::InvariantCulture).Replace('.','p')
    $runDir = Join-Path $WorkRoot "$($Fixture.fixture)_theta$thetaLabel"
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $arg = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot `
        -PqrPath (Join-Path $FixtureRoot $Fixture.pqr) -Backend tabi `
        -Mode treecode -FmmMode single-matvec-analytic -SurfacePoints 8 `
        -MaxIter 720 -Restart 128 -Tolerance 1e-6 -Threads 4 -SaltM 0.15 `
        -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
        -OpeningTheta $Theta -AuditSampleCount $auditTargets `
        -AuditMode sampled-detail -Preconditioner fmm-nearfield-ilu0 `
        -PreconditionerIluNeighborCap 64 -PreconditionerIluDropTol 0.0 `
        -NearIntegration analytic_singular -TaylorOrder 4 -FmmCartesianAudit on
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.timed_out) {
        return @([pscustomobject]@{fixture=$Fixture.fixture; theta=$Theta
            approximation_mode=''; status='timeout'; audit_rms=[double]::NaN
            k2_error_sq=[double]::NaN; total_error_sq=[double]::NaN
            dominant_channel=''; implemented_degree=0; configured_order=4
            recomposition_error=[double]::NaN; reference_gap=[double]::NaN
            moment_build_seconds=[double]::NaN; far_evaluation_seconds=[double]::NaN
            wall_seconds=$run.wall_seconds; peak_rss_mb=$run.peak_working_set_mb
            run_dir=$runDir})
    }
    if ($run.exit_code -ne 0) {
        throw "Z4J case failed: fixture=$($Fixture.fixture) theta=$Theta exit=$($run.exit_code)"
    }
    $case = Read-FmmCliCaseResult -RunDir $runDir `
        -ExpectedFmmMode single-matvec-analytic -RequireAudit $true
    $summaryFile = Join-Path $runDir 'pbe_tabi_operator_attribution_summary.csv'
    $levelFile = Join-Path $runDir 'pbe_tabi_operator_attribution_levels.csv'
    if (-not (Test-Path $summaryFile) -or -not (Test-Path $levelFile)) {
        throw "Z4J case missing attribution output: $runDir"
    }
    $summaryRows = @(Import-Csv $summaryFile)
    $modes = @($summaryRows.approximation_mode | Sort-Object -Unique)
    foreach ($required in @('legacy_degree1','cartesian_p4')) {
        if ($required -notin $modes) { throw "Z4J case missing mode $required" }
    }
    $output = New-Object System.Collections.Generic.List[object]
    foreach ($mode in @('legacy_degree1','cartesian_p4')) {
        $modeRows = @($summaryRows | Where-Object approximation_mode -eq $mode)
        $far = @($modeRows | Where-Object {
                $_.probe -eq 'solution' -and $_.region -eq 'far' -and
                $_.channel -in @('K1','K2','K3','K4') })
        $channelErrors = @{}
        foreach ($channel in @('K1','K2','K3','K4')) {
            $channelErrors[$channel] = Get-Z4jErrorSq @($far | Where-Object channel -eq $channel)
        }
        $dominant = @($channelErrors.GetEnumerator() | Sort-Object Value -Descending | Select-Object -First 1)[0]
        $aggregate = @($modeRows | Where-Object {
                $_.probe -eq 'solution' -and $_.row_component -eq 'both' -and
                $_.region -eq 'total' -and $_.channel -eq 'aggregate' }) |
            Select-Object -First 1
        $output.Add([pscustomobject]@{
            fixture=$Fixture.fixture; protein_name=$Fixture.protein_name
            pdb_id=$Fixture.pdb_id; theta=$Theta; approximation_mode=$mode
            status=$case.status; panels=$case.panels; true_residual=$case.true_residual
            audit_rms=Get-Z4jDouble $aggregate.aggregate_audit_rel
            k2_error_sq=[double]$channelErrors.K2; total_error_sq=Get-Z4jErrorSq $far
            dominant_channel=[string]$dominant.Key
            implemented_degree=[int]$aggregate.implemented_far_moment_degree
            configured_order=[int]$aggregate.configured_expansion_order
            recomposition_error=Get-Z4jDouble $aggregate.recomposition_rel_error
            reference_gap=Get-Z4jDouble $aggregate.reference_gap
            moment_build_seconds=Get-Z4jDouble $aggregate.moment_build_seconds 0.0
            far_evaluation_seconds=Get-Z4jDouble $aggregate.far_evaluation_seconds 0.0
            wall_seconds=$run.wall_seconds; peak_rss_mb=$run.peak_working_set_mb
            run_dir=$runDir
        }) | Out-Null
    }
    return $output.ToArray()
}

$fixtures = if ($Quick) { @('1sh1','1a63') } else { @('barnase','1a63','1sh1') }
$thetas = if ($Quick) { @(0.45) } else { @(0.6,0.45) }
$matrix = New-Object System.Collections.Generic.List[object]
$channelRows = New-Object System.Collections.Generic.List[object]
$levelRows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $fixtures) {
    $fixture = Get-Z4jFixture $fixtureName
    foreach ($theta in $thetas) {
        foreach ($row in @(Invoke-Z4jCase -Fixture $fixture -Theta $theta)) {
            $matrix.Add($row) | Out-Null
        }
        $runDir = $matrix[$matrix.Count - 1].run_dir
        $caseSummary = Join-Path $runDir 'pbe_tabi_operator_attribution_summary.csv'
        $caseLevels = Join-Path $runDir 'pbe_tabi_operator_attribution_levels.csv'
        foreach ($row in @(if (Test-Path $caseSummary) { Import-Csv $caseSummary })) {
            $row | Add-Member -NotePropertyName fixture -NotePropertyValue $fixtureName
            $row | Add-Member -NotePropertyName theta -NotePropertyValue $theta
            $channelRows.Add($row) | Out-Null
        }
        foreach ($row in @(if (Test-Path $caseLevels) { Import-Csv $caseLevels })) {
            $row | Add-Member -NotePropertyName fixture -NotePropertyValue $fixtureName
            $row | Add-Member -NotePropertyName theta -NotePropertyValue $theta
            $levelRows.Add($row) | Out-Null
        }
    }
}
$rows = $matrix.ToArray()
$rows | Export-Csv $matrixCsv -NoTypeInformation
$channelRows | Export-Csv $channelCsv -NoTypeInformation
$levelRows | Export-Csv $levelCsv -NoTypeInformation

$decision = 'z4j_cartesian_formula_gap'
if (@($rows | Where-Object status -eq 'timeout').Count) {
    $decision = 'z4j_runtime_limited'
} elseif (@($rows | Where-Object {
            $_.approximation_mode -eq 'legacy_degree1' -and
            [double]$_.reference_gap -gt 1e-10 }).Count) {
    $decision = 'z4j_reference_gap'
} else {
    $requiredTheta06 = @($rows | Where-Object {
            $_.approximation_mode -eq 'cartesian_p4' -and [double]$_.theta -eq 0.6 })
    $requiredTheta045 = @($rows | Where-Object {
            $_.approximation_mode -eq 'cartesian_p4' -and [double]$_.theta -eq 0.45 })
    if (-not $Quick -and $requiredTheta06.Count -eq 3 -and
        @($requiredTheta06 | Where-Object {[double]$_.audit_rms -gt 8e-3}).Count -eq 0) {
        $decision = 'z4j_p4_accuracy_pass'
    } elseif ($requiredTheta045.Count -eq $fixtures.Count -and
        @($requiredTheta045 | Where-Object {[double]$_.audit_rms -gt 8e-3}).Count -eq 0) {
        $decision = 'z4j_p4_tight_theta_pass'
    } else {
        $ratios = @()
        $k2Ratios = @()
        foreach ($fixtureName in $fixtures) {
            foreach ($theta in $thetas) {
                $legacy = @($rows | Where-Object { $_.fixture -eq $fixtureName -and
                    [double]$_.theta -eq $theta -and $_.approximation_mode -eq 'legacy_degree1' }) | Select-Object -First 1
                $p4 = @($rows | Where-Object { $_.fixture -eq $fixtureName -and
                    [double]$_.theta -eq $theta -and $_.approximation_mode -eq 'cartesian_p4' }) | Select-Object -First 1
                if ($null -ne $legacy -and $null -ne $p4) {
                    $ratios += [double]$legacy.total_error_sq / [Math]::Max(1e-300,[double]$p4.total_error_sq)
                    $k2Ratios += [double]$legacy.k2_error_sq / [Math]::Max(1e-300,[double]$p4.k2_error_sq)
                }
            }
        }
        if ($k2Ratios.Count -and ($k2Ratios | Measure-Object -Minimum).Minimum -ge 10) {
            if (@($rows | Where-Object { $_.approximation_mode -eq 'cartesian_p4' -and $_.dominant_channel -ne 'K2' }).Count) {
                $decision = 'z4j_channel_migration'
            } else { $decision = 'z4j_k2_improved_not_strict' }
        } elseif ($ratios.Count -and ($ratios | Measure-Object -Minimum).Minimum -ge 4) {
            $decision = 'z4j_k2_improved_not_strict'
        }
    }
}

[pscustomobject]@{decision=$decision; quick=[bool]$Quick; rows=$rows.Count; audit_targets=$auditTargets} |
    Export-Csv $decisionCsv -NoTypeInformation
@('scope: Z4J Cartesian p4 M2P accuracy gate',"decision: $decision",
  "quick: $([bool]$Quick)","rows: $($rows.Count)","audit_targets: $auditTargets",
  'production_fmm_operator: legacy_degree1 unchanged','production_default: grid remains default') |
    Set-Content $summaryPath -Encoding ASCII

$svg = @('<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="620" viewBox="0 0 1000 620">',
    '<rect width="1000" height="620" fill="white"/>',
    '<text x="35" y="38" font-family="Arial" font-size="20">Cartesian p4 M2P versus legacy degree-1 attribution</text>')
$y = 82
foreach ($fixtureName in $fixtures) {
    foreach ($theta in $thetas) {
        $legacy = @($rows | Where-Object {$_.fixture -eq $fixtureName -and [double]$_.theta -eq $theta -and $_.approximation_mode -eq 'legacy_degree1'}) | Select-Object -First 1
        $p4 = @($rows | Where-Object {$_.fixture -eq $fixtureName -and [double]$_.theta -eq $theta -and $_.approximation_mode -eq 'cartesian_p4'}) | Select-Object -First 1
        if ($null -eq $legacy -or $null -eq $p4) { continue }
        $maxAudit = [Math]::Max(1e-12,[Math]::Max([double]$legacy.audit_rms,[double]$p4.audit_rms))
        $svg += ('<text x="35" y="{0}" font-family="Arial" font-size="13">{1} theta={2}</text>' -f $y,$fixtureName,$theta)
        $svg += ('<rect x="260" y="{0}" width="{1}" height="14" fill="#d55e00"/>' -f ($y-13),(600*[double]$legacy.audit_rms/$maxAudit))
        $svg += ('<rect x="260" y="{0}" width="{1}" height="14" fill="#0072b2"/>' -f ($y+4),(600*[double]$p4.audit_rms/$maxAudit))
        $y += 55
    }
}
$svg += '<text x="260" y="595" font-family="Arial" font-size="12">orange: legacy degree-1; blue: Cartesian p4; each pair normalized locally</text>'
$svg += '</svg>'
$svg | Set-Content $plotPath -Encoding ASCII

$report = @('# TITAN FMM/M2L Scope Z4J Cartesian M2P Accuracy Report','',
    "Decision: $decision",'',
    'Z4J evaluates full degree-4 Cartesian source moments only inside sampled-detail attribution. Production GMRES and the production FMM operator remain degree-1.','')
foreach ($row in $rows) {
    $report += ('- {0} theta={1} {2}: audit={3}, K2_sq={4}, total_far_sq={5}, dominant={6}, degree={7}, moment_s={8}, far_eval_s={9}' -f
        $row.fixture,$row.theta,$row.approximation_mode,$row.audit_rms,$row.k2_error_sq,
        $row.total_error_sq,$row.dominant_channel,$row.implemented_degree,
        $row.moment_build_seconds,$row.far_evaluation_seconds)
}
$report += @('','Grid remains the production Poisson-Boltzmann default. A true local M2L integration is reserved for Z4K if this p4 gate passes scientific accuracy.')
$report | Set-Content $rootReport -Encoding ASCII

$badArgs = @(Get-ChildItem $WorkRoot -Recurse -Filter arg.in | ForEach-Object {
    $path=$_.FullName; Get-Content $path | Where-Object {$_ -match '^\s*-' -and $_ -notmatch '^\s*-[^=\s]+='} | ForEach-Object {"$path :: $_"}
})
if ($badArgs.Count) { throw "Generated arg.in contains bare options: $($badArgs -join '; ')" }
if (-not $AllowFailedDecision -and $decision -in @('z4j_reference_gap','z4j_runtime_limited','z4j_cartesian_formula_gap')) {
    throw "Scope Z4J decision: $decision"
}
Write-Host "PASS Scope Z4J decision=$decision rows=$($rows.Count)"
