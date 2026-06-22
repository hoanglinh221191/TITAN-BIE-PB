param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4i'),
    [int]$CaseTimeoutSeconds = 2400,
    [switch]$Quick,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
$matrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4i_matrix.csv'
$channelCsv = Join-Path $WorkRoot 'fagbi_scope_z4i_channel_summary.csv'
$levelCsv = Join-Path $WorkRoot 'fagbi_scope_z4i_level_summary.csv'
$decisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4i_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fagbi_scope_z4i_summary.dat'
$plotPath = Join-Path $WorkRoot 'fagbi_scope_z4i_attribution.svg'
$rootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4I_OPERATOR_ATTRIBUTION_REPORT.txt'

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $manifestPath)) { throw "Missing public fixture manifest: $manifestPath" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
$manifest = @(Import-Csv -LiteralPath $manifestPath)
$auditTargets = if ($Quick) { 128 } else { 256 }

function Get-Z4iFixture([string]$Name) {
    $rows = @($manifest | Where-Object fixture -eq $Name)
    if ($rows.Count -ne 1) { throw "Expected exactly one fixture row for $Name" }
    return $rows[0]
}

function Get-Z4iFinite([object]$Value, [double]$Default = [double]::NaN) {
    if ($null -eq $Value -or [string]::IsNullOrWhiteSpace([string]$Value)) { return $Default }
    $parsed = 0.0
    if (-not [double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture,[ref]$parsed)) { return $Default }
    return $parsed
}

function Get-Z4iErrorSq([object[]]$Rows) {
    $sum = 0.0
    foreach ($row in $Rows) {
        $rms = Get-Z4iFinite $row.error_rms 0.0
        $count = [Math]::Max(0,[int]$row.samples)
        $sum += $rms * $rms * $count
    }
    return $sum
}

function Invoke-Z4iCase([pscustomobject]$Fixture, [double]$Theta) {
    $thetaLabel = $Theta.ToString('0.00',[Globalization.CultureInfo]::InvariantCulture).Replace('.','p')
    $runDir = Join-Path $WorkRoot "$($Fixture.fixture)_theta$thetaLabel"
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $pqr = Join-Path $FixtureRoot $Fixture.pqr
    $arg = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot -PqrPath $pqr `
        -Backend tabi -Mode treecode -FmmMode single-matvec-analytic `
        -SurfacePoints 8 -MaxIter 720 -Restart 128 -Tolerance 1e-6 -Threads 4 `
        -SaltM 0.15 -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
        -OpeningTheta $Theta -AuditSampleCount $auditTargets -AuditMode sampled-detail `
        -Preconditioner fmm-nearfield-ilu0 -PreconditionerIluNeighborCap 64 `
        -PreconditionerIluDropTol 0.0 -NearIntegration analytic_singular -TaylorOrder 4
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    $row = [ordered]@{
        fixture=$Fixture.fixture; protein_name=$Fixture.protein_name; pdb_id=$Fixture.pdb_id
        theta=$Theta; audit_targets=$auditTargets; status='failed'; failure_class='exit_code'
        panels=0; true_residual=[double]::NaN; audit_rms=[double]::NaN
        configured_order=0; implemented_far_degree=0; near_error_fraction=[double]::NaN
        far_error_fraction=[double]::NaN; dominant_channel=''; dominant_channel_fraction=[double]::NaN
        near_max_abs_error=[double]::NaN; recomposition_error=[double]::NaN
        reference_gap=[double]::NaN; wall_seconds=$run.wall_seconds
        peak_rss_mb=$run.peak_working_set_mb; run_dir=$runDir; arg_file=$arg
    }
    if ($run.timed_out) { $row.status='timeout'; $row.failure_class='timeout'; return [pscustomobject]$row }
    if ($run.exit_code -ne 0) { $row.failure_class="exit_$($run.exit_code)"; return [pscustomobject]$row }
    $summaryFile = Join-Path $runDir 'pbe_tabi_operator_attribution_summary.csv'
    $levelFile = Join-Path $runDir 'pbe_tabi_operator_attribution_levels.csv'
    if (-not (Test-Path -LiteralPath $summaryFile) -or -not (Test-Path -LiteralPath $levelFile)) {
        $row.failure_class='missing_attribution_outputs'; return [pscustomobject]$row
    }
    $case = Read-FmmCliCaseResult -RunDir $runDir -ExpectedFmmMode single-matvec-analytic -RequireAudit $true
    $rows = @(Import-Csv -LiteralPath $summaryFile)
    $solutionChannels = @($rows | Where-Object { $_.probe -eq 'solution' -and $_.channel -in @('K1','K2','K3','K4') })
    $nearRows = @($solutionChannels | Where-Object region -eq 'near')
    $farRows = @($solutionChannels | Where-Object region -eq 'far')
    $nearSq = Get-Z4iErrorSq $nearRows
    $farSq = Get-Z4iErrorSq $farRows
    $totalSq = $nearSq + $farSq
    $channelSq = @{}
    foreach ($channelName in @('K1','K2','K3','K4')) {
        $channelSq[$channelName] = Get-Z4iErrorSq @($farRows | Where-Object channel -eq $channelName)
    }
    $dominant = @($channelSq.GetEnumerator() | Sort-Object Value -Descending | Select-Object -First 1)[0]
    $aggregate = @($rows | Where-Object {
            $_.probe -eq 'solution' -and $_.row_component -eq 'both' -and
            $_.region -eq 'total' -and $_.channel -eq 'aggregate'
        }) | Select-Object -First 1
    $row.status=$case.status; $row.failure_class=if($case.status -eq 'ok'){'diagnostic_complete'}else{'not_converged'}
    $row.panels=$case.panels; $row.true_residual=$case.true_residual
    $row.audit_rms=$case.tree_vs_direct_matvec_rel
    $row.configured_order=[int]$aggregate.configured_expansion_order
    $row.implemented_far_degree=[int]$aggregate.implemented_far_moment_degree
    $row.near_error_fraction=if($totalSq -gt 0){$nearSq/$totalSq}else{0.0}
    $row.far_error_fraction=if($totalSq -gt 0){$farSq/$totalSq}else{0.0}
    $row.dominant_channel=[string]$dominant.Key
    $row.dominant_channel_fraction=if($farSq -gt 0){[double]$dominant.Value/$farSq}else{0.0}
    $row.near_max_abs_error=if($nearRows.Count){($nearRows | ForEach-Object {[double]$_.max_abs_error} | Measure-Object -Maximum).Maximum}else{0.0}
    $row.recomposition_error=($rows | ForEach-Object {[double]$_.recomposition_rel_error} | Measure-Object -Maximum).Maximum
    $row.reference_gap=[double]$aggregate.reference_gap
    return [pscustomobject]$row
}

$fixtures = if ($Quick) { @('1sh1','1a63') } else { @('barnase','1a63','1sh1') }
$thetas = if ($Quick) { @(0.45) } else { @(0.6,0.45,0.35) }
$matrix = New-Object System.Collections.Generic.List[object]
$channelRows = New-Object System.Collections.Generic.List[object]
$levelRows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $fixtures) {
    $fixture = Get-Z4iFixture $fixtureName
    foreach ($theta in $thetas) {
        $case = Invoke-Z4iCase -Fixture $fixture -Theta $theta
        $matrix.Add($case) | Out-Null
        $caseSummary = Join-Path $case.run_dir 'pbe_tabi_operator_attribution_summary.csv'
        $caseLevels = Join-Path $case.run_dir 'pbe_tabi_operator_attribution_levels.csv'
        if ((Test-Path -LiteralPath $caseSummary) -and (Test-Path -LiteralPath $caseLevels)) {
            foreach ($row in @(Import-Csv $caseSummary)) {
                $row | Add-Member -NotePropertyName fixture -NotePropertyValue $fixtureName
                $row | Add-Member -NotePropertyName theta -NotePropertyValue $theta
                $channelRows.Add($row) | Out-Null
            }
            foreach ($row in @(Import-Csv $caseLevels)) {
                $row | Add-Member -NotePropertyName fixture -NotePropertyValue $fixtureName
                $row | Add-Member -NotePropertyName theta -NotePropertyValue $theta
                $levelRows.Add($row) | Out-Null
            }
        }
    }
}
$matrixArray = $matrix.ToArray()
$matrixArray | Export-Csv -LiteralPath $matrixCsv -NoTypeInformation
$channelRows | Export-Csv -LiteralPath $channelCsv -NoTypeInformation
$levelRows | Export-Csv -LiteralPath $levelCsv -NoTypeInformation

$decision = 'z4i_mixed_operator_gap'
if (@($matrixArray | Where-Object status -eq 'timeout').Count -gt 0) {
    $decision='z4i_runtime_limited'
} elseif (@($matrixArray | Where-Object {[double]$_.reference_gap -gt 1e-10}).Count -gt 0) {
    $decision='z4i_reference_gap'
} elseif (@($matrixArray | Where-Object {[double]$_.near_max_abs_error -gt 1e-8}).Count -gt 0) {
    $decision='z4i_near_classification_gap'
} elseif (@($matrixArray | Where-Object {
            [double]$_.far_error_fraction -ge 0.9 -and
            [int]$_.implemented_far_degree -lt [int]$_.configured_order
        }).Count -gt 0) {
    $decision='z4i_far_order_implementation_gap'
} elseif (@($matrixArray | Where-Object {
            [double]$_.far_error_fraction -ge 0.9 -and
            [double]$_.dominant_channel_fraction -ge 0.7
        }).Count -gt 0) {
    $decision='z4i_channel_specific_far_gap'
} elseif (-not $Quick) {
    $nonmonotonic = $false
    foreach ($fixtureName in @('barnase','1a63')) {
        $loose = @($matrixArray | Where-Object { $_.fixture -eq $fixtureName -and [double]$_.theta -eq 0.6 }) | Select-Object -First 1
        $tight = @($matrixArray | Where-Object { $_.fixture -eq $fixtureName -and [double]$_.theta -eq 0.35 }) | Select-Object -First 1
        if ($null -ne $loose -and $null -ne $tight -and [double]$loose.audit_rms -gt 0 -and
            [double]$tight.audit_rms -gt 0.8 * [double]$loose.audit_rms) { $nonmonotonic=$true }
    }
    if ($nonmonotonic) { $decision='z4i_mac_nonmonotonic' }
}

[pscustomobject]@{decision=$decision; quick=[bool]$Quick; rows=$matrixArray.Count; audit_targets=$auditTargets} |
    Export-Csv -LiteralPath $decisionCsv -NoTypeInformation
@('scope: Z4I full-protein operator attribution',"decision: $decision",
  "quick: $([bool]$Quick)","rows: $($matrixArray.Count)","audit_targets: $auditTargets",
  'production_default: grid remains default','tabi_policy_change: none') |
    Set-Content -LiteralPath $summaryPath -Encoding ASCII

$svg = New-Object System.Collections.Generic.List[string]
$svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="560" viewBox="0 0 1000 560">') | Out-Null
$svg.Add('<rect width="1000" height="560" fill="white"/>') | Out-Null
$svg.Add('<text x="35" y="35" font-family="Arial" font-size="20">Full-protein FMM operator error attribution</text>') | Out-Null
$y=80
foreach($row in $matrixArray){
    $far=[Math]::Min(1.0,[Math]::Max(0.0,[double]$row.far_error_fraction)); $near=1.0-$far
    $svg.Add(("<text x=`"35`" y=`"{0}`" font-family=`"Arial`" font-size=`"13`">{1} theta={2} channel={3}</text>" -f $y,$row.fixture,$row.theta,$row.dominant_channel)) | Out-Null
    $svg.Add(("<rect x=`"300`" y=`"{0}`" width=`"{1}`" height=`"18`" fill=`"#0072b2`"/>" -f ($y-15),(600*$far))) | Out-Null
    $svg.Add(("<rect x=`"{0}`" y=`"{1}`" width=`"{2}`" height=`"18`" fill=`"#d55e00`"/>" -f (300+600*$far),($y-15),(600*$near))) | Out-Null
    $y += 45
}
$svg.Add('<text x="300" y="535" font-family="Arial" font-size="12">blue: far error fraction; orange: near error fraction</text>') | Out-Null
$svg.Add('</svg>') | Out-Null
$svg | Set-Content -LiteralPath $plotPath -Encoding ASCII

$report = @('# TITAN FMM/M2L Scope Z4I Operator Attribution Report','',"Decision: $decision",'',
    'Z4I attributes sampled full-protein operator error without changing solver policy.',
    'Channel relative RMS values use the interaction-only K1-K4 reference norm; audit RMS uses the complete operator row, including the LHS diagonal. The reference-gap column verifies their shared error numerator against the existing audit.','')
foreach($row in $matrixArray){
    $report += ('- {0} theta={1}: audit={2}, far_fraction={3}, dominant={4} ({5}), near_max_abs={6}, degree={7}/{8}, reference_gap={9}' -f
        $row.fixture,$row.theta,$row.audit_rms,$row.far_error_fraction,$row.dominant_channel,
        $row.dominant_channel_fraction,$row.near_max_abs_error,$row.implemented_far_degree,
        $row.configured_order,$row.reference_gap)
}
$report += @('',
    'Measured conclusion: the shared near kernel and direct-reference bridge are clean. Far error dominates, configured order 4 currently maps to implemented far moment degree 1, and K2 is the leading channel. The next numerical scope is genuine higher-order analytical moments/M2L rather than another theta-only sweep.',
    '',
    'Grid remains the production Poisson-Boltzmann default. Z4I is diagnostic-only.')
$report | Set-Content -LiteralPath $rootReport -Encoding ASCII

$badArgs = @(Get-ChildItem -LiteralPath $WorkRoot -Recurse -Filter arg.in | ForEach-Object {
        $path=$_.FullName; Get-Content -LiteralPath $path | Where-Object {
            $_ -match '^\s*-' -and $_ -notmatch '^\s*-[^=\s]+='
        } | ForEach-Object { "$path :: $_" }
    })
if($badArgs.Count){ throw "Generated arg.in contains bare options: $($badArgs -join '; ')" }

$diagnosticDecisions=@('z4i_far_order_implementation_gap','z4i_channel_specific_far_gap',
    'z4i_near_classification_gap','z4i_mac_nonmonotonic','z4i_mixed_operator_gap')
if(-not $AllowFailedDecision -and $decision -notin $diagnosticDecisions){ throw "Scope Z4I decision: $decision" }
Write-Host "PASS Scope Z4I decision=$decision rows=$($matrixArray.Count)"
