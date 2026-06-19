param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4k'),
    [int]$CaseTimeoutSeconds = 2400,
    [switch]$Quick,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
$matrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4k_matrix.csv'
$channelCsv = Join-Path $WorkRoot 'fagbi_scope_z4k_channel_summary.csv'
$timingCsv = Join-Path $WorkRoot 'fagbi_scope_z4k_timing.csv'
$decisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4k_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fagbi_scope_z4k_summary.dat'
$plotPath = Join-Path $WorkRoot 'fagbi_scope_z4k_m2l_plot.svg'
$rootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4K_LOCAL_M2L_REPORT.txt'

if (-not (Test-Path $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path $manifestPath)) { throw "Missing public fixture manifest: $manifestPath" }
$TitanExe = (Resolve-Path $TitanExe).Path
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
$manifest = @(Import-Csv $manifestPath)
$auditTargets = if ($Quick) { 128 } else { 256 }

function Get-Z4kFixture([string]$Name) {
    $rows = @($manifest | Where-Object fixture -eq $Name)
    if ($rows.Count -ne 1) { throw "Expected exactly one fixture row for $Name" }
    return $rows[0]
}

function Get-Z4kDouble([object]$Value, [double]$Default = [double]::NaN) {
    $parsed = 0.0
    if ($null -eq $Value -or -not [double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture,[ref]$parsed)) {
        return $Default
    }
    return $parsed
}

function Get-Z4kModeParity([object[]]$Samples) {
    $totals = @{}
    foreach ($row in $Samples) {
        if ($row.approximation_mode -notin @('cartesian_p4','local_m2l_p4')) { continue }
        $key = "$($row.approximation_mode)|$($row.probe)|$($row.target_panel)|$($row.channel)"
        if (-not $totals.ContainsKey($key)) { $totals[$key] = 0.0 }
        $totals[$key] += Get-Z4kDouble $row.approximate 0.0
    }
    $diffSq=0.0; $refSq=0.0; $count=0
    foreach ($key in @($totals.Keys | Where-Object {$_ -like 'cartesian_p4|solution|*'})) {
        $localKey = $key -replace '^cartesian_p4\|','local_m2l_p4|'
        if (-not $totals.ContainsKey($localKey)) { continue }
        $delta=[double]$totals[$localKey]-[double]$totals[$key]
        $diffSq += $delta*$delta; $refSq += [double]$totals[$key]*[double]$totals[$key]
        ++$count
    }
    return [pscustomobject]@{relative_rms=[Math]::Sqrt($diffSq/[Math]::Max(1e-300,$refSq)); samples=$count}
}

function Invoke-Z4kCase([pscustomobject]$Fixture,[double]$Theta) {
    $thetaLabel=$Theta.ToString('0.00',[Globalization.CultureInfo]::InvariantCulture).Replace('.','p')
    $runDir=Join-Path $WorkRoot "$($Fixture.fixture)_theta$thetaLabel"
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $arg=Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot `
        -PqrPath (Join-Path $FixtureRoot $Fixture.pqr) -Backend tabi `
        -Mode treecode -FmmMode single-matvec-analytic -SurfacePoints 8 `
        -MaxIter 720 -Restart 128 -Tolerance 1e-6 -Threads 4 -SaltM 0.15 `
        -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
        -OpeningTheta $Theta -AuditSampleCount $auditTargets `
        -AuditMode sampled-detail -Preconditioner fmm-nearfield-ilu0 `
        -PreconditionerIluNeighborCap 64 -PreconditionerIluDropTol 0.0 `
        -NearIntegration analytic_singular -TaylorOrder 4 `
        -FmmCartesianAudit on -FmmLocal audit-p4
    $run=Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    if($run.timed_out){return [pscustomobject]@{fixture=$Fixture.fixture;theta=$Theta;status='timeout';run_dir=$runDir;wall_seconds=$run.wall_seconds}}
    if($run.exit_code -ne 0){throw "Z4K case failed fixture=$($Fixture.fixture) theta=$Theta exit=$($run.exit_code)"}
    $case=Read-FmmCliCaseResult -RunDir $runDir -ExpectedFmmMode single-matvec-analytic -RequireAudit $true
    $summary=@(Import-Csv (Join-Path $runDir 'pbe_tabi_operator_attribution_summary.csv'))
    $samples=@(Import-Csv (Join-Path $runDir 'pbe_tabi_operator_attribution_samples.csv'))
    $parity=Get-Z4kModeParity $samples
    $legacyAggregate=@($summary | Where-Object {$_.approximation_mode -eq 'legacy_degree1' -and $_.probe -eq 'solution' -and $_.row_component -eq 'both' -and $_.region -eq 'total' -and $_.channel -eq 'aggregate'}) | Select-Object -First 1
    if($null -eq $legacyAggregate){throw 'Missing legacy direct-reference bridge row'}
    $referenceBridgeGap=Get-Z4kDouble $legacyAggregate.reference_gap
    $rows=New-Object System.Collections.Generic.List[object]
    foreach($mode in @('cartesian_p4','local_m2l_p4')){
        $aggregate=@($summary | Where-Object {$_.approximation_mode -eq $mode -and $_.probe -eq 'solution' -and $_.row_component -eq 'both' -and $_.region -eq 'total' -and $_.channel -eq 'aggregate'}) | Select-Object -First 1
        if($null -eq $aggregate){throw "Missing Z4K aggregate for $mode"}
        $farTime=Get-Z4kDouble $aggregate.far_evaluation_seconds 0.0
        $m2p=@($summary | Where-Object {$_.approximation_mode -eq 'cartesian_p4' -and $_.probe -eq 'solution' -and $_.row_component -eq 'both' -and $_.region -eq 'total' -and $_.channel -eq 'aggregate'}) | Select-Object -First 1
        $speedup=(Get-Z4kDouble $m2p.far_evaluation_seconds 0.0)/[Math]::Max(1e-300,$farTime)
        $rows.Add([pscustomobject]@{
            fixture=$Fixture.fixture;protein_name=$Fixture.protein_name;pdb_id=$Fixture.pdb_id
            theta=$Theta;approximation_mode=$mode;status=$case.status;panels=$case.panels
            audit_rms=Get-Z4kDouble $aggregate.aggregate_audit_rel
            local_vs_m2p_rel=$parity.relative_rms;local_vs_m2p_samples=$parity.samples
            recomposition_error=Get-Z4kDouble $aggregate.recomposition_rel_error
            reference_gap=$referenceBridgeGap
            moment_seconds=Get-Z4kDouble $aggregate.moment_build_seconds 0.0
            far_seconds=$farTime;m2l_seconds=Get-Z4kDouble $aggregate.m2l_seconds 0.0
            l2p_seconds=Get-Z4kDouble $aggregate.l2p_seconds 0.0
            local_target_leaves=[int64]$aggregate.local_target_leaves
            local_m2l_pairs=[int64]$aggregate.local_m2l_pairs
            m2p_speedup=$speedup;true_residual=$case.true_residual
            wall_seconds=$run.wall_seconds;peak_rss_mb=$run.peak_working_set_mb;run_dir=$runDir
        }) | Out-Null
    }
    return $rows.ToArray()
}

$fixtures=if($Quick){@('1sh1')}else{@('barnase','1a63','1sh1')}
$thetas=if($Quick){@(0.6)}else{@(0.6,0.45)}
$matrix=New-Object System.Collections.Generic.List[object]
$channels=New-Object System.Collections.Generic.List[object]
foreach($fixtureName in $fixtures){
    $fixture=Get-Z4kFixture $fixtureName
    foreach($theta in $thetas){
        foreach($row in @(Invoke-Z4kCase $fixture $theta)){$matrix.Add($row)|Out-Null}
        $runDir=$matrix[$matrix.Count-1].run_dir
        $summaryFile=Join-Path $runDir 'pbe_tabi_operator_attribution_summary.csv'
        if(Test-Path $summaryFile){foreach($row in @(Import-Csv $summaryFile)){$row|Add-Member -NotePropertyName fixture -NotePropertyValue $fixtureName;$row|Add-Member -NotePropertyName theta -NotePropertyValue $theta;$channels.Add($row)|Out-Null}}
    }
}
$rows=$matrix.ToArray()
$rows|Export-Csv $matrixCsv -NoTypeInformation
$channels|Export-Csv $channelCsv -NoTypeInformation
$rows|Select fixture,theta,approximation_mode,moment_seconds,far_seconds,m2l_seconds,l2p_seconds,local_target_leaves,local_m2l_pairs,m2p_speedup,wall_seconds,peak_rss_mb|Export-Csv $timingCsv -NoTypeInformation

$decision='z4k_local_accuracy_speed_pass'
if(@($rows|Where-Object status -eq 'timeout').Count){$decision='z4k_runtime_limited'}
else{
    $local06=@($rows|Where-Object {$_.approximation_mode -eq 'local_m2l_p4' -and [double]$_.theta -eq 0.6})
    if(@($local06|Where-Object {[double]$_.reference_gap -gt 1e-10 -or [double]$_.recomposition_error -gt 1e-10 -or [double]$_.local_vs_m2p_rel -gt 2e-3}).Count){$decision='z4k_local_translation_gap'}
    elseif(-not $Quick -and @($local06|Where-Object {[double]$_.audit_rms -gt 8e-3}).Count){$decision='z4k_local_theta06_accuracy_gap'}
    elseif(@($local06|Where-Object {[double]$_.m2p_speedup -lt 4}).Count){$decision='z4k_local_speed_gap'}
}
[pscustomobject]@{decision=$decision;quick=[bool]$Quick;rows=$rows.Count;audit_targets=$auditTargets;gmres_enabled=($decision -eq 'z4k_local_accuracy_speed_pass' -and -not $Quick)}|Export-Csv $decisionCsv -NoTypeInformation
@('scope: Z4K target-leaf local M2L/L2P gate',"decision: $decision","quick: $([bool]$Quick)","rows: $($rows.Count)",'source_order: 4','target_order: 4','kernel_derivative_order: 8','production_default: grid remains default')|Set-Content $summaryPath -Encoding ASCII

$svg=@('<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="600" viewBox="0 0 1000 600">','<rect width="1000" height="600" fill="white"/>','<text x="35" y="38" font-family="Arial" font-size="20">Target-leaf M2L/L2P accuracy and speed</text>')
$y=80
foreach($row in @($rows|Where-Object approximation_mode -eq 'local_m2l_p4')){$svg+=('<text x="35" y="{0}" font-family="Arial" font-size="13">{1} theta={2}: audit={3}, parity={4}, speedup={5}x</text>' -f $y,$row.fixture,$row.theta,$row.audit_rms,$row.local_vs_m2p_rel,$row.m2p_speedup);$y+=42}
$svg+='</svg>';$svg|Set-Content $plotPath -Encoding ASCII

$report=@('# TITAN FMM/M2L Scope Z4K Local M2L/L2P Report','',"Decision: $decision",'','Z4K evaluates p4 source moments translated into p4 target-leaf local expansions using order-8 kernel derivatives. Production GMRES remains unchanged unless the full theta-0.6 gate passes.','')
foreach($row in $rows){$report+=('- {0} theta={1} {2}: audit={3}, local_vs_m2p={4}, far_s={5}, speedup={6}' -f $row.fixture,$row.theta,$row.approximation_mode,$row.audit_rms,$row.local_vs_m2p_rel,$row.far_seconds,$row.m2p_speedup)}
$conclusion=switch($decision){
    'z4k_local_theta06_accuracy_gap' {'The theta-0.6 gate failed on at least one required protein. GMRES-p4 remains disabled. Translation parity is diagnostic evidence only.'}
    'z4k_local_speed_gap' {'Local accuracy passed, but order-8 M2L construction did not reach the required 4x speedup. GMRES-p4 remains disabled.'}
    'z4k_local_accuracy_speed_pass' {'The candidate accuracy and speed gates passed; guarded GMRES verification may proceed.'}
    default {"The candidate remains fail-closed with decision $decision."}
}
$report+=@('',$conclusion,'Grid remains the production Poisson-Boltzmann default.')
$report|Set-Content $rootReport -Encoding ASCII

$bad=@(Get-ChildItem $WorkRoot -Recurse -Filter arg.in|ForEach-Object{$p=$_.FullName;Get-Content $p|Where-Object{$_ -match '^\s*-' -and $_ -notmatch '^\s*-[^=\s]+='}|ForEach-Object{"$p :: $_"}})
if($bad.Count){throw "Generated arg.in contains bare options: $($bad -join '; ')"}
if(-not $AllowFailedDecision -and $decision -ne 'z4k_local_accuracy_speed_pass'){throw "Scope Z4K decision: $decision"}
Write-Host "PASS Scope Z4K decision=$decision rows=$($rows.Count)"
