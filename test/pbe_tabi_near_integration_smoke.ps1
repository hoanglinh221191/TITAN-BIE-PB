param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_near_integration_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$InputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Convert-ToDouble([object]$Value) {
    return [double]::Parse([string]$Value, $InvariantCulture)
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI near integration root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Write-DummyPdb([string]$Path) {
    Set-Content -LiteralPath $Path -Value @(
        'CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   VAL A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  VAL A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   VAL A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   VAL A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    ) -Encoding ASCII
}

function New-ArgFile([string]$RunDir, [string]$NearMode, [string]$Geometry = 'triangle') {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $pqrPath = Join-Path $InputRoot 'sphere.pqr'
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $argLines = @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        '-1pdb',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe=tabi',
        "-pbe-test-pqr=$(Convert-ToTitanPath $pqrPath)",
        '-pbe-surface=vdw',
        '-pbe-timing=detail',
        '-pbe-tabi-sum=direct',
        '-pbe-tabi-formulation=jcp2013',
        '-pbe-tabi-surface-points=48',
        "-pbe-tabi-geometry=$Geometry",
        "-pbe-tabi-near-integration=$NearMode",
        '-nt=1'
    )
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    return $argPath
}

function Invoke-NearCase([string]$Name, [string]$NearMode) {
    $runDir = Join-Path $WorkRoot $Name
    New-CleanDir $runDir
    $argPath = New-ArgFile -RunDir $runDir -NearMode $NearMode
    Push-Location $runDir
    try {
        & $TitanExe -a $argPath | Out-Null
        if ($LASTEXITCODE -ne 0) {
            throw "TITAN exited with code $LASTEXITCODE for near=$NearMode"
        }
    } finally {
        Pop-Location
    }
    $audit = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_formulation_audit.csv'))[0]
    $timing = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_timing.csv'))[0]
    $energy = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))[0]
    if ($energy.status -ne 'ok') {
        throw "TABI near case did not finish ok: $($energy.status)"
    }
    return [PSCustomObject]@{
        case = $Name
        requested = $NearMode
        near_integration = $audit.near_integration
        near_kernel_policy = $audit.near_kernel_policy
        self_pairs = $audit.near_self_pairs
        one_ring_pairs = $audit.near_one_ring_pairs
        duffy_pairs = $audit.near_duffy_pairs
        centroid_pairs = $audit.near_centroid_pairs
        fallbacks = $audit.near_fallbacks
        near_k1_accepted = if ($audit.PSObject.Properties.Name -contains 'near_k1_accepted') { $audit.near_k1_accepted } else { '0' }
        near_k2_accepted = if ($audit.PSObject.Properties.Name -contains 'near_k2_accepted') { $audit.near_k2_accepted } else { '0' }
        near_k3_accepted = if ($audit.PSObject.Properties.Name -contains 'near_k3_accepted') { $audit.near_k3_accepted } else { '0' }
        near_k4_accepted = if ($audit.PSObject.Properties.Name -contains 'near_k4_accepted') { $audit.near_k4_accepted } else { '0' }
        near_k1_rejected = if ($audit.PSObject.Properties.Name -contains 'near_k1_rejected') { $audit.near_k1_rejected } else { '0' }
        near_k2_rejected = if ($audit.PSObject.Properties.Name -contains 'near_k2_rejected') { $audit.near_k2_rejected } else { '0' }
        near_k3_rejected = if ($audit.PSObject.Properties.Name -contains 'near_k3_rejected') { $audit.near_k3_rejected } else { '0' }
        near_k4_rejected = if ($audit.PSObject.Properties.Name -contains 'near_k4_rejected') { $audit.near_k4_rejected } else { '0' }
        near_one_ring_k1_accepted = if ($audit.PSObject.Properties.Name -contains 'near_one_ring_k1_accepted') { $audit.near_one_ring_k1_accepted } else { '0' }
        near_one_ring_k2_accepted = if ($audit.PSObject.Properties.Name -contains 'near_one_ring_k2_accepted') { $audit.near_one_ring_k2_accepted } else { '0' }
        near_one_ring_k3_accepted = if ($audit.PSObject.Properties.Name -contains 'near_one_ring_k3_accepted') { $audit.near_one_ring_k3_accepted } else { '0' }
        near_one_ring_k4_accepted = if ($audit.PSObject.Properties.Name -contains 'near_one_ring_k4_accepted') { $audit.near_one_ring_k4_accepted } else { '0' }
        near_one_ring_k1_rejected = if ($audit.PSObject.Properties.Name -contains 'near_one_ring_k1_rejected') { $audit.near_one_ring_k1_rejected } else { '0' }
        near_one_ring_k2_rejected = if ($audit.PSObject.Properties.Name -contains 'near_one_ring_k2_rejected') { $audit.near_one_ring_k2_rejected } else { '0' }
        near_one_ring_k3_rejected = if ($audit.PSObject.Properties.Name -contains 'near_one_ring_k3_rejected') { $audit.near_one_ring_k3_rejected } else { '0' }
        near_one_ring_k4_rejected = if ($audit.PSObject.Properties.Name -contains 'near_one_ring_k4_rejected') { $audit.near_one_ring_k4_rejected } else { '0' }
        near_k1_policy = if ($audit.PSObject.Properties.Name -contains 'near_k1_policy') { $audit.near_k1_policy } else { 'unknown' }
        near_k2_policy = if ($audit.PSObject.Properties.Name -contains 'near_k2_policy') { $audit.near_k2_policy } else { 'unknown' }
        near_k3_policy = if ($audit.PSObject.Properties.Name -contains 'near_k3_policy') { $audit.near_k3_policy } else { 'unknown' }
        near_k4_policy = if ($audit.PSObject.Properties.Name -contains 'near_k4_policy') { $audit.near_k4_policy } else { 'unknown' }
        near_seconds = $audit.near_integration_seconds
        timing_near = $timing.tabi_near_integration
        energy_kJmol = $energy.E_rxn_kJmol
    }
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
New-CleanDir $WorkRoot

$auto = Invoke-NearCase -Name 'auto' -NearMode 'auto'
$centroid = Invoke-NearCase -Name 'centroid' -NearMode 'centroid'
$duffy = Invoke-NearCase -Name 'duffy' -NearMode 'duffy_adaptive'

if ($auto.near_integration -ne 'centroid') {
    throw "Expected auto to remain centroid until analytic ASI gates pass, got $($auto.near_integration)"
}
if ($centroid.near_integration -ne 'centroid') {
    throw "Expected centroid mode, got $($centroid.near_integration)"
}
if ($duffy.near_integration -ne 'duffy_adaptive') {
    throw "Expected duffy_adaptive mode, got $($duffy.near_integration)"
}
if ($auto.near_kernel_policy -ne 'centroid' -or $centroid.near_kernel_policy -ne 'centroid') {
    throw "Unexpected centroid policy: auto=$($auto.near_kernel_policy) centroid=$($centroid.near_kernel_policy)"
}
if ($duffy.near_kernel_policy -ne 'one_ring_cpv') {
    throw "Expected duffy_adaptive to report one_ring_cpv policy, got $($duffy.near_kernel_policy)"
}
if ([int64]$duffy.self_pairs -le 0 -or [int64]$duffy.duffy_pairs -le 0) {
    throw "Duffy TABI did not record self/duffy near pairs."
}
if ([int64]$duffy.fallbacks -ne 0) {
    throw "Duffy TABI unexpectedly used near-integration fallback count=$($duffy.fallbacks)"
}
foreach ($block in @('k1','k2','k3','k4')) {
    $acceptedName = "near_${block}_accepted"
    $rejectedName = "near_${block}_rejected"
    if (([int64]$duffy.$acceptedName + [int64]$duffy.$rejectedName) -le 0) {
        throw "Duffy TABI did not record any $block near-kernel decisions."
    }
    $oneRingAcceptedName = "near_one_ring_${block}_accepted"
    $oneRingRejectedName = "near_one_ring_${block}_rejected"
    if (([int64]$duffy.$oneRingAcceptedName + [int64]$duffy.$oneRingRejectedName) -le 0) {
        throw "Duffy TABI did not record any one-ring $block near-kernel decisions."
    }
}
if ((Convert-ToDouble $duffy.near_seconds) -le 0.0) {
    throw "Duffy TABI did not record near integration time."
}

$rejectDir = Join-Path $WorkRoot 'reject_point_panel_duffy'
New-CleanDir $rejectDir
$rejectArg = New-ArgFile -RunDir $rejectDir -NearMode 'duffy_adaptive' -Geometry 'point-panel'
Push-Location $rejectDir
try {
    $oldErrorActionPreference = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    try {
        & $TitanExe -a $rejectArg *> (Join-Path $rejectDir 'reject.out')
    } finally {
        $ErrorActionPreference = $oldErrorActionPreference
    }
    if ($LASTEXITCODE -eq 0) {
        throw 'Explicit duffy_adaptive near integration unexpectedly accepted point-panel geometry.'
    }
} finally {
    Pop-Location
}

$rows = @($auto, $centroid, $duffy)
$outCsv = Join-Path $WorkRoot 'pbe_tabi_near_integration_smoke.csv'
$rows | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
$global:LASTEXITCODE = 0
Write-Host "PBE TABI near integration smoke PASS: $outCsv"
