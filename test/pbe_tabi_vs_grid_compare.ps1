param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_vs_grid_compare'),
    [string[]]$Fixtures = @('sphere','two-atoms','arg','mini-protein','charged-peptide','compact-mini-protein'),
    [ValidateSet('default','barnase','paper2013','all')]
    [string]$FixtureSet = 'default',
    [string]$GoldenCsv = (Join-Path $PSScriptRoot 'pbe_inputs_local\pbe_delphi_golden.csv'),
    [string]$PyDelPhiRoot = (Join-Path (Split-Path -Parent $PSScriptRoot) '.codex_tmp\pyDelPhi'),
    [int]$RepeatCount = 1,
    [ValidateSet('accuracy','perf')]
    [string]$Mode = 'accuracy',
    [switch]$UseStage3USafePreset,
    [switch]$UseStage3UPerfPreset,
    [switch]$UseStage3UPerf2Preset,
    [switch]$UseStage3AIPerfProbe,
    [switch]$UseStage3AJBudgetedPreconditioner,
    [switch]$UseStage3AKCoarseProbe,
    [switch]$ProfileMemory,
    [int]$BarnaseTabiSurfacePoints = 24,
    [int]$CaseTimeoutSeconds = 0,
    [switch]$AllowTimeout,
    [switch]$UseDeferredLeafPairWorklist,
    [switch]$UseTargetLeafBlockBatching,
    [switch]$UseSoaLeafPairKernel,
    [string[]]$TabiExtraArgs = @(),
    [int]$Threads = 1
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$InputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'
$EffectiveThreads = if (($UseStage3USafePreset -or $UseStage3UPerfPreset -or $UseStage3UPerf2Preset -or $UseStage3AIPerfProbe -or $UseStage3AJBudgetedPreconditioner -or $UseStage3AKCoarseProbe) -and -not $PSBoundParameters.ContainsKey('Threads')) { 8 } else { [math]::Max(1, $Threads) }
$NormalizedTabiExtraArgs = @()
foreach ($extraArg in $TabiExtraArgs) {
    foreach ($part in ([string]$extraArg -split ',')) {
        $trimmed = $part.Trim()
        if (-not [string]::IsNullOrWhiteSpace($trimmed)) {
            $NormalizedTabiExtraArgs += $trimmed
        }
    }
}
$TabiExtraArgs = $NormalizedTabiExtraArgs

$FixtureMap = @{
    'sphere' = 'sphere.pqr'
    'two-atoms' = 'two-atoms.pqr'
    'arg' = 'arg.pqr'
    'mini-protein' = 'mini-protein.pqr'
    'charged-peptide' = 'charged-peptide.pqr'
    'compact-mini-protein' = 'compact-mini-protein.pqr'
    'research-medium' = 'rec_027.pqr'
    'research-large' = 'rec_0.pqr'
    '1a63' = '1a63.pqr'
    'barnase_linear_salt0' = 'barnase.pqr'
    'barnase_linear_salt015' = 'barnase.pqr'
    'barnase_linear_salt2' = 'barnase.pqr'
}

$BarnaseFixtures = @('barnase_linear_salt0','barnase_linear_salt015','barnase_linear_salt2')
$PaperFixtures = @('sphere')

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Convert-ToDoubleOrNull([object]$Value) {
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') {
        return $null
    }
    $parsed = 0.0
    if ([double]::TryParse($text, [System.Globalization.NumberStyles]::Float, $InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $null
}

function Get-RowValue([object]$Row, [string]$Name, [string]$Default = 'NA') {
    if ($null -ne $Row -and $Row.PSObject.Properties.Name -contains $Name) {
        $value = [string]$Row.$Name
        if (-not [string]::IsNullOrWhiteSpace($value)) {
            return $value
        }
    }
    return $Default
}

function Format-NumberOrNA([object]$Value) {
    $parsed = Convert-ToDoubleOrNull $Value
    if ($null -eq $parsed -or [double]::IsNaN($parsed) -or [double]::IsInfinity($parsed)) {
        return 'NA'
    }
    return $parsed.ToString('G17', $InvariantCulture)
}

function Format-DivideOrNA([object]$Numerator, [object]$Denominator) {
    $num = Convert-ToDoubleOrNull $Numerator
    $den = Convert-ToDoubleOrNull $Denominator
    if ($null -eq $num -or $null -eq $den -or [math]::Abs($den) -le 1.0e-300) {
        return 'NA'
    }
    return Format-NumberOrNA ($num / $den)
}

function Get-FirstRowValue([object]$Row, [string[]]$Names, [string]$Default = 'NA') {
    foreach ($name in $Names) {
        $value = Get-RowValue $Row $name ''
        if (-not [string]::IsNullOrWhiteSpace($value) -and $value -ne 'NA') {
            return $value
        }
    }
    return $Default
}

function Get-FirstRowValueAcross([object[]]$Rows, [string[]]$Names, [string]$Default = 'NA') {
    foreach ($row in $Rows) {
        foreach ($name in $Names) {
            $value = Get-RowValue $row $name ''
            if (-not [string]::IsNullOrWhiteSpace($value) -and $value -ne 'NA') {
                return $value
            }
        }
    }
    return $Default
}

function Get-LastCsvRow([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) { return $null }
    $rows = @(Import-Csv -LiteralPath $Path)
    if ($rows.Count -eq 0) { return $null }
    return $rows[$rows.Count - 1]
}

function Sum-RowNumericValues([object]$Row, [string[]]$Names) {
    $sum = 0.0
    $seen = $false
    foreach ($name in $Names) {
        $value = Convert-ToDoubleOrNull (Get-RowValue $Row $name 'NA')
        if ($null -ne $value) {
            $sum += $value
            $seen = $true
        }
    }
    if (-not $seen) { return 'NA' }
    return Format-NumberOrNA $sum
}

function Get-GoldenRow([string]$Fixture) {
    if (-not (Test-Path -LiteralPath $GoldenCsv)) { return $null }
    $rows = @(Import-Csv -LiteralPath $GoldenCsv | Where-Object { $_.fixture -eq $Fixture } | Select-Object -First 1)
    if ($rows.Count -eq 0) { return $null }
    return $rows[0]
}

function Read-DelPhiValueTable([string]$Path) {
    $values = @{}
    foreach ($line in Get-Content -LiteralPath $Path) {
        $trim = $line.Trim()
        if ([string]::IsNullOrWhiteSpace($trim) -or $trim.StartsWith('!') -or $trim.StartsWith('#')) { continue }
        $parts = @($trim -split '\s+' | Where-Object { $_ })
        if ($parts.Count -lt 2) { continue }
        $key = $parts[0].ToUpperInvariant()
        $number = 0.0
        for ($i = $parts.Count - 1; $i -ge 1; --$i) {
            if ([double]::TryParse($parts[$i], [System.Globalization.NumberStyles]::Float, $InvariantCulture, [ref]$number)) {
                $values[$key] = $number
                break
            }
        }
    }
    return $values
}

function Convert-DelPhiPdbToPqr([string]$PdbPath, [string]$CrgPath, [string]$SizPath, [string]$OutPath) {
    $charges = Read-DelPhiValueTable $CrgPath
    $radii = Read-DelPhiValueTable $SizPath
    $out = New-Object System.Collections.Generic.List[string]
    foreach ($line in Get-Content -LiteralPath $PdbPath) {
        if (-not ($line.StartsWith('ATOM') -or $line.StartsWith('HETATM'))) { continue }
        $atomName = $line.Substring(12, [Math]::Min(4, $line.Length - 12)).Trim().ToUpperInvariant()
        $charge = if ($charges.ContainsKey($atomName)) { $charges[$atomName] } else { 0.0 }
        $radius = if ($radii.ContainsKey($atomName)) { $radii[$atomName] } else { 1.5 }
        $prefix = if ($line.Length -ge 54) { $line.Substring(0, 54) } else { $line.PadRight(54) }
        $element = if ($line.Length -ge 78) { $line.Substring(76, 2).Trim() } else { $atomName.Substring(0, 1) }
        $out.Add(("{0}{1,9:F4}{2,8:F4}       {3}" -f $prefix, $charge, $radius, $element))
    }
    $out.Add('END')
    Set-Content -LiteralPath $OutPath -Value $out -Encoding ASCII
}

function Convert-ResearchPqr([string]$SourcePath, [string]$RunDir, [string]$OutputName) {
    if (-not (Test-Path -LiteralPath $SourcePath)) { return $SourcePath }
    $outPath = Join-Path $RunDir $OutputName
    $out = New-Object System.Collections.Generic.List[string]
    foreach ($line in Get-Content -LiteralPath $SourcePath) {
        $trim = $line.Trim()
        if ($trim.StartsWith('ATOM') -or $trim.StartsWith('HETATM')) {
            $tokens = @($trim -split '\s+' | Where-Object { $_ })
            if ($tokens.Count -ge 10) {
                $radius = 0.0
                if ([double]::TryParse($tokens[$tokens.Count - 1],
                        [System.Globalization.NumberStyles]::Float,
                        $InvariantCulture, [ref]$radius) -and $radius -le 0.0) {
                    $tokens[$tokens.Count - 1] = '0.600'
                    $out.Add(($tokens -join ' '))
                    continue
                }
            }
        }
        $out.Add($line)
    }
    Set-Content -LiteralPath $outPath -Value $out -Encoding ASCII
    return $outPath
}

function Resolve-ComparePqr([string]$Fixture, [string]$RunDir) {
    $pqrName = $FixtureMap[$Fixture]
    if ([string]::IsNullOrWhiteSpace($pqrName)) { throw "Unknown fixture: $Fixture" }
    $localPqr = Join-Path $InputRoot $pqrName
    if (Test-Path -LiteralPath $localPqr) { return $localPqr }
    if ($Fixture -eq 'research-medium') {
        return Convert-ResearchPqr -SourcePath (Join-Path $RepoRoot 'example\1\rec_027.pqr') -RunDir $RunDir -OutputName 'research-medium_normalized.pqr'
    }
    if ($Fixture -eq 'research-large') {
        return Convert-ResearchPqr -SourcePath (Join-Path $RepoRoot 'example\4\rec_0.pqr') -RunDir $RunDir -OutputName 'research-large_normalized.pqr'
    }
    if ($Fixture -notlike 'barnase*') { return $localPqr }
    $fixtureDir = Join-Path $PyDelPhiRoot 'pydelphi\data\test_examples\barnase'
    $directPqr = Join-Path $fixtureDir $pqrName
    if (Test-Path -LiteralPath $directPqr) { return $directPqr }
    $pdbPath = Join-Path $fixtureDir 'barnase.pdb'
    $crgPath = Join-Path $fixtureDir 'amber.crg'
    $sizPath = Join-Path $fixtureDir 'amber.siz'
    if ((Test-Path -LiteralPath $pdbPath) -and
        (Test-Path -LiteralPath $crgPath) -and
        (Test-Path -LiteralPath $sizPath)) {
        $generated = Join-Path $RunDir $pqrName
        Convert-DelPhiPdbToPqr -PdbPath $pdbPath -CrgPath $crgPath -SizPath $sizPath -OutPath $generated
        return $generated
    }
    return $localPqr
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootFull = [System.IO.Path]::GetFullPath($Root).TrimEnd('\') + '\'
    if (-not $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to clean path outside TABI-vs-grid root: $full"
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
    $lines = @(
        'CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    )
    Set-Content -LiteralPath $Path -Value $lines -Encoding ASCII
}

function New-ArgFile([string]$RunDir, [string]$Backend, [string]$PqrPath, [string]$RunMode, [string]$Fixture) {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $gold = Get-GoldenRow $Fixture
    $salt = Get-RowValue $gold 'salt_M' '0.15'
    $isResearch = $Fixture -eq 'research-medium' -or $Fixture -eq 'research-large'
    $scale = if ($Fixture -like 'barnase*') { Get-RowValue $gold 'scale' '2.0' } else { '0.35' }
    $gridSize = if ($Fixture -like 'barnase*') { Get-RowValue $gold 'grid_size' '' } else { '' }
    $surface = if ($Fixture -like 'barnase*') { 'vdw' } else { 'vdw' }
    $margin = if ($Fixture -like 'barnase*') { '0.0' } else { '2.0' }
    $maxGrid = if ($Fixture -like 'barnase*' -or $isResearch) { '4000000' } else { '300000' }
    $maxIter = if ($Fixture -like 'barnase*' -or $isResearch) { '1200' } else { '800' }
    $rms = if ($Fixture -like 'barnase*' -or $isResearch) { '5e-4' } else { '1e-4' }
    $argLines = @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        '-1pdb',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        "-pbe=$Backend",
        "-pbe-test-pqr=$(Convert-ToTitanPath $PqrPath)",
        "-pbe-surface=$surface",
        '-pbe-probe=1.4',
        "-pbe-scale=$scale",
        "-pbe-margin=$margin",
        "-pbe-max-grid=$maxGrid",
        '-pbe-map=off',
        "-pbe-max-iter=$maxIter",
        "-pbe-rms=$rms",
        '-pbe-dphi=1e-4',
        '-pbe-timing=summary',
        "-pbe-salt=$salt",
        "-nt=$EffectiveThreads"
    )
    if (-not [string]::IsNullOrWhiteSpace($gridSize)) {
        $argLines += "-pbe-test-grid-size=$gridSize"
    }
    if ($Backend -eq 'tabi') {
        $argLines += '-pbe-tabi-derivatives=analytic'
        if ($UseStage3USafePreset -or $UseStage3UPerfPreset -or $UseStage3UPerf2Preset -or $UseStage3AIPerfProbe -or $UseStage3AJBudgetedPreconditioner -or $UseStage3AKCoarseProbe) {
            $presetName = if ($UseStage3AIPerfProbe -or $UseStage3AJBudgetedPreconditioner -or $UseStage3AKCoarseProbe -or $UseStage3UPerf2Preset) { 'stage3u-perf2' } elseif ($UseStage3UPerfPreset) { 'stage3u-perf' } else { 'stage3u-safe' }
            $argLines += "-pbe-tabi-tree-preset=$presetName"
            if ($Fixture -like 'barnase*') {
                $argLines += "-pbe-tabi-surface-points=$BarnaseTabiSurfacePoints"
            } else {
                $argLines += '-pbe-tabi-surface-points=48'
            }
            $argLines += '-pbe-tabi-energy-projection=boundary'
            if ($CaseTimeoutSeconds -gt 0 -or $Fixture -like 'barnase*') {
                $argLines += '-pbe-tabi-pre-solver-progress=on'
                $argLines += '-pbe-tabi-progress=auto'
                $argLines += '-pbe-tabi-progress-interval=5'
            }
        } elseif ($Fixture -like 'barnase*') {
            $argLines += '-pbe-tabi-sum=treecode'
            $argLines += '-pbe-tabi-formulation=jcp2013-treecode'
            $argLines += "-pbe-tabi-surface-points=$BarnaseTabiSurfacePoints"
            $argLines += '-pbe-tabi-energy-projection=boundary'
        } elseif ($RunMode -eq 'accuracy') {
            $argLines += '-pbe-tabi-sum=direct'
            $argLines += '-pbe-tabi-formulation=jcp2013'
            $argLines += '-pbe-tabi-energy-projection=boundary'
        } else {
            $argLines += '-pbe-tabi-sum=treecode'
            $argLines += '-pbe-tabi-formulation=jcp2013-treecode'
            $argLines += '-pbe-tabi-energy-projection=boundary'
        }
        if ($UseStage3AIPerfProbe -or $UseStage3AJBudgetedPreconditioner -or $UseStage3AKCoarseProbe) {
            $argLines += '-pbe-tabi-matvec-cost-profile=summary'
            $argLines += '-pbe-tabi-leaf-kernel-profile=summary'
            $argLines += '-pbe-tabi-leaf-pair-scheduler=target-leaf-block'
            $argLines += '-pbe-tabi-jcp-direct-kernel=fused-soa-avx2'
            $argLines += '-pbe-tabi-leaf-pair-min-source-panels=64'
        }
        if ($UseStage3AJBudgetedPreconditioner) {
            $argLines += '-pbe-tabi-precond=block-budgeted'
            $argLines += '-pbe-tabi-precond-block-near=buffer-aware'
            $argLines += '-pbe-tabi-precond-block-panel-cap=32'
            $argLines += '-pbe-tabi-precond-memory-budget-mb=256'
            $argLines += '-pbe-tabi-gmres-restart=64'
        }
        if ($UseStage3AKCoarseProbe) {
            $argLines += '-pbe-tabi-precond=coarse-galerkin-probe'
            $argLines += '-pbe-tabi-precond-coarse-aggregate-cap=16'
            $argLines += '-pbe-tabi-precond-coarse-base=diagonal'
            $argLines += '-pbe-tabi-precond-block-near=buffer-aware'
            $argLines += '-pbe-tabi-precond-memory-budget-mb=256'
            $argLines += '-pbe-tabi-gmres-restart=64'
        }
        foreach ($extraArg in $TabiExtraArgs) {
            if (-not [string]::IsNullOrWhiteSpace($extraArg)) {
                $argLines += $extraArg
            }
        }
        if ($UseDeferredLeafPairWorklist) {
            $argLines += '-pbe-tabi-leaf-pair-scheduler=deferred-worklist'
            $argLines += '-pbe-tabi-leaf-pair-worklist-max-active=2048'
            $argLines += '-pbe-tabi-leaf-pair-worklist-verify=off'
        }
        if ($UseTargetLeafBlockBatching) {
            $argLines += '-pbe-tabi-leaf-pair-scheduler=target-leaf-block'
            $argLines += '-pbe-tabi-leaf-block-verify=off'
        }
        if ($UseSoaLeafPairKernel) {
            $argLines += '-pbe-tabi-leaf-kernel-profile=summary'
            $argLines += '-pbe-tabi-jcp-direct-kernel=soa-avx2'
        }
    }
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    return $argPath
}

function Invoke-TitanWithTimeout([string]$ExePath, [string]$ArgPath, [string]$RunDir, [int]$TimeoutSeconds) {
    $stdoutPath = Join-Path $RunDir 'titan_stdout.txt'
    $stderrPath = Join-Path $RunDir 'titan_stderr.txt'
    $psi = New-Object System.Diagnostics.ProcessStartInfo
    $psi.FileName = $ExePath
    $psi.Arguments = ('-a "{0}"' -f $ArgPath)
    $psi.WorkingDirectory = $RunDir
    $psi.UseShellExecute = $false
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    $process = New-Object System.Diagnostics.Process
    $process.StartInfo = $psi
    $sw = [System.Diagnostics.Stopwatch]::StartNew()
    [void]$process.Start()
    $completed = if ($TimeoutSeconds -gt 0) {
        $process.WaitForExit($TimeoutSeconds * 1000)
    } else {
        $process.WaitForExit()
        $true
    }
    if (-not $completed) {
        try {
            $process.Kill($true)
        } catch {
            try { $process.Kill() } catch {}
        }
    }
    $stdout = $process.StandardOutput.ReadToEnd()
    $stderr = $process.StandardError.ReadToEnd()
    $sw.Stop()
    Set-Content -LiteralPath $stdoutPath -Value $stdout -Encoding ASCII
    Set-Content -LiteralPath $stderrPath -Value $stderr -Encoding ASCII
    $exitCode = if ($completed) { $process.ExitCode } else { 'timeout' }
    $process.Dispose()
    return [PSCustomObject]@{
        completed = $completed
        exit_code = $exitCode
        total_seconds = $sw.Elapsed.TotalSeconds
        stdout_path = $stdoutPath
        stderr_path = $stderrPath
    }
}

function Invoke-PbeCase([int]$Repeat, [string]$Fixture, [string]$Backend) {
    $runDir = Join-Path $WorkRoot ("repeat_{0}\{1}_{2}" -f $Repeat, $Fixture, $Backend)
    New-CleanDir $runDir
    $pqrPath = Resolve-ComparePqr -Fixture $Fixture -RunDir $runDir
    if (-not (Test-Path -LiteralPath $pqrPath)) {
        throw "Missing PQR fixture: $pqrPath"
    }
    $argPath = New-ArgFile -RunDir $runDir -Backend $Backend -PqrPath $pqrPath -RunMode $Mode -Fixture $Fixture
    $runResult = Invoke-TitanWithTimeout -ExePath $TitanExe -ArgPath $argPath -RunDir $runDir -TimeoutSeconds $CaseTimeoutSeconds
    if ((-not $runResult.completed) -and (-not $AllowTimeout)) {
        throw "TITAN timed out after $CaseTimeoutSeconds seconds for $Fixture/$Backend"
    }
    if ($runResult.completed -and $runResult.exit_code -ne 0) {
        throw "TITAN exited with code $($runResult.exit_code) for $Fixture/$Backend"
    }

    $energyPath = Join-Path $runDir 'pbe_energy.csv'
    $timingPath = Join-Path $runDir 'pbe_timing.csv'
    $energy = Get-LastCsvRow $energyPath
    $timing = Get-LastCsvRow $timingPath
    $longProgress = Get-LastCsvRow (Join-Path $runDir 'pbe_tabi_longrun_progress.csv')
    $prepareProgress = Get-LastCsvRow (Join-Path $runDir 'pbe_tabi_prepare_progress.csv')
    if (($null -eq $energy -or $null -eq $timing) -and (-not ($AllowTimeout -and -not $runResult.completed))) {
        foreach ($required in @($energyPath, $timingPath)) {
            if (-not (Test-Path -LiteralPath $required)) {
                throw "Missing PBE output for ${Fixture}/${Backend}: $required"
            }
        }
    }
    $memory = $null
    $memoryPath = Join-Path $runDir 'pbe_memory_profile.csv'
    if ($ProfileMemory -or (Test-Path -LiteralPath $memoryPath)) {
        if (-not (Test-Path -LiteralPath $memoryPath)) {
            throw "Missing PBE memory profile for ${Fixture}/${Backend}: $memoryPath"
        }
        $memory = Get-LastCsvRow $memoryPath
    }
    $audit = $null
    $auditPath = Join-Path $runDir 'pbe_tabi_formulation_audit.csv'
    if ($Backend -eq 'tabi' -and $runResult.completed) {
        if (-not (Test-Path -LiteralPath $auditPath)) {
            throw "Missing TABI formulation audit for ${Fixture}/${Backend}: $auditPath"
        }
        $audit = Get-LastCsvRow $auditPath
    }
    if ($null -ne $energy -and $energy.status -ne 'ok' -and $runResult.completed) {
        throw "PBE case did not finish ok for ${Fixture}/${Backend}: $($energy.status)"
    }

    $energyKJ = Convert-ToDoubleOrNull (Get-RowValue $energy 'E_rxn_kJmol' 'NA')
    if ($runResult.completed -and ($null -eq $energyKJ -or [double]::IsNaN($energyKJ) -or [double]::IsInfinity($energyKJ))) {
        throw "Non-finite PBE energy for ${Fixture}/${Backend}: $(Get-RowValue $energy 'E_rxn_kJmol' 'NA')"
    }
    $defaultTabiTreePreset = if ($UseStage3AIPerfProbe -or $UseStage3AJBudgetedPreconditioner -or $UseStage3AKCoarseProbe -or $UseStage3UPerf2Preset) { 'stage3u-perf2' } elseif ($UseStage3UPerfPreset) { 'stage3u-perf' } elseif ($UseStage3USafePreset) { 'stage3u-safe' } else { 'off' }
    $metricRows = @($timing, $longProgress) | Where-Object { $null -ne $_ }
    $caseStatus = if (-not $runResult.completed) {
        'timeout'
    } else {
        Get-RowValue $energy 'status' 'missing_output'
    }

    return [PSCustomObject]@{
        repeat = $Repeat
        fixture = $Fixture
        backend = $Backend
        completed = if ($runResult.completed) { 'true' } else { 'false' }
        exit_code = [string]$runResult.exit_code
        progress_phase = if ($Backend -eq 'tabi') { Get-RowValue $longProgress 'phase' 'NA' } else { 'grid' }
        tabi_prepare_progress_phase = if ($Backend -eq 'tabi') { Get-RowValue $prepareProgress 'phase' 'NA' } else { 'grid' }
        mode = Get-RowValue $timing 'mode' (Get-RowValue $energy 'solver_mode' 'NA')
        atom_count = Get-RowValue $energy 'atom_count' 'NA'
        grid_points = Get-RowValue $energy 'grid_points' 'NA'
        panels = if ($Backend -eq 'tabi') { Get-RowValue $timing 'tabi_surface_panels' '0' } else { '0' }
        tabi_tree_preset = if ($Backend -eq 'tabi') { Get-RowValue $timing 'tabi_tree_preset' $defaultTabiTreePreset } else { 'grid' }
        energy_kJmol = Get-RowValue $energy 'E_rxn_kJmol' 'NA'
        delta_vs_grid_abs = 'NA'
        delta_vs_grid_rel = 'NA'
        energy_ratio_vs_grid = 'NA'
        time_ratio_tabi_over_grid = 'NA'
        absolute_time_delta_vs_grid = 'NA'
        energy_delta_tabi_minus_grid = 'NA'
        energy_delta_rel_vs_grid = 'NA'
        peak_memory_ratio_tabi_over_grid = 'NA'
        rss_ratio_tabi_over_grid = 'NA'
        kirkwood_error_rel = 'NA'
        oracle_engine = 'none'
        oracle_status = 'not_requested'
        formulation = if ($Backend -eq 'tabi') { Get-RowValue $audit 'formulation' 'NA' } else { 'grid' }
        energy_projection = if ($Backend -eq 'tabi') { Get-RowValue $audit 'energy_projection' 'NA' } else { 'grid' }
        derivative_mode = if ($Backend -eq 'tabi') { Get-RowValue $audit 'derivative_mode' (Get-RowValue $timing 'tabi_derivative_mode' 'NA') } else { 'finite_difference' }
        near_integration = if ($Backend -eq 'tabi') { Get-RowValue $audit 'near_integration' (Get-RowValue $timing 'tabi_near_integration' 'NA') } else { 'grid' }
        near_kernel_policy = if ($Backend -eq 'tabi') { Get-RowValue $audit 'near_kernel_policy' (Get-RowValue $timing 'tabi_near_kernel_policy' 'NA') } else { 'grid' }
        near_k1_policy = if ($Backend -eq 'tabi') { Get-RowValue $audit 'near_k1_policy' (Get-RowValue $timing 'tabi_near_k1_policy' 'NA') } else { 'grid' }
        near_k2_policy = if ($Backend -eq 'tabi') { Get-RowValue $audit 'near_k2_policy' (Get-RowValue $timing 'tabi_near_k2_policy' 'NA') } else { 'grid' }
        near_k3_policy = if ($Backend -eq 'tabi') { Get-RowValue $audit 'near_k3_policy' (Get-RowValue $timing 'tabi_near_k3_policy' 'NA') } else { 'grid' }
        near_k4_policy = if ($Backend -eq 'tabi') { Get-RowValue $audit 'near_k4_policy' (Get-RowValue $timing 'tabi_near_k4_policy' 'NA') } else { 'grid' }
        fallback_used = if ($Backend -eq 'tabi') { Get-RowValue $audit 'fallback_used' 'NA' } else { 'false' }
        total_seconds = Get-RowValue $timing 'total_seconds' ([math]::Round($runResult.total_seconds, 6).ToString($InvariantCulture))
        setup_seconds = Get-RowValue $timing 'setup_seconds' 'NA'
        surface_seconds = Get-RowValue $timing 'surface_seconds' 'NA'
        solver_seconds = Get-RowValue $timing 'solver_seconds' 'NA'
        energy_seconds = Get-RowValue $timing 'energy_seconds' 'NA'
        tabi_surface_build_seconds = Get-RowValue $timing 'tabi_surface_build_seconds' 'NA'
        tabi_surface_neighbor_seconds = Get-RowValue $timing 'tabi_surface_neighbor_seconds' 'NA'
        tabi_surface_mesh_seconds = Get-RowValue $timing 'tabi_surface_mesh_seconds' 'NA'
        tabi_solid_angle_seconds = Get-RowValue $timing 'tabi_solid_angle_seconds' 'NA'
        tabi_rhs_seconds = Get-RowValue $timing 'tabi_rhs_seconds' 'NA'
        tabi_gmres_seconds = Get-RowValue $timing 'tabi_gmres_seconds' 'NA'
        tabi_matvec_seconds_total = Get-RowValue $timing 'tabi_matvec_seconds_total' 'NA'
        tabi_matvec_calls = Get-FirstRowValueAcross $metricRows @('tabi_matvec_calls','matvec_calls') 'NA'
        tabi_tree_build_seconds = Get-FirstRowValueAcross $metricRows @('tabi_tree_build_seconds') 'NA'
        tabi_tree_moment_seconds = Get-FirstRowValueAcross $metricRows @('tabi_tree_moment_seconds','tree_moment_seconds') 'NA'
        tabi_tree_traversal_seconds = Get-FirstRowValueAcross $metricRows @('tabi_tree_traversal_seconds','tree_traversal_seconds') 'NA'
        tabi_matvec_cost_profile = Get-FirstRowValueAcross $metricRows @('tabi_matvec_cost_profile','matvec_cost_profile') 'NA'
        tabi_leaf_panels_per_second = Get-FirstRowValueAcross $metricRows @('tabi_leaf_panels_per_second','leaf_panels_per_second') 'NA'
        tabi_clean_leaf_panels_per_second = Get-FirstRowValueAcross $metricRows @('tabi_clean_leaf_panels_per_second','clean_leaf_panels_per_second') 'NA'
        tabi_tree_traversal_fraction = Get-FirstRowValueAcross $metricRows @('tabi_tree_traversal_fraction','tree_traversal_fraction') 'NA'
        tabi_preconditioner_requested = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_requested','preconditioner_requested') 'NA' } else { 'grid' }
        tabi_preconditioner_effective = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_effective','preconditioner_effective') 'NA' } else { 'grid' }
        tabi_preconditioner_block_near = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_block_near','preconditioner_block_near') 'NA' } else { 'grid' }
        tabi_preconditioner_block_panel_cap = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_block_panel_cap','preconditioner_block_panel_cap') 'NA' } else { 'NA' }
        tabi_preconditioner_memory_budget_bytes = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_memory_budget_bytes','preconditioner_memory_budget_bytes') 'NA' } else { 'NA' }
        tabi_preconditioner_estimated_bytes = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_estimated_bytes','preconditioner_estimated_bytes') 'NA' } else { 'NA' }
        tabi_preconditioner_valid_blocks = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_valid_blocks','preconditioner_valid_blocks') 'NA' } else { 'NA' }
        tabi_preconditioner_valid_panels = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_valid_panels','preconditioner_valid_panels') 'NA' } else { 'NA' }
        tabi_preconditioner_ill_conditioned_blocks = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_ill_conditioned_blocks','preconditioner_ill_conditioned_blocks') 'NA' } else { 'NA' }
        tabi_preconditioner_setup_seconds = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_setup_seconds','preconditioner_setup_seconds') 'NA' } else { 'NA' }
        tabi_preconditioner_apply_seconds = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_apply_seconds','preconditioner_apply_seconds') 'NA' } else { 'NA' }
        tabi_preconditioner_apply_calls = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_apply_calls','preconditioner_apply_calls') 'NA' } else { 'NA' }
        tabi_preconditioner_coarse_aggregate_cap = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_coarse_aggregate_cap','preconditioner_coarse_aggregate_cap') 'NA' } else { 'NA' }
        tabi_preconditioner_coarse_aggregates = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_coarse_aggregates','preconditioner_coarse_aggregates') 'NA' } else { 'NA' }
        tabi_preconditioner_coarse_dimension = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_coarse_dimension','preconditioner_coarse_dimension') 'NA' } else { 'NA' }
        tabi_preconditioner_coarse_base = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_coarse_base','preconditioner_coarse_base') 'NA' } else { 'NA' }
        tabi_preconditioner_coarse_matrix_bytes = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_coarse_matrix_bytes','preconditioner_coarse_matrix_bytes') 'NA' } else { 'NA' }
        tabi_preconditioner_coarse_pivot_ratio = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_coarse_pivot_ratio','preconditioner_coarse_pivot_ratio') 'NA' } else { 'NA' }
        tabi_preconditioner_coarse_status = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_coarse_status','preconditioner_coarse_status') 'NA' } else { 'NA' }
        tabi_preconditioner_coarse_setup_matvec_calls = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_coarse_setup_matvec_calls','preconditioner_coarse_setup_matvec_calls') 'NA' } else { 'NA' }
        tabi_preconditioner_coarse_setup_seconds = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_coarse_setup_seconds','preconditioner_coarse_setup_seconds') 'NA' } else { 'NA' }
        tabi_preconditioner_coarse_apply_seconds = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_coarse_apply_seconds','preconditioner_coarse_apply_seconds') 'NA' } else { 'NA' }
        tabi_preconditioner_coarse_apply_calls = if ($Backend -eq 'tabi') { Get-FirstRowValueAcross $metricRows @('tabi_preconditioner_coarse_apply_calls','preconditioner_coarse_apply_calls') 'NA' } else { 'NA' }
        tabi_tree_nodes = Get-FirstRowValueAcross $metricRows @('tabi_tree_nodes','tree_nodes') 'NA'
        tabi_tree_leaf_count = Get-RowValue $timing 'tabi_tree_leaf_count' 'NA'
        tabi_near_interactions = Get-RowValue $timing 'tabi_near_interactions' 'NA'
        tabi_singular_near_interactions = Get-RowValue $timing 'tabi_singular_near_interactions' 'NA'
        tabi_leaf_near_interactions = Get-RowValue $timing 'tabi_leaf_near_interactions' 'NA'
        tabi_far_node_interactions = Get-FirstRowValueAcross $metricRows @('tabi_far_node_interactions','far_node_interactions') 'NA'
        tabi_leaf_panel_interactions = Get-FirstRowValueAcross $metricRows @('tabi_leaf_panel_interactions','leaf_panel_interactions') 'NA'
        tabi_leaf_tile_interactions = Get-RowValue $timing 'tabi_leaf_tile_interactions' 'NA'
        tabi_leaf_tile_panel_interactions = Get-RowValue $timing 'tabi_leaf_tile_panel_interactions' 'NA'
        tabi_leaf_residual_panel_interactions = Get-RowValue $timing 'tabi_leaf_residual_panel_interactions' 'NA'
        tabi_smooth_leaf_node_interactions = Get-RowValue $timing 'tabi_smooth_leaf_node_interactions' 'NA'
        tabi_leaf_pruned_interactions = Get-RowValue $timing 'tabi_leaf_pruned_interactions' 'NA'
        tabi_leaf_near_policy = Get-RowValue $timing 'tabi_leaf_near_policy' 'NA'
        tabi_leaf_panel_kernel = Get-RowValue $timing 'tabi_leaf_panel_kernel' 'NA'
        tabi_leaf_tile_status = Get-RowValue $timing 'tabi_leaf_tile_status' 'NA'
        tabi_leaf_panel_eval_seconds = Get-RowValue $timing 'tabi_leaf_panel_eval_seconds' 'NA'
        tabi_leaf_tile_eval_seconds = Get-RowValue $timing 'tabi_leaf_tile_eval_seconds' 'NA'
        tabi_leaf_residual_eval_seconds = Get-RowValue $timing 'tabi_leaf_residual_eval_seconds' 'NA'
        tabi_smooth_leaf_node_eval_seconds = Get-RowValue $timing 'tabi_smooth_leaf_node_eval_seconds' 'NA'
        tabi_leaf_panel_pairs_per_second = Get-RowValue $timing 'tabi_leaf_panel_pairs_per_second' 'NA'
        tabi_singular_near_eval_seconds = Get-RowValue $timing 'tabi_singular_near_eval_seconds' 'NA'
        tabi_leaf_near_eval_seconds = Get-RowValue $timing 'tabi_leaf_near_eval_seconds' 'NA'
        tabi_near_self_pairs = Get-RowValue $timing 'tabi_near_self_pairs' 'NA'
        tabi_near_one_ring_pairs = Get-RowValue $timing 'tabi_near_one_ring_pairs' 'NA'
        tabi_near_duffy_pairs = Get-RowValue $timing 'tabi_near_duffy_pairs' 'NA'
        tabi_near_centroid_pairs = Get-RowValue $timing 'tabi_near_centroid_pairs' 'NA'
        tabi_near_fallbacks = Get-RowValue $timing 'tabi_near_fallbacks' 'NA'
        tabi_near_k1_accepted = Get-RowValue $timing 'tabi_near_k1_accepted' 'NA'
        tabi_near_k1_rejected = Get-RowValue $timing 'tabi_near_k1_rejected' 'NA'
        tabi_near_k2_accepted = Get-RowValue $timing 'tabi_near_k2_accepted' 'NA'
        tabi_near_k2_rejected = Get-RowValue $timing 'tabi_near_k2_rejected' 'NA'
        tabi_near_k3_accepted = Get-RowValue $timing 'tabi_near_k3_accepted' 'NA'
        tabi_near_k3_rejected = Get-RowValue $timing 'tabi_near_k3_rejected' 'NA'
        tabi_near_k4_accepted = Get-RowValue $timing 'tabi_near_k4_accepted' 'NA'
        tabi_near_k4_rejected = Get-RowValue $timing 'tabi_near_k4_rejected' 'NA'
        tabi_near_one_ring_k1_accepted = Get-RowValue $timing 'tabi_near_one_ring_k1_accepted' 'NA'
        tabi_near_one_ring_k1_rejected = Get-RowValue $timing 'tabi_near_one_ring_k1_rejected' 'NA'
        tabi_near_one_ring_k2_accepted = Get-RowValue $timing 'tabi_near_one_ring_k2_accepted' 'NA'
        tabi_near_one_ring_k2_rejected = Get-RowValue $timing 'tabi_near_one_ring_k2_rejected' 'NA'
        tabi_near_one_ring_k3_accepted = Get-RowValue $timing 'tabi_near_one_ring_k3_accepted' 'NA'
        tabi_near_one_ring_k3_rejected = Get-RowValue $timing 'tabi_near_one_ring_k3_rejected' 'NA'
        tabi_near_one_ring_k4_accepted = Get-RowValue $timing 'tabi_near_one_ring_k4_accepted' 'NA'
        tabi_near_one_ring_k4_rejected = Get-RowValue $timing 'tabi_near_one_ring_k4_rejected' 'NA'
        tabi_near_integration_seconds = Get-RowValue $timing 'tabi_near_integration_seconds' 'NA'
        tabi_far_interactions = Get-RowValue $timing 'tabi_far_interactions' 'NA'
        tabi_leaf_pair_batch_interactions = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_batch_interactions','leaf_pair_batch_interactions') 'NA'
        tabi_leaf_pair_batch_targets = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_batch_targets','leaf_pair_batch_targets') 'NA'
        tabi_leaf_pair_batch_source_panels = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_batch_source_panels','leaf_pair_batch_source_panels') 'NA'
        tabi_leaf_pair_batch_residual_panels = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_batch_residual_panels','leaf_pair_batch_residual_panels') 'NA'
        tabi_leaf_pair_batch_fallbacks = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_batch_fallbacks','leaf_pair_batch_fallbacks') 'NA'
        tabi_leaf_pair_batch_seconds = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_batch_seconds','leaf_pair_batch_seconds') 'NA'
        tabi_leaf_pair_batch_kernel = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_batch_kernel','leaf_pair_batch_kernel') 'NA'
        tabi_leaf_pair_batch_vectorized_sources = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_batch_vectorized_sources','leaf_pair_batch_vectorized_sources') 'NA'
        tabi_leaf_pair_batch_vector_lanes = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_batch_vector_lanes','leaf_pair_batch_vector_lanes') 'NA'
        tabi_leaf_pair_batch_pairs_per_second = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_batch_pairs_per_second','leaf_pair_batch_pairs_per_second') 'NA'
        tabi_leaf_pair_selection_candidate_scans = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_selection_candidate_scans','leaf_pair_selection_candidate_scans') 'NA'
        tabi_leaf_pair_selection_clean_rejects = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_selection_clean_rejects','leaf_pair_selection_clean_rejects') 'NA'
        tabi_leaf_pair_selection_probe_calls = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_selection_probe_calls','leaf_pair_selection_probe_calls') 'NA'
        tabi_leaf_pair_selection_probe_node_visits = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_selection_probe_node_visits','leaf_pair_selection_probe_node_visits') 'NA'
        tabi_leaf_pair_selection_probe_seconds = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_selection_probe_seconds','leaf_pair_selection_probe_seconds') 'NA'
        tabi_leaf_pair_selection_seconds = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_selection_seconds','leaf_pair_selection_seconds') 'NA'
        tabi_leaf_pair_selection_probe_cache_hits = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_selection_probe_cache_hits','leaf_pair_selection_probe_cache_hits') 'NA'
        tabi_leaf_pair_selection_probe_cache_misses = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_selection_probe_cache_misses','leaf_pair_selection_probe_cache_misses') 'NA'
        tabi_leaf_pair_selection_insufficient_targets = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_selection_insufficient_targets','leaf_pair_selection_insufficient_targets') 'NA'
        tabi_leaf_pair_selection_same_target_leaf_rejects = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_selection_same_target_leaf_rejects','leaf_pair_selection_same_target_leaf_rejects') 'NA'
        tabi_leaf_pair_directness_probe = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_directness_probe','leaf_pair_directness_probe') 'NA'
        tabi_leaf_pair_directness_cache = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_directness_cache','leaf_pair_directness_cache') 'NA'
        tabi_leaf_pair_directness_verify = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_directness_verify','leaf_pair_directness_verify') 'NA'
        tabi_leaf_pair_directness_ancestor_probe_calls = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_directness_ancestor_probe_calls','leaf_pair_directness_ancestor_probe_calls') 'NA'
        tabi_leaf_pair_directness_ancestor_probe_node_visits = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_directness_ancestor_probe_node_visits','leaf_pair_directness_ancestor_probe_node_visits') 'NA'
        tabi_leaf_pair_directness_ancestor_probe_seconds = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_directness_ancestor_probe_seconds','leaf_pair_directness_ancestor_probe_seconds') 'NA'
        tabi_leaf_pair_directness_pair_cache_hits = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_directness_pair_cache_hits','leaf_pair_directness_pair_cache_hits') 'NA'
        tabi_leaf_pair_directness_pair_cache_misses = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_directness_pair_cache_misses','leaf_pair_directness_pair_cache_misses') 'NA'
        tabi_leaf_pair_directness_all_direct_proofs = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_directness_all_direct_proofs','leaf_pair_directness_all_direct_proofs') 'NA'
        tabi_leaf_pair_directness_mixed_fallbacks = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_directness_mixed_fallbacks','leaf_pair_directness_mixed_fallbacks') 'NA'
        tabi_leaf_pair_directness_verify_checks = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_directness_verify_checks','leaf_pair_directness_verify_checks') 'NA'
        tabi_leaf_pair_directness_verify_mismatches = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_directness_verify_mismatches','leaf_pair_directness_verify_mismatches') 'NA'
        tabi_leaf_pair_scheduler = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_scheduler','leaf_pair_scheduler') 'NA'
        tabi_leaf_pair_worklist_max_active = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_worklist_max_active','leaf_pair_worklist_max_active') 'NA'
        tabi_leaf_pair_worklist_verify = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_worklist_verify','leaf_pair_worklist_verify') 'NA'
        tabi_leaf_pair_worklist_tasks = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_worklist_tasks','leaf_pair_worklist_tasks') 'NA'
        tabi_leaf_pair_worklist_groups = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_worklist_groups','leaf_pair_worklist_groups') 'NA'
        tabi_leaf_pair_worklist_batched_groups = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_worklist_batched_groups','leaf_pair_worklist_batched_groups') 'NA'
        tabi_leaf_pair_worklist_scalar_flush_tasks = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_worklist_scalar_flush_tasks','leaf_pair_worklist_scalar_flush_tasks') 'NA'
        tabi_leaf_pair_worklist_forced_flushes = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_worklist_forced_flushes','leaf_pair_worklist_forced_flushes') 'NA'
        tabi_leaf_pair_worklist_active_groups_max = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_worklist_active_groups_max','leaf_pair_worklist_active_groups_max') 'NA'
        tabi_leaf_pair_worklist_seconds = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_worklist_seconds','leaf_pair_worklist_seconds') 'NA'
        tabi_leaf_pair_worklist_verify_checks = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_worklist_verify_checks','leaf_pair_worklist_verify_checks') 'NA'
        tabi_leaf_pair_worklist_verify_mismatches = Get-FirstRowValueAcross $metricRows @('tabi_leaf_pair_worklist_verify_mismatches','leaf_pair_worklist_verify_mismatches') 'NA'
        tabi_leaf_block_verify = Get-FirstRowValueAcross $metricRows @('tabi_leaf_block_verify','leaf_block_verify') 'NA'
        tabi_leaf_block_target_leaves = Get-FirstRowValueAcross $metricRows @('tabi_leaf_block_target_leaves','leaf_block_target_leaves') 'NA'
        tabi_leaf_block_targets = Get-FirstRowValueAcross $metricRows @('tabi_leaf_block_targets','leaf_block_targets') 'NA'
        tabi_leaf_block_source_groups = Get-FirstRowValueAcross $metricRows @('tabi_leaf_block_source_groups','leaf_block_source_groups') 'NA'
        tabi_leaf_block_batched_groups = Get-FirstRowValueAcross $metricRows @('tabi_leaf_block_batched_groups','leaf_block_batched_groups') 'NA'
        tabi_leaf_block_scalar_tail_tasks = Get-FirstRowValueAcross $metricRows @('tabi_leaf_block_scalar_tail_tasks','leaf_block_scalar_tail_tasks') 'NA'
        tabi_leaf_block_group_seconds = Get-FirstRowValueAcross $metricRows @('tabi_leaf_block_group_seconds','leaf_block_group_seconds') 'NA'
        tabi_leaf_block_flush_seconds = Get-FirstRowValueAcross $metricRows @('tabi_leaf_block_flush_seconds','leaf_block_flush_seconds') 'NA'
        tabi_leaf_block_writeback_seconds = Get-FirstRowValueAcross $metricRows @('tabi_leaf_block_writeback_seconds','leaf_block_writeback_seconds') 'NA'
        tabi_leaf_block_verify_checks = Get-FirstRowValueAcross $metricRows @('tabi_leaf_block_verify_checks','leaf_block_verify_checks') 'NA'
        tabi_leaf_block_verify_mismatches = Get-FirstRowValueAcross $metricRows @('tabi_leaf_block_verify_mismatches','leaf_block_verify_mismatches') 'NA'
        tabi_status = Get-RowValue $timing 'tabi_status' 'NA'
        tabi_convergence_reason = Get-RowValue $timing 'tabi_convergence_reason' 'NA'
        tabi_residual = if ($Backend -eq 'tabi') { Get-FirstRowValue $timing @('tabi_true_residual_rel','tabi_sampled_true_residual_rel','tabi_final_residual') } else { 'NA' }
        tabi_far_expansion_mean = if ($Backend -eq 'tabi') { Get-RowValue $timing 'tabi_far_expansion_rel_error_mean' 'NA' } else { 'NA' }
        tabi_accepted_missed_total = if ($Backend -eq 'tabi') { Sum-RowNumericValues $timing @('tabi_convergence_k1_high_error_missed','tabi_convergence_k2_high_error_missed','tabi_convergence_k3_high_error_missed','tabi_convergence_k4_high_error_missed') } else { 'NA' }
        reaction_energy_source = Get-RowValue $energy 'reaction_energy_source' 'NA'
        memory_estimated_bytes_total = Get-RowValue $memory 'estimated_bytes_total' 'NA'
        memory_estimated_bytes_surface = Get-RowValue $memory 'estimated_bytes_surface' 'NA'
        memory_estimated_bytes_solver = Get-RowValue $memory 'estimated_bytes_solver' 'NA'
        memory_estimated_bytes_operator = Get-RowValue $memory 'estimated_bytes_operator' 'NA'
        memory_estimated_bytes_tree = Get-RowValue $memory 'estimated_bytes_tree' 'NA'
        memory_rss_before_bytes = Get-RowValue $memory 'rss_before_bytes' 'NA'
        memory_rss_after_bytes = Get-RowValue $memory 'rss_after_bytes' 'NA'
        memory_rss_delta_bytes = Get-RowValue $memory 'rss_delta_bytes' 'NA'
        memory_rss_delta_signed_bytes = Get-RowValue $memory 'rss_delta_signed_bytes' 'NA'
        memory_peak_working_set_bytes = Get-RowValue $memory 'peak_working_set_bytes' 'NA'
        status = $caseStatus
    }
}

function Get-MedianText([object[]]$Values) {
    $numbers = New-Object System.Collections.Generic.List[double]
    foreach ($value in $Values) {
        $parsed = Convert-ToDoubleOrNull $value
        if ($null -ne $parsed) {
            $numbers.Add($parsed)
        }
    }
    if ($numbers.Count -eq 0) {
        return 'NA'
    }
    $sorted = @($numbers | Sort-Object)
    $mid = [int]($sorted.Count / 2)
    $median = if (($sorted.Count % 2) -eq 1) {
        $sorted[$mid]
    } else {
        0.5 * ($sorted[$mid - 1] + $sorted[$mid])
    }
    return ([math]::Round($median, 9)).ToString($InvariantCulture)
}

function New-MedianRows([object[]]$Rows) {
    $numericColumns = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::OrdinalIgnoreCase)
    foreach ($name in @(
        'atom_count','grid_points','panels','energy_kJmol',
        'delta_vs_grid_abs','delta_vs_grid_rel','energy_ratio_vs_grid',
        'time_ratio_tabi_over_grid','absolute_time_delta_vs_grid',
        'energy_delta_tabi_minus_grid',
        'energy_delta_rel_vs_grid','peak_memory_ratio_tabi_over_grid',
        'rss_ratio_tabi_over_grid','kirkwood_error_rel','total_seconds',
        'setup_seconds','surface_seconds','solver_seconds','energy_seconds',
        'tabi_surface_build_seconds','tabi_surface_neighbor_seconds',
        'tabi_surface_mesh_seconds','tabi_solid_angle_seconds',
        'tabi_rhs_seconds','tabi_gmres_seconds','tabi_matvec_seconds_total',
        'tabi_matvec_calls','tabi_tree_build_seconds',
        'tabi_tree_moment_seconds','tabi_tree_traversal_seconds',
        'tabi_leaf_panels_per_second',
        'tabi_clean_leaf_panels_per_second',
        'tabi_tree_traversal_fraction',
        'tabi_preconditioner_block_panel_cap',
        'tabi_preconditioner_memory_budget_bytes',
        'tabi_preconditioner_estimated_bytes',
        'tabi_preconditioner_valid_blocks',
        'tabi_preconditioner_valid_panels',
        'tabi_preconditioner_ill_conditioned_blocks',
        'tabi_preconditioner_setup_seconds',
        'tabi_preconditioner_apply_seconds',
        'tabi_preconditioner_apply_calls',
        'tabi_preconditioner_coarse_aggregate_cap',
        'tabi_preconditioner_coarse_aggregates',
        'tabi_preconditioner_coarse_dimension',
        'tabi_preconditioner_coarse_matrix_bytes',
        'tabi_preconditioner_coarse_pivot_ratio',
        'tabi_preconditioner_coarse_setup_matvec_calls',
        'tabi_preconditioner_coarse_setup_seconds',
        'tabi_preconditioner_coarse_apply_seconds',
        'tabi_preconditioner_coarse_apply_calls',
        'tabi_tree_nodes','tabi_tree_leaf_count','tabi_near_interactions',
        'tabi_singular_near_interactions','tabi_leaf_near_interactions',
        'tabi_far_node_interactions','tabi_leaf_panel_interactions',
        'tabi_leaf_tile_interactions','tabi_leaf_tile_panel_interactions',
        'tabi_leaf_residual_panel_interactions',
        'tabi_smooth_leaf_node_interactions','tabi_leaf_pruned_interactions',
        'tabi_leaf_panel_eval_seconds',
        'tabi_leaf_tile_eval_seconds','tabi_leaf_residual_eval_seconds',
        'tabi_smooth_leaf_node_eval_seconds',
        'tabi_leaf_panel_pairs_per_second',
        'tabi_singular_near_eval_seconds',
        'tabi_leaf_near_eval_seconds',
        'tabi_near_self_pairs','tabi_near_one_ring_pairs',
        'tabi_near_duffy_pairs','tabi_near_centroid_pairs',
        'tabi_near_fallbacks','tabi_near_k1_accepted',
        'tabi_near_k1_rejected','tabi_near_k2_accepted',
        'tabi_near_k2_rejected','tabi_near_k3_accepted',
        'tabi_near_k3_rejected','tabi_near_k4_accepted',
        'tabi_near_k4_rejected','tabi_near_one_ring_k1_accepted',
        'tabi_near_one_ring_k1_rejected',
        'tabi_near_one_ring_k2_accepted',
        'tabi_near_one_ring_k2_rejected',
        'tabi_near_one_ring_k3_accepted',
        'tabi_near_one_ring_k3_rejected',
        'tabi_near_one_ring_k4_accepted',
        'tabi_near_one_ring_k4_rejected',
        'tabi_near_integration_seconds',
        'tabi_far_interactions','tabi_residual','tabi_far_expansion_mean',
        'tabi_accepted_missed_total','tabi_leaf_pair_batch_interactions',
        'tabi_leaf_pair_batch_targets','tabi_leaf_pair_batch_source_panels',
        'tabi_leaf_pair_batch_residual_panels',
        'tabi_leaf_pair_batch_fallbacks','tabi_leaf_pair_batch_seconds',
        'tabi_leaf_pair_batch_vectorized_sources',
        'tabi_leaf_pair_batch_vector_lanes',
        'tabi_leaf_pair_batch_pairs_per_second',
        'tabi_leaf_pair_selection_candidate_scans',
        'tabi_leaf_pair_selection_clean_rejects',
        'tabi_leaf_pair_selection_probe_calls',
        'tabi_leaf_pair_selection_probe_node_visits',
        'tabi_leaf_pair_selection_probe_seconds',
        'tabi_leaf_pair_selection_seconds',
        'tabi_leaf_pair_selection_probe_cache_hits',
        'tabi_leaf_pair_selection_probe_cache_misses',
        'tabi_leaf_pair_selection_insufficient_targets',
        'tabi_leaf_pair_selection_same_target_leaf_rejects',
        'tabi_leaf_pair_directness_ancestor_probe_calls',
        'tabi_leaf_pair_directness_ancestor_probe_node_visits',
        'tabi_leaf_pair_directness_ancestor_probe_seconds',
        'tabi_leaf_pair_directness_pair_cache_hits',
        'tabi_leaf_pair_directness_pair_cache_misses',
        'tabi_leaf_pair_directness_all_direct_proofs',
        'tabi_leaf_pair_directness_mixed_fallbacks',
        'tabi_leaf_pair_directness_verify_checks',
        'tabi_leaf_pair_directness_verify_mismatches',
        'memory_estimated_bytes_total',
        'memory_estimated_bytes_surface','memory_estimated_bytes_solver',
        'memory_estimated_bytes_operator','memory_estimated_bytes_tree',
        'memory_rss_before_bytes','memory_rss_after_bytes',
        'memory_rss_delta_bytes','memory_rss_delta_signed_bytes',
        'memory_peak_working_set_bytes')) {
        [void]$numericColumns.Add($name)
    }
    $out = New-Object System.Collections.Generic.List[object]
    foreach ($group in ($Rows | Group-Object -Property fixture,backend)) {
        $first = $group.Group[0]
        $ordered = [ordered]@{}
        foreach ($property in $first.PSObject.Properties.Name) {
            if ($property -eq 'repeat') {
                continue
            }
            if ($numericColumns.Contains($property)) {
                $ordered[$property] = Get-MedianText @($group.Group | ForEach-Object { $_.$property })
            } else {
                $ordered[$property] = $first.$property
            }
        }
        $ordered['timing_stat'] = 'median'
        $ordered['repeat_count'] = $group.Count
        $out.Add([PSCustomObject]$ordered)
    }
    return $out
}

$presetSwitches = @($UseStage3USafePreset, $UseStage3UPerfPreset, $UseStage3UPerf2Preset, $UseStage3AIPerfProbe, $UseStage3AJBudgetedPreconditioner, $UseStage3AKCoarseProbe) | Where-Object { $_ }
if (@($presetSwitches).Count -gt 1) {
    throw 'Use only one of -UseStage3USafePreset, -UseStage3UPerfPreset, -UseStage3UPerf2Preset, -UseStage3AIPerfProbe, -UseStage3AJBudgetedPreconditioner, or -UseStage3AKCoarseProbe.'
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
if (-not (Test-Path -LiteralPath $InputRoot)) {
    throw "Missing PBE fixture directory: $InputRoot"
}

if ($FixtureSet -eq 'barnase') {
    $Fixtures = $BarnaseFixtures
} elseif ($FixtureSet -eq 'paper2013') {
    $Fixtures = $PaperFixtures
    if (Test-Path -LiteralPath (Join-Path $InputRoot '1a63.pqr')) {
        $Fixtures += '1a63'
    } else {
        Write-Warning "Skipping paper2013 1A63 fixture: missing test\pbe_inputs_local\1a63.pqr"
    }
} elseif ($FixtureSet -eq 'all') {
    $Fixtures = @('sphere','two-atoms','arg','mini-protein','charged-peptide','compact-mini-protein') + $BarnaseFixtures + $PaperFixtures
    if (Test-Path -LiteralPath (Join-Path $InputRoot '1a63.pqr')) {
        $Fixtures += '1a63'
    }
} elseif (($FixtureSet -eq 'default') -and (-not $PSBoundParameters.ContainsKey('Fixtures')) -and
          (Test-Path -LiteralPath (Join-Path $InputRoot '1a63.pqr'))) {
    $Fixtures += '1a63'
}
$Fixtures = @($Fixtures | ForEach-Object {
    ([string]$_ -split ',') | ForEach-Object { $_.Trim() } | Where-Object { $_ }
} | Select-Object -Unique)
if ($Fixtures.Count -eq 0) {
    throw 'No fixtures selected.'
}
$RepeatCount = [math]::Max(1, $RepeatCount)
$WorkRoot = (New-Item -ItemType Directory -Path $WorkRoot -Force).FullName

$rawRows = New-Object System.Collections.Generic.List[object]
for ($repeat = 1; $repeat -le $RepeatCount; ++$repeat) {
    foreach ($fixture in $Fixtures) {
        $grid = Invoke-PbeCase -Repeat $repeat -Fixture $fixture -Backend 'grid'
        $tabi = Invoke-PbeCase -Repeat $repeat -Fixture $fixture -Backend 'tabi'
        $gridEnergy = Convert-ToDoubleOrNull $grid.energy_kJmol
        $tabiEnergy = Convert-ToDoubleOrNull $tabi.energy_kJmol
        $grid.delta_vs_grid_abs = '0'
        $grid.delta_vs_grid_rel = '0'
        $grid.energy_ratio_vs_grid = '1'
        $grid.time_ratio_tabi_over_grid = '1'
        $grid.absolute_time_delta_vs_grid = '0'
        $grid.energy_delta_tabi_minus_grid = '0'
        $grid.energy_delta_rel_vs_grid = '0'
        $grid.peak_memory_ratio_tabi_over_grid = '1'
        $grid.rss_ratio_tabi_over_grid = '1'
        if ($null -ne $gridEnergy -and $null -ne $tabiEnergy) {
            $delta = [math]::Abs($tabiEnergy - $gridEnergy)
            $tabi.delta_vs_grid_abs = $delta.ToString($InvariantCulture)
            $tabi.delta_vs_grid_rel = ($delta / [math]::Max(1.0, [math]::Abs($gridEnergy))).ToString($InvariantCulture)
            $ratioDenom = if ([math]::Abs($gridEnergy) -gt 1.0e-12) { $gridEnergy } else { 1.0e-12 }
            $tabi.energy_ratio_vs_grid = ($tabiEnergy / $ratioDenom).ToString($InvariantCulture)
            $tabi.energy_delta_tabi_minus_grid = Format-NumberOrNA ($tabiEnergy - $gridEnergy)
            $tabi.energy_delta_rel_vs_grid = Format-DivideOrNA ([math]::Abs($tabiEnergy - $gridEnergy)) ([math]::Max(1.0, [math]::Abs($gridEnergy)))
        }
        $tabi.time_ratio_tabi_over_grid = Format-DivideOrNA $tabi.total_seconds $grid.total_seconds
        $tabi.absolute_time_delta_vs_grid = Format-NumberOrNA (([double]$tabi.total_seconds) - ([double]$grid.total_seconds))
        $tabi.peak_memory_ratio_tabi_over_grid = Format-DivideOrNA $tabi.memory_peak_working_set_bytes $grid.memory_peak_working_set_bytes
        $tabi.rss_ratio_tabi_over_grid = Format-DivideOrNA $tabi.memory_rss_after_bytes $grid.memory_rss_after_bytes
        if ($fixture -eq 'sphere' -and $tabi.status -eq 'ok') {
            $tabiAudit = @(Import-Csv -LiteralPath (Join-Path $WorkRoot ("repeat_{0}\{1}_{2}\pbe_tabi_formulation_audit.csv" -f $repeat, $fixture, 'tabi')))[0]
            $refEnergy = Convert-ToDoubleOrNull $tabiAudit.reference_energy_kJmol
            if ($null -ne $refEnergy) {
                $tabi.kirkwood_error_rel = ([math]::Abs($tabiEnergy - $refEnergy) / [math]::Max(1.0e-12, [math]::Abs($refEnergy))).ToString($InvariantCulture)
                $grid.kirkwood_error_rel = ([math]::Abs($gridEnergy - $refEnergy) / [math]::Max(1.0e-12, [math]::Abs($refEnergy))).ToString($InvariantCulture)
            }
        }
        $rawRows.Add($grid)
        $rawRows.Add($tabi)
        Write-Host ("PBE TABI-vs-grid {0} mode={1} repeat={2}: grid={3} kJ/mol tabi={4} kJ/mol abs_delta={5}" -f
            $fixture, $Mode, $repeat, $grid.energy_kJmol, $tabi.energy_kJmol, $tabi.delta_vs_grid_abs)
    }
}

$rawCsv = Join-Path $WorkRoot 'pbe_tabi_vs_grid_compare_repeats.csv'
$rawRows | Export-Csv -LiteralPath $rawCsv -NoTypeInformation -Encoding ASCII
$medianRows = New-MedianRows $rawRows.ToArray()
$outCsv = Join-Path $WorkRoot 'pbe_tabi_vs_grid_compare.csv'
$medianRows | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII

Write-Host "PBE TABI-vs-grid compare COMPLETE: $outCsv"
Write-Host "PBE TABI-vs-grid raw repeats: $rawCsv"
