$script:FmmCliInvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

function Convert-FmmCliToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Convert-FmmCliDouble([object]$Value, [double]$Default = [double]::NaN) {
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') {
        return $Default
    }
    return [double]::Parse($text, $script:FmmCliInvariantCulture)
}

function Convert-FmmCliInt64([object]$Value, [int64]$Default = 0) {
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') {
        return $Default
    }
    return [int64]::Parse($text, $script:FmmCliInvariantCulture)
}

function Test-FmmCliFinite([double]$Value) {
    return (-not [double]::IsNaN($Value)) -and
        (-not [double]::IsInfinity($Value))
}

function Convert-FmmCliResidualHistory([string]$History) {
    $rows = New-Object System.Collections.Generic.List[object]
    if ([string]::IsNullOrWhiteSpace($History)) {
        return @()
    }
    $tokens = @($History -split '[;| ]+' | Where-Object {
        -not [string]::IsNullOrWhiteSpace($_)
    })
    $iteration = 0
    foreach ($token in $tokens) {
        $valueText = $token
        if ($token -match '^[0-9]+[:=](.+)$') {
            $iteration = [int]($token -replace '[:=].+$', '')
            $valueText = $Matches[1]
        }
        $value = Convert-FmmCliDouble $valueText
        if (Test-FmmCliFinite $value) {
            $rows.Add([pscustomobject]@{
                iteration = $iteration
                residual = $value
            }) | Out-Null
            $iteration += 1
        }
    }
    return $rows.ToArray()
}

function Get-FmmCliSafeName([string]$Text) {
    return ($Text -replace '[^A-Za-z0-9_.-]', '_')
}

function Assert-FmmCliUnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootFull = [System.IO.Path]::GetFullPath($Root).TrimEnd('\') + '\'
    if (-not $full.StartsWith($rootFull,
            [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to clean path outside FMM/M2L work root: $full"
    }
}

function New-FmmCliCleanDir([string]$Path, [string]$WorkRoot) {
    Assert-FmmCliUnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Read-FmmCliDelPhiValueTable([string]$Path) {
    $values = @{}
    foreach ($line in Get-Content -LiteralPath $Path) {
        $trim = $line.Trim()
        if ([string]::IsNullOrWhiteSpace($trim) -or $trim.StartsWith('!') -or
            $trim.StartsWith('#')) {
            continue
        }
        $parts = @($trim -split '\s+' | Where-Object { $_ })
        if ($parts.Count -lt 2) {
            continue
        }
        $number = 0.0
        for ($i = $parts.Count - 1; $i -ge 1; --$i) {
            if ([double]::TryParse(
                    $parts[$i],
                    [System.Globalization.NumberStyles]::Float,
                    $script:FmmCliInvariantCulture,
                    [ref]$number)) {
                $values[$parts[0].ToUpperInvariant()] = $number
                break
            }
        }
    }
    return $values
}

function Convert-FmmCliDelPhiPdbToPqr(
    [string]$PdbPath,
    [string]$CrgPath,
    [string]$SizPath,
    [string]$OutPath
) {
    $charges = Read-FmmCliDelPhiValueTable $CrgPath
    $radii = Read-FmmCliDelPhiValueTable $SizPath
    $out = New-Object System.Collections.Generic.List[string]
    foreach ($line in Get-Content -LiteralPath $PdbPath) {
        if (-not ($line.StartsWith('ATOM') -or $line.StartsWith('HETATM'))) {
            continue
        }
        $atomName = $line.Substring(12, [Math]::Min(4, $line.Length - 12)).
            Trim().ToUpperInvariant()
        $charge = if ($charges.ContainsKey($atomName)) {
            $charges[$atomName]
        } else {
            0.0
        }
        $radius = if ($radii.ContainsKey($atomName)) {
            $radii[$atomName]
        } else {
            1.5
        }
        $prefix = if ($line.Length -ge 54) {
            $line.Substring(0, 54)
        } else {
            $line.PadRight(54)
        }
        $element = if ($line.Length -ge 78) {
            $line.Substring(76, 2).Trim()
        } else {
            $atomName.Substring(0, 1)
        }
        $out.Add(("{0}{1,9:F4}{2,8:F4}       {3}" -f
                $prefix, $charge, $radius, $element))
    }
    $out.Add('END')
    Set-Content -LiteralPath $OutPath -Value $out -Encoding ASCII
}

function Resolve-FmmCliFixture(
    [string]$Fixture,
    [string]$RepoRoot,
    [string]$RunDir = ''
) {
    switch ($Fixture) {
        'research-medium' { return Join-Path $RepoRoot 'example\1\rec_027.pqr' }
        'research-large' { return Join-Path $RepoRoot 'example\4\rec_0.pqr' }
        'research-alt' { return Join-Path $RepoRoot 'example\3\rec_00.pqr' }
        'barnase' {
            $localPqr = Join-Path $RepoRoot 'test\input\barnase.pqr'
            if (Test-Path -LiteralPath $localPqr) {
                return $localPqr
            }
            $fixtureDir = Join-Path $RepoRoot '.codex_tmp\pyDelPhi\pydelphi\data\test_examples\barnase'
            $directPqr = Join-Path $fixtureDir 'barnase.pqr'
            if (Test-Path -LiteralPath $directPqr) {
                return $directPqr
            }
            $pdbPath = Join-Path $fixtureDir 'barnase.pdb'
            $crgPath = Join-Path $fixtureDir 'amber.crg'
            $sizPath = Join-Path $fixtureDir 'amber.siz'
            if ((Test-Path -LiteralPath $pdbPath) -and
                (Test-Path -LiteralPath $crgPath) -and
                (Test-Path -LiteralPath $sizPath) -and
                -not [string]::IsNullOrWhiteSpace($RunDir)) {
                New-Item -ItemType Directory -Path $RunDir -Force | Out-Null
                $generated = Join-Path $RunDir 'barnase.pqr'
                Convert-FmmCliDelPhiPdbToPqr -PdbPath $pdbPath -CrgPath $crgPath `
                    -SizPath $sizPath -OutPath $generated
                return $generated
            }
            return $null
        }
        default {
            if (Test-Path -LiteralPath $Fixture) {
                return [System.IO.Path]::GetFullPath($Fixture)
            }
            throw "Unknown FMM/M2L fixture: $Fixture"
        }
    }
}

function Write-FmmCliDummyPdb([string]$Path) {
    $lines = @(
        'CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       0.000   1.450   0.000  1.00  0.00           C',
        'END'
    )
    Set-Content -LiteralPath $Path -Value $lines -Encoding ASCII
}

function Write-FmmCliSampledPqr(
    [string]$InputPath,
    [string]$OutputPath,
    [int]$Cap
) {
    $atomLines = @(Get-Content -LiteralPath $InputPath | Where-Object {
        $_.StartsWith('ATOM') -or $_.StartsWith('HETATM')
    })
    if ($atomLines.Count -eq 0) {
        throw "PQR fixture has no atoms: $InputPath"
    }
    $n = if ($Cap -le 0) {
        $atomLines.Count
    } else {
        [Math]::Min($Cap, $atomLines.Count)
    }
    $sampled = New-Object System.Collections.Generic.List[string]
    for ($i = 0; $i -lt $n; ++$i) {
        $idx = [Math]::Floor(($i * $atomLines.Count) / [double]$n)
        $fields = @($atomLines[[int]$idx] -split '\s+' | Where-Object {
            $_ -ne ''
        })
        if ($fields.Count -lt 10) {
            throw "Malformed PQR row selected from ${InputPath}: $($atomLines[[int]$idx])"
        }
        $numeric = New-Object System.Collections.Generic.List[double]
        foreach ($field in $fields) {
            $number = 0.0
            if ([double]::TryParse(
                    $field,
                    [System.Globalization.NumberStyles]::Float,
                    $script:FmmCliInvariantCulture,
                    [ref]$number)) {
                $numeric.Add($number)
            }
        }
        if ($numeric.Count -lt 5) {
            throw "Malformed PQR numeric fields selected from ${InputPath}: $($atomLines[[int]$idx])"
        }
        $last = $numeric.Count
        $x = $numeric[$last - 5]
        $y = $numeric[$last - 4]
        $z = $numeric[$last - 3]
        $q = $numeric[$last - 2]
        $radius = [Math]::Max(0.6, $numeric[$last - 1])
        $atomName = $fields[2]
        $resName = $fields[3]
        $chain = if ($fields.Count -ge 11) { $fields[4] } else { 'A' }
        $resnr = if ($fields.Count -ge 11) { $fields[5] } else { $fields[4] }
        $line = [string]::Format(
            $script:FmmCliInvariantCulture,
            'ATOM {0,6} {1,-4} {2,-3} {3,1} {4,5} {5,10:F3} {6,8:F3} {7,8:F3} {8,8:F3} {9,6:F3}',
            ($i + 1), $atomName, $resName, $chain, $resnr,
            $x, $y, $z, $q, $radius)
        $sampled.Add($line)
    }
    $sampled.Add('END')
    Set-Content -LiteralPath $OutputPath -Value $sampled -Encoding ASCII
    return [pscustomobject]@{
        input_atoms = $atomLines.Count
        sampled_atoms = $n
    }
}

function Write-FmmCliArgFile(
    [string]$RunDir,
    [string]$RepoRoot,
    [string]$PqrPath,
    [ValidateSet('tabi','grid')][string]$Backend = 'tabi',
    [ValidateSet('direct','treecode')][string]$Mode,
    [ValidateSet('off','single-matvec-analytic')][string]$FmmMode,
    [int]$SurfacePoints,
    [int]$MaxIter,
    [int]$Restart,
    [double]$Tolerance,
    [int]$Threads,
    [double]$SaltM,
    [double]$InternalDielectric,
    [double]$ExternalDielectric,
    [double]$TemperatureK,
    [double]$OpeningTheta = 0.8,
    [int]$AuditSampleCount = 0,
    [ValidateSet('auto','off','sampled','sampled-detail')]
    [string]$AuditMode = 'auto',
    [string]$Preconditioner = 'diagonal',
    [int]$PreconditionerBlockPanelCap = 32,
    [int]$PreconditionerOverlapPanelCap = 64,
    [int]$PreconditionerIluNeighborCap = 192,
    [double]$PreconditionerIluDropTol = 0.0,
    [ValidateSet('centroid','analytic_singular','duffy_adaptive')]
    [string]$NearIntegration = 'centroid',
    [ValidateSet(3,4)][int]$TaylorOrder = 4,
    [ValidateSet('off','on')][string]$FmmCartesianAudit = 'off',
    [ValidateSet('off','audit-p4','single-matvec-p4','gmres-p4')][string]$FmmLocal = 'off',
    [ValidateSet('off','p3-p4-tail','operator-tail-static')][string]$FmmLocalRefine = 'off',
    [double]$FmmLocalRefineTolerance = 4.0e-3,
    [ValidateSet('off','kernel-jet','green-jet-compressed')]
    [string]$FmmLocalCache = 'off',
    [ValidateRange(64,1024)][int]$FmmLocalCacheBudgetMb = 256,
    [ValidateSet('scalar','packed-avx2')][string]$FmmLocalFallback = 'scalar',
    [ValidateSet('legacy-rescan','metric-cache')]
    [string]$FmmLocalPlanBuilder = 'legacy-rescan',
    [ValidateSet(4,8,16)][int]$FmmLocalRefineProbeCount = 4,
    [ValidateSet('rademacher','solution','mixed')]
    [string]$FmmLocalRefineProbeMode = 'rademacher'
) {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-FmmCliDummyPdb $pdbPath
    $libPath = Convert-FmmCliToTitanPath(
        (Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $format = {
        param([double]$Value)
        return $Value.ToString('G15', $script:FmmCliInvariantCulture)
    }
    $argLines = @(
        "-f=$(Convert-FmmCliToTitanPath $pdbPath)",
        '-1pdb=true',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        "-pbe=$Backend",
        "-pbe-test-pqr=$(Convert-FmmCliToTitanPath $PqrPath)",
        '-pbe-surface=vdw',
        "-pbe-max-iter=$MaxIter",
        "-pbe-rms=$(& $format $Tolerance)",
        "-pbe-dphi=$(& $format $Tolerance)",
        "-pbe-int-diel=$(& $format $InternalDielectric)",
        "-pbe-ext-diel=$(& $format $ExternalDielectric)",
        "-pbe-salt=$(& $format $SaltM)",
        "-pbe-temp=$(& $format $TemperatureK)",
        '-pbe-timing=detail',
        "-nt=$Threads"
    )
    if ($Backend -eq 'tabi') {
        $argLines += @(
        "-pbe-tabi-sum=$Mode",
        '-pbe-tabi-formulation=jcp2013-treecode',
        "-pbe-tabi-fmm-m2l=$FmmMode",
        "-pbe-tabi-fmm-cartesian-audit=$FmmCartesianAudit",
        "-pbe-tabi-fmm-local=$FmmLocal",
        "-pbe-tabi-fmm-local-refine=$FmmLocalRefine",
        "-pbe-tabi-fmm-local-refine-tol=$(& $format $FmmLocalRefineTolerance)",
        "-pbe-tabi-fmm-local-cache=$FmmLocalCache",
        "-pbe-tabi-fmm-local-cache-budget-mb=$FmmLocalCacheBudgetMb",
        "-pbe-tabi-fmm-local-fallback=$FmmLocalFallback",
        "-pbe-tabi-fmm-local-plan-builder=$FmmLocalPlanBuilder",
        "-pbe-tabi-fmm-local-refine-probe-count=$FmmLocalRefineProbeCount",
        "-pbe-tabi-fmm-local-refine-probe-mode=$FmmLocalRefineProbeMode",
        "-pbe-tabi-surface-points=$SurfacePoints",
        '-pbe-tabi-tree-leaf-size=64',
        "-pbe-tabi-opening-theta=$(& $format $OpeningTheta)",
        "-pbe-tabi-near-integration=$NearIntegration",
        '-pbe-tabi-quad-buffer=off',
        "-pbe-tabi-taylor-order=$TaylorOrder",
        "-pbe-tabi-precond=$Preconditioner",
        "-pbe-tabi-precond-block-panel-cap=$PreconditionerBlockPanelCap",
        "-pbe-tabi-precond-overlap-panel-cap=$PreconditionerOverlapPanelCap",
        "-pbe-tabi-precond-ilu-neighbor-cap=$PreconditionerIluNeighborCap",
        "-pbe-tabi-precond-ilu-drop-tol=$(& $format $PreconditionerIluDropTol)",
        "-pbe-tabi-gmres-restart=$Restart",
        '-pbe-tabi-true-residual-probe=restart',
        '-pbe-tabi-true-residual-probe-scope=sampled',
        '-pbe-tabi-true-residual-sample-count=512'
        )
        $resolvedAuditMode = if ($AuditMode -eq 'auto') {
            if ($AuditSampleCount -gt 0) { 'sampled' } else { 'off' }
        } else { $AuditMode }
        if ($resolvedAuditMode -ne 'off') {
            if ($AuditSampleCount -le 0) {
                throw "Audit mode $resolvedAuditMode requires AuditSampleCount > 0"
            }
            $argLines += "-pbe-tabi-matvec-audit=$resolvedAuditMode"
            $argLines += "-pbe-tabi-matvec-audit-sample-count=$AuditSampleCount"
        } else {
            $argLines += '-pbe-tabi-matvec-audit=off'
        }
    }
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    $badOptions = @($argLines | Where-Object {
        -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*='
    })
    if ($badOptions.Count -ne 0) {
        throw "Generated arg.in contains non equals-syntax options: $($badOptions -join '; ')"
    }
    return $argPath
}

function Invoke-FmmCliTitanCase(
    [string]$TitanExe,
    [string]$RunDir,
    [string]$ArgPath,
    [int]$CaseTimeoutSeconds = 0
) {
    $stdoutPath = Join-Path $RunDir 'stdout.log'
    $stderrPath = Join-Path $RunDir 'stderr.log'
    $psi = New-Object System.Diagnostics.ProcessStartInfo
    $psi.FileName = $TitanExe
    $psi.WorkingDirectory = $RunDir
    $psi.UseShellExecute = $false
    $psi.CreateNoWindow = $true
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    $psi.Arguments = '-a "' + ($ArgPath -replace '"', '\"') + '"'
    $process = New-Object System.Diagnostics.Process
    $process.StartInfo = $psi
    $stopwatch = [System.Diagnostics.Stopwatch]::StartNew()
    [void]$process.Start()
    $stdoutTask = $process.StandardOutput.ReadToEndAsync()
    $stderrTask = $process.StandardError.ReadToEndAsync()
    $deadline = if ($CaseTimeoutSeconds -gt 0) {
        [DateTime]::UtcNow.AddSeconds($CaseTimeoutSeconds)
    } else {
        [DateTime]::MaxValue
    }
    $peakWorkingSetBytes = 0L
    $finished = $false
    while (-not $process.HasExited) {
        try {
            $process.Refresh()
            $peakWorkingSetBytes = [Math]::Max(
                $peakWorkingSetBytes,
                [int64]$process.WorkingSet64)
            $peakWorkingSetBytes = [Math]::Max(
                $peakWorkingSetBytes,
                [int64]$process.PeakWorkingSet64)
        } catch {}
        if ([DateTime]::UtcNow -ge $deadline) {
            break
        }
        Start-Sleep -Milliseconds 100
    }
    if ($process.HasExited) {
        $finished = $true
        try {
            $process.Refresh()
            $peakWorkingSetBytes = [Math]::Max(
                $peakWorkingSetBytes,
                [int64]$process.PeakWorkingSet64)
        } catch {}
    }
    $timedOut = -not $finished
    if ($timedOut) {
        try { $process.Kill() } catch {}
        $process.WaitForExit()
        try {
            $process.Refresh()
            $peakWorkingSetBytes = [Math]::Max(
                $peakWorkingSetBytes,
                [int64]$process.PeakWorkingSet64)
        } catch {}
    }
    $stdoutTask.Wait()
    $stderrTask.Wait()
    $stopwatch.Stop()
    Set-Content -LiteralPath $stdoutPath -Value $stdoutTask.Result -Encoding UTF8
    Set-Content -LiteralPath $stderrPath -Value $stderrTask.Result -Encoding UTF8
    return [pscustomobject]@{
        exit_code = if ($timedOut) { 124 } else { $process.ExitCode }
        timed_out = $timedOut
        wall_seconds = $stopwatch.Elapsed.TotalSeconds
        peak_working_set_bytes = $peakWorkingSetBytes
        peak_working_set_mb = $peakWorkingSetBytes / 1MB
    }
}

function Read-FmmCliCaseResult(
    [string]$RunDir,
    [string]$ExpectedFmmMode,
    [bool]$RequireAudit
) {
    $energyPath = Join-Path $RunDir 'pbe_energy.csv'
    $timingPath = Join-Path $RunDir 'pbe_timing.csv'
    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    foreach ($required in @($energyPath, $timingPath, $profilePath)) {
        if (-not (Test-Path -LiteralPath $required)) {
            throw "Missing FMM/M2L CLI output: $required"
        }
    }
    $energy = @(Import-Csv -LiteralPath $energyPath)[0]
    $timing = @(Import-Csv -LiteralPath $timingPath)[0]
    $profile = @(Import-Csv -LiteralPath $profilePath)[0]
    foreach ($column in @(
        'fmm_m2l_mode','fmm_m2l_near_pairs','fmm_m2l_far_pairs',
        'fmm_m2l_setup_seconds','fmm_m2l_apply_seconds',
        'fmm_m2l_finite_outputs','tree_vs_direct_matvec_rel',
        'matvec_audit_rows')) {
        if (-not ($profile.PSObject.Properties.Name -contains $column)) {
            throw "pbe_tabi_profile.csv missing FMM/M2L column: $column"
        }
    }
    $actualMode = [string]$profile.fmm_m2l_mode
    if ($ExpectedFmmMode -eq 'single-matvec-analytic') {
        if ($actualMode -ne 'single-matvec-analytic') {
            throw "FMM CLI row did not report FMM/M2L mode: $actualMode"
        }
        if ((Convert-FmmCliInt64 $profile.fmm_m2l_far_pairs) -le 0) {
            throw 'FMM CLI row reported zero far pairs.'
        }
        if ([string]$profile.fmm_m2l_finite_outputs -ne 'true') {
            throw "FMM CLI row reported non-finite outputs: $($profile.fmm_m2l_finite_outputs)"
        }
    } elseif (-not ([string]::IsNullOrWhiteSpace($actualMode) -or
            $actualMode -eq 'off')) {
        throw "Direct CLI row unexpectedly reported FMM/M2L mode: $actualMode"
    }
    if ([string]$energy.reaction_energy_source -ne 'tabi_boundary_integral') {
        throw "FMM/M2L case used fallback/source $($energy.reaction_energy_source)"
    }
    $auditRows = Convert-FmmCliInt64 $profile.matvec_audit_rows
    $auditRel = Convert-FmmCliDouble $profile.tree_vs_direct_matvec_rel
    if ($RequireAudit -and ($auditRows -le 0 -or -not (Test-FmmCliFinite $auditRel))) {
        throw "FMM/M2L sampled audit missing: rows=$auditRows rel=$auditRel"
    }
    return [pscustomobject]@{
        status = if ($energy.status -eq 'ok' -and $profile.status -eq 'ok') {
            'ok'
        } else {
            "energy=$($energy.status);profile=$($profile.status)"
        }
        panels = Convert-FmmCliInt64 $profile.panels
        iterations = Convert-FmmCliInt64 $profile.iterations
        final_residual = Convert-FmmCliDouble $profile.final_residual
        true_residual = Convert-FmmCliDouble $profile.true_residual_rel
        residual_history = [string]$profile.residual_history
        energy_kj_mol = Convert-FmmCliDouble $energy.E_rxn_kJmol
        matvec_seconds_avg = Convert-FmmCliDouble $profile.matvec_seconds_avg
        total_seconds = Convert-FmmCliDouble $profile.total_seconds
        fmm_m2l_mode = $actualMode
        fmm_near_pairs = Convert-FmmCliInt64 $profile.fmm_m2l_near_pairs
        fmm_far_pairs = Convert-FmmCliInt64 $profile.fmm_m2l_far_pairs
        fmm_setup_seconds = Convert-FmmCliDouble $profile.fmm_m2l_setup_seconds 0.0
        fmm_apply_seconds = Convert-FmmCliDouble $profile.fmm_m2l_apply_seconds 0.0
        fmm_finite_outputs = [string]$profile.fmm_m2l_finite_outputs
        tree_vs_direct_matvec_rel = $auditRel
        matvec_audit_rows = $auditRows
        matvec_audit_seconds = Convert-FmmCliDouble $profile.matvec_audit_seconds 0.0
        reaction_energy_source = [string]$energy.reaction_energy_source
        timing_backend = [string]$timing.pbe_backend
        memory_bytes = Convert-FmmCliDouble $profile.memory_bytes
    }
}

function Get-FmmCliFailureClass(
    [int]$DirectExit,
    [int]$FmmExit,
    [AllowNull()][object]$Direct,
    [AllowNull()][object]$Fmm,
    [string]$Fixture,
    [bool]$RequireSpeed,
    [bool]$RequireAudit = $true,
    [double]$MinimumSpeedup = 1.2,
    [double]$AuditGate = 1.0e-2
) {
    if ($DirectExit -eq 124 -or $FmmExit -eq 124) {
        return 'timeout_not_verified'
    }
    if ($DirectExit -ne 0 -or $FmmExit -ne 0 -or $null -eq $Direct -or
        $null -eq $Fmm) {
        return 'status_fail'
    }
    if ($Direct.status -ne 'ok' -or $Fmm.status -ne 'ok') {
        return 'status_fail'
    }
    if ($Fmm.fmm_m2l_mode -ne 'single-matvec-analytic' -or
        $Fmm.fmm_far_pairs -le 0 -or $Fmm.fmm_finite_outputs -ne 'true') {
        return 'telemetry_fail'
    }
    foreach ($value in @(
        $Direct.energy_kj_mol, $Fmm.energy_kj_mol,
        $Direct.true_residual, $Fmm.true_residual,
        $Direct.matvec_seconds_avg, $Fmm.matvec_seconds_avg)) {
        if (-not (Test-FmmCliFinite $value)) {
            return 'nonfinite_fail'
        }
    }
    $residualGate = [Math]::Max(1.0e-4,
        10.0 * [Math]::Max(1.0e-12, [Math]::Abs($Direct.true_residual)))
    if ([Math]::Abs($Fmm.true_residual) -gt $residualGate) {
        return 'convergence_fail'
    }
    if ($RequireAudit -and ($Fmm.matvec_audit_rows -le 0 -or
        -not (Test-FmmCliFinite $Fmm.tree_vs_direct_matvec_rel) -or
        $Fmm.tree_vs_direct_matvec_rel -gt $AuditGate)) {
        return 'accuracy_fail'
    }
    $directMagnitude = [Math]::Max(1.0, [Math]::Abs($Direct.energy_kj_mol))
    $energyGate = if ($Fixture -eq 'research-large') {
        [Math]::Max(5.0, 0.08 * $directMagnitude)
    } else {
        [Math]::Max(2.0, 0.05 * $directMagnitude)
    }
    if ([Math]::Abs($Fmm.energy_kj_mol - $Direct.energy_kj_mol) -gt
        $energyGate) {
        return 'energy_fail'
    }
    $speedup = $Direct.matvec_seconds_avg / $Fmm.matvec_seconds_avg
    if ($RequireSpeed -and $speedup -lt $MinimumSpeedup) {
        return 'speed_fail'
    }
    return 'pass'
}
