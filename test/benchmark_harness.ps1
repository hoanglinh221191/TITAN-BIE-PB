param(
    [string]$SessionName = "",
    [string[]]$Cases = @(),
    [switch]$SkipInvalidCases,
    [ValidateSet("all", "small", "medium", "large")]
    [string]$Suite = "all",
    [ValidateSet("unknown", "native", "portable")]
    [string]$BuildProfile = "unknown",
    [string]$TitanExePath = ""
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Convert-ToTitanPath {
    param([Parameter(Mandatory = $true)][string]$Path)
    return ($Path -replace "\\", "/")
}

function Get-TitanRelativePath {
    param(
        [Parameter(Mandatory = $true)][string]$FromDir,
        [Parameter(Mandatory = $true)][string]$ToPath
    )

    $fromFull = (Resolve-Path -LiteralPath $FromDir).Path
    $toFull = (Resolve-Path -LiteralPath $ToPath).Path
    $fromUri = [System.Uri]::new(($fromFull.TrimEnd("\") + "\"))
    $toUri = [System.Uri]::new($toFull)
    $relative = $fromUri.MakeRelativeUri($toUri).ToString()
    return [System.Uri]::UnescapeDataString($relative)
}

function Format-Seconds {
    param([double]$Value)
    return ("{0:N3}" -f $Value)
}

function Format-MiB {
    param([Int64]$Value)
    return ("{0:N2}" -f ($Value / 1MB))
}

function New-SessionRoot {
    param(
        [Parameter(Mandatory = $true)][string]$TestRoot,
        [string]$PreferredName
    )

    if ([string]::IsNullOrWhiteSpace($PreferredName)) {
        $PreferredName = "benchmark_run_{0}" -f (Get-Date -Format "yyyyMMdd_HHmmss")
    }

    $root = Join-Path $TestRoot $PreferredName
    if (Test-Path $root) {
        throw "Benchmark session directory already exists: $root"
    }

    New-Item -ItemType Directory -Path $root | Out-Null
    return $root
}

function Copy-BenchmarkAssets {
    param(
        [Parameter(Mandatory = $true)][string]$SourceDir,
        [Parameter(Mandatory = $true)][string]$DestinationDir,
        [Parameter(Mandatory = $true)][string[]]$Assets
    )

    foreach ($asset in $Assets) {
        $sourcePath = Join-Path $SourceDir $asset
        $destinationPath = Join-Path $DestinationDir $asset

        if (-not (Test-Path $sourcePath)) {
            throw "Missing benchmark asset: $sourcePath"
        }

        if ((Get-Item -LiteralPath $sourcePath).PSIsContainer) {
            Copy-Item -LiteralPath $sourcePath -Destination $destinationPath -Recurse -Force
        } else {
            Copy-Item -LiteralPath $sourcePath -Destination $destinationPath -Force
        }
    }
}

function Write-BenchmarkTextFiles {
    param(
        [Parameter(Mandatory = $true)][string]$DestinationDir,
        [Parameter(Mandatory = $true)][object[]]$Files,
        [Parameter(Mandatory = $true)][hashtable]$PlaceholderMap
    )

    foreach ($file in $Files) {
        $relativePath = [string]$file.Name
        $targetPath = Join-Path $DestinationDir $relativePath
        $targetDir = Split-Path -Parent $targetPath
        if (-not [string]::IsNullOrWhiteSpace($targetDir) -and
            -not (Test-Path -LiteralPath $targetDir)) {
            New-Item -ItemType Directory -Path $targetDir -Force | Out-Null
        }

        $content = [string]$file.Content
        foreach ($key in $PlaceholderMap.Keys) {
            $content = $content.Replace($key, $PlaceholderMap[$key])
        }
        Set-Content -LiteralPath $targetPath -Value $content -NoNewline -Encoding ASCII
    }
}

function Test-CaseExpectations {
    param(
        [Parameter(Mandatory = $true)][hashtable]$Case,
        [Parameter(Mandatory = $true)][string]$RunDir
    )

    $notes = @()
    $passed = $true

    if ($Case.ContainsKey("RequiredFiles")) {
        foreach ($name in @($Case.RequiredFiles)) {
            $candidate = Join-Path $RunDir ([string]$name)
            if (-not (Test-Path -LiteralPath $candidate)) {
                $passed = $false
                $notes += "Missing expected output: $name"
            }
        }
    }

    if ($Case.ContainsKey("RequiredPaths")) {
        foreach ($path in @($Case.RequiredPaths)) {
            $candidate = [string]$path
            if (-not (Test-Path -LiteralPath $candidate)) {
                $passed = $false
                $notes += "Missing expected output path: $candidate"
            }
        }
    }

    if ($Case.ContainsKey("ForbiddenFiles")) {
        foreach ($name in @($Case.ForbiddenFiles)) {
            $candidate = Join-Path $RunDir ([string]$name)
            if (Test-Path -LiteralPath $candidate) {
                $passed = $false
                $notes += "Unexpected output present: $name"
            }
        }
    }

    if ($Case.ContainsKey("TextChecks")) {
        foreach ($check in @($Case.TextChecks)) {
            $relativePath = [string]$check.Path
            $expectedText = [string]$check.Contains
            $candidate = Join-Path $RunDir $relativePath
            if (-not (Test-Path -LiteralPath $candidate)) {
                $passed = $false
                $notes += "Missing text-check file: $relativePath"
                continue
            }

            $text = Get-Content -LiteralPath $candidate -Raw
            if (-not $text.Contains($expectedText)) {
                $passed = $false
                $notes += "Missing expected text in ${relativePath}: $expectedText"
            }
        }
    }

    if ($Case.ContainsKey("RegexChecks")) {
        foreach ($check in @($Case.RegexChecks)) {
            $relativePath = [string]$check.Path
            $pattern = [string]$check.Pattern
            $candidate = Join-Path $RunDir $relativePath
            if (-not (Test-Path -LiteralPath $candidate)) {
                $passed = $false
                $notes += "Missing regex-check file: $relativePath"
                continue
            }

            $text = Get-Content -LiteralPath $candidate -Raw
            if (-not [regex]::IsMatch($text, $pattern, [System.Text.RegularExpressions.RegexOptions]::Multiline)) {
                $passed = $false
                $notes += "Missing expected regex in ${relativePath}: $pattern"
            }
        }
    }

    if ($Case.ContainsKey("PathTextChecks")) {
        foreach ($check in @($Case.PathTextChecks)) {
            $candidate = [string]$check.Path
            $expectedText = [string]$check.Contains
            if (-not (Test-Path -LiteralPath $candidate)) {
                $passed = $false
                $notes += "Missing text-check file: $candidate"
                continue
            }

            $text = Get-Content -LiteralPath $candidate -Raw
            if (-not $text.Contains($expectedText)) {
                $passed = $false
                $notes += "Missing expected text in ${candidate}: $expectedText"
            }
        }
    }

    return [PSCustomObject]@{
        Passed = $passed
        Notes = $notes
    }
}

function Get-LogTimings {
    param([Parameter(Mandatory = $true)][string]$LogPath)

    $timings = @{}
    if (-not (Test-Path $LogPath)) {
        return $timings
    }

    $pattern = '^\s*-\s+(?<name>[^:]+?)\s*:\s*(?<value>[0-9.]+\s+s)\s*$'
    foreach ($line in Get-Content $LogPath) {
        $match = [regex]::Match($line, $pattern)
        if ($match.Success) {
            $timings[$match.Groups["name"].Value.Trim()] = $match.Groups["value"].Value.Trim()
        }
    }

    return $timings
}

function Get-NanInfMatches {
    param(
        [Parameter(Mandatory = $true)][string]$RunDir,
        [string[]]$AllowedPatterns = @()
    )

    $textExtensions = @(".txt", ".dat", ".csv", ".xvg", ".pdb", ".dx")
    $skipNames = @("meta.json", "summary.md", "log.txt", "run.log", "run.stdout.log", "run.stderr.log")
    $allowedInfPatterns = @(
        'eps_rf=inf',
        '^MBAR Residual Inf\s*:'
    )
    $allowedInfPatterns += @($AllowedPatterns)
    $results = @()

    $files = Get-ChildItem -Path $RunDir -File -Recurse | Where-Object {
        $textExtensions -contains $_.Extension.ToLowerInvariant() -and
        -not ($skipNames -contains $_.Name)
    }

    foreach ($file in $files) {
        $hits = Select-String -Path $file.FullName -Pattern '(?i)(?<![A-Za-z])(?:nan|inf)(?![A-Za-z])' -ErrorAction SilentlyContinue
        foreach ($hit in $hits) {
            $lineText = $hit.Line.Trim()
            $allowed = $false
            foreach ($pattern in $allowedInfPatterns) {
                if ($lineText -match $pattern) {
                    $allowed = $true
                    break
                }
            }
            if (-not $allowed -and $file.Name -eq "umap_equilibrium_summary.dat" -and $lineText -match '(?i)\s+nan$') {
                $allowed = $true
            }
            if ($allowed) {
                continue
            }
            $results += [PSCustomObject]@{
                File = $file.FullName
                Line = $hit.LineNumber
                Text = $lineText
            }
        }
    }

    return $results
}

function Invoke-TitanCase {
    param(
        [Parameter(Mandatory = $true)][hashtable]$Case,
        [Parameter(Mandatory = $true)][string]$SessionRoot,
        [Parameter(Mandatory = $true)][string]$SessionTag,
        [Parameter(Mandatory = $true)][string]$TitanExe,
        [Parameter(Mandatory = $true)][string]$RepoRoot
    )

    $runDir = Join-Path $SessionRoot ("{0}_{1}" -f $SessionTag, $Case.Name)
    New-Item -ItemType Directory -Path $runDir | Out-Null
    $repoRootLocal = $RepoRoot
    $placeholderMap = @{
        "__REL_EXAMPLE2__" = Convert-ToTitanPath (Get-TitanRelativePath -FromDir $runDir -ToPath (Join-Path $repoRootLocal "example\\2"))
        "__REL_EXAMPLE3__" = Convert-ToTitanPath (Get-TitanRelativePath -FromDir $runDir -ToPath (Join-Path $repoRootLocal "example\\3"))
        "__REL_EXAMPLE5__" = Convert-ToTitanPath (Get-TitanRelativePath -FromDir $runDir -ToPath (Join-Path $repoRootLocal "example\\5"))
        "__REL_EXAMPLE6__" = Convert-ToTitanPath (Get-TitanRelativePath -FromDir $runDir -ToPath (Join-Path $repoRootLocal "example\\6_xtc"))
        "__REL_LIB_FF__" = Convert-ToTitanPath (Get-TitanRelativePath -FromDir $runDir -ToPath (Join-Path $repoRootLocal "lib_ff"))
    }

    if ($Case.ContainsKey("ArgContent")) {
        $argPath = Join-Path $runDir $Case.ArgFileName
        $argContent = $Case.ArgContent
        foreach ($key in $placeholderMap.Keys) {
            $argContent = $argContent.Replace($key, $placeholderMap[$key])
        }
        Set-Content -Path $argPath -Value $argContent -NoNewline -Encoding ASCII
        $argPathForProcess = $Case.ArgFileName
    } else {
        $argPath = $Case.ArgPath
        $argPathForProcess = Convert-ToTitanPath (Get-TitanRelativePath -FromDir $runDir -ToPath $argPath)
    }

    if ($Case.ContainsKey("Assets") -and (@($Case.Assets).Count -gt 0)) {
        $assetSourceDir = if ($Case.ContainsKey("AssetSourceDir")) { $Case.AssetSourceDir } else { Split-Path -Parent $argPath }
        Copy-BenchmarkAssets -SourceDir $assetSourceDir -DestinationDir $runDir -Assets $Case.Assets
    }
    if ($Case.ContainsKey("ExtraFiles") -and (@($Case.ExtraFiles).Count -gt 0)) {
        Write-BenchmarkTextFiles -DestinationDir $runDir -Files @($Case.ExtraFiles) -PlaceholderMap $placeholderMap
    }
    if ($Case.ContainsKey("AssetCopies") -and (@($Case.AssetCopies).Count -gt 0)) {
        foreach ($assetCopy in $Case.AssetCopies) {
            Copy-BenchmarkAssets -SourceDir $assetCopy.SourceDir -DestinationDir $runDir -Assets $assetCopy.Assets
        }
    }
    if ($Case.ContainsKey("CleanupPaths") -and (@($Case.CleanupPaths).Count -gt 0)) {
        foreach ($cleanupPath in @($Case.CleanupPaths)) {
            $targetPath = [string]$cleanupPath
            if (Test-Path -LiteralPath $targetPath) {
                Remove-Item -LiteralPath $targetPath -Force
            }
        }
    }

    $stdoutPath = Join-Path $runDir "run.stdout.log"
    $stderrPath = Join-Path $runDir "run.stderr.log"
    $runLogPath = Join-Path $runDir "run.log"

    $psi = New-Object System.Diagnostics.ProcessStartInfo
    $psi.FileName = $TitanExe
    $psi.Arguments = ('-a "{0}"' -f $argPathForProcess)
    $psi.WorkingDirectory = $runDir
    $psi.UseShellExecute = $false
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    $psi.CreateNoWindow = $true

    $process = New-Object System.Diagnostics.Process
    $process.StartInfo = $psi

    $stopwatch = [System.Diagnostics.Stopwatch]::StartNew()
    [void]$process.Start()
    $stdoutStream = [System.IO.File]::Open($stdoutPath, [System.IO.FileMode]::Create, [System.IO.FileAccess]::Write, [System.IO.FileShare]::Read)
    $stderrStream = [System.IO.File]::Open($stderrPath, [System.IO.FileMode]::Create, [System.IO.FileAccess]::Write, [System.IO.FileShare]::Read)
    $stdoutTask = $process.StandardOutput.BaseStream.CopyToAsync($stdoutStream)
    $stderrTask = $process.StandardError.BaseStream.CopyToAsync($stderrStream)

    $peakWorkingSet = [Int64]0
    $peakPrivateBytes = [Int64]0

    while (-not $process.HasExited) {
        Start-Sleep -Milliseconds 200
        try { $process.Refresh() } catch { }

        try {
            $snapshot = Get-Process -Id $process.Id -ErrorAction Stop
            $peakWorkingSet = [Math]::Max($peakWorkingSet, [Int64]$snapshot.WorkingSet64)
            $peakWorkingSet = [Math]::Max($peakWorkingSet, [Int64]$snapshot.PeakWorkingSet64)
            $peakPrivateBytes = [Math]::Max($peakPrivateBytes, [Int64]$snapshot.PrivateMemorySize64)
        } catch {
        }
    }

    $process.WaitForExit()
    $stdoutTask.Wait()
    $stderrTask.Wait()
    $stopwatch.Stop()
    try { $process.WaitForExit(1000) | Out-Null } catch { }
    try {
        $process.Refresh()
        $peakWorkingSet = [Math]::Max($peakWorkingSet, [Int64]$process.PeakWorkingSet64)
        $peakPrivateBytes = [Math]::Max($peakPrivateBytes, [Int64]$process.PrivateMemorySize64)
    } catch {
    }

    if ($stdoutStream) {
        $stdoutStream.Flush()
        $stdoutStream.Dispose()
    }
    if ($stderrStream) {
        $stderrStream.Flush()
        $stderrStream.Dispose()
    }

    $stdoutText = if (Test-Path $stdoutPath) { [string](Get-Content $stdoutPath -Raw) } else { "" }
    $stderrText = if (Test-Path $stderrPath) { [string](Get-Content $stderrPath -Raw) } else { "" }
    [System.IO.File]::WriteAllText($stdoutPath, $stdoutText)
    [System.IO.File]::WriteAllText($stderrPath, $stderrText)

    $combinedLog = [string]$stdoutText
    if (-not [string]::IsNullOrEmpty($stderrText)) {
        if (-not [string]::IsNullOrEmpty($combinedLog)) {
            $combinedLog += [Environment]::NewLine
        }
        $combinedLog += $stderrText
    }
    [System.IO.File]::WriteAllText($runLogPath, $combinedLog)

    $logPath = Join-Path $runDir "log.txt"
    $timings = Get-LogTimings -LogPath $logPath
    $nanInfMatches = @()
    if ($Case.ExpectedOutcome -eq "success") {
        $allowedNanInfPatterns = if ($Case.ContainsKey("AllowedNanInfPatterns")) { @($Case.AllowedNanInfPatterns) } else { @() }
        $nanInfMatches = @(Get-NanInfMatches -RunDir $runDir -AllowedPatterns $allowedNanInfPatterns)
    }

    $combinedOutput = ""
    if (Test-Path $stdoutPath) {
        $combinedOutput += (Get-Content $stdoutPath -Raw)
    }
    if (Test-Path $stderrPath) {
    if (-not [string]::IsNullOrEmpty($combinedOutput)) {
        $combinedOutput += "`n"
    }
        $combinedOutput += (Get-Content $stderrPath -Raw)
    }

    $exitCode = $process.ExitCode
    $passed = $false
    $notes = @()

    if ($Case.ExpectedOutcome -eq "success") {
        $passed = ($exitCode -eq 0 -and $nanInfMatches.Count -eq 0)
        if ($exitCode -ne 0) {
            $notes += "Unexpected exit code $exitCode"
        }
        if ($nanInfMatches.Count -gt 0) {
            $notes += "Detected nan/inf in text outputs"
        }
        $expectationResult = Test-CaseExpectations -Case $Case -RunDir $runDir
        if (-not $expectationResult.Passed) {
            $passed = $false
            $notes += @($expectationResult.Notes)
        }
    } else {
        $expectedExit = $Case.ExpectedExitCode
        $expectedMessage = $Case.ExpectedMessage
        $passed = ($exitCode -eq $expectedExit -and $combinedOutput.Contains($expectedMessage))
        if ($exitCode -ne $expectedExit) {
            $notes += "Expected exit $expectedExit but got $exitCode"
        }
        if (-not $combinedOutput.Contains($expectedMessage)) {
            $notes += "Expected message not found: $expectedMessage"
        }
    }

    $metaPath = Join-Path $runDir "meta.json"

    $result = [ordered]@{
        name = $Case.Name
        category = $Case.Category
        scale = $Case.Scale
        expected_outcome = $Case.ExpectedOutcome
        exit_code = $exitCode
        passed = $passed
        wall_seconds = [Math]::Round($stopwatch.Elapsed.TotalSeconds, 3)
        peak_working_set_bytes = $peakWorkingSet
        peak_working_set_mib = [Math]::Round(($peakWorkingSet / 1MB), 3)
        peak_private_bytes = $peakPrivateBytes
        peak_private_mib = [Math]::Round(($peakPrivateBytes / 1MB), 3)
        run_directory = $runDir
        arg_path = $argPath
        log_path = $logPath
        stdout_path = $stdoutPath
        stderr_path = $stderrPath
        meta_path = $metaPath
        timings = $timings
        nan_inf_match_count = $nanInfMatches.Count
        nan_inf_matches = $nanInfMatches
        notes = $notes
    }

    $timingEntries = @()
    foreach ($timingEntry in $timings.GetEnumerator()) {
        $timingEntries += [PSCustomObject]@{
            name = [string]$timingEntry.Key
            value = [string]$timingEntry.Value
        }
    }
    $metaResult = [PSCustomObject]@{
        name = [string]$Case.Name
        category = [string]$Case.Category
        scale = [string]$Case.Scale
        expected_outcome = [string]$Case.ExpectedOutcome
        exit_code = [int]$exitCode
        passed = [bool]$passed
        wall_seconds = [Math]::Round($stopwatch.Elapsed.TotalSeconds, 3)
        peak_working_set_bytes = [Int64]$peakWorkingSet
        peak_private_bytes = [Int64]$peakPrivateBytes
        run_directory = [string]$runDir
        arg_path = [string]$argPath
        log_path = [string]$logPath
        stdout_path = [string]$stdoutPath
        stderr_path = [string]$stderrPath
        timings = $timingEntries
        nan_inf_match_count = [int]$nanInfMatches.Count
        notes = @($notes)
    }
    $metaResult | ConvertTo-Json -Depth 6 | Set-Content -Path $metaPath

    $result.meta_path = $metaPath
    return [PSCustomObject]$result
}

function Update-ResultMeta {
    param([Parameter(Mandatory = $true)][pscustomobject]$Result)

    if (-not $Result.PSObject.Properties.Name.Contains("meta_path")) {
        return
    }
    $metaPath = [string]$Result.meta_path
    if ([string]::IsNullOrWhiteSpace($metaPath)) {
        return
    }

    $timingEntries = @()
    foreach ($timingEntry in $Result.timings.GetEnumerator()) {
        $timingEntries += [PSCustomObject]@{
            name = [string]$timingEntry.Key
            value = [string]$timingEntry.Value
        }
    }

    $metaResult = [PSCustomObject]@{
        name = [string]$Result.name
        category = [string]$Result.category
        scale = [string]$Result.scale
        expected_outcome = [string]$Result.expected_outcome
        exit_code = [int]$Result.exit_code
        passed = [bool]$Result.passed
        wall_seconds = [double]$Result.wall_seconds
        peak_working_set_bytes = [Int64]$Result.peak_working_set_bytes
        peak_private_bytes = [Int64]$Result.peak_private_bytes
        run_directory = [string]$Result.run_directory
        arg_path = [string]$Result.arg_path
        log_path = [string]$Result.log_path
        stdout_path = [string]$Result.stdout_path
        stderr_path = [string]$Result.stderr_path
        timings = $timingEntries
        nan_inf_match_count = [int]$Result.nan_inf_match_count
        notes = @($Result.notes)
    }
    $metaResult | ConvertTo-Json -Depth 6 | Set-Content -Path $metaPath
}

function Add-ResultFailureNote {
    param(
        [Parameter(Mandatory = $true)][pscustomobject]$Result,
        [Parameter(Mandatory = $true)][string]$Note
    )

    $notes = @($Result.notes)
    if (-not ($notes -contains $Note)) {
        $Result.notes = @($notes + $Note)
    }
    $Result.passed = $false
    Update-ResultMeta -Result $Result
}

function Invoke-CrossCaseHashValidation {
    param(
        [Parameter(Mandatory = $true)][pscustomobject[]]$Results,
        [Parameter(Mandatory = $true)][object[]]$Rules
    )

    $comparisons = @()
    foreach ($rule in @($Rules)) {
        $label = [string]$rule.Label
        $caseAName = [string]$rule.CaseA
        $caseBName = [string]$rule.CaseB
        $fileA = [string]$rule.FileA
        $fileB = if ($rule.ContainsKey("FileB")) { [string]$rule.FileB } else { $fileA }
        $caseA = $Results | Where-Object { $_.name -eq $caseAName } | Select-Object -First 1
        $caseB = $Results | Where-Object { $_.name -eq $caseBName } | Select-Object -First 1

        $entry = [ordered]@{
            label = $label
            case_a = $caseAName
            file_a = $fileA
            case_b = $caseBName
            file_b = $fileB
            identical = $false
            skipped = $false
            note = ""
        }

        if (-not $caseA -or -not $caseB) {
            $entry.skipped = $true
            $entry.note = "One or both cases were not part of this run."
            $comparisons += [PSCustomObject]$entry
            continue
        }
        if (-not $caseA.passed -or -not $caseB.passed) {
            $entry.skipped = $true
            $entry.note = "One or both prerequisite cases did not pass."
            $comparisons += [PSCustomObject]$entry
            continue
        }

        $pathA = Join-Path $caseA.run_directory $fileA
        $pathB = Join-Path $caseB.run_directory $fileB
        if (-not (Test-Path -LiteralPath $pathA) -or -not (Test-Path -LiteralPath $pathB)) {
            $entry.skipped = $true
            $entry.note = "One or both comparison files are missing."
            $comparisons += [PSCustomObject]$entry
            continue
        }

        $hashA = (Get-FileHash -Algorithm SHA256 -LiteralPath $pathA).Hash
        $hashB = (Get-FileHash -Algorithm SHA256 -LiteralPath $pathB).Hash
        $entry.identical = ($hashA -eq $hashB)
        if (-not $entry.identical) {
            $note = "Cross-case hash mismatch for ${label}: $caseAName/$fileA vs $caseBName/$fileB"
            Add-ResultFailureNote -Result $caseA -Note $note
            Add-ResultFailureNote -Result $caseB -Note $note
            $entry.note = $note
        }
        $comparisons += [PSCustomObject]$entry
    }

    return $comparisons
}

function Get-GistComparison {
    param(
        [Parameter(Mandatory = $true)][pscustomobject]$Nt1Case,
        [Parameter(Mandatory = $true)][pscustomobject]$Nt4Case
    )

    $comparison = @()
    $nt1Files = @(Get-ChildItem -Path $Nt1Case.run_directory -Filter "gist_*.dx" -File | Sort-Object Name)
    $nt4Files = @(Get-ChildItem -Path $Nt4Case.run_directory -Filter "gist_*.dx" -File | Sort-Object Name)

    if ($nt1Files.Count -ne $nt4Files.Count) {
        return [PSCustomObject]@{
            identical = $false
            files = @()
            note = "Different number of gist_*.dx files"
        }
    }

    for ($i = 0; $i -lt $nt1Files.Count; $i++) {
        $file1 = $nt1Files[$i]
        $file2 = $nt4Files[$i]
        $hash1 = (Get-FileHash -Algorithm SHA256 -LiteralPath $file1.FullName).Hash
        $hash2 = (Get-FileHash -Algorithm SHA256 -LiteralPath $file2.FullName).Hash
        $comparison += [PSCustomObject]@{
            name = $file1.Name
            identical = ($hash1 -eq $hash2)
            nt1_hash = $hash1
            nt4_hash = $hash2
        }
    }

    return [PSCustomObject]@{
        identical = (@($comparison).Count -gt 0 -and (@($comparison | Where-Object { -not $_.identical }).Count -eq 0))
        files = $comparison
        note = ""
    }
}

function Write-Summary {
    param(
        [Parameter(Mandatory = $true)][string]$SummaryPath,
        [Parameter(Mandatory = $true)][string]$SessionRoot,
        [Parameter(Mandatory = $true)][pscustomobject[]]$Results,
        [Parameter()][pscustomobject]$GistComparison,
        [Parameter()][pscustomobject[]]$CrossCaseComparisons,
        [Parameter(Mandatory = $true)][string]$Suite,
        [Parameter(Mandatory = $true)][string]$BuildProfile
    )

    $generatedAt = Get-Date -Format "yyyy-MM-dd HH:mm:ss zzz"
    $suiteLabel = if ($Suite -eq "all") { "all-default (small+medium)" } else { $Suite }
    $lines = @()
    $lines += "# Benchmark Harness Summary"
    $lines += ""
    $lines += "- Generated: $generatedAt"
    $lines += "- Session root: ``$SessionRoot``"
    $lines += "- Suite: ``$suiteLabel``"
    $lines += "- Build profile: ``$BuildProfile``"
    $lines += ""
    $lines += "## Results"
    $lines += ""
    $lines += "| Case | Scale | Outcome | Exit | Wall (s) | Peak WS (MiB) | Peak Private (MiB) | TOTAL TRAJECTORY | GLOBAL TOTAL | NaN/Inf |"
    $lines += "| --- | --- | --- | --- | ---: | ---: | ---: | --- | --- | ---: |"

    foreach ($result in $Results) {
        $traj = if ($result.timings.ContainsKey("TOTAL TRAJECTORY TIME")) { $result.timings["TOTAL TRAJECTORY TIME"] } else { "-" }
        $global = if ($result.timings.ContainsKey("GLOBAL TOTAL RUN TIME")) { $result.timings["GLOBAL TOTAL RUN TIME"] } else { "-" }
        $outcome = if ($result.passed) { "pass" } else { "fail" }
        $lines += "| $($result.name) | $($result.scale) | $outcome | $($result.exit_code) | $(Format-Seconds $result.wall_seconds) | $(Format-MiB $result.peak_working_set_bytes) | $(Format-MiB $result.peak_private_bytes) | $traj | $global | $($result.nan_inf_match_count) |"
    }

    $lines += ""
    $lines += "## Notes"
    $lines += ""

    foreach ($result in $Results) {
        if ($result.notes.Count -eq 0) {
            continue
        }
        $lines += "- $($result.name): $($result.notes -join '; ')"
    }

    if ($Results | Where-Object { $_.expected_outcome -eq "expected_failure" }) {
        foreach ($result in $Results | Where-Object { $_.expected_outcome -eq "expected_failure" }) {
            $stderr = if (Test-Path $result.stderr_path) { (Get-Content $result.stderr_path -Raw).Trim() } else { "" }
            $stdout = if (Test-Path $result.stdout_path) { (Get-Content $result.stdout_path -Raw).Trim() } else { "" }
            $message = if (-not [string]::IsNullOrWhiteSpace($stderr)) { $stderr } else { $stdout }
            if (-not [string]::IsNullOrWhiteSpace($message)) {
                $lines += "- $($result.name) message: $message"
            }
        }
    }

    if ($GistComparison) {
        $lines += ""
        $lines += "## GIST DX Comparison"
        $lines += ""
        $lines += "- Hash-identical across nt=1 and nt=4: $($GistComparison.identical)"
        if (@($GistComparison.files).Count -gt 0) {
            $lines += ""
            $lines += "| File | Identical |"
            $lines += "| --- | --- |"
            foreach ($entry in $GistComparison.files) {
                $lines += "| $($entry.name) | $($entry.identical) |"
            }
        }
        if (-not [string]::IsNullOrWhiteSpace($GistComparison.note)) {
            $lines += "- $($GistComparison.note)"
        }
    }

    if ($CrossCaseComparisons -and @($CrossCaseComparisons).Count -gt 0) {
        $lines += ""
        $lines += "## Cross-Case Hash Checks"
        $lines += ""
        $lines += "| Rule | Identical | Status |"
        $lines += "| --- | --- | --- |"
        foreach ($entry in $CrossCaseComparisons) {
            $status = if ($entry.skipped) { $entry.note } elseif ([string]::IsNullOrWhiteSpace($entry.note)) { "ok" } else { $entry.note }
            $lines += "| $($entry.label) | $($entry.identical) | $status |"
        }
    }

    $lines += ""
    $lines += "## Follow-up"
    $lines += ""
    $lines += "- Small suite phu hop cho smoke/regression nhanh sau khi sua hot path."
    $lines += "- Medium suite phu hop cho data-layout refactor va postproc regression."
    $lines += "- Large suite phu hop cho memory-sensitive work va so sanh native/portable."
    $lines += "- Sau Phase 0, uu tien code nen la DCCM/allostery pruning va MSM spectral cost."

    Set-Content -Path $SummaryPath -Value ($lines -join "`r`n")
}

$scriptRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Split-Path -Parent $scriptRoot
$testRoot = $scriptRoot

if ([string]::IsNullOrWhiteSpace($TitanExePath)) {
    $titanExeCandidates = @()
    $titanExeCandidates += Join-Path $repoRoot "titan.exe"
    if ($BuildProfile -ne "unknown") {
        $titanExeCandidates += Join-Path $testRoot ("bin\\{0}\\titan.exe" -f $BuildProfile)
    }
    $titanExe = $null
    foreach ($candidate in $titanExeCandidates) {
        if (Test-Path $candidate) {
            $titanExe = (Resolve-Path -LiteralPath $candidate).Path
            break
        }
    }
    if ($null -eq $titanExe) {
        throw "Could not find titan.exe for build profile '$BuildProfile'. Tried: $($titanExeCandidates -join ', ')"
    }
} else {
    $titanExe = (Resolve-Path -LiteralPath $TitanExePath).Path
}

$example2Arg = Join-Path $repoRoot "example\\2\\arg.in"
$example2Dir = Join-Path $repoRoot "example\\2"
$example3Dir = Join-Path $repoRoot "example\\3"
$example5Arg = Join-Path $repoRoot "example\\5\\arg.in"
$example5SqArg = Join-Path $repoRoot "example\\5\\arg_sq_channels.in"
$sqLocalArg = Join-Path $repoRoot "test\\sq_inputs_local\\arg_sq_local.in"
$corrKnnInputDir = Join-Path $repoRoot "test\\corr_knn_inputs_local"
$hbkInputDir = Join-Path $repoRoot "test\\hbond_kinetics_inputs_local"
$graphInputDir = Join-Path $repoRoot "test\\graph_inputs_local"
$corrKnnLocalArg = Join-Path $corrKnnInputDir "arg_corr_knn_local.in"
$corrKnnNoTeArg = Join-Path $corrKnnInputDir "arg_corr_knn_no_te.in"
$allosteryHistLocalArg = Join-Path $corrKnnInputDir "arg_allostery_hist_local.in"
$allosteryKdeLocalArg = Join-Path $corrKnnInputDir "arg_allostery_kde_local.in"
$allosteryKnnLocalArg = Join-Path $corrKnnInputDir "arg_allostery_knn_local.in"
$corrKnnAllosteryReuseArg = Join-Path $corrKnnInputDir "arg_corr_knn_allostery_reuse.in"
$hbkLocalArg = Join-Path $hbkInputDir "arg_hbk_local.in"
$hbkWaterMissingArg = Join-Path $hbkInputDir "arg_hbk_water_missing.in"
$graphDccmLocalArg = Join-Path $graphInputDir "arg_graph_dccm_local.in"
$graphLspmdLocalArg = Join-Path $graphInputDir "arg_graph_lspmd_local.in"
$graphAllosteryGcLocalArg = Join-Path $graphInputDir "arg_graph_allostery_gc_local.in"
$graphAllosteryTeLocalArg = Join-Path $graphInputDir "arg_graph_allostery_te_local.in"
$graphTptLocalArg = Join-Path $graphInputDir "arg_graph_tpt_local.in"
$graphAutoLocalArg = Join-Path $graphInputDir "arg_graph_auto_local.in"
$corrKnnLocalArgText = (Get-Content -LiteralPath $corrKnnLocalArg -Raw).Replace("../../lib_ff", "__REL_LIB_FF__")
$corrKnnNoTeArgText = (Get-Content -LiteralPath $corrKnnNoTeArg -Raw).Replace("../../lib_ff", "__REL_LIB_FF__")
$allosteryHistLocalArgText = (Get-Content -LiteralPath $allosteryHistLocalArg -Raw).Replace("../../lib_ff", "__REL_LIB_FF__")
$allosteryKdeLocalArgText = (Get-Content -LiteralPath $allosteryKdeLocalArg -Raw).Replace("../../lib_ff", "__REL_LIB_FF__")
$allosteryKnnLocalArgText = (Get-Content -LiteralPath $allosteryKnnLocalArg -Raw).Replace("../../lib_ff", "__REL_LIB_FF__")
$corrKnnAllosteryReuseArgText = (Get-Content -LiteralPath $corrKnnAllosteryReuseArg -Raw).Replace("../../lib_ff", "__REL_LIB_FF__")
$hbkLocalArgText = (Get-Content -LiteralPath $hbkLocalArg -Raw).Replace("../../lib_ff", "__REL_LIB_FF__")
$hbkWaterMissingArgText = (Get-Content -LiteralPath $hbkWaterMissingArg -Raw).Replace("../../lib_ff", "__REL_LIB_FF__")
$graphDccmLocalArgText = (Get-Content -LiteralPath $graphDccmLocalArg -Raw).Replace("../../lib_ff", "__REL_LIB_FF__")
$graphLspmdLocalArgText = (Get-Content -LiteralPath $graphLspmdLocalArg -Raw).Replace("../../lib_ff", "__REL_LIB_FF__")
$graphAllosteryGcLocalArgText = (Get-Content -LiteralPath $graphAllosteryGcLocalArg -Raw).Replace("../../lib_ff", "__REL_LIB_FF__")
$graphAllosteryTeLocalArgText = (Get-Content -LiteralPath $graphAllosteryTeLocalArg -Raw).Replace("../../lib_ff", "__REL_LIB_FF__")
$graphTptLocalArgText = (Get-Content -LiteralPath $graphTptLocalArg -Raw).Replace("../../lib_ff", "__REL_LIB_FF__")
$graphAutoLocalArgText = (Get-Content -LiteralPath $graphAutoLocalArg -Raw).Replace("../../lib_ff", "__REL_LIB_FF__")
$example6Dir = Join-Path $repoRoot "example\\6_xtc"
$xdrMiniArg = Join-Path $repoRoot "example\\6_xtc\\arg_mini.in"
$xdrTicaUmapArg = Join-Path $repoRoot "example\\6_xtc\\arg_tica_umap.in"
$example8Arg = Join-Path $repoRoot "example\\8cbf_ba1\\arg.in"
$gistNt1Arg = Join-Path $repoRoot "example\\6_xtc\\arg_xdr_water_gist_nt1.in"
$gistNt4Arg = Join-Path $repoRoot "example\\6_xtc\\arg_xdr_water_gist.in"
$gistDir = Join-Path $repoRoot "example\\6_xtc"
$mbarInputDir = Join-Path $repoRoot "test\\mbar_inputs_local"
$example5SqOutputDir = Join-Path $repoRoot "example\\5\\md1"
$example5SqRequiredPaths = @(
    (Join-Path $example5SqOutputDir "sq_summary.dat"),
    (Join-Path $example5SqOutputDir "sq_rdf_all.dat"),
    (Join-Path $example5SqOutputDir "sq_rdf_backbone.dat"),
    (Join-Path $example5SqOutputDir "sq_rdf_sticker.dat"),
    (Join-Path $example5SqOutputDir "sq_rdf_backbone__sticker.dat"),
    (Join-Path $example5SqOutputDir "sq_structure_factor_all.dat"),
    (Join-Path $example5SqOutputDir "sq_structure_factor_backbone.dat"),
    (Join-Path $example5SqOutputDir "sq_structure_factor_sticker.dat"),
    (Join-Path $example5SqOutputDir "sq_structure_factor_backbone__sticker.dat")
)

$screeningValidArg = @"
-ff=AMBER14SB
-pbci=0
-g=4-3
-nt=4
-entropy
-eu=kCal/mol
-energy-res
-split
-1pdb
-hb
-emap-level=all
-b=100000
-e=110000
-n=__REL_EXAMPLE2__/index.ndx
-common_topo=__REL_EXAMPLE2__/topol.top
-noprint-side_chain_contact
-f=__REL_EXAMPLE2__/md1
-lspmd
-elec-screening
-saltcon=0.15
"@

$screeningInvalidSaltArg = @"
-ff=AMBER14SB
-pbci=0
-g=4-3
-nt=4
-entropy
-eu=kCal/mol
-energy-res
-split
-1pdb
-hb
-emap-level=all
-b=100000
-e=110000
-n=__REL_EXAMPLE2__/index.ndx
-common_topo=__REL_EXAMPLE2__/topol.top
-noprint-side_chain_contact
-f=__REL_EXAMPLE2__/md1
-lspmd
-elec-screening
-saltcon=0
"@

$screeningInvalidEpsilonArg = @"
-ff=AMBER14SB
-pbci=0
-g=4-3
-nt=4
-entropy
-eu=kCal/mol
-energy-res
-split
-1pdb
-hb
-emap-level=all
-b=100000
-e=110000
-n=__REL_EXAMPLE2__/index.ndx
-common_topo=__REL_EXAMPLE2__/topol.top
-noprint-side_chain_contact
-f=__REL_EXAMPLE2__/md1
-lspmd
-elec-screening
-saltcon=0.15
-epsilon_in=0
"@

$pcaArg = @"
-ff=amber99sb-ildn
-path-a=__REL_LIB_FF__
-path-n=__REL_LIB_FF__
-topo=__REL_EXAMPLE3__/OLA_GMX.top
-nt=4
-split
-eu=kcal/mol
-1pdb
-n=__REL_EXAMPLE3__/index.ndx
-g=1-13
-n_wat=17
-pbcn
-e=110000
-pca
-f=__REL_EXAMPLE3__/md1.pdb
"@

$sqExample5Arg = @"
-ff=martini
-g=1-12
-nt=4
-entropy
-eu=kCal/mol
-energy-res
-split
-1pdb
-emap-level=all
-e=1010000
-n=__REL_EXAMPLE5__/index.ndx
-common_topo=__REL_EXAMPLE5__/topol.top
-noprint-side_chain_contact
-f=__REL_EXAMPLE5__/md1.pdb
-llps
-llps-dimer-size=2
-llps-rcut=0.7
-llps-perc=0.5
-llps-contacts
-llps-sticker=BB,SC1,SC2,SC3,SC4,BB1,BB2,BB3
-sq
-sq-channel=all
-sq-channel=backbone
-sq-channel=sticker
-sq-rdf-bin=2.0
-sq-rdf-max=30.0
-sq-q-bin=0.1
-sq-q-max=5.0
"@

$mbarFixedArg = @"
-s=__REL_EXAMPLE6__/md1_snap0.pdb
-ff=charmm36
-g=1-20
-common_topo=__REL_EXAMPLE6__/topol.top
-i=
./xtc_list_twice.txt
-dt=400
-e=1200
-n=__REL_EXAMPLE6__/index1.ndx
-pbcn
-eu=kcal/mol
-tica
-tica-lag=1
-tica-components=2
-tica-rank=8
-msm
-msm-clusters=4
-msm-lag=1
-nt=2
-mbar
-mbar-temps=
./temps_valid.txt
-mbar-ener=
./energy_list_valid.txt
-mbar-target-temp=300
-mbar-domain=tica
-mbar-window-match=auto
-mbar-scan-points=16
"@

$mbarPerFrameArg = @"
-s=__REL_EXAMPLE6__/md1_snap0.pdb
-ff=charmm36
-g=1-20
-common_topo=__REL_EXAMPLE6__/topol.top
-i=
./xtc_list_twice.txt
-dt=400
-e=1200
-n=__REL_EXAMPLE6__/index1.ndx
-pbcn
-eu=kcal/mol
-tica
-tica-lag=1
-tica-components=2
-tica-rank=8
-msm
-msm-clusters=4
-msm-lag=1
-nt=2
-mbar
-mbar-temps=
./temps_tremd_valid.txt
-mbar-ener=
./energy_list_tremd_nvt.txt
-mbar-target-temp=300
-mbar-domain=tica
-mbar-ensemble=nvt
-mbar-temp-mode=per_frame
-mbar-window-match=auto
-mbar-scan-points=8
"@

$mbarDemuxArg = @"
-s=__REL_EXAMPLE6__/md1_snap0.pdb
-ff=charmm36
-g=1-20
-common_topo=__REL_EXAMPLE6__/topol.top
-i=
./xtc_list_twice.txt
-dt=400
-e=1200
-n=__REL_EXAMPLE6__/index1.ndx
-pbcn
-eu=kcal/mol
-tica
-tica-lag=1
-tica-components=2
-tica-rank=8
-msm
-msm-clusters=4
-msm-lag=1
-nt=2
-mbar
-mbar-temps=
./temps_tremd_valid.txt
-mbar-ener=
./energy_list_tremd_nvt.txt
-mbar-target-temp=300
-mbar-domain=tica
-mbar-ensemble=nvt
-mbar-temp-mode=per_frame
-mbar-demux-map=
./demux_map_tremd.txt
-mbar-window-match=auto
-mbar-scan-points=8
"@

$mbarNptArg = @"
-s=__REL_EXAMPLE6__/md1_snap0.pdb
-ff=charmm36
-g=1-20
-common_topo=__REL_EXAMPLE6__/topol.top
-i=
./xtc_list_twice.txt
-dt=400
-e=1200
-n=__REL_EXAMPLE6__/index1.ndx
-pbcn
-eu=kcal/mol
-pca
-nt=2
-mbar
-mbar-temps=
./temps_tremd_valid.txt
-mbar-ener=
./energy_list_tremd_npt.txt
-mbar-target-temp=300
-mbar-domain=pca
-mbar-ensemble=npt
-mbar-temp-mode=per_frame
-mbar-window-match=auto
-mbar-scan-points=8
"@

$mbarUmapArg = @"
-s=__REL_EXAMPLE6__/md1_snap0.pdb
-ff=charmm36
-g=1-20
-common_topo=__REL_EXAMPLE6__/topol.top
-i=
./xtc_list_twice.txt
-dt=400
-e=1200
-n=__REL_EXAMPLE6__/index1.ndx
-pbcn
-eu=kcal/mol
-tica
-tica-lag=1
-tica-components=2
-tica-rank=8
-msm
-msm-clusters=4
-msm-lag=1
-umap
-nt=2
-mbar
-mbar-temps=
./temps_valid.txt
-mbar-ener=
./energy_list_valid.txt
-mbar-target-temp=300
-mbar-domain=tica
-mbar-window-match=auto
-mbar-scan-points=16
-eq-weight-source=mbar
"@

$mbarDccmContextArg = @"
-s=__REL_EXAMPLE6__/md1_snap0.pdb
-ff=charmm36
-g=1-20
-common_topo=__REL_EXAMPLE6__/topol.top
-i=
./xtc_list_twice.txt
-dt=400
-e=1200
-n=__REL_EXAMPLE6__/index1.ndx
-pbcn
-eu=kcal/mol
-tica
-tica-lag=1
-tica-components=2
-tica-rank=8
-msm
-msm-clusters=4
-msm-lag=1
-allostery
-allostery-no-te
-allostery-tau=1
-lspmd
-nt=2
-mbar
-mbar-temps=
./temps_valid.txt
-mbar-ener=
./energy_list_valid.txt
-mbar-target-temp=300
-mbar-domain=tica
-mbar-window-match=auto
-mbar-scan-points=16
-eq-weight-source=mbar
"@

$mbarPcaSelectorArg = @"
-s=__REL_EXAMPLE6__/md1_snap0.pdb
-ff=charmm36
-g=1-20
-common_topo=__REL_EXAMPLE6__/topol.top
-i=
./xtc_list_twice.txt
-dt=400
-e=1200
-n=__REL_EXAMPLE6__/index1.ndx
-pbcn
-eu=kcal/mol
-pca
-fes-grid=40
-nt=2
-mbar
-mbar-temps=
./temps_valid.txt
-mbar-ener=
./energy_list_valid.txt
-mbar-target-temp=300
-mbar-domain=pca
-mbar-window-match=auto
-mbar-scan-points=16
-eq-weight-source=mbar
"@

$mbarTptOverlayArg = @"
-s=__REL_EXAMPLE6__/md1_snap0.pdb
-ff=charmm36
-g=1-20
-common_topo=__REL_EXAMPLE6__/topol.top
-i=
./xtc_list_twice.txt
-dt=400
-e=1200
-n=__REL_EXAMPLE6__/index1.ndx
-pbcn
-eu=kcal/mol
-tica
-tica-lag=1
-tica-components=2
-tica-rank=8
-msm
-msm-clusters=4
-msm-lag=1
-nt=2
-mbar
-mbar-temps=
./temps_valid.txt
-mbar-ener=
./energy_list_valid.txt
-mbar-target-temp=300
-mbar-domain=tica
-mbar-window-match=auto
-mbar-scan-points=16
-eq-weight-source=mbar
-tpt
-tpt-source=micro:0
-tpt-sink=micro:1
"@

$eqAutoGuardArg = @"
-s=__REL_EXAMPLE6__/md1_snap0.pdb
-ff=charmm36
-g=1-20
-common_topo=__REL_EXAMPLE6__/topol.top
-i=
./xtc_list_twice.txt
-dt=400
-e=1200
-n=__REL_EXAMPLE6__/index1.ndx
-pbcn
-eu=kcal/mol
-tica
-tica-lag=1
-tica-components=2
-tica-rank=8
-msm
-msm-clusters=4
-msm-lag=1
-umap
-tpt
-tpt-source=micro:0
-tpt-sink=micro:1
-allostery
-allostery-no-te
-allostery-tau=1
-lspmd
-nt=2
-mbar
-mbar-temps=
./temps_valid.txt
-mbar-ener=
./energy_list_valid.txt
-mbar-target-temp=300
-mbar-domain=tica
-mbar-window-match=auto
-mbar-scan-points=16
-eq-weight-source=auto
"@

$eqUniformGuardArg = @"
-s=__REL_EXAMPLE6__/md1_snap0.pdb
-ff=charmm36
-g=1-20
-common_topo=__REL_EXAMPLE6__/topol.top
-i=
./xtc_list_twice.txt
-dt=400
-e=1200
-n=__REL_EXAMPLE6__/index1.ndx
-pbcn
-eu=kcal/mol
-tica
-tica-lag=1
-tica-components=2
-tica-rank=8
-msm
-msm-clusters=4
-msm-lag=1
-umap
-tpt
-tpt-source=micro:0
-tpt-sink=micro:1
-allostery
-allostery-no-te
-allostery-tau=1
-lspmd
-nt=2
-mbar
-mbar-temps=
./temps_valid.txt
-mbar-ener=
./energy_list_valid.txt
-mbar-target-temp=300
-mbar-domain=tica
-mbar-window-match=auto
-mbar-scan-points=16
-eq-weight-source=uniform
"@

$eqDhamGuardArg = @"
-s=__REL_EXAMPLE6__/md1_snap0.pdb
-ff=charmm36
-g=1-20
-common_topo=__REL_EXAMPLE6__/topol.top
-i=
./xtc_list_twice.txt
-dt=400
-e=1200
-n=__REL_EXAMPLE6__/index1.ndx
-pbcn
-eu=kcal/mol
-tica
-tica-lag=1
-tica-components=2
-tica-rank=8
-msm
-msm-clusters=4
-msm-lag=1
-umap
-tpt
-tpt-source=micro:0
-tpt-sink=micro:1
-allostery
-allostery-no-te
-allostery-tau=1
-lspmd
-dham
-dham-temps=
./temps_valid.txt
-dham-ener=
./energy_list_valid.txt
-dham-target-temp=300
-dham-window-match=auto
-eq-weight-source=dham
-nt=2
"@

$mbarMismatchGuardArg = @"
-s=__REL_EXAMPLE6__/md1_snap0.pdb
-ff=charmm36
-g=1-20
-common_topo=__REL_EXAMPLE6__/topol.top
-i=
./xtc_list_twice.txt
-dt=400
-e=1200
-n=__REL_EXAMPLE6__/index1.ndx
-pbcn
-eu=kcal/mol
-tica
-tica-lag=1
-tica-components=2
-tica-rank=8
-msm
-msm-clusters=4
-msm-lag=1
-umap
-tpt
-tpt-source=micro:0
-tpt-sink=micro:1
-allostery
-allostery-no-te
-allostery-tau=1
-lspmd
-nt=2
-mbar
-mbar-temps=
./temps_valid.txt
-mbar-ener=
./energy_list_mismatch.txt
-mbar-target-temp=300
-mbar-domain=tica
-mbar-window-match=auto
-eq-weight-source=mbar
"@

$mbarXtcList = @"
__REL_EXAMPLE6__/md1_1_mini.xtc
__REL_EXAMPLE6__/md1_1_mini.xtc
"@

$mbarTempsValid = @"
290.0
310.0
"@

$mbarEnergyListFixed = @"
./energy_rep1.xvg
./energy_rep2.xvg
"@

$mbarEnergyListMismatch = @"
./energy_rep1.xvg
./energy_rep2_shifted.xvg
./energy_rep1.xvg
"@

$mbarEnergyListTremdNvt = @"
./energy_rep1_tremd_nvt.xvg
./energy_rep2_tremd_nvt.xvg
"@

$mbarEnergyListTremdNpt = @"
./energy_rep1_tremd_npt.xvg
./energy_rep2_tremd_npt.xvg
"@

$mbarDemuxMap = @"
0 1
0 1
1 0
1 0
"@

$caseTable = @(
    @{
        Name = "chunk_streamer_full"
        Category = "chunk_streamer"
        Scale = "medium"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgPath = $example2Arg
    },
    @{
        Name = "screening_valid"
        Category = "screening"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_valid.in"
        ArgContent = $screeningValidArg
        ForbiddenFiles = @(
            "dccm_equilibrium_context.dat",
            "allostery_equilibrium_ranking.csv",
            "lspmd_equilibrium_ranking.csv"
        )
    },
    @{
        Name = "screening_invalid_salt"
        Category = "screening"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "expected_failure"
        ArgFileName = "arg_invalid_salt.in"
        ArgContent = $screeningInvalidSaltArg
        ExpectedExitCode = 1
        ExpectedMessage = "salt concentration must be > 0"
    },
    @{
        Name = "screening_invalid_epsilon"
        Category = "screening"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "expected_failure"
        ArgFileName = "arg_invalid_epsilon.in"
        ArgContent = $screeningInvalidEpsilonArg
        ExpectedExitCode = 1
        ExpectedMessage = "-epsilon_in must be > 0"
    },
    @{
        Name = "xdr_gbsa_mini"
        Category = "xdr"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgPath = $xdrMiniArg
    },
    @{
        Name = "mbar_fixed_tica"
        Category = "mbar"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_mbar_fixed.in"
        ArgContent = $mbarFixedArg
        AssetSourceDir = $mbarInputDir
        Assets = @(
            "energy_rep1.xvg",
            "energy_rep2.xvg"
        )
        ExtraFiles = @(
            @{ Name = "xtc_list_twice.txt"; Content = $mbarXtcList },
            @{ Name = "temps_valid.txt"; Content = $mbarTempsValid },
            @{ Name = "energy_list_valid.txt"; Content = $mbarEnergyListFixed }
        )
        RequiredFiles = @(
            "mbar_summary.dat",
            "mbar_window_audit.dat",
            "mbar_frame_weights.dat",
            "mbar_cv.dat",
            "mbar_tica_fes_2d.dat",
            "mbar_tica_cluster_representatives.csv"
        )
        TextChecks = @(
            @{ Path = "mbar_summary.dat"; Contains = "MBAR Converged                 : yes" },
            @{ Path = "mbar_summary.dat"; Contains = "MBAR Temp Mode                 : fixed" }
        )
    },
    @{
        Name = "mbar_pca_selector"
        Category = "mbar"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_mbar_pca_selector.in"
        ArgContent = $mbarPcaSelectorArg
        AssetSourceDir = $mbarInputDir
        Assets = @(
            "energy_rep1.xvg",
            "energy_rep2.xvg"
        )
        ExtraFiles = @(
            @{ Name = "xtc_list_twice.txt"; Content = $mbarXtcList },
            @{ Name = "temps_valid.txt"; Content = $mbarTempsValid },
            @{ Name = "energy_list_valid.txt"; Content = $mbarEnergyListFixed }
        )
        RequiredFiles = @(
            "mbar_summary.dat",
            "mbar_pca_fes_2d.dat",
            "mbar_pca_cluster_representatives.csv",
            "pca_fes_2d.dat",
            "cluster_representatives.dat"
        )
        TextChecks = @(
            @{ Path = "log.txt"; Contains = "Equilibrium Weight Source      : mbar" },
            @{ Path = "log.txt"; Contains = "PCA Weight Mode                : MBAR-weighted" }
        )
    },
    @{
        Name = "mbar_umap_sidecar"
        Category = "mbar"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_mbar_umap.in"
        ArgContent = $mbarUmapArg
        AssetSourceDir = $mbarInputDir
        Assets = @(
            "energy_rep1.xvg",
            "energy_rep2.xvg"
        )
        ExtraFiles = @(
            @{ Name = "xtc_list_twice.txt"; Content = $mbarXtcList },
            @{ Name = "temps_valid.txt"; Content = $mbarTempsValid },
            @{ Name = "energy_list_valid.txt"; Content = $mbarEnergyListFixed }
        )
        RequiredFiles = @(
            "mbar_summary.dat",
            "msm_equilibrium_summary.dat",
            "tica_fes_2d.dat",
            "titan_umap_embedding.csv",
            "umap_fes_2d.dat",
            "umap_cluster_representatives.csv",
            "umap_equilibrium_summary.dat"
        )
        TextChecks = @(
            @{ Path = "log.txt"; Contains = "Equilibrium Weight Source      : mbar" },
            @{ Path = "log.txt"; Contains = "TICA Weight Mode               : MBAR-weighted" },
            @{ Path = "umap_equilibrium_summary.dat"; Contains = "# source mbar" },
            @{ Path = "umap_equilibrium_summary.dat"; Contains = "# embedding_dims 2" }
        )
    },
    @{
        Name = "mbar_per_frame_tica"
        Category = "mbar"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_mbar_per_frame.in"
        ArgContent = $mbarPerFrameArg
        AssetSourceDir = $mbarInputDir
        Assets = @(
            "energy_rep1_tremd_nvt.xvg",
            "energy_rep2_tremd_nvt.xvg"
        )
        ExtraFiles = @(
            @{ Name = "xtc_list_twice.txt"; Content = $mbarXtcList },
            @{ Name = "temps_tremd_valid.txt"; Content = $mbarTempsValid },
            @{ Name = "energy_list_tremd_nvt.txt"; Content = $mbarEnergyListTremdNvt }
        )
        RequiredFiles = @(
            "mbar_summary.dat",
            "mbar_window_audit.dat",
            "mbar_frame_weights.dat",
            "mbar_cv.dat",
            "mbar_tica_fes_2d.dat",
            "mbar_tica_cluster_representatives.csv"
        )
        TextChecks = @(
            @{ Path = "mbar_summary.dat"; Contains = "MBAR Converged                 : yes" },
            @{ Path = "mbar_summary.dat"; Contains = "MBAR Temp Mode                 : per_frame" }
        )
    },
    @{
        Name = "sq_local_fixture"
        Category = "sq"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgPath = $sqLocalArg
        RequiredFiles = @(
            "sq_summary.dat",
            "sq_rdf_all.dat",
            "sq_structure_factor_all.dat"
        )
        TextChecks = @(
            @{ Path = "sq_summary.dat"; Contains = "module: titan_sq" },
            @{ Path = "sq_summary.dat"; Contains = "first_valid_box_type: cubic" },
            @{ Path = "sq_structure_factor_all.dat"; Contains = "# method sq_from_rdf_lorch" },
            @{ Path = "log.txt"; Contains = "SQ Analysis" }
        )
    },
    @{
        Name = "hbk_local_fixture"
        Category = "hbk"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_hbk_local.in"
        ArgContent = $hbkLocalArgText
        AssetSourceDir = $corrKnnInputDir
        Assets = @("mini_corr_knn.pdb", "index.ndx")
        RequiredFiles = @(
            "hbond_kinetics_summary.dat",
            "hbond_kinetics_hbond_pairs.csv",
            "hbond_kinetics_salt_pairs.csv",
            "hbond_kinetics_hbond_ct.dat",
            "hbond_kinetics_salt_ct.dat"
        )
        TextChecks = @(
            @{ Path = "hbond_kinetics_summary.dat"; Contains = "module: titan_hbond_kinetics" },
            @{ Path = "log.txt"; Contains = "HBond Kinetics Analysis" }
        )
    },
    @{
        Name = "hbk_water_missing_fixture"
        Category = "hbk"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_hbk_water_missing.in"
        ArgContent = $hbkWaterMissingArgText
        AssetSourceDir = $corrKnnInputDir
        Assets = @("mini_corr_knn.pdb", "index.ndx")
        RequiredFiles = @(
            "hbond_kinetics_summary.dat",
            "hbond_kinetics_water_pairs.csv",
            "hbond_kinetics_water_ct.dat",
            "hbond_kinetics_water_lifetimes.dat"
        )
        TextChecks = @(
            @{ Path = "hbond_kinetics_summary.dat"; Contains = "water_stream_available: no" },
            @{ Path = "hbond_kinetics_summary.dat"; Contains = "water.intra_mol1.status: skipped_no_water_stream" }
        )
        RegexChecks = @(
            @{ Path = "hbond_kinetics_water_pairs.csv"; Pattern = "\AScope,WaterID,AtomAMol,AtomAChain,AtomAResID,AtomAResName,AtomAAtom,AtomBMol,AtomBChain,AtomBResID,AtomBResName,AtomBAtom,Occupancy_%,Events,MeanLifetime_frames,MeanLifetime_ps,MaxLifetime_frames,MaxLifetime_ps\r?\n?\z" }
        )
    },
    @{
        Name = "corr_knn_local_fixture"
        Category = "corr_knn"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_corr_knn_local.in"
        ArgContent = $corrKnnLocalArgText
        AssetSourceDir = $corrKnnInputDir
        Assets = @("mini_corr_knn.pdb", "index.ndx")
        RequiredFiles = @(
            "corr_knn_summary.dat",
            "corr_knn_mi_pairs.csv",
            "corr_knn_gc_pairs.csv",
            "corr_knn_te_pairs.csv"
        )
        TextChecks = @(
            @{ Path = "corr_knn_summary.dat"; Contains = "# method exact_ksg_mi + exact_frenzel_pompe_te" },
            @{ Path = "corr_knn_summary.dat"; Contains = "candidate_pair_count " },
            @{ Path = "log.txt"; Contains = "Corr-kNN" }
        )
    },
    @{
        Name = "corr_knn_local_no_te"
        Category = "corr_knn"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_corr_knn_no_te.in"
        ArgContent = $corrKnnNoTeArgText
        AssetSourceDir = $corrKnnInputDir
        Assets = @("mini_corr_knn.pdb", "index.ndx")
        RequiredFiles = @(
            "corr_knn_summary.dat",
            "corr_knn_mi_pairs.csv",
            "corr_knn_gc_pairs.csv",
            "corr_knn_te_pairs.csv"
        )
        RegexChecks = @(
            @{ Path = "corr_knn_te_pairs.csv"; Pattern = "\ASource,Target,TE\r?\n?\z" }
        )
        TextChecks = @(
            @{ Path = "log.txt"; Contains = "Corr-kNN TE Stage" }
        )
    },
    @{
        Name = "allostery_hist_local"
        Category = "allostery"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_allostery_hist_local.in"
        ArgContent = $allosteryHistLocalArgText
        AssetSourceDir = $corrKnnInputDir
        Assets = @("mini_corr_knn.pdb", "index.ndx")
        RequiredFiles = @(
            "raw_mi_edges.csv",
            "network_gc_edges.csv",
            "network_te_edges.csv"
        )
        TextChecks = @(
            @{ Path = "log.txt"; Contains = "Allostery Estimator" },
            @{ Path = "log.txt"; Contains = "hist" }
        )
        RegexChecks = @(
            @{ Path = "network_te_edges.csv"; Pattern = "\ASource,Target,Net_TE\r?\n?\z" }
        )
    },
    @{
        Name = "allostery_kde_local"
        Category = "allostery"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_allostery_kde_local.in"
        ArgContent = $allosteryKdeLocalArgText
        AssetSourceDir = $corrKnnInputDir
        Assets = @("mini_corr_knn.pdb", "index.ndx")
        RequiredFiles = @(
            "raw_mi_edges.csv",
            "network_gc_edges.csv",
            "network_te_edges.csv"
        )
        TextChecks = @(
            @{ Path = "log.txt"; Contains = "kde" }
        )
    },
    @{
        Name = "allostery_knn_local"
        Category = "allostery"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_allostery_knn_local.in"
        ArgContent = $allosteryKnnLocalArgText
        AssetSourceDir = $corrKnnInputDir
        Assets = @("mini_corr_knn.pdb", "index.ndx")
        RequiredFiles = @(
            "raw_mi_edges.csv",
            "network_gc_edges.csv",
            "network_te_edges.csv"
        )
        RegexChecks = @(
            @{ Path = "network_te_edges.csv"; Pattern = "\ASource,Target,Net_TE\r?\n?\z" }
        )
        TextChecks = @(
            @{ Path = "log.txt"; Contains = "titan_corr_knn_core" }
        )
    },
    @{
        Name = "corr_knn_allostery_reuse_local"
        Category = "corr_knn"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_corr_knn_allostery_reuse.in"
        ArgContent = $corrKnnAllosteryReuseArgText
        AssetSourceDir = $corrKnnInputDir
        Assets = @("mini_corr_knn.pdb", "index.ndx")
        RequiredFiles = @(
            "corr_knn_summary.dat",
            "corr_knn_mi_pairs.csv",
            "corr_knn_gc_pairs.csv",
            "corr_knn_te_pairs.csv",
            "raw_mi_edges.csv",
            "network_gc_edges.csv",
            "network_te_edges.csv"
        )
        TextChecks = @(
            @{ Path = "log.txt"; Contains = "Allostery kNN Reuse" },
            @{ Path = "log.txt"; Contains = "reused corr-kNN result" }
        )
        RegexChecks = @(
            @{ Path = "corr_knn_te_pairs.csv"; Pattern = "\ASource,Target,TE\r?\n?\z" },
            @{ Path = "network_te_edges.csv"; Pattern = "\ASource,Target,Net_TE\r?\n?\z" }
        )
    },
    @{
        Name = "graph_dccm_local"
        Category = "graph"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_graph_dccm_local.in"
        ArgContent = $graphDccmLocalArgText
        AssetSourceDir = $corrKnnInputDir
        Assets = @("mini_corr_knn.pdb", "index.ndx")
        RequiredFiles = @(
            "graph_summary.dat",
            "graph_dccm_edges.csv",
            "graph_dccm_nodes.csv",
            "graph_dccm_communities.csv",
            "graph_dccm_community_summary.csv",
            "graph_dccm_supergraph.csv",
            "graph_dccm_bottleneck_edges.csv",
            "graph_dccm_bottleneck_nodes.csv",
            "network_communities.csv",
            "network_community_summary.csv",
            "network_supergraph.csv"
        )
        TextChecks = @(
            @{ Path = "graph_summary.dat"; Contains = "[source.dccm]" },
            @{ Path = "graph_summary.dat"; Contains = "community_available: yes" },
            @{ Path = "log.txt"; Contains = "Graph Summary File" }
        )
    },
    @{
        Name = "graph_lspmd_local"
        Category = "graph"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_graph_lspmd_local.in"
        ArgContent = $graphLspmdLocalArgText
        AssetSourceDir = $corrKnnInputDir
        Assets = @("mini_corr_knn.pdb", "index.ndx")
        RequiredFiles = @(
            "graph_summary.dat",
            "graph_lspmd_edges.csv",
            "graph_lspmd_nodes.csv",
            "graph_lspmd_communities.csv",
            "graph_lspmd_bottleneck_edges.csv",
            "graph_lspmd_bottleneck_nodes.csv",
            "lspmd_degree_centrality.csv",
            "lspmd_betweenness_centrality.csv"
        )
        TextChecks = @(
            @{ Path = "graph_summary.dat"; Contains = "[source.lspmd]" },
            @{ Path = "log.txt"; Contains = "LSP-MD Analysis" }
        )
    },
    @{
        Name = "graph_allostery_gc_local"
        Category = "graph"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_graph_allostery_gc_local.in"
        ArgContent = $graphAllosteryGcLocalArgText
        AssetSourceDir = $corrKnnInputDir
        Assets = @("mini_corr_knn.pdb", "index.ndx")
        RequiredFiles = @(
            "graph_summary.dat",
            "graph_allostery-gc_edges.csv",
            "graph_allostery-gc_nodes.csv",
            "graph_allostery-gc_communities.csv",
            "graph_allostery-gc_bottleneck_edges.csv",
            "graph_allostery-gc_bottleneck_nodes.csv",
            "network_gc_edges.csv"
        )
        TextChecks = @(
            @{ Path = "graph_summary.dat"; Contains = "[source.allostery-gc]" },
            @{ Path = "graph_summary.dat"; Contains = "community_available:" }
        )
    },
    @{
        Name = "graph_allostery_te_local"
        Category = "graph"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_graph_allostery_te_local.in"
        ArgContent = $graphAllosteryTeLocalArgText
        AssetSourceDir = $corrKnnInputDir
        Assets = @("mini_corr_knn.pdb", "index.ndx")
        RequiredFiles = @(
            "graph_summary.dat",
            "graph_allostery-te_edges.csv",
            "graph_allostery-te_nodes.csv",
            "graph_allostery-te_bottleneck_edges.csv",
            "graph_allostery-te_bottleneck_nodes.csv",
            "network_te_edges.csv"
        )
        TextChecks = @(
            @{ Path = "graph_summary.dat"; Contains = "[source.allostery-te]" },
            @{ Path = "graph_summary.dat"; Contains = "community_available: no" }
        )
    },
    @{
        Name = "graph_tpt_local"
        Category = "graph"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_graph_tpt_local.in"
        ArgContent = $graphTptLocalArgText
        AssetSourceDir = $corrKnnInputDir
        Assets = @("mini_corr_knn.pdb", "index.ndx")
        RequiredFiles = @(
            "graph_summary.dat",
            "graph_tpt_edges.csv",
            "graph_tpt_nodes.csv",
            "graph_tpt_bottleneck_edges.csv",
            "graph_tpt_bottleneck_nodes.csv",
            "tpt_bottleneck_edges.dat",
            "tpt_bottleneck_states.dat"
        )
        TextChecks = @(
            @{ Path = "graph_summary.dat"; Contains = "[source.tpt]" },
            @{ Path = "graph_summary.dat"; Contains = "community_available: no" }
        )
    },
    @{
        Name = "graph_auto_local"
        Category = "graph"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_graph_auto_local.in"
        ArgContent = $graphAutoLocalArgText
        AssetSourceDir = $corrKnnInputDir
        Assets = @("mini_corr_knn.pdb", "index.ndx")
        RequiredFiles = @(
            "graph_summary.dat",
            "graph_dccm_edges.csv",
            "graph_lspmd_edges.csv",
            "graph_allostery-gc_edges.csv",
            "graph_allostery-te_edges.csv",
            "graph_tpt_edges.csv"
        )
        TextChecks = @(
            @{ Path = "graph_summary.dat"; Contains = "requested_sources: auto" },
            @{ Path = "graph_summary.dat"; Contains = "[source.dccm]" },
            @{ Path = "graph_summary.dat"; Contains = "[source.tpt]" }
        )
    },
    @{
        Name = "pca_example3"
        Category = "pca"
        Scale = "medium"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_pca_only.in"
        ArgContent = $pcaArg
    },
    @{
        Name = "llps_example5"
        Category = "llps"
        Scale = "medium"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgPath = $example5Arg
    },
    @{
        Name = "sq_example5_channels"
        Category = "sq"
        Scale = "medium"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_sq_channels.in"
        ArgContent = $sqExample5Arg
        CleanupPaths = $example5SqRequiredPaths
        RequiredPaths = $example5SqRequiredPaths
        TextChecks = @(
            @{ Path = "log.txt"; Contains = "SQ Analysis" }
        )
        PathTextChecks = @(
            @{ Path = (Join-Path $example5SqOutputDir "sq_summary.dat"); Contains = "channel.backbone.status: active" },
            @{ Path = (Join-Path $example5SqOutputDir "sq_summary.dat"); Contains = "channel.sticker.status: active" },
            @{ Path = (Join-Path $example5SqOutputDir "sq_summary.dat"); Contains = "pair.backbone__sticker.ordered_pairs_per_frame:" }
        )
    },
    @{
        Name = "gist_nt1"
        Category = "gist"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgPath = $gistNt1Arg
    },
    @{
        Name = "gist_nt4"
        Category = "gist"
        Scale = "medium"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgPath = $gistNt4Arg
    },
    @{
        Name = "mbar_demux_tica"
        Category = "mbar"
        Scale = "small"
        RunByDefault = $false
        ExpectedOutcome = "success"
        ArgFileName = "arg_mbar_demux.in"
        ArgContent = $mbarDemuxArg
        AssetSourceDir = $mbarInputDir
        Assets = @(
            "energy_rep1_tremd_nvt.xvg",
            "energy_rep2_tremd_nvt.xvg"
        )
        ExtraFiles = @(
            @{ Name = "xtc_list_twice.txt"; Content = $mbarXtcList },
            @{ Name = "temps_tremd_valid.txt"; Content = $mbarTempsValid },
            @{ Name = "energy_list_tremd_nvt.txt"; Content = $mbarEnergyListTremdNvt },
            @{ Name = "demux_map_tremd.txt"; Content = $mbarDemuxMap }
        )
        RequiredFiles = @(
            "mbar_summary.dat",
            "mbar_window_audit.dat",
            "mbar_frame_weights.dat",
            "mbar_cv.dat",
            "mbar_tica_fes_2d.dat",
            "mbar_tica_cluster_representatives.csv"
        )
        TextChecks = @(
            @{ Path = "mbar_summary.dat"; Contains = "MBAR Converged                 : yes" },
            @{ Path = "mbar_summary.dat"; Contains = "MBAR Temp Mode                 : per_frame" },
            @{ Path = "mbar_summary.dat"; Contains = "MBAR Demux Map                 : ./demux_map_tremd.txt" }
        )
    },
    @{
        Name = "xdr_tica_umap"
        Category = "xdr"
        Scale = "large"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgPath = $xdrTicaUmapArg
        ForbiddenFiles = @(
            "umap_fes_2d.dat",
            "umap_cluster_representatives.csv",
            "umap_equilibrium_summary.dat"
        )
    },
    @{
        Name = "mbar_dccm_lspmd_context"
        Category = "mbar"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_mbar_dccm_lspmd.in"
        ArgContent = $mbarDccmContextArg
        AssetSourceDir = $mbarInputDir
        Assets = @(
            "energy_rep1.xvg",
            "energy_rep2.xvg"
        )
        ExtraFiles = @(
            @{ Name = "xtc_list_twice.txt"; Content = $mbarXtcList },
            @{ Name = "temps_valid.txt"; Content = $mbarTempsValid },
            @{ Name = "energy_list_valid.txt"; Content = $mbarEnergyListFixed }
        )
        RequiredFiles = @(
            "mbar_summary.dat",
            "dccm_equilibrium_context.dat",
            "allostery_equilibrium_ranking.csv",
            "lspmd_equilibrium_ranking.csv",
            "raw_mi_edges.csv",
            "network_gc_edges.csv",
            "lspmd_edges.csv",
            "lspmd_degree_centrality.csv",
            "lspmd_betweenness_centrality.csv"
        )
        TextChecks = @(
            @{ Path = "log.txt"; Contains = "Equilibrium Weight Source      : mbar" },
            @{ Path = "dccm_equilibrium_context.dat"; Contains = "# source mbar" },
            @{ Path = "allostery_equilibrium_ranking.csv"; Contains = "# note contextual_only" },
            @{ Path = "lspmd_equilibrium_ranking.csv"; Contains = "# note contextual_only" }
        )
    },
    @{
        Name = "mbar_tpt_overlay"
        Category = "mbar"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_mbar_tpt_overlay.in"
        ArgContent = $mbarTptOverlayArg
        AssetSourceDir = $mbarInputDir
        Assets = @(
            "energy_rep1.xvg",
            "energy_rep2.xvg"
        )
        ExtraFiles = @(
            @{ Name = "xtc_list_twice.txt"; Content = $mbarXtcList },
            @{ Name = "temps_valid.txt"; Content = $mbarTempsValid },
            @{ Name = "energy_list_valid.txt"; Content = $mbarEnergyListFixed }
        )
        RequiredFiles = @(
            "mbar_summary.dat",
            "msm_equilibrium_summary.dat",
            "tpt_summary.dat",
            "tpt_committor.dat",
            "tpt_flux_edges.dat",
            "tpt_paths.dat",
            "tpt_state_labels.dat",
            "tpt_all_pairs_kinetics.dat",
            "tpt_equilibrium_overlay.dat"
        )
        TextChecks = @(
            @{ Path = "log.txt"; Contains = "Equilibrium Weight Source      : mbar" },
            @{ Path = "log.txt"; Contains = "TPT Model Source               : msm_built_model" },
            @{ Path = "tpt_equilibrium_overlay.dat"; Contains = "# source mbar" }
        )
    },
    @{
        Name = "mbar_npt_pca"
        Category = "mbar"
        Scale = "medium"
        RunByDefault = $false
        ExpectedOutcome = "success"
        ArgFileName = "arg_mbar_npt_pca.in"
        ArgContent = $mbarNptArg
        AssetSourceDir = $mbarInputDir
        Assets = @(
            "energy_rep1_tremd_npt.xvg",
            "energy_rep2_tremd_npt.xvg"
        )
        ExtraFiles = @(
            @{ Name = "xtc_list_twice.txt"; Content = $mbarXtcList },
            @{ Name = "temps_tremd_valid.txt"; Content = $mbarTempsValid },
            @{ Name = "energy_list_tremd_npt.txt"; Content = $mbarEnergyListTremdNpt }
        )
        RequiredFiles = @(
            "mbar_summary.dat",
            "mbar_window_audit.dat",
            "mbar_frame_weights.dat",
            "mbar_cv.dat",
            "mbar_pca_fes_2d.dat",
            "mbar_pca_cluster_representatives.csv"
        )
        TextChecks = @(
            @{ Path = "mbar_summary.dat"; Contains = "MBAR Converged                 : yes" },
            @{ Path = "mbar_summary.dat"; Contains = "MBAR Ensemble                  : npt" },
            @{ Path = "mbar_summary.dat"; Contains = "MBAR Domain                    : pca" }
        )
    },
    @{
        Name = "mbar_auto_guard_matrix"
        Category = "mbar"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_mbar_auto_guard.in"
        ArgContent = $eqAutoGuardArg
        AssetSourceDir = $mbarInputDir
        Assets = @(
            "energy_rep1.xvg",
            "energy_rep2.xvg"
        )
        ExtraFiles = @(
            @{ Name = "xtc_list_twice.txt"; Content = $mbarXtcList },
            @{ Name = "temps_valid.txt"; Content = $mbarTempsValid },
            @{ Name = "energy_list_valid.txt"; Content = $mbarEnergyListFixed }
        )
        RequiredFiles = @(
            "mbar_summary.dat",
            "tica_fes_2d.dat",
            "titan_umap_embedding.csv",
            "tpt_summary.dat",
            "tpt_committor.dat",
            "raw_mi_edges.csv",
            "network_gc_edges.csv",
            "lspmd_edges.csv"
        )
        ForbiddenFiles = @(
            "msm_equilibrium_summary.dat",
            "tpt_equilibrium_overlay.dat",
            "umap_fes_2d.dat",
            "umap_cluster_representatives.csv",
            "umap_equilibrium_summary.dat",
            "dccm_equilibrium_context.dat",
            "allostery_equilibrium_ranking.csv",
            "lspmd_equilibrium_ranking.csv"
        )
        TextChecks = @(
            @{ Path = "log.txt"; Contains = "Equilibrium Weight Source      : uniform" },
            @{ Path = "log.txt"; Contains = "TICA Weight Mode               : uniform" }
        )
    },
    @{
        Name = "mbar_uniform_guard_matrix"
        Category = "mbar"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_mbar_uniform_guard.in"
        ArgContent = $eqUniformGuardArg
        AssetSourceDir = $mbarInputDir
        Assets = @(
            "energy_rep1.xvg",
            "energy_rep2.xvg"
        )
        ExtraFiles = @(
            @{ Name = "xtc_list_twice.txt"; Content = $mbarXtcList },
            @{ Name = "temps_valid.txt"; Content = $mbarTempsValid },
            @{ Name = "energy_list_valid.txt"; Content = $mbarEnergyListFixed }
        )
        RequiredFiles = @(
            "mbar_summary.dat",
            "tica_fes_2d.dat",
            "titan_umap_embedding.csv",
            "tpt_summary.dat",
            "tpt_committor.dat",
            "raw_mi_edges.csv",
            "network_gc_edges.csv",
            "lspmd_edges.csv"
        )
        ForbiddenFiles = @(
            "msm_equilibrium_summary.dat",
            "tpt_equilibrium_overlay.dat",
            "umap_fes_2d.dat",
            "umap_cluster_representatives.csv",
            "umap_equilibrium_summary.dat",
            "dccm_equilibrium_context.dat",
            "allostery_equilibrium_ranking.csv",
            "lspmd_equilibrium_ranking.csv"
        )
        TextChecks = @(
            @{ Path = "log.txt"; Contains = "Equilibrium Weight Source      : uniform" },
            @{ Path = "log.txt"; Contains = "TICA Weight Mode               : uniform" }
        )
    },
    @{
        Name = "dham_guard_matrix"
        Category = "dham"
        Scale = "small"
        RunByDefault = $false
        ExpectedOutcome = "success"
        ArgFileName = "arg_dham_guard.in"
        ArgContent = $eqDhamGuardArg
        AssetSourceDir = $mbarInputDir
        Assets = @(
            "energy_rep1.xvg",
            "energy_rep2.xvg"
        )
        ExtraFiles = @(
            @{ Name = "xtc_list_twice.txt"; Content = $mbarXtcList },
            @{ Name = "temps_valid.txt"; Content = $mbarTempsValid },
            @{ Name = "energy_list_valid.txt"; Content = $mbarEnergyListFixed }
        )
        RequiredFiles = @(
            "dham_stationary_distribution.dat",
            "dham_free_energy_profile.dat",
            "dham_unbiased_tpm.dat",
            "dham_window_audit.dat",
            "tica_fes_2d.dat",
            "titan_umap_embedding.csv",
            "tpt_summary.dat",
            "raw_mi_edges.csv",
            "network_gc_edges.csv",
            "lspmd_edges.csv"
        )
        ForbiddenFiles = @(
            "mbar_summary.dat",
            "msm_equilibrium_summary.dat",
            "tpt_equilibrium_overlay.dat",
            "umap_fes_2d.dat",
            "umap_cluster_representatives.csv",
            "umap_equilibrium_summary.dat",
            "dccm_equilibrium_context.dat",
            "allostery_equilibrium_ranking.csv",
            "lspmd_equilibrium_ranking.csv"
        )
        TextChecks = @(
            @{ Path = "log.txt"; Contains = "Equilibrium Weight Source      : dham" },
            @{ Path = "log.txt"; Contains = "TICA Weight Mode               : DHAM-weighted" }
        )
        AllowedNanInfPatterns = @(
            '^\d+\s+inf$',
            '^mfpt_source_to_sink_lags inf$',
            '^\d+\s+\d+\s+.*\sinf\s+nan(\s+.*)?$'
        )
    },
    @{
        Name = "mbar_unresolved_guard"
        Category = "mbar"
        Scale = "small"
        RunByDefault = $true
        ExpectedOutcome = "success"
        ArgFileName = "arg_mbar_unresolved_guard.in"
        ArgContent = $mbarMismatchGuardArg
        AssetSourceDir = $mbarInputDir
        Assets = @(
            "energy_rep1.xvg",
            "energy_rep2_shifted.xvg"
        )
        ExtraFiles = @(
            @{ Name = "xtc_list_twice.txt"; Content = $mbarXtcList },
            @{ Name = "temps_valid.txt"; Content = $mbarTempsValid },
            @{ Name = "energy_list_mismatch.txt"; Content = $mbarEnergyListMismatch }
        )
        RequiredFiles = @(
            "mbar_summary.dat",
            "mbar_window_audit.dat",
            "mbar_frame_weights.dat",
            "tica_fes_2d.dat",
            "titan_umap_embedding.csv",
            "tpt_summary.dat",
            "raw_mi_edges.csv",
            "network_gc_edges.csv",
            "lspmd_edges.csv"
        )
        ForbiddenFiles = @(
            "msm_equilibrium_summary.dat",
            "tpt_equilibrium_overlay.dat",
            "umap_fes_2d.dat",
            "umap_cluster_representatives.csv",
            "umap_equilibrium_summary.dat",
            "dccm_equilibrium_context.dat",
            "allostery_equilibrium_ranking.csv",
            "lspmd_equilibrium_ranking.csv"
        )
        TextChecks = @(
            @{ Path = "mbar_summary.dat"; Contains = "MBAR Abort Reason              : temps/energy count mismatch" },
            @{ Path = "log.txt"; Contains = "Equilibrium Weight Source      : uniform" }
        )
    },
    @{
        Name = "allatom_8cbf_ba1"
        Category = "allatom"
        Scale = "large"
        RunByDefault = $false
        ExpectedOutcome = "success"
        ArgPath = $example8Arg
        AssetSourceDir = (Join-Path $repoRoot "example\\8cbf_ba1")
        Assets = @(
            "input.txt",
            "md1.pdb",
            "md2.pdb",
            "topol.top",
            "index.ndx",
            "toppar"
        )
    }
)

if ($SkipInvalidCases) {
    $caseTable = @($caseTable | Where-Object { $_.ExpectedOutcome -eq "success" })
}

if ($Cases.Count -gt 0) {
    $requested = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::OrdinalIgnoreCase)
    foreach ($name in $Cases) {
        [void]$requested.Add($name)
    }
    $caseTable = @($caseTable | Where-Object { $requested.Contains($_.Name) })
    if ($caseTable.Count -eq 0) {
        throw "No benchmark cases matched the requested filter."
    }
} else {
    $caseTable = @($caseTable | Where-Object {
        -not $_.ContainsKey("RunByDefault") -or $_.RunByDefault
    })
    if ($Suite -ne "all") {
        $caseTable = @($caseTable | Where-Object { $_.Scale -eq $Suite })
        if ($caseTable.Count -eq 0) {
            throw "No benchmark cases matched the requested suite."
        }
    }
}

$sessionTag = if ([string]::IsNullOrWhiteSpace($SessionName)) { "benchmark_run_{0}" -f (Get-Date -Format "yyyyMMdd_HHmmss") } else { $SessionName }
$sessionRoot = $testRoot
$results = @()

foreach ($case in $caseTable) {
    $results += Invoke-TitanCase -Case $case -SessionRoot $sessionRoot -SessionTag $sessionTag -TitanExe $titanExe -RepoRoot $repoRoot
}

$gistComparison = $null
$gistNt1 = $results | Where-Object { $_.name -eq "gist_nt1" } | Select-Object -First 1
$gistNt4 = $results | Where-Object { $_.name -eq "gist_nt4" } | Select-Object -First 1
if ($gistNt1 -and $gistNt4 -and $gistNt1.passed -and $gistNt4.passed) {
    $gistComparison = Get-GistComparison -Nt1Case $gistNt1 -Nt4Case $gistNt4
}

$crossCaseHashRules = @(
    @{
        Label = "auto_vs_uniform_tica_fes"
        CaseA = "mbar_auto_guard_matrix"
        FileA = "tica_fes_2d.dat"
        CaseB = "mbar_uniform_guard_matrix"
        FileB = "tica_fes_2d.dat"
    },
    @{
        Label = "auto_vs_uniform_umap_embedding"
        CaseA = "mbar_auto_guard_matrix"
        FileA = "titan_umap_embedding.csv"
        CaseB = "mbar_uniform_guard_matrix"
        FileB = "titan_umap_embedding.csv"
    },
    @{
        Label = "auto_vs_uniform_tpt_summary"
        CaseA = "mbar_auto_guard_matrix"
        FileA = "tpt_summary.dat"
        CaseB = "mbar_uniform_guard_matrix"
        FileB = "tpt_summary.dat"
    },
    @{
        Label = "auto_vs_uniform_raw_mi"
        CaseA = "mbar_auto_guard_matrix"
        FileA = "raw_mi_edges.csv"
        CaseB = "mbar_uniform_guard_matrix"
        FileB = "raw_mi_edges.csv"
    }
)
$crossCaseComparisons = Invoke-CrossCaseHashValidation -Results $results -Rules $crossCaseHashRules

$summaryPath = Join-Path $testRoot ("benchmark_summary_{0}.md" -f (Get-Date -Format "yyyyMMdd_HHmmss"))
Write-Summary -SummaryPath $summaryPath -SessionRoot $sessionRoot -Results $results -GistComparison $gistComparison -CrossCaseComparisons $crossCaseComparisons -Suite $Suite -BuildProfile $BuildProfile

$latestSummaryPath = Join-Path $testRoot "benchmark_latest.md"
Copy-Item -LiteralPath $summaryPath -Destination $latestSummaryPath -Force

Write-Host "Benchmark session complete."
Write-Host "Session root: $sessionRoot"
Write-Host "Summary: $summaryPath"
