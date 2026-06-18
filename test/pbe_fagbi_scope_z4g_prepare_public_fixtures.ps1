param(
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4g_fixture_prep'),
    [string]$PythonExe = 'C:\Users\ADMIN\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe',
    [string]$Pdb2PqrTarget = '',
    [switch]$RefreshFixtures
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$pdb2pqrVersion = '3.7.1'
$pdb2pqrArguments = '--ff=AMBER --ffout=AMBER --keep-chain --titration-state-method=propka --with-ph=7.0'
$radiusFloorAngstrom = 0.6
$preparationPolicy = "$pdb2pqrArguments; postprocess-radius-floor=0.6000A"
$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
if ([string]::IsNullOrWhiteSpace($Pdb2PqrTarget)) {
    $Pdb2PqrTarget = Join-Path $RepoRoot '.codex_tmp\pdb2pqr-3.7.1'
}

function Get-Z4gSha256([string]$Path) {
    return (Get-FileHash -Algorithm SHA256 -LiteralPath $Path).Hash.ToLowerInvariant()
}

function Get-Z4gAtomCount([string]$Path) {
    return @(Get-Content -LiteralPath $Path | Where-Object {
            $_.StartsWith('ATOM') -or $_.StartsWith('HETATM')
        }).Count
}

function Write-Z4gModelChainPdb(
    [string]$InputPath,
    [string]$OutputPath,
    [string]$Chain = 'A',
    [int]$Model = 1
) {
    $hasModels = @(Select-String -LiteralPath $InputPath -Pattern '^MODEL\s').Count -gt 0
    $activeModel = if ($hasModels) { 0 } else { 1 }
    $out = New-Object System.Collections.Generic.List[string]
    foreach ($line in Get-Content -LiteralPath $InputPath) {
        if ($line.StartsWith('MODEL')) {
            $number = 0
            if ([int]::TryParse($line.Substring(5).Trim(), [ref]$number)) {
                $activeModel = $number
            }
            continue
        }
        if ($line.StartsWith('ENDMDL')) {
            if ($activeModel -eq $Model) { break }
            $activeModel = 0
            continue
        }
        if ($activeModel -ne $Model -or -not $line.StartsWith('ATOM')) {
            continue
        }
        $lineChain = if ($line.Length -gt 21) { [string]$line[21] } else { '' }
        if ($lineChain -eq $Chain) {
            $out.Add($line) | Out-Null
        }
    }
    if ($out.Count -eq 0) {
        throw "No ATOM rows found for model $Model chain $Chain in $InputPath"
    }
    $out.Add('END') | Out-Null
    Set-Content -LiteralPath $OutputPath -Value $out -Encoding ASCII
}

function Install-Z4gPdb2Pqr {
    if (-not (Test-Path -LiteralPath $PythonExe)) {
        throw "Bundled Python is missing: $PythonExe"
    }
    New-Item -ItemType Directory -Path $Pdb2PqrTarget -Force | Out-Null
    $marker = Join-Path $Pdb2PqrTarget 'pdb2pqr\__init__.py'
    if (-not (Test-Path -LiteralPath $marker)) {
        & $PythonExe -m pip install --disable-pip-version-check --target $Pdb2PqrTarget `
            ("pdb2pqr=={0}" -f $pdb2pqrVersion)
        if ($LASTEXITCODE -ne 0) { throw 'Unable to install pinned PDB2PQR.' }
    }
}

function Invoke-Z4gPdb2Pqr([string]$InputPath, [string]$OutputPath) {
    $oldPythonPath = $env:PYTHONPATH
    try {
        $env:PYTHONPATH = if ([string]::IsNullOrWhiteSpace($oldPythonPath)) {
            $Pdb2PqrTarget
        } else {
            "$Pdb2PqrTarget;$oldPythonPath"
        }
        & $PythonExe -m pdb2pqr --ff=AMBER --ffout=AMBER --keep-chain `
            --titration-state-method=propka --with-ph=7.0 $InputPath $OutputPath
        if ($LASTEXITCODE -ne 0) { throw "PDB2PQR failed for $InputPath" }
    } finally {
        $env:PYTHONPATH = $oldPythonPath
    }
}

function Normalize-Z4gPqrRadii([string]$Path, [double]$MinimumRadius) {
    $culture = [System.Globalization.CultureInfo]::InvariantCulture
    $output = New-Object System.Collections.Generic.List[string]
    $clamped = 0
    foreach ($line in Get-Content -LiteralPath $Path) {
        if (-not ($line.StartsWith('ATOM') -or $line.StartsWith('HETATM'))) {
            $output.Add($line) | Out-Null
            continue
        }
        $fields = @($line -split '\s+' | Where-Object { $_ -ne '' })
        if ($fields.Count -lt 10) { throw "Malformed PDB2PQR atom row: $line" }
        $radius = 0.0
        if (-not [double]::TryParse($fields[-1], [System.Globalization.NumberStyles]::Float, $culture, [ref]$radius)) {
            throw "Malformed PDB2PQR radius: $line"
        }
        if ($radius -le 0.0) {
            $tokenStart = $line.LastIndexOf($fields[-1], [System.StringComparison]::Ordinal)
            if ($tokenStart -lt 0) { throw "Unable to locate PDB2PQR radius token: $line" }
            $line = $line.Substring(0, $tokenStart) + $MinimumRadius.ToString('F4', $culture)
            ++$clamped
        }
        $output.Add($line) | Out-Null
    }
    Set-Content -LiteralPath $Path -Value $output -Encoding ASCII
    return $clamped
}

$specs = @(
    [pscustomobject]@{
        fixture='barnase'; pdb_id='PyDelPhi benchmark';
        protein_name='Barnase'; model=1; chain='A';
        source_url='https://github.com/delphi001/pyDelPhi/tree/master/pydelphi/data/test_examples/barnase';
        download_url=''; raw_name='barnase_source.pdb'; pdb_name='barnase_model1_chainA.pdb';
        pqr_name='barnase_full_amber_ph7.pqr'
    },
    [pscustomobject]@{
        fixture='1a63'; pdb_id='1A63';
        protein_name='RNA-binding domain of E. coli Rho factor'; model=1; chain='A';
        source_url='https://www.rcsb.org/structure/1A63';
        download_url='https://files.rcsb.org/download/1A63.pdb'; raw_name='1a63_source.pdb';
        pdb_name='1a63_model1_chainA.pdb'; pqr_name='1a63_full_amber_ph7.pqr'
    },
    [pscustomobject]@{
        fixture='1sh1'; pdb_id='1SH1';
        protein_name='Stichodactyla helianthus neurotoxin I'; model=1; chain='A';
        source_url='https://www.rcsb.org/structure/1SH1';
        download_url='https://files.rcsb.org/download/1SH1.pdb'; raw_name='1sh1_source.pdb';
        pdb_name='1sh1_model1_chainA.pdb'; pqr_name='1sh1_full_amber_ph7.pqr'
    }
)

if ($RefreshFixtures) {
    New-Item -ItemType Directory -Path $FixtureRoot -Force | Out-Null
    New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
    Install-Z4gPdb2Pqr
    $rows = New-Object System.Collections.Generic.List[object]
    foreach ($spec in $specs) {
        $rawPath = Join-Path $FixtureRoot $spec.raw_name
        if ($spec.fixture -eq 'barnase') {
            $candidates = @(
                (Join-Path $RepoRoot '.codex_tmp\pyDelPhi\pydelphi\data\test_examples\barnase\barnase.pdb'),
                'C:\Users\ADMIN\OneDrive\Documents\antigravity\elec_vdw\.codex_tmp\pyDelPhi\pydelphi\data\test_examples\barnase\barnase.pdb'
            )
            $source = @($candidates | Where-Object { Test-Path -LiteralPath $_ }) | Select-Object -First 1
            if ([string]::IsNullOrWhiteSpace([string]$source)) {
                throw 'Missing public PyDelPhi Barnase source fixture.'
            }
            Copy-Item -LiteralPath $source -Destination $rawPath -Force
        } else {
            Invoke-WebRequest -UseBasicParsing -Uri $spec.download_url -OutFile $rawPath
        }
        $pdbPath = Join-Path $FixtureRoot $spec.pdb_name
        $pqrPath = Join-Path $FixtureRoot $spec.pqr_name
        Write-Z4gModelChainPdb -InputPath $rawPath -OutputPath $pdbPath `
            -Chain $spec.chain -Model $spec.model
        Invoke-Z4gPdb2Pqr -InputPath $pdbPath -OutputPath $pqrPath
        $radiusClampedAtoms = Normalize-Z4gPqrRadii -Path $pqrPath -MinimumRadius $radiusFloorAngstrom
        $atomCount = Get-Z4gAtomCount $pqrPath
        if ($atomCount -lt 100) { throw "Unexpectedly small full PQR: $pqrPath" }
        if (@(Get-Content -LiteralPath $pqrPath | Where-Object { $_.StartsWith('HETATM') }).Count -ne 0) {
            throw "Unsupported HETATM rows remain in $pqrPath"
        }
        $rows.Add([pscustomobject]@{
                fixture=$spec.fixture; pdb_id=$spec.pdb_id; protein_name=$spec.protein_name;
                model=$spec.model; chain=$spec.chain; source_url=$spec.source_url;
                raw_pdb=$spec.raw_name; prepared_pdb=$spec.pdb_name; pqr=$spec.pqr_name;
                raw_pdb_sha256=(Get-Z4gSha256 $rawPath);
                prepared_pdb_sha256=(Get-Z4gSha256 $pdbPath);
                pdb_sha256=(Get-Z4gSha256 $pdbPath); pqr_sha256=(Get-Z4gSha256 $pqrPath);
                atom_count=$atomCount; pdb2pqr_version=$pdb2pqrVersion;
                preparation=$preparationPolicy; radius_floor_angstrom=$radiusFloorAngstrom;
                radius_clamped_atoms=$radiusClampedAtoms
            }) | Out-Null
    }
    $rows.ToArray() | Export-Csv -LiteralPath $manifestPath -NoTypeInformation
}

if (-not (Test-Path -LiteralPath $manifestPath)) {
    throw "Missing cached Z4G fixture manifest. Run with -RefreshFixtures: $manifestPath"
}
$manifest = @(Import-Csv -LiteralPath $manifestPath)
if ($manifest.Count -ne 3) { throw "Expected three Z4G fixtures, found $($manifest.Count)" }
foreach ($row in $manifest) {
    $pdbPath = Join-Path $FixtureRoot $row.prepared_pdb
    $pqrPath = Join-Path $FixtureRoot $row.pqr
    $rawPath = Join-Path $FixtureRoot $row.raw_pdb
    if (-not (Test-Path -LiteralPath $pdbPath) -or -not (Test-Path -LiteralPath $pqrPath)) {
        throw "Cached Z4G fixture is incomplete for $($row.fixture). Run -RefreshFixtures."
    }
    if (-not (Test-Path -LiteralPath $rawPath) -or
        (Get-Z4gSha256 $rawPath) -ne $row.raw_pdb_sha256 -or
        (Get-Z4gSha256 $pdbPath) -ne $row.prepared_pdb_sha256 -or
        (Get-Z4gSha256 $pdbPath) -ne $row.pdb_sha256 -or
        (Get-Z4gSha256 $pqrPath) -ne $row.pqr_sha256) {
        throw "Checksum drift for Z4G fixture $($row.fixture). Run -RefreshFixtures deliberately."
    }
    if ((Get-Z4gAtomCount $pqrPath) -ne [int]$row.atom_count) {
        throw "Atom-count drift for Z4G fixture $($row.fixture)."
    }
    if ($row.pdb2pqr_version -ne $pdb2pqrVersion -or $row.preparation -ne $preparationPolicy -or
        $row.PSObject.Properties.Name -notcontains 'radius_floor_angstrom' -or
        [double]$row.radius_floor_angstrom -ne $radiusFloorAngstrom) {
        throw "Preparation policy drift for Z4G fixture $($row.fixture)."
    }
}

Write-Host "PASS Scope Z4G public fixtures: $manifestPath"
