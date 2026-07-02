param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_cli_config_parser')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

function Invoke-PbctopoParserCase(
    [string]$Name,
    [string[]]$CaseArgs,
    [string]$ExpectedMessage
) {
    $stdout = Join-Path $WorkRoot "$Name.stdout.log"
    $stderr = Join-Path $WorkRoot "$Name.stderr.log"
    $process = Start-Process -FilePath $TitanExe `
        -ArgumentList $CaseArgs `
        -RedirectStandardOutput $stdout `
        -RedirectStandardError $stderr `
        -NoNewWindow `
        -Wait `
        -PassThru
    $exitCode = $process.ExitCode
    $log = (Get-Content -Raw $stdout) + (Get-Content -Raw $stderr)
    if ($exitCode -eq 0) {
        throw "$Name unexpectedly succeeded"
    }
    if ($log -notmatch $ExpectedMessage) {
        throw "$Name error was unclear: $log"
    }
}

Invoke-PbctopoParserCase 'bad_bond_cutoff' `
    @('-pbctopo', '-pbctopo-bond-cutoff=-1') `
    'pbctopo-bond-cutoff.*> 0'
Invoke-PbctopoParserCase 'bad_clash_cutoff' `
    @('-pbctopo', '-pbctopo-clash-cutoff=0') `
    'pbctopo-clash-cutoff.*> 0'
Invoke-PbctopoParserCase 'bad_contact_cutoff' `
    @('-pbctopo', '-pbctopo-contact-cutoff=-0.5') `
    'pbctopo-contact-cutoff.*> 0'
Invoke-PbctopoParserCase 'bad_continuity_cutoff' `
    @('-pbctopo', '-pbctopo-continuity-cutoff=0') `
    'pbctopo-continuity-cutoff.*> 0'
Invoke-PbctopoParserCase 'bad_exact_limit' `
    @('-pbctopo', '-pbctopo-exact-component-limit=0') `
    'pbctopo-exact-component-limit.*> 0'
Invoke-PbctopoParserCase 'bad_exact_limit_too_high' `
    @('-pbctopo', '-pbctopo-exact-component-limit=17') `
    'pbctopo-exact-component-limit.*<= 16'
Invoke-PbctopoParserCase 'bad_beam_width' `
    @('-pbctopo', '-pbctopo-beam-width=0') `
    'pbctopo-beam-width.*> 0'
Invoke-PbctopoParserCase 'bad_beam_width_too_high' `
    @('-pbctopo', '-pbctopo-beam-width=4097') `
    'pbctopo-beam-width.*<= 4096'
Invoke-PbctopoParserCase 'bad_anchor_backend' `
    @('-pbctopo', '-pbctopo-anchor-backend=bad') `
    'pbctopo-anchor-backend.*auto.*prim.*spatial-all'
Invoke-PbctopoParserCase 'bad_pbctopo_mode' `
    @('-pbctopo', '-pbctopo-mode=bad') `
    'pbctopo-mode.*vibe.*llps-local'
Invoke-PbctopoParserCase 'missing_pbctopo_mode' `
    @('-pbctopo', '-pbctopo-mode') `
    'pbctopo-mode.*requires a value'
Invoke-PbctopoParserCase 'missing_bond_cutoff' `
    @('-pbctopo', '-pbctopo-bond-cutoff') `
    'pbctopo-bond-cutoff.*requires a value'

Write-Host 'PASS pbctopo_cli_config_parser_smoke'
