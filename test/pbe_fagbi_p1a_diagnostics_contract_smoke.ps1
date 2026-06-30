param()

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$Root = Split-Path -Parent $PSScriptRoot
Set-Location $Root

function Assert-Contains {
  param(
    [Parameter(Mandatory=$true)][string]$Path,
    [Parameter(Mandatory=$true)][string]$Pattern,
    [Parameter(Mandatory=$true)][string]$Message
  )
  $text = Get-Content -Raw -LiteralPath $Path
  if ($text -notmatch $Pattern) {
    throw $Message
  }
}

$compactHeader =
  'trajectory_id,frame_index,global_index,fixture,atom_count,panels,energy_kJmol,iterations,internal_residual,sampled_direct_residual,audit_rows,wall_seconds,rss_mb,p_order,theta,guard,tail_tol,ilu_cap,far_pairs,fallback_pairs,status'

Assert-Contains -Path 'src/arguments.h' `
  -Pattern 'tabi_diagnostics\{"off"\}' `
  -Message 'P1A missing PBE TABI diagnostics option default.'

Assert-Contains -Path 'src/parse_arguments.cpp' `
  -Pattern 'pbe-tabi-diagnostics' `
  -Message 'P1A missing -pbe-tabi-diagnostics parser wiring.'

Assert-Contains -Path 'src/parse_arguments.cpp' `
  -Pattern '--diagnostics' `
  -Message 'P1A missing public --diagnostics alias.'

Assert-Contains -Path 'src/arguments.cpp' `
  -Pattern 'tabi_diagnostics.*off.*audit.*full' `
  -Message 'P1A missing diagnostics validation for off/audit/full.'

Assert-Contains -Path 'src/tabi/analyzer_impl.cpp' `
  -Pattern 'pbe_tabi_audit_profile\.csv' `
  -Message 'P1A missing compact paper-facing TABI audit profile output.'

Assert-Contains -Path 'src/tabi/analyzer_impl.cpp' `
  -Pattern ([regex]::Escape($compactHeader)) `
  -Message 'P1A compact audit profile schema drifted or is missing.'

Assert-Contains -Path 'docs/developer_diagnostics.md' `
  -Pattern 'diagnostics.*off.*audit.*full' `
  -Message 'P1A missing developer diagnostics tier documentation.'

'p1a diagnostics contract ok'
