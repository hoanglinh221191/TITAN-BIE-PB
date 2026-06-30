# TITAN TABI Diagnostics

This note separates release-facing TABI output from developer-only telemetry.

## Tiers

Use `-pbe-tabi-diagnostics=<off|audit|full>` or the public alias
`--diagnostics=<off|audit|full>`.

- `off`: default. TITAN writes the normal production outputs such as
  `pbe_energy.csv` and `pbe_summary.dat`; no compact TABI audit profile is
  produced.
- `audit`: writes `pbe_tabi_audit_profile.csv`, a compact paper-facing table
  intended for archived benchmark rows and reviewer reproduction.
- `full`: writes the same compact audit profile and permits the existing
  developer telemetry requested through older switches such as
  `-pbe-timing=detail`, `-pbe-tabi-matvec-audit=sampled-detail`, and the
  individual FMM/M2L profile flags.

The compact audit profile has one row per TABI solve with:

```text
fixture, panels, energy, iterations, internal_residual,
sampled_direct_residual, audit_rows, wall_seconds, rss_mb,
p_order, theta, guard, tail_tol, ilu_cap, far_pairs, fallback_pairs
```

The long developer CSVs retain detailed estimator, over-rejection, expansion,
leaf-work, fallback, and memory counters. They are useful for development and
forensics, but they are not the manuscript-facing contract.
