# FMM/M2L Scope G Real-Surface Offline Comparison Design

## Goal

Validate the research FMM/M2L prototype against TITAN's production TABI surface
and direct JCP operator on real PQR geometry without wiring FMM into GMRES or
changing any production default.

## Architecture

The Scope G probe reads a PQR fixture, builds panels through
`TitanPBE::Tabi::build_surface`, and deterministically samples at most a bounded
number of panels for the quadratic direct comparison. The sampled production
panels are converted to `PrototypePanel` values with deterministic `phi` and
`phin` unknowns.

The comparison has two independent layers:

1. Run production `apply_jcp_treecode_operator` in direct-fallback mode with
   zero LHS coefficients and negate its output to obtain the off-diagonal JCP
   integral reference.
2. Compare that reference to the prototype analytic K1-K4 direct formulas,
   then run the existing Scope E K1 and Scope F K2-K4 M2L gates on the same
   real-surface panel set.

This avoids validating the prototype only against itself. No panel-export file,
new production CLI, solver hook, or GMRES path is introduced.

## Inputs And Outputs

The PowerShell harness accepts `research-medium`, `research-large`, or an
explicit PQR path. It writes:

- `test/tmp/fmm_m2l_scope_g/fmm_m2l_scope_g_matrix.csv`
- `test/tmp/fmm_m2l_scope_g/fmm_m2l_scope_g_summary.dat`
- `TITAN_FMM_M2L_SCOPE_G_REPORT.txt`

The CSV records atom count, full surface panel count, sampled panel count,
production direct parity, K1 accuracy, derivative-channel accuracy, near/far
work, timings, and the final decision.

## Gate

- finite outputs;
- production direct parity aggregate relative RMS `<= 1e-10`;
- K1 relative RMS `<= 5e-2`;
- derivative aggregate relative RMS `<= 1e-1`;
- derivative maximum relative RMS `<= 2e-1` for channels whose reference RMS
  is at least `0.1%` of the strongest derivative channel;
- nonzero far work;
- combined M2L apply time lower than combined direct time.

`research-large` is gated behind a passing `research-medium` run. Passing
Scope G permits planning a hidden solver-integration experiment, but does not
permit changing `-pbe=grid` or exposing FMM as a production backend.

## Error Handling

Malformed PQR rows are rejected with the line number. Non-positive radii use
the existing research-fixture normalization policy of `0.6 A`. Degenerate or
empty surfaces fail clearly. Sampling is deterministic and preserves panel
order.

## Verification

The C++ probe first runs on a small generated production surface. The
PowerShell harness then runs `research-medium`; `research-large` runs only after
the medium gate passes. Native TITAN build and existing Scope A-F smoke remain
required regressions.
