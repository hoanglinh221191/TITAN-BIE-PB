# Five-test-system PBE preconditioner validation

This directory archives the manuscript-facing preconditioner cohort generated
with the hardened analytical patch-ray SES implementation. The cohort contains
five PBE test geometries, three surface representations, and six requested
preconditioners (90 rows total).

The two source matrices are:

- `data/grid_sphere_preconditioner_matrix.csv` for Grid-SES and the
  sphere-patch algebraic stress representation;
- `data/analytic_ses_preconditioner_matrix.csv` for the analytical patch-ray
  SES.

Every row records the executable SHA-256, solver-contract SHA-256, prepared-PQR
SHA-256, requested and effective preconditioner, residual history, true
residual, sampled operator audit, timing, and memory. A regularized ILU(0)
label denotes diagonal pivot regularization inside the requested zero-fill
factorization; it is not a fallback to diagonal preconditioning.

## Validate and regenerate the tables

From the repository root:

```powershell
python .\test\pbe_precond_cohort_n5_tables.py `
  --matrix .\validation\n5_preconditioner_cohort\data\grid_sphere_preconditioner_matrix.csv `
  --matrix .\validation\n5_preconditioner_cohort\data\analytic_ses_preconditioner_matrix.csv `
  --manifest .\test\pbe_inputs_public\public_fixture_manifest.csv `
  --out .\validation\n5_preconditioner_cohort\tables_n5_regenerated.md
```

For a full rerun after building `titan.exe`:

```powershell
.\test\pbe_precond_cohort_n5_gate.ps1
```

The full gate caches a row only when the executable, solver contract, numerical
protocol, and PQR hashes match. Use `-Force` only for a deliberate complete
rerun.

`SHA256SUMS.csv` covers the archived matrices and generated table. The matrix
rows were produced by executable SHA-256
`da4c4f5d2224cd20e0c617f46061ec5dabda2436aeadf58ba2da3601288f21da`.
