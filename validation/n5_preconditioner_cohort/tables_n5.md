# PBE preconditioner cohort: 5 test systems (BARNASE, 1A63, 1SH1, 1A7M, 1UBQ)

Table 3: PBE test systems and panel counts

| PBE test system | Geometry source | Atoms | Surface | Panels | Use |
|---|---|---|---|---|---|
| Barnase | 1BRS | 1,700 | Grid-SES | 25,548 | Independent surface-generator control |
| 1A63 Rho domain | 1A63 | 2,065 | Grid-SES | 36,688 | Independent surface-generator control |
| 1SH1 neurotoxin | 1SH1 | 696 | Grid-SES | 14,728 | Independent surface-generator control |
| 1A7M | 1A7M | 2,804 | Grid-SES | 39,794 | Independent surface-generator control |
| 1UBQ ubiquitin | 1UBQ | 1,231 | Grid-SES | 19,800 | Independent surface-generator control |
| Barnase | 1BRS | 1,700 | Analytic patch-ray SES | 25,422 | Primary physical protocol |
| 1A63 Rho domain | 1A63 | 2,065 | Analytic patch-ray SES | 36,486 | Primary physical protocol |
| 1SH1 neurotoxin | 1SH1 | 696 | Analytic patch-ray SES | 14,682 | Primary physical protocol |
| 1A7M | 1A7M | 2,804 | Analytic patch-ray SES | 39,538 | Primary physical protocol |
| 1UBQ ubiquitin | 1UBQ | 1,231 | Analytic patch-ray SES | 19,682 | Primary physical protocol |
| Barnase | 1BRS | 1,700 | Sphere patch | 14,942 | Algebraic stress test only |
| 1A63 Rho domain | 1A63 | 2,065 | Sphere patch | 17,250 | Algebraic stress test only |
| 1SH1 neurotoxin | 1SH1 | 696 | Sphere patch | 6,673 | Algebraic stress test only |
| 1A7M | 1A7M | 2,804 | Sphere patch | 21,042 | Algebraic stress test only |
| 1UBQ ubiquitin | 1UBQ | 1,231 | Sphere patch | 11,242 | Algebraic stress test only |

Table 4: GMRES iterations for the six preconditioners (a value of 720 denotes exhaustion of the iteration budget)

| PBE test system | Surface | Panels | Diag. | Leaf block | FMM RAS | kNN RAS | kNN ILU | FMM ILU |
|---|---|---|---|---|---|---|---|---|
| Barnase | Watertight Grid-SES | 25,548 | 720 | 33 | 11 | 11 | 10 | 10 |
| 1A63 | Watertight Grid-SES | 36,688 | 720 | 76 | 13 | 13 | 11 | 11 |
| 1SH1 | Watertight Grid-SES | 14,728 | 720 | 34 | 14 | 14 | 10 | 10 |
| 1A7M | Watertight Grid-SES | 39,794 | 720 | 53 | 23 | 23 | 13 | 13 |
| 1UBQ | Watertight Grid-SES | 19,800 | 720 | 31 | 12 | 12 | 9 | 9 |
| Barnase | Analytic patch-ray SES | 25,422 | 512 | 37 | 17 | 17 | 11 | 11 |
| 1A63 | Analytic patch-ray SES | 36,486 | 720 | 68 | 13 | 13 | 11 | 11 |
| 1SH1 | Analytic patch-ray SES | 14,682 | 720 | 25 | 16 | 16 | 10 | 10 |
| 1A7M | Analytic patch-ray SES | 39,538 | 720 | 56 | 27 | 27 | 15 | 15 |
| 1UBQ | Analytic patch-ray SES | 19,682 | 720 | 25 | 15 | 15 | 9 | 9 |
| Barnase | Sphere-patch stress | 14,942 | 720 | 720 | 720 | 720 | 34 | 34 |
| 1A63 | Sphere-patch stress | 17,250 | 720 | 720 | 720 | 720 | 48 | 48 |
| 1SH1 | Sphere-patch stress | 6,673 | 720 | 720 | 720 | 720 | 15 | 15 |
| 1A7M | Sphere-patch stress | 21,042 | 720 | 720 | 720 | 720 | 51 | 51 |
| 1UBQ | Sphere-patch stress | 11,242 | 720 | 720 | 720 | 720 | 41 | 41 |

Table 5: bounded graph and factor diagnostics for the two ILU(0) routes

| PBE test system | Surface | Graph | Panels | Blocks/N | Setup (s) | Apply (s) | Iter. | RSS (MB) |
|---|---|---|---|---|---|---|---|---|
| Barnase | Watertight Grid-SES | distance-kNN | 25548 | 64 | 2.86 | 0.19 | 10 | 178.6 |
| Barnase | Watertight Grid-SES | FMM-near | 25548 | 64 | 5.43 | 0.19 | 10 | 178.9 |
| 1A63 | Watertight Grid-SES | distance-kNN | 36688 | 64 | 6.12 | 0.31 | 11 | 266.9 |
| 1A63 | Watertight Grid-SES | FMM-near | 36688 | 64 | 8.11 | 0.30 | 11 | 266.7 |
| 1SH1 | Watertight Grid-SES | distance-kNN | 14728 | 64 | 1.58 | 0.10 | 10 | 105.4 |
| 1SH1 | Watertight Grid-SES | FMM-near | 14728 | 64 | 2.78 | 0.11 | 10 | 105.9 |
| 1A7M | Watertight Grid-SES | distance-kNN | 39794 | 64 | 8.31 | 0.43 | 13 | 284.0 |
| 1A7M | Watertight Grid-SES | FMM-near | 39794 | 64 | 11.71 | 0.40 | 13 | 284.9 |
| 1UBQ | Watertight Grid-SES | distance-kNN | 19800 | 64 | 2.35 | 0.18 | 9 | 146.4 |
| 1UBQ | Watertight Grid-SES | FMM-near | 19800 | 64 | 4.66 | 0.17 | 9 | 146.9 |
| Barnase | Analytic patch-ray SES | distance-kNN | 25422 | 64 | 3.71 | 0.22 | 11 | 179.6 |
| Barnase | Analytic patch-ray SES | FMM-near | 25422 | 64 | 5.43 | 0.22 | 11 | 179.9 |
| 1A63 | Analytic patch-ray SES | distance-kNN | 36486 | 64 | 5.31 | 0.30 | 11 | 266.9 |
| 1A63 | Analytic patch-ray SES | FMM-near | 36486 | 64 | 4.64 | 0.23 | 11 | 267.2 |
| 1SH1 | Analytic patch-ray SES | distance-kNN | 14682 | 64 | 1.18 | 0.09 | 10 | 106.4 |
| 1SH1 | Analytic patch-ray SES | FMM-near | 14682 | 64 | 2.63 | 0.10 | 10 | 106.4 |
| 1A7M | Analytic patch-ray SES | distance-kNN | 39538 | 64 | 6.58 | 0.37 | 15 | 284.2 |
| 1A7M | Analytic patch-ray SES | FMM-near | 39538 | 64 | 10.67 | 0.47 | 15 | 283.9 |
| 1UBQ | Analytic patch-ray SES | distance-kNN | 19682 | 64 | 2.39 | 0.14 | 9 | 147.3 |
| 1UBQ | Analytic patch-ray SES | FMM-near | 19682 | 64 | 4.58 | 0.17 | 9 | 147.1 |
| Barnase | Sphere-patch stress | distance-kNN | 14942 | 64 | 1.52 | 0.26 | 34 | 119.9 |
| Barnase | Sphere-patch stress | FMM-near | 14942 | 64 | 2.97 | 0.26 | 34 | 120.3 |
| 1A63 | Sphere-patch stress | distance-kNN | 17250 | 64 | 1.52 | 0.49 | 48 | 158.6 |
| 1A63 | Sphere-patch stress | FMM-near | 17250 | 64 | 1.90 | 0.53 | 48 | 158.8 |
| 1SH1 | Sphere-patch stress | distance-kNN | 6673 | 64 | 0.40 | 0.06 | 15 | 60.9 |
| 1SH1 | Sphere-patch stress | FMM-near | 6673 | 64 | 1.02 | 0.07 | 15 | 61.1 |
| 1A7M | Sphere-patch stress | distance-kNN | 21042 | 64 | 4.34 | 0.80 | 51 | 179.1 |
| 1A7M | Sphere-patch stress | FMM-near | 21042 | 64 | 6.09 | 0.76 | 51 | 179.7 |
| 1UBQ | Sphere-patch stress | distance-kNN | 11242 | 64 | 1.02 | 0.31 | 41 | 99.7 |
| 1UBQ | Sphere-patch stress | FMM-near | 11242 | 64 | 0.95 | 0.29 | 41 | 99.8 |

Derived quantities for the Abstract, 3.2 and 3.3

- cohort size: 5 PBE test systems x 3 surface representations
- graph ILU(0) iteration range (all surfaces): 9-51
- panel range: 6,673-39,794
- normalized memory: 7.0-9.2 MB per 1,000 panels
- sampled operator audit RMS, ILU rows (physical surfaces): 3.501e-03 to 7.882e-03
- true residual, ILU rows (physical surfaces): 7.374e-06 to 9.129e-05
- sampled operator audit RMS, ILU rows (sphere-patch stress): 6.823e-03 to 1.002e+00
- true residual, ILU rows (sphere-patch stress): 2.225e-05 to 6.871e-05

- Watertight Grid-SES: ILU(0) 9-13 iterations; non-factorized controls 11-720 (5/20 exhausted the 720-iteration budget)
- Analytic patch-ray SES: ILU(0) 9-15 iterations; non-factorized controls 13-720 (4/20 exhausted the 720-iteration budget)
- Sphere-patch stress: ILU(0) 15-51 iterations; non-factorized controls 720-720 (20/20 exhausted the 720-iteration budget)

Consistency checks

- non-factorized control rows that exhausted the budget (expected, not an error): 29
- ILU rows with status != ok: 0
- ILU rows failing the preconditioner quality gate: 0
- kNN vs FMM-near ILU iteration mismatches: none (both routes agree on every row)

