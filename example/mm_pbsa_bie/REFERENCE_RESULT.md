# Reference full-trajectory result

The bundled protein-ligand trajectory contains 52 coordinate models spanning
89.8-100.0 ns. TITAN uses the first model as the coordinate/topology seed and
evaluates 51 PBC-corrected frames from 90.0 to 100.0 ns with four CPU threads
and the settings in `arg_bie_mm_pbsa.in`.

All 153 BIE endpoint solves (complex, receptor, and ligand for every frame)
completed successfully. The complex solve used sparse near-field ILU(0) in
all 51 profiles, without fallback. The sampled accelerated-versus-direct
operator discrepancy was at most `3.38e-3`, and the final true residual was at
most `9.93e-5`.

The trajectory-mean energy cycle was:

| Component | Mean (kcal/mol) | Standard deviation (kcal/mol) |
|---|---:|---:|
| Molecular-mechanics electrostatics | -24.6874 | 5.1356 |
| Molecular-mechanics van der Waals | -35.4776 | 2.9098 |
| Molecular mechanics | -60.1650 | 7.2341 |
| Polar solvation difference | 56.6212 | 17.7772 |
| Nonpolar SASA difference | -6.1280 | 0.1604 |
| Total, without entropy | -9.6717 | 19.1605 |
| Interaction-entropy correction, `-T Delta S` | 12.0562 | n/a |
| Final total, with interaction entropy | 2.3846 | n/a |

The no-entropy MM/PBSA estimate was negative for 39 of 51 frames and positive
for 12 frames. Its range was -113.4759 to 29.2355 kcal/mol. After adding the
single trajectory-level interaction-entropy correction, 16 frame-associated
values were negative and 35 were positive. Thus the PBC-corrected trajectory
is favorable at the mean enthalpic MM/PBSA level and nearly neutral, but mildly
unfavorable, after the interaction-entropy estimate.

This differs from the historical MM/GBSA result for the same trajectory
(`-48.3711` kcal/mol without entropy; `-43.1954` kcal/mol with interaction
entropy). The molecular-mechanics components agree closely between runs:
electrostatics is `-24.6874` versus `-24.6861` kcal/mol, van der Waals is
`-35.4776` versus `-35.4776` kcal/mol, and the nonpolar term is `-6.1280`
versus `-6.1231` kcal/mol. The main change is the PB polar-solvation penalty
(`+56.6212` kcal/mol) versus the GB value (`+17.9157` kcal/mol), followed by a
larger interaction-entropy correction. The different final sign therefore
comes from the continuum-solvation and entropy models, not from
periodic-boundary handling or the ligand ILU cap.

The three BIE endpoint solves required 69.61 s per analyzed frame on average.
The ligand endpoint averaged only 0.26 s, so retaining the common ILU neighbor
cap of 64 has negligible practical cost for this small ligand. Peak reported
memory was 204.9 MiB. Hybrid decomposition produced 434 protein-residue and
ligand-atom units, and the maximum additive closure error was `1e-10` kJ/mol.

These values verify the executable workflow and additive decomposition over
the complete sampled trajectory. They are not an experimentally calibrated
binding-affinity prediction.
