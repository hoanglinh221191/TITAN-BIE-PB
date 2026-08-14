# TITAN public source package

This staging directory contains the current TITAN C++ source, build files,
manual, license, citation metadata, and runnable examples.

- `src/`: solver and analysis source, including the BIE PB and MM/PBSA paths.
- `test/`: bounded N=5 preconditioner drivers and public prepared inputs.
- `validation/n5_preconditioner_cohort/`: two raw matrix CSV files, generated
  manuscript tables, validation instructions, and a SHA-256 manifest.
- `example/pbe/`: compact single-system solvation, surface-potential, and
  binding-cycle examples.
- `example/mm_pbsa_bie/`: complete trajectory-native BIE MM/PBSA example
  with reference outputs.
- `example/compact_mm_pbsa/`: two-model compact submission companion.

No executable, build object, temporary benchmark directory, or private raw
simulation archive is included.
