# TITAN Guarded p4 FAGBI Paper Package

This package is a tracked, paper-facing bundle for the guarded p4 TABI treecode evidence.

## Scope

- Production PBE default remains `-pbe=grid`.
- TABI treecode default remains the degree-one engine.
- Guarded p4 remains an advanced explicit TABI treecode path.
- Z4Q supplement rows cover full public proteins 1UBQ, 2CI2, and 1VII.

## Contents

- `manuscript/`: refreshed Z4P DOCX/PDF render and source Markdown.
- `submission/`: Z4T submission-polished DOCX/PDF, locked references, BibTeX, claim trace, citation-gap closure, and page-render QA.
- `tables/`: Z4O snapshot tables, Z4P decisions, and full Z4Q supplement matrix/repeatability.
- `figures/`: supplement repeatability figure.
- `provenance/`: claim trace, DOCX QA, public fixture manifest, and package checksums.
- `reports/`: root Z4O/Z4P/Z4Q evidence reports copied at packaging time.
- `release_payload/`: manuscript-release payload with reference outputs for Tables 3-6 and Supplementary Tables S1-S6, APBS input/log/grid-setting files for Supplementary Figure S1, traceability tables, checksums, and build instructions for the tagged release `fagbi-submission-v0.1`.

## Public Release Preparation

The release manifest has been prepared for the MIT-licensed tagged release `fagbi-submission-v0.1` at `https://github.com/hoanglinh221191/TITAN-fagbi-submission`. The Zenodo DOI remains pending until the public tagged release is deposited. The top-level audit files are:

- `release_artifact_manifest.csv`
- `release_artifact_manifest.json`
- `release_readiness_summary.dat`
- `release_payload/README_RELEASE_ARCHIVE.md`
- `release_payload/table_figure_traceability.csv`
- `release_payload/checksums_sha256.csv`

## Decision

Package decision: `z4s_paper_package_ready`.

The package was generated from tracked evidence artifacts and keeps generated scratch files out of `test/tmp`.

## Submission Polish

Scope Z4T adds a journal-neutral full manuscript draft in `submission/` with eight locked references,
closed citation gaps, a rendered PDF, fifteen page PNGs for layout QA, and the decision
`z4t_submission_polish_ready`. It does not add new numerical evidence or change solver policy:
`-pbe=grid` remains the production default, TABI treecode remains degree-one by default, and guarded
p4 remains an advanced explicit path.
