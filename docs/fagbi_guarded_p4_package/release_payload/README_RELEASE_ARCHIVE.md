# TITAN FAGBI manuscript release payload

This directory prepares the manuscript-associated release payload for the tagged repository archive fagbi-submission-v0.1.

## Code and data availability statement

The TITAN source code corresponding to this manuscript has been prepared for release under the MIT License as the tagged release fagbi-submission-v0.1 at https://github.com/hoanglinh221191/TITAN-fagbi-submission. The persistent Zenodo DOI will be inserted after deposition of the tagged release and before journal submission. The archived release will include source code, build instructions, benchmark manifests, input metadata, scripts, figures, APBS input files, APBS logs, APBS grid-size settings used for Supplementary Figure S1, and reference outputs for Tables 3-6 and Supplementary Tables S1-S6. Public protein structures are identified by their RCSB PDB accession codes; derived surface metadata, run manifests, and tabulated outputs corresponding to the manuscript tables are included in the release manifest.

## Payload contents

- reference_outputs/: CSV outputs copied from the verified local benchmark runs for Tables 3-6 and Supplementary Tables S1-S6.
- apbs_context/: APBS input files, APBS stdout/stderr logs, APBS grid-size settings, and the PQR inputs used for Supplementary Figure S1.
- table_figure_traceability.csv: table/figure to source-output/script mapping.
- reference_output_copy_manifest.csv: source-to-archive copy manifest for tabulated reference outputs.
- checksums_sha256.csv: SHA256 checksums for release-facing source, scripts, figures, inputs, and outputs.

The release manifest at ../release_artifact_manifest.csv is the top-level audit table.
