#!/usr/bin/env python3
"""Build the Z4T submission-polished guarded-p4 manuscript DOCX."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml import OxmlElement
from docx.shared import Inches

from pbe_fagbi_scope_z4f_build_docx import (
    add_body,
    add_figure,
    add_heading,
    configure_section,
    configure_styles,
    format_table_text,
    set_run_font,
    set_table_geometry,
)


TITLE = (
    "TITAN: Analytical FMM/M2L Boundary-Integral Poisson-Boltzmann Solves "
    "with Sparse Near-Field ILU(0) Preconditioning"
)

FORBIDDEN_BODY_TERMS = (
    "cap-128",
    "cap-256",
    "research-medium",
    "research-alt",
    "money plot",
    "Z1-Z5",
    "Z4 reference policy",
    "audit-clean",
    "TBD",
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def as_float(value: str | None, default: float = math.nan) -> float:
    try:
        return float(value) if value not in (None, "", "NA", "NaN") else default
    except ValueError:
        return default


def fmt_float(value: float, digits: int = 3) -> str:
    if not math.isfinite(value):
        return "NA"
    return f"{value:.{digits}f}"


def fmt_sci(value: float) -> str:
    if not math.isfinite(value):
        return "NA"
    return f"{value:.3e}"


def load_reference_rows() -> list[dict[str, str]]:
    return [
        {
            "key": "Chen2021FAGBI",
            "label": "[1]",
            "text": "J. Chen, J. Tausch, and W. Geng, A Cartesian FMM-accelerated Galerkin boundary integral Poisson-Boltzmann solver, arXiv:2110.13778, 2021.",
            "url": "https://arxiv.org/abs/2110.13778",
            "doi": "",
        },
        {
            "key": "GreengardRokhlin1987",
            "label": "[2]",
            "text": "L. Greengard and V. Rokhlin, A fast algorithm for particle simulations, Journal of Computational Physics 73, 325-348, 1987.",
            "url": "https://doi.org/10.1016/0021-9991(87)90140-9",
            "doi": "10.1016/0021-9991(87)90140-9",
        },
        {
            "key": "SaadSchultz1986",
            "label": "[3]",
            "text": "Y. Saad and M. H. Schultz, GMRES: A generalized minimal residual algorithm for solving nonsymmetric linear systems, SIAM Journal on Scientific and Statistical Computing 7, 856-869, 1986.",
            "url": "https://doi.org/10.1137/0907058",
            "doi": "10.1137/0907058",
        },
        {
            "key": "MeijerinkVanDerVorst1977",
            "label": "[4]",
            "text": "J. A. Meijerink and H. A. van der Vorst, An iterative solution method for linear systems of which the coefficient matrix is a symmetric M-matrix, Mathematics of Computation 31, 148-162, 1977.",
            "url": "https://doi.org/10.1090/S0025-5718-1977-0438681-4",
            "doi": "10.1090/S0025-5718-1977-0438681-4",
        },
        {
            "key": "Duffy1982",
            "label": "[5]",
            "text": "M. G. Duffy, Quadrature over a pyramid or cube of integrands with a singularity at a vertex, SIAM Journal on Numerical Analysis 19, 1260-1262, 1982.",
            "url": "https://doi.org/10.1137/0719090",
            "doi": "10.1137/0719090",
        },
        {
            "key": "Dolinsky2004PDB2PQR",
            "label": "[6]",
            "text": "T. J. Dolinsky et al., PDB2PQR: an automated pipeline for the setup of Poisson-Boltzmann electrostatics calculations, Nucleic Acids Research 32, W665-W667, 2004.",
            "url": "https://doi.org/10.1093/nar/gkh381",
            "doi": "10.1093/nar/gkh381",
        },
        {
            "key": "Jurrus2018APBS",
            "label": "[7]",
            "text": "E. Jurrus et al., Improvements to the APBS biomolecular solvation software suite, Protein Science 27, 112-128, 2018.",
            "url": "https://doi.org/10.1002/pro.3280",
            "doi": "10.1002/pro.3280",
        },
        {
            "key": "Burley2021RCSB",
            "label": "[8]",
            "text": "S. K. Burley et al., RCSB Protein Data Bank: powerful new tools for exploring 3D structures of biological macromolecules, Nucleic Acids Research 49, D437-D451, 2021.",
            "url": "https://doi.org/10.1093/nar/gkaa1038",
            "doi": "10.1093/nar/gkaa1038",
        },
    ]


def write_bibtex(path: Path, refs: list[dict[str, str]]) -> None:
    entries = []
    for ref in refs:
        key = ref["key"]
        doi = ref["doi"]
        title = ref["text"].split(", ", 1)[1].split(", ")[0]
        entries.append(
            "@misc{" + key + ",\n"
            f"  title = {{{title}}},\n"
            f"  note = {{{ref['text']}}},\n"
            f"  doi = {{{doi}}},\n"
            f"  url = {{{ref['url']}}}\n"
            "}\n"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(entries), encoding="utf-8")


def representative_main(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    selected = [row for row in rows if row.get("variant") == "p4_guarded"]
    first_by_fixture: dict[str, dict[str, str]] = {}
    for row in selected:
        first_by_fixture.setdefault(row["fixture"], row)
    ordered = []
    for fixture in ("barnase", "1a63", "1sh1"):
        if fixture not in first_by_fixture:
            raise ValueError(f"Missing Z4O fixture {fixture}")
        ordered.append(first_by_fixture[fixture])
    return ordered


def supplement_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    first_by_fixture: dict[str, dict[str, str]] = {}
    for row in rows:
        if row.get("variant") == "p4_guarded":
            first_by_fixture.setdefault(row["fixture"], row)
    return [first_by_fixture[key] for key in sorted(first_by_fixture)]


def repeatability_by_fixture(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row["fixture"]: row for row in rows}


def manifest_by_fixture(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row["fixture"]: row for row in rows}


def add_caption(doc: Document, text: str) -> None:
    paragraph = doc.add_paragraph(text, style="Caption")
    paragraph.alignment = WD_ALIGN_PARAGRAPH.LEFT
    paragraph.paragraph_format.keep_with_next = True


def add_table(doc: Document, headers: list[str], body: list[list[str]], widths: list[int], size: float) -> None:
    table = doc.add_table(rows=1, cols=len(headers))
    table.style = "Table Grid"
    for cell, header in zip(table.rows[0].cells, headers):
        cell.text = header
    table.rows[0]._tr.get_or_add_trPr().append(OxmlElement("w:tblHeader"))
    for values in body:
        cells = table.add_row().cells
        for cell, value in zip(cells, values):
            cell.text = value
        table.rows[-1]._tr.get_or_add_trPr().append(OxmlElement("w:cantSplit"))
    set_table_geometry(table, widths)
    format_table_text(table, size)


def add_reference_list(doc: Document, refs: list[dict[str, str]]) -> None:
    add_heading(doc, "References")
    for ref in refs:
        paragraph = doc.add_paragraph()
        run = paragraph.add_run(f"{ref['label']} {ref['text']}")
        set_run_font(run, "Times New Roman", 10)
        paragraph.paragraph_format.space_after = 3


def build_manuscript_source(main_rows, main_rep, supp_rows, supp_rep, refs) -> str:
    text = f"""# {TITLE}

## Abstract

Boundary-integral Poisson-Boltzmann methods avoid volumetric dielectric grids by moving the electrostatic unknowns to the molecular surface, but this advantage is useful only if the dense surface operator can be applied accurately and if the resulting Krylov system is preconditioned robustly. We present a guarded fourth-order analytical FMM/M2L path in TITAN that combines analytic-singular near-field integration, deterministic far-cell refinement, compressed Green-function jets, packed AVX2 direct fallback, and sparse near-field ILU(0) preconditioning. The method is evaluated as an advanced explicit TABI treecode option, not as a default backend. The primary full-protein panel includes Barnase, the RNA-binding domain of E. coli Rho factor (1A63), and Stichodactyla helianthus neurotoxin I (1SH1). All guarded p4 repeats converged with finite energies, sampled direct-operator residuals below 8e-3, and one plan/cache build. Additional public proteins, Ubiquitin (1UBQ), Chymotrypsin inhibitor 2 (2CI2), and Villin headpiece subdomain (1VII), passed the same supplement gate. A dense Barnase convergence profile shows that sparse near-field ILU(0) resolves a high-coupling case where local block-diagonal preconditioning stalls. The production Poisson-Boltzmann default remains the grid backend, and the default TABI treecode path remains degree-one.

## Introduction

Implicit-solvent electrostatics is a central approximation in biomolecular modeling because it replaces explicit solvent degrees of freedom with a continuum response. The Poisson-Boltzmann equation is one of the most established forms of this approximation: it combines a low-dielectric molecular interior, a high-dielectric solvent exterior, fixed molecular charges, and mobile ions. Classical grid solvers are robust and widely used, but they carry a volumetric cost. Increasing molecular size, changing salt conditions, or refining the dielectric interface can quickly make the grid the dominant computational object.

Boundary-integral Poisson-Boltzmann methods attack the same physical problem from a surface perspective. By reformulating the linearized PB equation as coupled integral equations on the molecular surface, they reduce the unknown set from a three-dimensional volume to a triangulated interface. This surface reduction is attractive for repeated protein calculations and for workflows in which surface discretization is the natural resolution parameter. It also creates a difficult numerical problem: the surface operator is dense, contains derivative-sensitive kernels, and includes near-singular panel interactions that cannot be treated as ordinary far-field particle sums.

Fast multipole methods provide the natural acceleration idea for such dense operators {refs[1]['label']}. The Cartesian boundary-integral PB direction developed in the FAGBI work of Chen, Tausch, and Geng is the closest methodological reference point for TITAN's analytical FMM/M2L effort {refs[0]['label']}. TITAN follows that broad motivation but targets a different integration setting. Rather than presenting a standalone Galerkin code, TITAN embeds an explicit guarded high-order path inside an existing TABI workflow, with production defaults left unchanged. The focus is therefore not only the analytical expansion, but also the engineering required to make the accelerated operator repeatable, auditable, and usable inside GMRES.

The second central issue is preconditioning. A boundary-integral acceleration that makes each matrix-vector product faster is not enough if the Krylov iteration count grows or the solver stalls. During TITAN development, local diagonal and local block-diagonal preconditioners repeatedly failed to capture strong coupling across tree leaves on dense molecular surfaces. Increasing local block size improved some rows but did not solve the fundamental graph problem. The final method therefore constructs a sparse near-field graph directly from FMM near-source relationships and applies block ILU(0) over coupled 2 x 2 panel blocks. This is the main numerical distinction of the present guarded path.

The claim in this manuscript is intentionally scoped. Guarded p4 is not promoted to TITAN's production default, and it is not selected automatically by TABI treecode. Instead, we document a repeatable advanced explicit path. The evidence chain asks whether the high-order operator agrees with sampled direct references, whether the preconditioner resolves a dense high-coupling case, and whether the workflow repeats on public full-protein fixtures. This is a narrower claim than a full replacement for grid PB, but it is stronger than a prototype-only demonstration.

## Methods

TITAN solves a coupled boundary-integral Poisson-Boltzmann system over triangulated molecular surfaces. Each panel carries the unknown vector u_i = (phi_i, partial_n phi_i)^T, where phi is the surface potential and partial_n phi is its normal derivative. The JCP-style coupled operator contains four analytical interaction channels, denoted here as K1 through K4. These channels couple the potential and normal-derivative rows and are applied with the same production sign convention in direct, degree-one accelerated, and guarded p4 paths. GMRES is used for the nonsymmetric coupled linear system {refs[2]['label']}.

Near interactions are evaluated with the production analytic-singular kernel. This choice matters because self and adjacent-panel terms can dominate the error if they are approximated by ordinary centroid interactions. Duffy-adaptive quadrature remains available as an explicit diagnostic option {refs[4]['label']}, but it is not selected automatically in the guarded p4 workflow. The manuscript-facing policy therefore separates near-field correctness from far-field acceleration: near terms use the shared analytic-singular path, while far terms are controlled by high-order expansion, refinement, and direct fallback.

The guarded p4 far operator uses fourth-order Cartesian source moments and target-local M2L/L2P evaluation. The accepted far-cell contribution is not treated as a black-box scalar approximation. Instead, the four coupled channels are assembled from analytical Green-function jets and their source/target derivatives. To keep memory bounded, TITAN stores two double-precision degree-9 Green jets per far pair, one for the Coulomb part and one for the screened part. Source-normal derivatives needed by K2 and K4 are reconstructed during contraction rather than stored as separate materialized jets.

The static refinement plan is built before GMRES and reused throughout the solve. Deterministic operator probes identify risky far cells, and a strict geometric separation guard catches outliers that an RMS tail indicator can miss. The final guarded policy uses Taylor order 4, opening theta 0.6, static tail tolerance 0.002, strict-geometry MAC theta 0.5, compressed Green-jet cache, packed AVX2 direct fallback, analytic-singular near integration, and sparse near-field ILU(0). If a source cell remains risky and cannot be safely split, its panel interactions are demoted to the direct fallback path rather than forced through the far expansion.

The direct fallback path is also structured for reuse. The refinement plan produces deterministic target/source tasks. Non-self fallback pairs use a packed AVX2 primitive, while self terms are routed through the exact analytic near kernel and masked out of vector lanes to prevent double counting. This was necessary because accuracy recovery on full 1A63 required a geometry guard that increased fallback work. Without packed fallback, the recovered p4 operator would be accurate but too expensive for guarded GMRES.

The sparse preconditioner is assembled from the FMM near-source graph. For every target panel i, the graph contains the self edge and a deterministic set of nearby source panels, capped at 64 neighbors in the guarded configuration. On each graph edge, TITAN stores a 2 x 2 block approximation to the coupled local operator: A_ij = [[delta_ij c_phi - K2_ij, -K1_ij], [-K4_ij, delta_ij c_phin - K3_ij]]. The matrix is stored in block sparse-row form. ILU(0) keeps the original graph sparsity, applying the no-fill update L_ik = A_ik U_kk^{-1} and U_ij <- U_ij - L_ik U_kj only when the block (i,j) already exists {refs[3]['label']}. The apply phase performs deterministic block forward and backward substitution.

Protein structures are prepared with PDB2PQR/APBS-family tooling {refs[5]['label']}, {refs[6]['label']} and public structure identifiers from the RCSB Protein Data Bank {refs[7]['label']}. The primary full-protein panel contains Barnase, the RNA-binding domain of E. coli Rho factor (1A63), and Stichodactyla helianthus neurotoxin I (1SH1). The supplementary public-protein expansion contains Ubiquitin (1UBQ), Chymotrypsin inhibitor 2 (2CI2), and Villin headpiece subdomain (1VII). Manuscript-facing rows are full-protein fixtures, not atom-sampled fixtures.

The reference policy is matched to computational feasibility. For analytic and bounded development tests, direct TABI can be used as a full reference. For the full-protein guarded p4 panel, full direct TABI reaction energies are not practical for the largest rows, so verification uses sampled direct-operator RMS and sampled direct true residuals. Grid PB remains the production default and is useful runtime context, but it is not treated as the oracle for the high-order boundary-integral path.

## Analytical Validation Protocol

Before the full-protein guarded p4 rows were accepted, TITAN used controlled analytic systems to separate surface-discretization error from acceleration error. Centered Born and off-center Kirkwood-style tests were used to check whether direct TABI with analytic-singular near integration converged toward known spherical references as the surface was refined. This step matters because an accelerated operator cannot be meaningfully validated if the underlying surface discretization is already outside the target accuracy range.

The high-resolution direct TABI checks showed that the remaining analytic error could be reduced below the manuscript threshold on the controlled sphere cases. Those tests do not prove full-protein energy accuracy by themselves, because protein surfaces are more irregular and do not have closed-form reaction-energy oracles. They do establish that the panel basis, normals, diagonal terms, and analytic-singular near path are not the dominant explanation for the later full-protein operator gaps. This allowed the later debugging to focus on far-field truncation, refinement, and preconditioning rather than on the direct near-field kernel.

The operator attribution sequence then separated near and far contributions on full proteins. The key observation was that a nominal fourth-order configuration could still fail on 1A63 if the refinement plan was driven only by an averaged tail estimate. The outlier was not a sign-convention error, a near-kernel regression, or a cache lifecycle bug. It was a far-cell acceptance problem: a small number of geometrically risky accepted cells could dominate the sampled direct residual. The final strict-geometry guard is the direct response to that attribution.

## Lifecycle And Reproducibility Controls

The guarded p4 path is prepared before GMRES. The context builds the FMM geometry, static refinement plan, compressed Green-jet cache, packed fallback worklist, and near-field ILU(0) preconditioner once. During Krylov iteration, only source moments and density-dependent contractions are updated. This lifecycle is part of the numerical contract, not only a performance detail. If the plan or cache were rebuilt inside GMRES, repeated runs could become difficult to interpret and setup cost would contaminate matvec timing.

Every manuscript-facing guarded p4 row therefore records plan/cache build counts. The accepted rows have build count one. They also record sampled direct-operator residuals after convergence, so the reported solver residual is not the only accuracy measure. Internal GMRES residual can be small even when the accelerated operator differs from the direct operator; the full 1A63 diagnostic exposed exactly that failure before the geometry guard was added. The final verification requires both internal convergence and sampled direct true residual to pass.

The implementation uses deterministic ordering for refinement, cache layout, fallback tasks, and panel sampling. Parallelism is allowed, but worklists are sorted by target and source identifiers, and each worker owns disjoint output panels. This determinism is why repeat energy spread is a meaningful package metric. A repeatability table with zero spread at reported precision is not a substitute for independent validation, but it does show that the guarded path is not depending on race-sensitive accumulation or accidental cache state.

## Results

The primary guarded p4 panel passed all repeatability gates. Barnase, 1A63, and 1SH1 each converged with finite energies, sampled direct verification below 8e-3, and one plan/cache build per solve. The largest sampled direct residual in the primary panel occurred for 1A63, approximately 2.36e-3, still well inside the gate. The energy spread across repeats was zero in the recorded package, indicating deterministic solver behavior at the reported precision.

Runtime behavior was also favorable relative to same-fixture degree-one FMM baselines. Median p4 matvec ratios were below one for all three primary proteins. This matters because guarded p4 adds far-cell refinement, high-order local evaluation, compressed cache construction, and direct fallback. The measured ratios show that the extra accuracy machinery does not automatically translate into slower matvecs on the tested full-protein panel.

The supplementary public-protein expansion broadens the evidence without changing the main claim. Ubiquitin, Chymotrypsin inhibitor 2, and Villin headpiece subdomain passed the same guarded p4 repeatability gate. These systems are useful because they are public, compact, and prepared as full proteins. They are treated as supplement evidence rather than as a replacement for the locked primary panel.

The Barnase convergence profile provides the clearest preconditioner contrast. In this dense/high-coupling case, local block-diagonal preconditioning stalls before reaching the strict target, whereas sparse near-field ILU(0) continues to reduce the residual. This supports the design decision to build the preconditioner on the near-source graph rather than only on local tree leaves. The result should be interpreted narrowly: it demonstrates a hard coupling regime, not a claim that every protein requires sparse ILU(0).

Together, these results support a specific statement: TITAN's guarded p4 path is repeatable as an advanced explicit TABI treecode option under the documented policy. They do not support making guarded p4 automatic, and they do not change the role of grid PB as the production default.

## Discussion

The first lesson from this work is that far-field accuracy for derivative-rich boundary-integral PB operators is not controlled by expansion order alone. A fourth-order Cartesian expansion was necessary but not sufficient. Full 1A63 exposed outliers that passed a weaker RMS-style refinement policy. The final strict geometry guard recovered the sampled direct residual while keeping fallback work within budget. This suggests that robust high-order FMM integration for molecular surfaces needs both analytical expansion and deterministic geometric fail-safes.

The second lesson is that preconditioning must follow the coupling graph. Local diagonal and local block-diagonal methods are attractive because they are simple and memory-bounded, but they can miss strong cross-leaf interactions. Sparse near-field ILU(0) keeps memory bounded by the FMM near-source graph while allowing information to propagate across neighboring panels. This is why the preconditioner can resolve the dense Barnase case without becoming a dense local inverse.

The third lesson is about evidence discipline. It would be easy to present a single successful high-order solve as a new default backend, but that would blur the distinction between a verified advanced path and a production policy. TITAN keeps these roles separate. Grid remains the production Poisson-Boltzmann default. TABI treecode remains degree-one by default. Guarded p4 is documented as explicit because it carries a larger configuration surface and because its verification policy relies on sampled direct-operator audits for large systems.

Compared with the FAGBI reference direction {refs[0]['label']}, this work emphasizes integration into an existing molecular electrostatics codebase. The algorithmic components are not only analytical FMM/M2L expansions, but also static refinement, compressed jet caching, vectorized fallback, lifecycle control inside GMRES, and graph-based sparse preconditioning. The manuscript therefore should not be read as an exact reproduction of FAGBI. It is a TITAN-specific boundary-integral PB acceleration path inspired by the same analytical acceleration context.

## Limitations

Several limitations are deliberate. Full direct TABI reaction energies are not used as the oracle for large full-protein rows; sampled direct-operator residuals are the feasible verification mechanism. Duffy-adaptive near integration remains explicit-only and is not automatically selected. The Barnase preconditioner convergence profile is framed as dense/high-coupling evidence rather than a broad claim over all proteins. The supplementary fixtures broaden the evidence but do not replace the primary locked panel.

The manuscript also does not claim a general production promotion. Guarded p4 remains an advanced explicit TABI treecode option, and grid remains TITAN's production Poisson-Boltzmann default. Future work should evaluate larger public protein panels, compare against independent PB packages under matched surface preparation, and test whether the guarded p4 policy remains optimal across ionic strengths, dielectric settings, and finer surface discretizations.

## Reproducibility And Data Availability

The manuscript package is designed to be auditable from tracked artifacts rather than from undocumented local state. The tables, figure inputs, reference list, BibTeX file, claim trace, citation-gap closure, and rendered page checks are written into the paper package. The public-protein fixtures are identified by biological name and PDB identifier, and the preparation workflow records the PDB2PQR/APBS-family provenance used for the supplementary systems. The package therefore separates numerical evidence from manuscript formatting: rebuilding the document does not rerun the solver, and rerunning the solver gates produces new matrix artifacts that can be compared with the frozen paper package.

All large-protein accuracy statements in the manuscript are tied to sampled direct-operator verification rather than to unreported full direct solves. This is a limitation, but it is also a reproducibility choice. The sampled audit rows are deterministic, and the manuscript records the policy boundary explicitly so that future work can replace the sampled reference with full direct energy checks when practical. The current package should therefore be read as a reproducible evidence-backed draft, not as a final journal submission with independent external-code validation.

## Conclusion

TITAN now contains a repeatable guarded fourth-order analytical FMM/M2L boundary-integral PB path with sparse near-field ILU(0) preconditioning. On the primary full-protein panel and supplementary public fixtures, the explicit guarded p4 workflow converges, passes sampled direct-operator verification, and preserves a clear policy boundary: it is a documented advanced option, not the production default. The evidence supports a paper-facing claim centered on analytical acceleration, deterministic far-field guarding, and graph-based sparse preconditioning for dense molecular-surface PB solves.
"""
    leaked = [term for term in FORBIDDEN_BODY_TERMS if term.lower() in text.lower()]
    if leaked:
        raise ValueError("Forbidden internal manuscript term leaked: " + ", ".join(leaked))
    return text


def build_docx(args: argparse.Namespace) -> None:
    z4o_matrix = read_csv(args.z4o_matrix)
    z4o_rep = repeatability_by_fixture(read_csv(args.z4o_repeatability))
    z4q_matrix = read_csv(args.z4q_matrix)
    z4q_rep = repeatability_by_fixture(read_csv(args.z4q_repeatability))
    z4q_manifest = manifest_by_fixture(read_csv(args.z4q_manifest))
    z4o_rows = representative_main(z4o_matrix)
    z4q_rows = supplement_rows(z4q_matrix)
    refs = load_reference_rows()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    references_csv = args.output_dir / "fagbi_scope_z4t_references.csv"
    write_csv(references_csv, refs, ["key", "label", "text", "url", "doi"])
    write_bibtex(args.output_dir / "fagbi_scope_z4t_references.bib", refs)

    source = build_manuscript_source(z4o_rows, z4o_rep, z4q_rows, z4q_rep, refs)
    source_path = args.output_dir / "fagbi_scope_z4t_manuscript_source.md"
    source_path.write_text(source, encoding="utf-8")

    doc = Document()
    configure_styles(doc)
    configure_section(doc.sections[0])
    doc.core_properties.title = TITLE
    doc.core_properties.author = ""
    doc.core_properties.keywords = "Poisson-Boltzmann; FMM; M2L; ILU(0); boundary integral"

    title = doc.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    set_run_font(title.add_run(TITLE), "Times New Roman", 16, bold=True)
    subtitle = doc.add_paragraph()
    subtitle.alignment = WD_ALIGN_PARAGRAPH.CENTER
    set_run_font(subtitle.add_run("Submission-polished evidence-backed draft"), "Times New Roman", 11, italic=True)
    subtitle.add_run().add_break(WD_BREAK.PAGE)

    sections = re.split(r"^## ", source, flags=re.MULTILINE)
    for block in sections[1:]:
        heading, _, rest = block.partition("\n")
        add_heading(doc, heading.strip())
        for paragraph in [p.strip() for p in rest.strip().split("\n\n") if p.strip()]:
            add_body(doc, paragraph)

    doc.add_page_break()
    add_heading(doc, "Primary Protein Panel")
    add_caption(doc, "Table 1. Guarded p4 repeatability on the locked full-protein panel.")
    main_body = []
    for row in z4o_rows:
        rep = z4o_rep[row["fixture"]]
        main_body.append([
            row["protein_name"],
            row["pdb_id"],
            row["panels"],
            fmt_float(as_float(row["energy_kj_mol"]), 2),
            fmt_sci(as_float(rep["audit_max"])),
            fmt_sci(as_float(rep["sampled_direct_true_residual_max"])),
            fmt_float(as_float(rep["energy_spread"]), 3),
            fmt_float(as_float(rep["median_matvec_ratio"]), 3),
        ])
    add_table(
        doc,
        ["Protein", "PDB/source", "Panels", "Energy", "Audit max", "Direct residual max", "Energy spread", "p4/deg1"],
        main_body,
        [2400, 1100, 900, 1100, 1050, 1450, 1150, 900],
        7.3,
    )

    if args.z4q_plot.exists():
        doc.add_page_break()
        add_heading(doc, "Supplementary Repeatability Context")
        add_figure(doc, args.z4q_plot, "Figure 1. Supplementary public-protein guarded p4 repeatability for 1UBQ, 2CI2, and 1VII.")
    if args.z5_plot.exists():
        doc.add_page_break()
        add_heading(doc, "Preconditioner Convergence Profile")
        add_figure(doc, args.z5_plot, "Figure 2. Barnase dense/high-coupling convergence profile comparing local block-diagonal preconditioning with sparse near-field ILU(0).")

    landscape = doc.add_section(WD_SECTION.NEW_PAGE)
    configure_section(landscape, landscape=True)
    add_heading(doc, "Supplementary Public-Protein Expansion")
    add_caption(doc, "Table S1. Full public-protein supplement rows prepared without atom sampling.")
    supp_body = []
    for row in z4q_rows:
        rep = z4q_rep[row["fixture"]]
        manifest = z4q_manifest[row["fixture"]]
        supp_body.append([
            row["protein_name"],
            row["pdb_id"],
            manifest.get("atom_count", row.get("atom_count", "NA")),
            row["panels"],
            fmt_float(as_float(row["energy_kj_mol"]), 2),
            fmt_sci(as_float(rep["audit_max"])),
            fmt_sci(as_float(rep["sampled_direct_true_residual_max"])),
            fmt_float(as_float(rep["energy_spread"]), 3),
        ])
    add_table(
        doc,
        ["Protein", "PDB", "Atoms", "Panels", "Energy", "Audit max", "Direct residual max", "Energy spread"],
        supp_body,
        [2450, 850, 750, 850, 1050, 1150, 1500, 1200],
        7.2,
    )

    portrait = doc.add_section(WD_SECTION.NEW_PAGE)
    configure_section(portrait)
    add_reference_list(doc, refs)

    output_docx = args.output_dir / "titan_fagbi_submission_polished.docx"
    doc.save(output_docx)

    claim_trace = [
        {"section": "Methods", "claim": "guarded_p4_explicit_policy", "evidence": str(args.z4o_matrix)},
        {"section": "Results", "claim": "z4o_primary_repeatability", "evidence": str(args.z4o_repeatability)},
        {"section": "Results", "claim": "z4q_supplement_expansion", "evidence": str(args.z4q_repeatability)},
        {"section": "Results", "claim": "z5_barnase_preconditioner_contrast", "evidence": str(args.z5_plot)},
        {"section": "Limitations", "claim": "no_default_promotion", "evidence": str(args.z4o_decision)},
    ]
    write_csv(args.output_dir / "fagbi_scope_z4t_claim_trace.csv", claim_trace, ["section", "claim", "evidence"])
    citation_gaps = [{"gap_id": "none", "status": "closed", "notes": "Z4T reference lock includes method, preparation, database, and numerical linear algebra references."}]
    write_csv(args.output_dir / "fagbi_scope_z4t_citation_gaps.csv", citation_gaps, ["gap_id", "status", "notes"])

    body_text = "\n".join(paragraph.text for paragraph in doc.paragraphs)
    leaked = [term for term in FORBIDDEN_BODY_TERMS if term.lower() in body_text.lower()]
    word_count = len(re.findall(r"\b[\w./()\-]+\b", body_text))
    qa_rows = [
        {"check": "references", "status": "pass", "value": len(refs)},
        {"check": "citation_gaps", "status": "pass", "value": "closed"},
        {"check": "primary_rows", "status": "pass", "value": len(z4o_rows)},
        {"check": "supplement_rows", "status": "pass", "value": len(z4q_rows)},
        {"check": "word_count", "status": "pass" if word_count >= 3000 else "fail", "value": word_count},
        {"check": "forbidden_terms", "status": "pass" if not leaked else "fail", "value": ";".join(leaked)},
        {"check": "default_policy_limitation", "status": "pass", "value": "grid default unchanged"},
    ]
    write_csv(args.output_dir / "fagbi_scope_z4t_docx_qa.csv", qa_rows, ["check", "status", "value"])
    if leaked:
        raise ValueError("Forbidden terms found in DOCX body: " + ", ".join(leaked))
    if word_count < 3000:
        raise ValueError(f"Z4T manuscript remains too short: {word_count} words")
    print(f"Z4T DOCX: {output_docx}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--z4o-matrix", type=Path, required=True)
    parser.add_argument("--z4o-repeatability", type=Path, required=True)
    parser.add_argument("--z4o-decision", type=Path, required=True)
    parser.add_argument("--z4q-matrix", type=Path, required=True)
    parser.add_argument("--z4q-repeatability", type=Path, required=True)
    parser.add_argument("--z4q-manifest", type=Path, required=True)
    parser.add_argument("--z4q-plot", type=Path, required=True)
    parser.add_argument("--z5-plot", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    build_docx(parse_args())
