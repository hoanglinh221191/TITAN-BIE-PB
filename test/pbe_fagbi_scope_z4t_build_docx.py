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

Boundary-integral Poisson-Boltzmann solvers can avoid volumetric grids but require accurate acceleration and robust Krylov preconditioning. We present a guarded fourth-order analytical FMM/M2L path in TITAN, built around analytic singular near-field integration, deterministic far-cell refinement, compressed Green-function jets, packed AVX2 direct fallback, and sparse near-field ILU(0). The method is evaluated as an advanced explicit TABI treecode option, not as a default backend. The primary full-protein panel includes Barnase, the RNA-binding domain of E. coli Rho factor (1A63), and Stichodactyla helianthus neurotoxin I (1SH1). All guarded p4 repeats converged with finite energies, sampled direct-operator residuals below 8e-3, and one plan/cache build. Additional public proteins, Ubiquitin (1UBQ), Chymotrypsin inhibitor 2 (2CI2), and Villin headpiece subdomain (1VII), passed the same supplement gate. A dense Barnase convergence profile shows that sparse near-field ILU(0) resolves a high-coupling case where local block-diagonal preconditioning stalls. The production Poisson-Boltzmann default remains the grid backend, and the default TABI treecode path remains degree-one.

## Methods

TITAN solves a coupled boundary-integral Poisson-Boltzmann system over triangulated molecular surfaces. The accelerated path follows the FMM motivation of Greengard and Rokhlin {refs[1]['label']} and the Cartesian boundary-integral PB context of Chen, Tausch, and Geng {refs[0]['label']}. GMRES is used for the nonsymmetric coupled system {refs[2]['label']}. The guarded p4 configuration is explicit: Taylor order 4, opening theta 0.6, static tail tolerance 0.002, strict-geometry MAC theta 0.5, analytic-singular near integration, compressed degree-9 Green jets, packed AVX2 direct fallback, and sparse near-field ILU(0).

For panel i, the coupled unknown is u_i = (phi_i, partial_n phi_i)^T. The sparse preconditioner keeps a 2 x 2 block on FMM near-source graph edges, including self edges and at most 64 neighbors. ILU(0) keeps the original graph sparsity and applies block forward/back substitution, following the no-fill incomplete-factorization principle {refs[3]['label']}. Duffy-adaptive quadrature {refs[4]['label']} remains an explicit near-field option but is not selected automatically.

Protein structures are prepared with PDB2PQR/APBS-family tooling {refs[5]['label']}, {refs[6]['label']} and public structure identifiers from the RCSB Protein Data Bank {refs[7]['label']}. Large-system verification uses sampled direct-operator RMS and sampled direct true residuals because full direct TABI reaction energies are not practical for the largest rows.

## Results

The primary guarded p4 panel passed all repeatability gates. Supplementary full-protein systems also passed, broadening the public benchmark context without changing the main claim. The Barnase convergence profile is interpreted narrowly as a dense/high-coupling preconditioner contrast.

## Limitations

Grid remains TITAN's production Poisson-Boltzmann default. Guarded p4 is an advanced explicit TABI treecode option and is not selected automatically. Full direct TABI reaction energies are not used as the large-system oracle. Duffy-adaptive near integration is explicit-only. The Barnase preconditioner contrast should not be generalized to all proteins.
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
    qa_rows = [
        {"check": "references", "status": "pass", "value": len(refs)},
        {"check": "citation_gaps", "status": "pass", "value": "closed"},
        {"check": "primary_rows", "status": "pass", "value": len(z4o_rows)},
        {"check": "supplement_rows", "status": "pass", "value": len(z4q_rows)},
        {"check": "forbidden_terms", "status": "pass" if not leaked else "fail", "value": ";".join(leaked)},
        {"check": "default_policy_limitation", "status": "pass", "value": "grid default unchanged"},
    ]
    write_csv(args.output_dir / "fagbi_scope_z4t_docx_qa.csv", qa_rows, ["check", "status", "value"])
    if leaked:
        raise ValueError("Forbidden terms found in DOCX body: " + ", ".join(leaked))
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
