#!/usr/bin/env python3
"""
pbe_precond_cohort_n5_tables.py

Aggregates the six-preconditioner ablation matrices produced by
pbe_precond_cohort_n5_gate.ps1 and prints, in manuscript order:

    Table 3  cohort and panel counts per surface representation
    Table 4  GMRES iterations for the six preconditioners
    Table 5  bounded graph and factor diagnostics for the two ILU(0) routes
             (all three surfaces, including analytic patch-ray SES)
    Derived  the numeric ranges quoted in the Abstract, 3.2 and 3.3

Usage:
    python test\\pbe_precond_cohort_n5_tables.py
    python test\\pbe_precond_cohort_n5_tables.py --out tables_n5.md
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from itertools import product

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DEFAULT_MATRICES = [
    os.path.join(REPO_ROOT, "test", "tmp", "fagbi_graph_preconditioner_ablation",
                 "fagbi_graph_preconditioner_matrix.csv"),
    os.path.join(REPO_ROOT, "test", "tmp", "fagbi_nodepatch_analytic_ses_preconditioner",
                 "analytic_ses_ablation", "fagbi_graph_preconditioner_matrix.csv"),
]
DEFAULT_MANIFEST = os.path.join(REPO_ROOT, "test", "pbe_inputs_public",
                                "public_fixture_manifest.csv")

# manuscript ordering
FIXTURE_ORDER = ["barnase", "1a63", "1sh1", "1a7m", "1ubq"]
SURFACE_ORDER = ["grid-ses", "analytic-ses", "sphere-patch"]

SURFACE_LABEL = {
    "grid-ses": "Watertight Grid-SES",
    "analytic-ses": "Analytic patch-ray SES",
    "sphere-patch": "Sphere-patch stress",
}
SURFACE_SHORT = {
    "grid-ses": "Grid-SES",
    "analytic-ses": "Analytic patch-ray SES",
    "sphere-patch": "Sphere patch",
}
SURFACE_USE = {
    "grid-ses": "Independent surface-generator control",
    "analytic-ses": "Primary physical protocol",
    "sphere-patch": "Algebraic stress test only",
}

# column order of Table 4, mapped to ablation stage labels
STAGES = [
    ("diagonal", "Diag."),
    ("leaf_block_jacobi", "Leaf block"),
    ("fmm_near_graph_ras", "FMM RAS"),
    ("distance_knn_graph_ras", "kNN RAS"),
    ("distance_knn_ilu0", "kNN ILU"),
    ("fmm_near_graph_ilu0", "FMM ILU"),
]
ILU_STAGES = ["distance_knn_ilu0", "fmm_near_graph_ilu0"]
NON_FACTORIZED = ["diagonal", "leaf_block_jacobi", "fmm_near_graph_ras", "distance_knn_graph_ras"]
GRAPH_LABEL = {"distance_knn_ilu0": "distance-kNN", "fmm_near_graph_ilu0": "FMM-near"}

SYSTEM_LABEL = {
    "barnase": ("Barnase", "1BRS"),
    "1a63": ("1A63 Rho domain", "1A63"),
    "1sh1": ("1SH1 neurotoxin", "1SH1"),
    "1a7m": ("1A7M", "1A7M"),
    "1ubq": ("1UBQ ubiquitin", "1UBQ"),
}

ITER_BUDGET = 720
EXPECTED_ROW_COUNT = len(FIXTURE_ORDER) * len(SURFACE_ORDER) * len(STAGES)
EXPECTED_EFFECTIVE = {
    "diagonal": {"diagonal"},
    "leaf_block_jacobi": {"fmm_block_jacobi"},
    "fmm_near_graph_ras": {"fmm_overlap_schwarz"},
    "distance_knn_graph_ras": {"distance_knn_overlap_schwarz"},
    # Regularization changes only the numerical pivot policy of the requested
    # ILU(0) factorization.  It is not a fallback to diagonal scaling.
    "distance_knn_ilu0": {"distance_knn_ilu0", "distance_knn_ilu0_regularized"},
    "fmm_near_graph_ilu0": {"fmm_nearfield_ilu0", "fmm_nearfield_ilu0_regularized"},
}


def _f(value, default=float("nan")):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _i(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def load_matrices(paths):
    rows = []
    seen = set()
    for path in paths:
        if not os.path.exists(path):
            sys.exit(f"error: missing matrix CSV: {path}")
        with open(path, newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                key = (row.get("fixture", ""), row.get("surface_protocol", ""),
                       row.get("preconditioner_label", ""))
                if key in seen:
                    sys.exit(f"error: duplicate ablation row: {'/'.join(key)}")
                seen.add(key)
                rows.append(row)
    if not rows:
        sys.exit("error: no ablation rows found. Run pbe_precond_cohort_n5_gate.ps1 first.")
    return rows


def load_manifest(path):
    metadata = {}
    if not os.path.exists(path):
        sys.exit(f"error: missing fixture manifest: {path}")
    with open(path, newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            metadata[row["fixture"].strip().lower()] = row
    return metadata


def validate_cohort(rows, manifest):
    errors = []
    stage_names = [stage for stage, _ in STAGES]
    expected = set(product(FIXTURE_ORDER, SURFACE_ORDER, stage_names))
    observed = {
        (r.get("fixture", "").strip().lower(),
         r.get("surface_protocol", "").strip().lower(),
         r.get("preconditioner_label", "").strip().lower())
        for r in rows
    }
    missing = sorted(expected - observed)
    extra = sorted(observed - expected)
    if len(rows) != EXPECTED_ROW_COUNT:
        errors.append(f"row count {len(rows)} != {EXPECTED_ROW_COUNT}")
    if missing:
        errors.append("missing rows: " + ", ".join("/".join(k) for k in missing[:12]))
    if extra:
        errors.append("unexpected rows: " + ", ".join("/".join(k) for k in extra[:12]))

    for fixture in FIXTURE_ORDER:
        if fixture not in manifest:
            errors.append(f"manifest row missing for {fixture}")

    titan_hashes = {r.get("titan_sha256", "").strip().lower() for r in rows}
    contract_hashes = {r.get("solver_contract_sha256", "").strip().lower() for r in rows}
    if len(titan_hashes) != 1 or "" in titan_hashes:
        errors.append(f"mixed or missing titan hashes: {sorted(titan_hashes)}")
    if len(contract_hashes) != 1 or "" in contract_hashes:
        errors.append(f"mixed or missing solver-contract hashes: {sorted(contract_hashes)}")

    panel_groups = {}
    for row in rows:
        fixture = row.get("fixture", "").strip().lower()
        surface = row.get("surface_protocol", "").strip().lower()
        stage = row.get("preconditioner_label", "").strip().lower()
        row_id = row.get("row_id") or f"{fixture}/{surface}/{stage}"
        panels = _i(row.get("panels"))
        iterations = _i(row.get("iterations"))
        if panels <= 0:
            errors.append(f"{row_id}: invalid panel count {panels}")
        panel_groups.setdefault((fixture, surface), set()).add(panels)
        if iterations <= 0 or iterations > ITER_BUDGET:
            errors.append(f"{row_id}: invalid iteration count {iterations}")
        status = row.get("status", "").strip().lower()
        residual = _f(row.get("true_residual"))
        nonconverged_control = (
            stage in NON_FACTORIZED
            and status == "failed"
            and math.isfinite(residual)
        )
        if status != "ok" and not nonconverged_control:
            errors.append(f"{row_id}: status={row.get('status', '')}")
        expected_effective = EXPECTED_EFFECTIVE.get(stage, set())
        observed_effective = row.get("preconditioner_effective", "").strip().lower()
        if observed_effective not in expected_effective:
            errors.append(
                f"{row_id}: effective={row.get('preconditioner_effective', '')}, "
                f"expected one of {sorted(expected_effective)}"
            )
        manifest_row = manifest.get(fixture)
        if manifest_row and row.get("pqr_sha256", "").strip().lower() != \
                manifest_row.get("pqr_sha256", "").strip().lower():
            errors.append(f"{row_id}: PQR hash differs from manifest")

        if stage in ILU_STAGES:
            audit = _f(row.get("audit_rms"))
            quality = row.get("quality_status", "").strip().lower()
            if not math.isfinite(residual) or residual > 1e-4:
                errors.append(f"{row_id}: true residual {residual!r} exceeds 1e-4")
            if quality not in ("pass", ""):
                errors.append(f"{row_id}: quality_status={quality}")
            if surface != "sphere-patch" and (not math.isfinite(audit) or audit > 1e-2):
                errors.append(f"{row_id}: physical-surface audit {audit!r} exceeds 1e-2")

    for key, values in sorted(panel_groups.items()):
        if len(values) != 1:
            errors.append(f"{key[0]}/{key[1]}: inconsistent panel counts {sorted(values)}")

    if errors:
        detail = "\n  - ".join(errors[:40])
        sys.exit(f"error: n=5 cohort validation failed:\n  - {detail}")


def order_key(seq, value):
    return seq.index(value) if value in seq else len(seq)


def build_index(rows):
    idx = {}
    for row in rows:
        fixture = row.get("fixture", "").strip().lower()
        surface = row.get("surface_protocol", "").strip().lower()
        stage = row.get("preconditioner_label", "").strip().lower()
        idx[(fixture, surface, stage)] = row
    return idx


def present(rows, key):
    return sorted({r.get(key, "").strip().lower() for r in rows if r.get(key, "").strip()},
                  key=lambda v: order_key(FIXTURE_ORDER if key == "fixture" else SURFACE_ORDER, v))


def fmt_int(value):
    return f"{value:,}" if value else "-"


def table3(idx, fixtures, surfaces, manifest, out):
    out.append("Table 3: PBE test systems and panel counts")
    out.append("")
    out.append("| PBE test system | Geometry source | Atoms | Surface | Panels | Use |")
    out.append("|---|---|---|---|---|---|")
    for surface in surfaces:
        for fixture in fixtures:
            row = next((idx.get((fixture, surface, s)) for s in ILU_STAGES
                        if idx.get((fixture, surface, s))), None)
            if row is None:
                continue
            name, pdb = SYSTEM_LABEL.get(fixture, (fixture.upper(), fixture.upper()))
            out.append("| {} | {} | {} | {} | {} | {} |".format(
                name, pdb, fmt_int(_i(manifest.get(fixture, {}).get("atom_count"))), SURFACE_SHORT[surface],
                fmt_int(_i(row.get("panels"))), SURFACE_USE[surface]))
    out.append("")


def table4(idx, fixtures, surfaces, out):
    out.append("Table 4: GMRES iterations for the six preconditioners "
               f"(a value of {ITER_BUDGET} denotes exhaustion of the iteration budget)")
    out.append("")
    out.append("| PBE test system | Surface | Panels | " + " | ".join(h for _, h in STAGES) + " |")
    out.append("|---|---|---|" + "---|" * len(STAGES))
    for surface in surfaces:
        for fixture in fixtures:
            cells, panels = [], 0
            any_row = False
            for stage, _ in STAGES:
                row = idx.get((fixture, surface, stage))
                if row is None:
                    cells.append("-")
                    continue
                any_row = True
                panels = panels or _i(row.get("panels"))
                cells.append(str(_i(row.get("iterations"))))
            if not any_row:
                continue
            name = SYSTEM_LABEL.get(fixture, (fixture.upper(), ""))[0].split()[0]
            out.append("| {} | {} | {} | {} |".format(
                name, SURFACE_LABEL[surface], fmt_int(panels), " | ".join(cells)))
    out.append("")


def table5(idx, fixtures, surfaces, out):
    out.append("Table 5: bounded graph and factor diagnostics for the two ILU(0) routes")
    out.append("")
    out.append("| PBE test system | Surface | Graph | Panels | Blocks/N | Setup (s) | Apply (s) | Iter. | RSS (MB) |")
    out.append("|---|---|---|---|---|---|---|---|---|")
    for surface in surfaces:
        for fixture in fixtures:
            for stage in ILU_STAGES:
                row = idx.get((fixture, surface, stage))
                if row is None:
                    continue
                name = SYSTEM_LABEL.get(fixture, (fixture.upper(), ""))[0].split()[0]
                out.append("| {} | {} | {} | {} | {} | {:.2f} | {:.2f} | {} | {:.1f} |".format(
                    name, SURFACE_LABEL[surface], GRAPH_LABEL[stage],
                    _i(row.get("panels")), _i(row.get("graph_cap"), 64),
                    _f(row.get("setup_seconds"), 0.0), _f(row.get("apply_seconds"), 0.0),
                    _i(row.get("iterations")), _f(row.get("peak_rss_mb"), 0.0)))
    out.append("")


def derived(rows, idx, fixtures, surfaces, out):
    out.append("Derived quantities for the Abstract, 3.2 and 3.3")
    out.append("")

    ilu = [r for r in rows if r.get("preconditioner_label") in ILU_STAGES]
    ilu_iters = [_i(r.get("iterations")) for r in ilu if _i(r.get("iterations"))]
    panels = [_i(r.get("panels")) for r in ilu if _i(r.get("panels"))]
    n_prot = len({r.get("fixture") for r in ilu})
    out.append(f"- cohort size: {n_prot} PBE test systems x {len(surfaces)} surface representations")
    if ilu_iters:
        out.append(f"- graph ILU(0) iteration range (all surfaces): {min(ilu_iters)}-{max(ilu_iters)}")
    if panels:
        out.append(f"- panel range: {min(panels):,}-{max(panels):,}")

    # normalized memory, MB per 1000 panels
    norm = [_f(r.get("peak_rss_mb")) / (_i(r.get("panels")) / 1000.0)
            for r in ilu if _i(r.get("panels")) and not math.isnan(_f(r.get("peak_rss_mb")))]
    if norm:
        out.append(f"- normalized memory: {min(norm):.1f}-{max(norm):.1f} MB per 1,000 panels")

    # audit / residual envelopes.
    # Sphere-patch rows are reported separately: that representation is an
    # algebraic stress test, and degree-one p1 is not expected to clear its
    # sampled-operator audit (see 3.1).
    phys = [r for r in ilu if r.get("surface_protocol", "").lower() != "sphere-patch"]
    stress = [r for r in ilu if r.get("surface_protocol", "").lower() == "sphere-patch"]
    for label, subset in (("physical surfaces", phys), ("sphere-patch stress", stress)):
        audits = [_f(r.get("audit_rms")) for r in subset if not math.isnan(_f(r.get("audit_rms")))]
        resids = [_f(r.get("true_residual")) for r in subset
                  if not math.isnan(_f(r.get("true_residual")))]
        if audits:
            out.append(f"- sampled operator audit RMS, ILU rows ({label}): "
                       f"{min(audits):.3e} to {max(audits):.3e}")
        if resids:
            out.append(f"- true residual, ILU rows ({label}): "
                       f"{min(resids):.3e} to {max(resids):.3e}")

    out.append("")
    for surface in surfaces:
        sub = [r for r in rows if r.get("surface_protocol", "").lower() == surface]
        if not sub:
            continue
        ilu_s = [_i(r.get("iterations")) for r in sub
                 if r.get("preconditioner_label") in ILU_STAGES and _i(r.get("iterations"))]
        non_s = [_i(r.get("iterations")) for r in sub
                 if r.get("preconditioner_label") in NON_FACTORIZED and _i(r.get("iterations"))]
        stalled = sum(1 for v in non_s if v >= ITER_BUDGET)
        line = f"- {SURFACE_LABEL[surface]}: ILU(0) {min(ilu_s)}-{max(ilu_s)} iterations" if ilu_s else f"- {SURFACE_LABEL[surface]}:"
        if non_s:
            line += (f"; non-factorized controls {min(non_s)}-{max(non_s)}"
                     f" ({stalled}/{len(non_s)} exhausted the {ITER_BUDGET}-iteration budget)")
        out.append(line)

    out.append("")
    # integrity checks a reviewer would run
    out.append("Consistency checks")
    out.append("")
    # Non-convergence of the non-factorized controls is the expected result and
    # is not an error; only ILU rows must all be clean.
    bad_status = [r["row_id"] for r in rows
                  if r.get("preconditioner_label") in ILU_STAGES
                  and r.get("status", "ok") != "ok"]
    stalled_controls = [r["row_id"] for r in rows
                        if r.get("preconditioner_label") in NON_FACTORIZED
                        and _i(r.get("iterations")) >= ITER_BUDGET]
    fallback = [r["row_id"] for r in rows
                if r.get("preconditioner_label") in ILU_STAGES
                and r.get("quality_status", "pass") not in ("pass", "")]
    out.append(f"- non-factorized control rows that exhausted the budget "
               f"(expected, not an error): {len(stalled_controls)}")
    mismatch = []
    for fixture in fixtures:
        for surface in surfaces:
            a = idx.get((fixture, surface, "distance_knn_ilu0"))
            b = idx.get((fixture, surface, "fmm_near_graph_ilu0"))
            if a and b and _i(a.get("iterations")) != _i(b.get("iterations")):
                mismatch.append(f"{fixture}/{surface}: kNN {_i(a.get('iterations'))} vs "
                                f"FMM {_i(b.get('iterations'))}")
    out.append(f"- ILU rows with status != ok: {len(bad_status)}" +
               (" -> " + ", ".join(bad_status[:8]) if bad_status else ""))
    out.append(f"- ILU rows failing the preconditioner quality gate: {len(fallback)}" +
               (" -> " + ", ".join(fallback[:8]) if fallback else ""))
    out.append("- kNN vs FMM-near ILU iteration mismatches: " +
               (", ".join(mismatch) if mismatch else "none (both routes agree on every row)"))
    out.append("")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--matrix", action="append", default=None,
                        help="ablation matrix CSV (repeatable); defaults to both standard paths")
    parser.add_argument("--manifest", default=DEFAULT_MANIFEST,
                        help="public fixture manifest, for atom counts")
    parser.add_argument("--out", default=None, help="write markdown to this file as well as stdout")
    parser.add_argument("--validate-only", action="store_true",
                        help="validate the complete n=5 cohort without rendering tables")
    args = parser.parse_args()

    rows = load_matrices(args.matrix or DEFAULT_MATRICES)
    manifest = load_manifest(args.manifest)
    validate_cohort(rows, manifest)
    if args.validate_only:
        print(f"validated: {EXPECTED_ROW_COUNT} rows, one executable hash, one solver contract")
        return
    idx = build_index(rows)
    fixtures = present(rows, "fixture")
    surfaces = present(rows, "surface_protocol")

    out = []
    out.append(f"# PBE preconditioner cohort: {len(fixtures)} test systems "
               f"({', '.join(f.upper() for f in fixtures)})")
    out.append("")
    table3(idx, fixtures, surfaces, manifest, out)
    table4(idx, fixtures, surfaces, out)
    table5(idx, fixtures, surfaces, out)
    derived(rows, idx, fixtures, surfaces, out)

    text = "\n".join(out)
    print(text)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")
        print(f"\nwritten: {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
