#!/usr/bin/env python
r"""Report amylin PMF deltaG for monomer and fibril systems.

The reported deltaG is plateau/tail mean minus the PMF minimum. Bootstrap
uncertainty is computed from bootstrap profiles when available; for fibril
4631, only WHAM bsResult mean/stddev is present, so the reported uncertainty is
propagated from the bootstrap pointwise stddev over the same min/tail regions.
"""

from __future__ import annotations

import csv
import math
import re
import shutil
import statistics
import tempfile
from dataclasses import dataclass
from pathlib import Path


KJ_TO_KCAL = 1.0 / 4.184


@dataclass(frozen=True)
class Case:
    state: str
    ligand: str
    pmf: Path
    bootstrap_profiles: Path | None
    bootstrap_result: Path | None
    plateau_range: tuple[float, float]
    scale_bootstrap_sigma: float = 1.0
    note: str = ""


@dataclass
class Result:
    state: str
    ligand: str
    pmf_path: Path
    bootstrap_source: str
    unit: str
    plateau_range: tuple[float, float]
    n_plateau_points: int
    min_x: float
    min_y: float
    plateau_mean: float
    delta_g: float
    bootstrap_n: int | None
    bootstrap_mean: float | None
    bootstrap_sd: float | None
    note: str


def detect_unit(path: Path) -> str:
    header: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(("#", "@")):
                header.append(line.lower())
                continue
            if line and not line.startswith("&"):
                break
    text = "\n".join(header)
    if "kcal" in text or "-unit kcal" in text:
        return "kcal/mol"
    if "kj" in text or "-unit kj" in text:
        return "kJ/mol"
    return "unknown"


def parse_xvg_xy(path: Path) -> tuple[list[float], list[float], str]:
    unit = detect_unit(path)
    scale = KJ_TO_KCAL if unit == "kJ/mol" else 1.0
    xs: list[float] = []
    ys: list[float] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith(("#", "@", "&")):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]) * scale)
            except ValueError:
                continue
    if not xs:
        raise ValueError(f"No numeric data found in {path}")
    return xs, ys, unit


def parse_xvg_xydy(path: Path) -> tuple[list[float], list[float], list[float], str]:
    unit = detect_unit(path)
    scale = KJ_TO_KCAL if unit == "kJ/mol" else 1.0
    xs: list[float] = []
    ys: list[float] = []
    dys: list[float] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith(("#", "@", "&")):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]) * scale)
                dys.append(float(parts[2]) * scale)
            except ValueError:
                continue
    if not xs:
        raise ValueError(f"No xydy numeric data found in {path}")
    return xs, ys, dys, unit


def parse_bootstrap_profiles(path: Path) -> tuple[list[tuple[list[float], list[float]]], str]:
    unit = detect_unit(path)
    scale = KJ_TO_KCAL if unit == "kJ/mol" else 1.0
    profiles: list[tuple[list[float], list[float]]] = []
    xs: list[float] = []
    ys: list[float] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if line == "&":
                if xs:
                    profiles.append((xs, ys))
                    xs, ys = [], []
                continue
            if not line or line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]) * scale)
            except ValueError:
                continue
    if xs:
        profiles.append((xs, ys))
    if not profiles:
        raise ValueError(f"No bootstrap profiles found in {path}")
    return profiles, unit


def delta_from_xy(xs: list[float], ys: list[float], plateau_range: tuple[float, float]) -> tuple[float, float, float, float, int]:
    min_idx = min(range(len(ys)), key=lambda i: ys[i])
    lo, hi = plateau_range
    plateau = [y for x, y in zip(xs, ys) if lo <= x <= hi and math.isfinite(y)]
    if not plateau:
        raise ValueError(f"No PMF points inside plateau range {lo}-{hi}")
    plateau_mean = statistics.fmean(plateau)
    return xs[min_idx], ys[min_idx], plateau_mean, plateau_mean - ys[min_idx], len(plateau)


def propagated_xydy_delta(
    xs: list[float],
    ys: list[float],
    dys: list[float],
    plateau_range: tuple[float, float],
    sigma_scale: float,
) -> tuple[float, float, float, float, int]:
    min_idx = min(range(len(ys)), key=lambda i: ys[i])
    lo, hi = plateau_range
    plateau_sigmas = [dy * sigma_scale for x, dy in zip(xs, dys) if lo <= x <= hi and math.isfinite(dy)]
    if not plateau_sigmas:
        raise ValueError(f"No bootstrap stddev points inside plateau range {lo}-{hi}")
    # SD of a mean from independent pointwise bootstrap SDs; approximate because
    # bsResult does not include the profile covariance.
    plateau_mean_sigma = math.sqrt(sum(s * s for s in plateau_sigmas)) / len(plateau_sigmas)
    min_sigma = dys[min_idx] * sigma_scale
    delta_sigma = math.sqrt(plateau_mean_sigma * plateau_mean_sigma + min_sigma * min_sigma)
    return delta_sigma, min_sigma, plateau_mean_sigma, statistics.fmean(plateau_sigmas), len(plateau_sigmas)


def write_outputs(results: list[Result], output_prefix: Path) -> None:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = output_prefix.with_suffix(".csv")
    txt_path = output_prefix.with_suffix(".txt")
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "state",
                "ligand",
                "deltaG_pmf_kcal_mol",
                "bootstrap_mean_kcal_mol",
                "bootstrap_sd_kcal_mol",
                "bootstrap_n",
                "pmf_min_x_nm",
                "pmf_min_y_kcal_mol",
                "plateau_range_nm",
                "plateau_mean_kcal_mol",
                "plateau_points",
                "pmf_unit",
                "pmf_path",
                "bootstrap_source",
                "note",
            ]
        )
        for r in results:
            writer.writerow(
                [
                    r.state,
                    r.ligand,
                    f"{r.delta_g:.6f}",
                    "" if r.bootstrap_mean is None else f"{r.bootstrap_mean:.6f}",
                    "" if r.bootstrap_sd is None else f"{r.bootstrap_sd:.6f}",
                    "" if r.bootstrap_n is None else r.bootstrap_n,
                    f"{r.min_x:.6f}",
                    f"{r.min_y:.6f}",
                    f"{r.plateau_range[0]:.3f}-{r.plateau_range[1]:.3f}",
                    f"{r.plateau_mean:.6f}",
                    r.n_plateau_points,
                    r.unit,
                    r.pmf_path,
                    r.bootstrap_source,
                    r.note,
                ]
            )
    lines = ["Amylin PMF deltaG report (kcal/mol)", ""]
    for r in results:
        err = "NA" if r.bootstrap_sd is None else f"{r.bootstrap_sd:.3f}"
        boot = "" if r.bootstrap_mean is None else f"; bootstrap mean {r.bootstrap_mean:.3f}"
        lines.append(
            f"{r.state:7s} {r.ligand:9s}: deltaG = {r.delta_g:.3f} +/- {err}{boot}; "
            f"min {r.min_x:.3f} nm; plateau {r.plateau_range[0]:.2f}-{r.plateau_range[1]:.2f} nm"
        )
        if r.note:
            lines.append(f"  note: {r.note}")
    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def copy_outputs_to_d(src_prefix: Path, dst_prefix: Path) -> None:
    for ext in (".csv", ".txt"):
        shutil.copyfile(src_prefix.with_suffix(ext), dst_prefix.with_suffix(ext))


def main() -> int:
    cases = [
        Case(
            "monomer",
            "4497",
            Path(r"C:\tmp\amylin_deltaG_cache\monomer_4497_pmf.xvg"),
            Path(r"C:\tmp\amylin_deltaG_cache\monomer_4497_bsProfs.xvg"),
            None,
            (2.33, 3.79),
            note="scaled flattened PMF, factor 1.75",
        ),
        Case(
            "monomer",
            "4631",
            Path(r"C:\tmp\amylin_deltaG_cache\monomer_4631_pmf.xvg"),
            Path(r"C:\tmp\amylin_deltaG_cache\monomer_4631_bsProfs.xvg"),
            None,
            (3.42, 4.66),
            note="scaled flattened PMF, factor 1.75",
        ),
        Case(
            "fibril",
            "4497",
            Path(r"C:\Users\ADMIN\OneDrive\Documents\DTU_MsLinh\new_conf\US-4479\hist.xvg"),
            Path(r"C:\Users\ADMIN\OneDrive\Documents\DTU_MsLinh\new_conf\US-4479\bsProfs.xvg"),
            Path(r"C:\Users\ADMIN\OneDrive\Documents\DTU_MsLinh\new_conf\US-4479\bsResult.xvg"),
            (1.50, 2.00),
            note="tail range matches fibril plot cutoff; bootstrap n=5",
        ),
        Case(
            "fibril",
            "4631",
            Path(r"C:\Users\ADMIN\OneDrive\Documents\DTU_MsLinh\new_conf\US-4631\scaled_hist.xvg"),
            None,
            Path(r"C:\Users\ADMIN\OneDrive\Documents\DTU_MsLinh\new_conf\US-4631\bsResult.xvg"),
            (1.50, 1.895),
            scale_bootstrap_sigma=1.6,
            note="scaled_hist uses scale 1.6; bsResult has pointwise stddev only",
        ),
    ]

    results: list[Result] = []
    for case in cases:
        xs, ys, unit = parse_xvg_xy(case.pmf)
        min_x, min_y, plateau_mean, delta_g, n_plateau = delta_from_xy(xs, ys, case.plateau_range)
        bootstrap_source = "none"
        bootstrap_n: int | None = None
        bootstrap_mean: float | None = None
        bootstrap_sd: float | None = None

        if case.bootstrap_profiles is not None and case.bootstrap_profiles.exists():
            profiles, _ = parse_bootstrap_profiles(case.bootstrap_profiles)
            deltas = [delta_from_xy(px, py, case.plateau_range)[3] for px, py in profiles]
            bootstrap_source = str(case.bootstrap_profiles)
            bootstrap_n = len(deltas)
            bootstrap_mean = statistics.fmean(deltas)
            bootstrap_sd = statistics.stdev(deltas) if len(deltas) > 1 else 0.0
        elif case.bootstrap_result is not None and case.bootstrap_result.exists():
            bx, by, bdy, _ = parse_xvg_xydy(case.bootstrap_result)
            bootstrap_sd, _, _, _, bootstrap_n = propagated_xydy_delta(
                bx, by, bdy, case.plateau_range, case.scale_bootstrap_sigma
            )
            bootstrap_mean = delta_g
            bootstrap_source = str(case.bootstrap_result)

        results.append(
            Result(
                state=case.state,
                ligand=case.ligand,
                pmf_path=case.pmf,
                bootstrap_source=bootstrap_source,
                unit=unit,
                plateau_range=case.plateau_range,
                n_plateau_points=n_plateau,
                min_x=min_x,
                min_y=min_y,
                plateau_mean=plateau_mean,
                delta_g=delta_g,
                bootstrap_n=bootstrap_n,
                bootstrap_mean=bootstrap_mean,
                bootstrap_sd=bootstrap_sd,
                note=case.note,
            )
        )

    tmp_prefix = Path(tempfile.gettempdir()) / "amylin_monomer_fibril_deltaG_pmf_report"
    dst_prefix = Path(r"D:\Data_MD\amylin\amylin_monomer_fibril_deltaG_pmf_report")
    write_outputs(results, tmp_prefix)
    try:
        copy_outputs_to_d(tmp_prefix, dst_prefix)
        out_prefix = dst_prefix
    except OSError:
        out_prefix = tmp_prefix
    print((out_prefix.with_suffix(".txt")).read_text(encoding="utf-8"))
    print(f"CSV: {out_prefix.with_suffix('.csv')}")
    print(f"TXT: {out_prefix.with_suffix('.txt')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
