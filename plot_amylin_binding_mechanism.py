#!/usr/bin/env python
r"""Plot PMF and TITAN elec/vdW binding contributions for two amylin ligands.

Default inputs target:
  D:\Data_MD\amylin\md_2ligands\4497-2l86
  D:\Data_MD\amylin\md_2ligands\4631-2l86

Outputs are publication-oriented PNG/TIFF at 600 dpi plus a vector PDF.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


KJ_TO_KCAL = 1.0 / 4.184


@dataclass
class LigandData:
    label: str
    system_dir: Path
    pmf_x: list[float]
    pmf_kcal: list[float]
    pmf_min_x: float
    pmf_min_y: float
    elec_kcal: float
    elec_ci_kcal: float | None
    vdw_kcal: float
    vdw_ci_kcal: float | None
    total_kcal: float
    total_ci_kcal: float | None


def parse_xvg(path: Path) -> tuple[list[float], list[float]]:
    xs: list[float] = []
    ys: list[float] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
            except ValueError:
                continue
    if not xs:
        raise ValueError(f"No numeric XVG data found in {path}")
    return xs, ys


def parse_average_line(text: str, label_regex: str) -> tuple[float, float | None]:
    line_pattern = re.compile(rf"Average\s+{label_regex}\s+energy:", re.IGNORECASE)
    target_line = next((line for line in text.splitlines() if line_pattern.search(line)), None)
    if target_line is None:
        raise ValueError(f"Could not parse Average {label_regex} energy from cal.txt")
    mean_match = re.search(r"energy:\s*([-+0-9.eE]+)\s+kJ/mol", target_line, re.IGNORECASE)
    if not mean_match:
        raise ValueError(f"Could not parse mean for Average {label_regex} energy")
    ci_match = re.search(r"Confidence interval\s*\(95%\):\s*([-+0-9.eE]+)\s+kJ/mol", target_line, re.IGNORECASE)
    mean = float(mean_match.group(1)) * KJ_TO_KCAL
    ci = float(ci_match.group(1)) * KJ_TO_KCAL if ci_match else None
    return mean, ci


def parse_titan_cal(path: Path) -> tuple[float, float | None, float, float | None, float, float | None]:
    text = path.read_text(encoding="utf-8", errors="replace")
    elec, elec_ci = parse_average_line(text, r"Electrostatic")
    vdw, vdw_ci = parse_average_line(text, r"Van\s+der\s+waals")
    total, total_ci = parse_average_line(text, r"total")
    return elec, elec_ci, vdw, vdw_ci, total, total_ci


def read_optional_window_centers(system_dir: Path, windows: Iterable[int]) -> tuple[float, float] | None:
    """Best-effort read of umbrella window center range for PMF shading.

    The script does not depend on this because CHARMM/GROMACS helper CSV
    schemas vary. It looks for a window/id column and a center/coord/distance
    column; if unavailable, panel A still marks the PMF minima directly.
    """

    csv_path = system_dir / "umbrella_windows.csv"
    if not csv_path.exists():
        return None
    wanted = set(windows)
    with csv_path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return None
        fields = reader.fieldnames

        def pick(names: tuple[str, ...]) -> str | None:
            for field in fields:
                low = field.lower()
                if any(token in low for token in names):
                    return field
            return None

        window_col = pick(("window", "win", "idx", "index"))
        center_col = pick(("center", "coord", "distance", "pull", "rc"))
        if window_col is None or center_col is None:
            return None
        centers: list[float] = []
        for row in reader:
            try:
                window_id = int(float(row[window_col]))
                center = float(row[center_col])
            except (TypeError, ValueError):
                continue
            if window_id in wanted:
                centers.append(center)
    if not centers:
        return None
    return min(centers), max(centers)


def load_ligand_data(
    root: Path,
    system_name: str,
    label: str,
    pmf_file: str,
    energy_dir_name: str,
) -> LigandData:
    system_dir = root / system_name
    pmf_path = system_dir / pmf_file
    cal_path = system_dir / "data_MD" / energy_dir_name / "cal.txt"
    if not pmf_path.exists():
        raise FileNotFoundError(f"Missing PMF file: {pmf_path}")
    if not cal_path.exists():
        raise FileNotFoundError(f"Missing TITAN cal.txt: {cal_path}")

    # gmx wham was run with "-unit kCal"; these PMF files are already kcal/mol.
    pmf_x, pmf_kcal = parse_xvg(pmf_path)
    min_idx = min(range(len(pmf_kcal)), key=lambda i: pmf_kcal[i])
    elec, elec_ci, vdw, vdw_ci, total, total_ci = parse_titan_cal(cal_path)
    return LigandData(
        label=label,
        system_dir=system_dir,
        pmf_x=pmf_x,
        pmf_kcal=pmf_kcal,
        pmf_min_x=pmf_x[min_idx],
        pmf_min_y=pmf_kcal[min_idx],
        elec_kcal=elec,
        elec_ci_kcal=elec_ci,
        vdw_kcal=vdw,
        vdw_ci_kcal=vdw_ci,
        total_kcal=total,
        total_ci_kcal=total_ci,
    )


def write_summary_csv(path: Path, rows: list[LigandData]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "ligand",
                "system_dir",
                "pmf_min_rc",
                "pmf_min_kcal_mol",
                "elec_kcal_mol",
                "elec_ci95_kcal_mol",
                "vdw_kcal_mol",
                "vdw_ci95_kcal_mol",
                "total_kcal_mol",
                "total_ci95_kcal_mol",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.label,
                    row.system_dir,
                    f"{row.pmf_min_x:.6g}",
                    f"{row.pmf_min_y:.6g}",
                    f"{row.elec_kcal:.6g}",
                    "" if row.elec_ci_kcal is None else f"{row.elec_ci_kcal:.6g}",
                    f"{row.vdw_kcal:.6g}",
                    "" if row.vdw_ci_kcal is None else f"{row.vdw_ci_kcal:.6g}",
                    f"{row.total_kcal:.6g}",
                    "" if row.total_ci_kcal is None else f"{row.total_ci_kcal:.6g}",
                ]
            )


def add_panel_label(ax, label: str) -> None:
    ax.text(
        -0.12,
        1.06,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=13,
        fontweight="bold",
    )


def savefig_via_temp(fig, destination: Path, **kwargs) -> None:
    """Render locally first, then copy to destination.

    This avoids intermittent Pillow/Matplotlib streaming failures on removable
    drives while preserving the requested final output path.
    """

    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="amylin_figure_") as tmpdir:
        tmp_path = Path(tmpdir) / destination.name
        fig.savefig(tmp_path, **kwargs)
        shutil.copyfile(tmp_path, destination)


def make_figure(
    rows: list[LigandData],
    output_prefix: Path,
    dpi: int,
    zero_pmf: bool,
    x_unit: str,
    energy_window_label: str,
    highlight_windows: range,
    align_second_minimum: bool,
    pmf_xmax: float,
) -> None:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    try:
        import seaborn as sns

        sns.set_theme(style="whitegrid")
    except ImportError:
        pass

    tokens = {
        "surface": "#FCFCFD",
        "panel": "#FFFFFF",
        "ink": "#1F2430",
        "muted": "#6F768A",
        "grid": "#E6E8F0",
        "axis": "#D7DBE7",
    }
    blue = {"base": "#A3BEFA", "mid": "#5477C4", "dark": "#2E4780"}
    orange = {"base": "#F0986E", "mid": "#CC6F47", "dark": "#804126"}
    gold = {"base": "#FFE15B", "mid": "#B8A037", "dark": "#736422"}

    plt.rcParams.update(
        {
            "font.family": ["Segoe UI", "DejaVu Sans", "Arial", "sans-serif"],
            "axes.edgecolor": tokens["axis"],
            "axes.labelcolor": tokens["ink"],
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.color": tokens["ink"],
            "ytick.color": tokens["ink"],
            "figure.facecolor": tokens["surface"],
            "axes.facecolor": tokens["panel"],
        }
    )

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(7.2, 3.55),
        gridspec_kw={"width_ratios": [1.35, 1.0], "wspace": 0.34},
        constrained_layout=False,
    )
    fig.patch.set_facecolor(tokens["surface"])
    ax_a, ax_b = axes

    line_styles = ["-", "--"]
    line_colors = [blue["mid"], orange["mid"]]
    x_shifts = [0.0 for _ in rows]
    if align_second_minimum and len(rows) >= 2:
        reference_min_x = rows[0].pmf_min_x
        x_shifts[1] = reference_min_x - rows[1].pmf_min_x

    for idx, row in enumerate(rows):
        x_values = [value + x_shifts[idx] for value in row.pmf_x]
        y_values = row.pmf_kcal
        min_y = min(y_values)
        if zero_pmf:
            y_values = [value - min_y for value in y_values]
            min_marker_y = row.pmf_min_y - min_y
        else:
            min_marker_y = row.pmf_min_y
        ax_a.plot(
            x_values,
            y_values,
            color=line_colors[idx % len(line_colors)],
            linestyle=line_styles[idx % len(line_styles)],
            linewidth=1.8,
            label=row.label,
        )
        ax_a.scatter(
            [row.pmf_min_x + x_shifts[idx]],
            [min_marker_y],
            s=24,
            color=line_colors[idx % len(line_colors)],
            edgecolor=tokens["ink"],
            linewidth=0.6,
            zorder=5,
        )

    for idx, row in enumerate(rows):
        span = read_optional_window_centers(row.system_dir, highlight_windows)
        if span is None:
            continue
        ax_a.axvspan(
            span[0] + x_shifts[idx],
            span[1] + x_shifts[idx],
            color=gold["base"],
            alpha=0.12,
            linewidth=0,
        )

    ax_a.set_title("PMF profile", loc="left", fontweight="bold", color=tokens["ink"], pad=8)
    ax_a.set_xlabel(f"Reaction coordinate ({x_unit})", fontweight="bold")
    ax_a.set_ylabel("PMF (kcal/mol)", fontweight="bold")
    ax_a.set_xlim(0.0, pmf_xmax)
    ax_a.legend(frameon=False, loc="upper right", fontsize=8.5)
    ax_a.grid(False)
    ax_a.grid(False, axis="x")
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)
    ax_a.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))
    add_panel_label(ax_a, "A")

    x_positions = list(range(len(rows)))
    width = 0.34
    elec_values = [row.elec_kcal for row in rows]
    vdw_values = [row.vdw_kcal for row in rows]
    elec_err = [0.0 if row.elec_ci_kcal is None else row.elec_ci_kcal for row in rows]
    vdw_err = [0.0 if row.vdw_ci_kcal is None else row.vdw_ci_kcal for row in rows]

    bars_elec = ax_b.bar(
        [x - width / 2 for x in x_positions],
        elec_values,
        width=width,
        yerr=elec_err if any(elec_err) else None,
        capsize=2.8,
        color=blue["base"],
        edgecolor=blue["dark"],
        linewidth=0.8,
        label="Elec",
    )
    bars_vdw = ax_b.bar(
        [x + width / 2 for x in x_positions],
        vdw_values,
        width=width,
        yerr=vdw_err if any(vdw_err) else None,
        capsize=2.8,
        color=orange["base"],
        edgecolor=orange["dark"],
        linewidth=0.8,
        label="VdW",
    )
    ax_b.axhline(0, color=tokens["ink"], linewidth=0.8)
    ax_b.set_title("Minimum-region interaction energy", loc="left", fontweight="bold", color=tokens["ink"], pad=8)
    ax_b.set_ylabel("Energy contribution (kcal/mol)", fontweight="bold")
    ax_b.set_xticks(x_positions)
    ax_b.set_xticklabels([row.label for row in rows], rotation=0)
    ax_b.legend(
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.02, 0.98),
        borderaxespad=0.0,
        fontsize=8.5,
    )
    ax_b.grid(True, axis="y", color=tokens["grid"], linewidth=0.8)
    ax_b.grid(False, axis="x")
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)
    ax_b.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))

    for bars in (bars_elec, bars_vdw):
        for bar in bars:
            value = bar.get_height()
            ax_b.text(
                bar.get_x() + bar.get_width() / 2,
                value - 0.9 if value < 0 else value + 0.9,
                f"{value:.1f}",
                ha="center",
                va="top" if value < 0 else "bottom",
                fontsize=7.5,
                color=tokens["ink"],
            )

    add_panel_label(ax_b, "B")

    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.16, top=0.86)

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    png = output_prefix.with_suffix(".png")
    tif = output_prefix.with_suffix(".tiff")
    pdf = output_prefix.with_suffix(".pdf")
    savefig_via_temp(fig, png, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
    savefig_via_temp(fig, tif, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
    savefig_via_temp(fig, pdf, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create a two-panel PMF and elec/vdW contribution figure for amylin ligands.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(r"D:\Data_MD\amylin\md_2ligands"),
        help="Root containing the two ligand system folders.",
    )
    parser.add_argument(
        "--systems",
        nargs=2,
        default=["4497-2l86", "4631-2l86"],
        help="Two system folder names under --root.",
    )
    parser.add_argument(
        "--labels",
        nargs=2,
        default=["4497-2l86", "4631-2l86"],
        help="Labels shown in the figure.",
    )
    parser.add_argument(
        "--pmf-file",
        default="pmf_flat_scaled_1p75.xvg",
        help="PMF XVG filename inside each system folder.",
    )
    parser.add_argument(
        "--energy-dir",
        default=r"titan_energy_windows2_9_50_100ns_multi_xdr",
        help="Directory inside data_MD containing TITAN cal.txt.",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path(r"D:\Data_MD\amylin\md_2ligands\amylin_2ligands_binding_mechanism_600dpi"),
        help="Output path without extension.",
    )
    parser.add_argument("--dpi", type=int, default=600, help="Raster export DPI for PNG/TIFF.")
    parser.add_argument("--x-unit", default="nm", help="Reaction-coordinate unit for panel A.")
    parser.add_argument(
        "--keep-pmf-offset",
        action="store_true",
        help="Do not shift each PMF curve to its own minimum.",
    )
    parser.add_argument(
        "--pmf-xmax",
        type=float,
        default=3.5,
        help="Maximum x-axis value shown in PMF panel.",
    )
    parser.add_argument(
        "--no-align-second-minimum",
        action="store_true",
        help="Do not shift the second ligand PMF x-axis minimum onto the first ligand minimum.",
    )
    parser.add_argument(
        "--energy-window-label",
        default="windows 2-9, 50-100 ns",
        help="Subtitle text describing the TITAN energy window.",
    )
    parser.add_argument(
        "--highlight-windows",
        default="2-9",
        help="Umbrella window range used for optional PMF shading.",
    )
    return parser


def parse_window_range(text: str) -> range:
    match = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", text)
    if match:
        start, end = int(match.group(1)), int(match.group(2))
        if end < start:
            start, end = end, start
        return range(start, end + 1)
    values = [int(part.strip()) for part in text.split(",") if part.strip()]
    if not values:
        return range(0)
    return range(min(values), max(values) + 1)


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    rows = [
        load_ligand_data(
            root=args.root,
            system_name=system,
            label=label,
            pmf_file=args.pmf_file,
            energy_dir_name=args.energy_dir,
        )
        for system, label in zip(args.systems, args.labels)
    ]
    highlight_windows = parse_window_range(args.highlight_windows)
    make_figure(
        rows=rows,
        output_prefix=args.output_prefix,
        dpi=args.dpi,
        zero_pmf=not args.keep_pmf_offset,
        x_unit=args.x_unit,
        energy_window_label=args.energy_window_label,
        highlight_windows=highlight_windows,
        align_second_minimum=not args.no_align_second_minimum,
        pmf_xmax=args.pmf_xmax,
    )
    summary_path = args.output_prefix.with_name(args.output_prefix.name + "_summary_kcal.csv")
    write_summary_csv(summary_path, rows)

    print(f"Wrote: {args.output_prefix.with_suffix('.png')}")
    print(f"Wrote: {args.output_prefix.with_suffix('.tiff')}")
    print(f"Wrote: {args.output_prefix.with_suffix('.pdf')}")
    print(f"Wrote: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
