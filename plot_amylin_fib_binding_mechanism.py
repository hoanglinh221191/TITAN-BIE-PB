#!/usr/bin/env python
r"""Plot PMF and TITAN elec/vdW contributions for amylin fib ligands.

PMF inputs:
  4497_fib: C:\Users\ADMIN\OneDrive\Documents\DTU_MsLinh\new_conf\US-4479\hist.xvg
  4631_fib: C:\Users\ADMIN\OneDrive\Documents\DTU_MsLinh\new_conf\US-4631\scaled_hist.xvg

TITAN inputs:
  D:\Data_MD\amylin\4497_fib\data_MD\titan_energy_md2_9_multi_xdr_amber14sb_ligtopo\cal.txt
  D:\Data_MD\amylin\4631_fib\data_MD\titan_energy_md2_9_multi_xdr_amber14sb_ligtopo\cal.txt
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path


KJ_TO_KCAL = 1.0 / 4.184


@dataclass
class Curve:
    label: str
    pmf_path: Path
    cal_path: Path
    pmf_unit: str
    x: list[float]
    y_kcal: list[float]
    min_x: float
    min_y: float
    elec_kcal: float
    elec_ci_kcal: float | None
    vdw_kcal: float
    vdw_ci_kcal: float | None
    total_kcal: float
    total_ci_kcal: float | None


def detect_xvg_energy_unit(path: Path) -> str:
    header = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(("#", "@")):
                header.append(line.lower())
                continue
            if line:
                break
    text = "\n".join(header)
    if "kcal" in text or "kcal mol" in text or "-unit kcal" in text:
        return "kcal/mol"
    if "kj" in text or "kj mol" in text or "-unit kj" in text:
        return "kJ/mol"
    return "unknown"


def parse_xvg(path: Path) -> tuple[list[float], list[float], str]:
    unit = detect_xvg_energy_unit(path)
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
        raise ValueError(f"No numeric XVG data found in {path}")
    return xs, ys, unit


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


def load_curve(label: str, pmf_path: Path, cal_path: Path) -> Curve:
    x, y, unit = parse_xvg(pmf_path)
    min_idx = min(range(len(y)), key=lambda i: y[i])
    elec, elec_ci, vdw, vdw_ci, total, total_ci = parse_titan_cal(cal_path)
    return Curve(
        label=label,
        pmf_path=pmf_path,
        cal_path=cal_path,
        pmf_unit=unit,
        x=x,
        y_kcal=y,
        min_x=x[min_idx],
        min_y=y[min_idx],
        elec_kcal=elec,
        elec_ci_kcal=elec_ci,
        vdw_kcal=vdw,
        vdw_ci_kcal=vdw_ci,
        total_kcal=total,
        total_ci_kcal=total_ci,
    )


def savefig_via_temp(fig, destination: Path, **kwargs) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="amylin_fib_figure_") as tmpdir:
        tmp_path = Path(tmpdir) / destination.name
        fig.savefig(tmp_path, **kwargs)
        shutil.copyfile(tmp_path, destination)


def write_summary(path: Path, curves: list[Curve], x_shifts: list[float]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "ligand",
                "pmf_path",
                "pmf_detected_unit",
                "pmf_min_x_raw_nm",
                "pmf_x_shift_nm",
                "pmf_min_kcal_mol",
                "elec_kcal_mol",
                "elec_ci95_kcal_mol",
                "vdw_kcal_mol",
                "vdw_ci95_kcal_mol",
                "total_kcal_mol",
                "total_ci95_kcal_mol",
                "cal_path",
            ]
        )
        for curve, shift in zip(curves, x_shifts):
            writer.writerow(
                [
                    curve.label,
                    curve.pmf_path,
                    curve.pmf_unit,
                    f"{curve.min_x:.6g}",
                    f"{shift:.6g}",
                    f"{curve.min_y:.6g}",
                    f"{curve.elec_kcal:.6g}",
                    "" if curve.elec_ci_kcal is None else f"{curve.elec_ci_kcal:.6g}",
                    f"{curve.vdw_kcal:.6g}",
                    "" if curve.vdw_ci_kcal is None else f"{curve.vdw_ci_kcal:.6g}",
                    f"{curve.total_kcal:.6g}",
                    "" if curve.total_ci_kcal is None else f"{curve.total_ci_kcal:.6g}",
                    curve.cal_path,
                ]
            )


def add_panel_label(ax, text: str) -> None:
    ax.text(
        -0.12,
        1.04,
        text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=13,
        fontweight="bold",
        color="#111111",
    )


def make_figure(curves: list[Curve], output_prefix: Path, dpi: int, pmf_xmax: float) -> None:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    colors = ["#2F6DB3", "#C85A3A"]
    line_styles = ["-", "--"]
    elec_color = "#8EB6E5"
    vdw_color = "#F0A36E"
    ink = "#111111"
    grid = "#E5E7EB"

    plt.rcParams.update(
        {
            "font.family": ["Segoe UI", "DejaVu Sans", "Arial", "sans-serif"],
            "axes.labelsize": 10,
            "axes.labelcolor": ink,
            "axes.edgecolor": "#A8AFBD",
            "xtick.color": ink,
            "ytick.color": ink,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "figure.facecolor": "#FFFFFF",
            "axes.facecolor": "#FFFFFF",
        }
    )

    fig, (ax_a, ax_b) = plt.subplots(
        1,
        2,
        figsize=(7.25, 3.45),
        gridspec_kw={"width_ratios": [1.35, 1.0], "wspace": 0.36},
    )

    x_shifts = [0.0 for _ in curves]
    if len(curves) >= 2:
        x_shifts[1] = curves[0].min_x - curves[1].min_x

    for idx, curve in enumerate(curves):
        shifted_x = [value + x_shifts[idx] for value in curve.x]
        min_y = min(curve.y_kcal)
        shifted_y = [value - min_y for value in curve.y_kcal]
        ax_a.plot(
            shifted_x,
            shifted_y,
            color=colors[idx % len(colors)],
            linestyle=line_styles[idx % len(line_styles)],
            linewidth=1.9,
            label=curve.label,
        )
        ax_a.scatter(
            [curve.min_x + x_shifts[idx]],
            [curve.min_y - min_y],
            s=24,
            color=colors[idx % len(colors)],
            edgecolor=ink,
            linewidth=0.6,
            zorder=5,
        )

    ax_a.set_xlim(0.0, pmf_xmax)
    ax_a.set_xlabel("Reaction coordinate (nm)", fontweight="bold")
    ax_a.set_ylabel("PMF (kcal/mol)", fontweight="bold")
    ax_a.legend(frameon=False, loc="upper left", fontsize=8.5)
    ax_a.grid(False)
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)
    ax_a.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))
    ax_a.tick_params(axis="both", colors=ink)
    add_panel_label(ax_a, "A")

    x_positions = list(range(len(curves)))
    width = 0.34
    elec = [curve.elec_kcal for curve in curves]
    vdw = [curve.vdw_kcal for curve in curves]
    elec_err = [0.0 if curve.elec_ci_kcal is None else curve.elec_ci_kcal for curve in curves]
    vdw_err = [0.0 if curve.vdw_ci_kcal is None else curve.vdw_ci_kcal for curve in curves]

    bars_elec = ax_b.bar(
        [x - width / 2 for x in x_positions],
        elec,
        width=width,
        yerr=elec_err if any(elec_err) else None,
        capsize=2.8,
        color=elec_color,
        edgecolor="#245487",
        linewidth=0.8,
        label="Elec",
    )
    bars_vdw = ax_b.bar(
        [x + width / 2 for x in x_positions],
        vdw,
        width=width,
        yerr=vdw_err if any(vdw_err) else None,
        capsize=2.8,
        color=vdw_color,
        edgecolor="#874327",
        linewidth=0.8,
        label="VdW",
    )
    ax_b.axhline(0, color=ink, linewidth=0.8)
    ax_b.set_ylabel("Energy contribution (kcal/mol)", fontweight="bold")
    ax_b.set_xticks(x_positions)
    ax_b.set_xticklabels([curve.label for curve in curves])
    ax_b.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.04, 0.98), borderaxespad=0.0, fontsize=8.5)
    ax_b.grid(True, axis="y", color=grid, linewidth=0.8)
    ax_b.grid(False, axis="x")
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)
    ax_b.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))
    ax_b.tick_params(axis="both", colors=ink)

    for bars in (bars_elec, bars_vdw):
        for bar in bars:
            value = bar.get_height()
            ax_b.text(
                bar.get_x() + bar.get_width() / 2,
                value - 0.75 if value < 0 else value + 0.75,
                f"{value:.1f}",
                ha="center",
                va="top" if value < 0 else "bottom",
                fontsize=7.6,
                color=ink,
            )
    add_panel_label(ax_b, "B")

    fig.subplots_adjust(left=0.08, right=0.86, bottom=0.17, top=0.92)
    savefig_via_temp(fig, output_prefix.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    savefig_via_temp(fig, output_prefix.with_suffix(".tiff"), dpi=dpi, bbox_inches="tight")
    savefig_via_temp(fig, output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    write_summary(output_prefix.with_name(output_prefix.name + "_summary_kcal.csv"), curves, x_shifts)
    plt.close(fig)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create publication-style PMF/elec-vdW figure for amylin fib systems.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--pmf-xmax", type=float, default=2.0)
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path(r"D:\Data_MD\amylin\amylin_fib_binding_mechanism_600dpi"),
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    root = Path(r"D:\Data_MD\amylin")
    energy_dir = "titan_energy_md2_9_multi_xdr_amber14sb_ligtopo"
    curves = [
        load_curve(
            "4497_fib",
            Path(r"C:\Users\ADMIN\OneDrive\Documents\DTU_MsLinh\new_conf\US-4479\hist.xvg"),
            root / "4497_fib" / "data_MD" / energy_dir / "cal.txt",
        ),
        load_curve(
            "4631_fib",
            Path(r"C:\Users\ADMIN\OneDrive\Documents\DTU_MsLinh\new_conf\US-4631\scaled_hist.xvg"),
            root / "4631_fib" / "data_MD" / energy_dir / "cal.txt",
        ),
    ]
    make_figure(curves, args.output_prefix, args.dpi, args.pmf_xmax)
    for curve in curves:
        print(
            f"{curve.label}: PMF unit={curve.pmf_unit}, min_x={curve.min_x:.6g}, "
            f"elec={curve.elec_kcal:.4f} kcal/mol, vdw={curve.vdw_kcal:.4f} kcal/mol"
        )
    print(f"Wrote {args.output_prefix}.png/.tiff/.pdf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
