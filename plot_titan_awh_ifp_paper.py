#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
LOCAL_PYDEPS = SCRIPT_DIR / ".pydeps"
if LOCAL_PYDEPS.exists():
    sys.path.insert(0, str(LOCAL_PYDEPS))

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import colors
    from matplotlib.lines import Line2D
except ModuleNotFoundError as exc:
    raise SystemExit(
        "Missing plotting dependency. Please install matplotlib before running this script:\n"
        "  pip install matplotlib"
    ) from exc


AA_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

NA_ONE = {
    "A": "A",
    "ADE": "A",
    "DA": "A",
    "G": "G",
    "GUA": "G",
    "DG": "G",
    "C": "C",
    "CYT": "C",
    "DC": "C",
    "U": "U",
    "URA": "U",
    "DU": "U",
    "T": "T",
    "THY": "T",
    "DT": "T",
}

INTERACTION_DISPLAY = {
    "HBond": "H-bond",
    "SaltBridge": "Salt bridge",
    "PiPi": "Pi-pi",
    "CationPi": "Cation-pi",
    "Hydrophobic": "Hydrophobic",
    "Halogen": "Halogen",
}

INTERACTION_COLORS = {
    "HBond": "#1f78b4",
    "SaltBridge": "#d95f02",
    "PiPi": "#6a3d9a",
    "CationPi": "#b15928",
    "Hydrophobic": "#1b9e77",
    "Halogen": "#7570b3",
    "Total": "#111111",
}

WALKER_COLORS = {
    0: "#1b9e77",
    1: "#d95f02",
    2: "#7570b3",
    3: "#e7298a",
}


@dataclass
class BinSummary:
    bin_id: int
    rc_left: float
    rc_center: float
    rc_right: float
    frames: int
    replica_id: int
    time_ps: float | None
    rc_value: float | None
    delta_to_center: float | None


@dataclass
class ContactTurnover:
    machine_label: str
    display_label: str
    interaction_type: str
    values: np.ndarray


def read_text_lines(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8", errors="ignore").splitlines()


def load_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        return list(csv.DictReader(handle))


def set_publication_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 180,
            "savefig.dpi": 300,
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
            "font.size": 12,
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "axes.labelweight": "bold",
            "axes.titleweight": "bold",
            "axes.linewidth": 1.1,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "xtick.major.size": 5,
            "ytick.major.size": 5,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "legend.frameon": False,
            "legend.fontsize": 10,
            "grid.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def aa_code(resname: str) -> str:
    return AA_ONE.get(resname.upper(), resname[:1].upper())


def na_code(resname: str) -> str:
    return NA_ONE.get(resname.upper(), resname[:1].upper())


def pretty_residue_pair(rec_name: str, rec_id: int, lig_name: str, lig_id: int) -> str:
    return f"{aa_code(rec_name)}{rec_id}-{na_code(lig_name)}{lig_id}"


LABEL_RE = re.compile(
    r"(?P<chain_rec>[^:]+):(?P<resname_rec>[A-Za-z]+)(?P<resid_rec>\d+)-"
    r"(?P<chain_lig>[^:]+):(?P<resname_lig>[A-Za-z]+)(?P<resid_lig>\d+):"
    r"(?P<interaction>[A-Za-z0-9]+)"
)


def parse_contact_machine_label(label: str) -> tuple[str, str]:
    match = LABEL_RE.fullmatch(label.strip())
    if not match:
        return label, "Unknown"
    interaction = match.group("interaction")
    display = (
        f"{pretty_residue_pair(match.group('resname_rec'), int(match.group('resid_rec')), match.group('resname_lig'), int(match.group('resid_lig')))} "
        f"{INTERACTION_DISPLAY.get(interaction, interaction.lower())}"
    )
    return display, interaction


def load_bin_summary(path: Path) -> tuple[dict[str, float], list[BinSummary]]:
    metadata: dict[str, float] = {}
    entries: list[BinSummary] = []
    for raw in read_text_lines(path):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("# "):
            parts = line[2:].split(maxsplit=1)
            if len(parts) == 2:
                key, value = parts
                try:
                    metadata[key] = float(value)
                except ValueError:
                    continue
            continue
        if line.startswith("bin_id"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            continue

        def parse_optional(text: str) -> float | None:
            return None if text == "NA" else float(text)

        entries.append(
            BinSummary(
                bin_id=int(fields[0]),
                rc_left=float(fields[1]),
                rc_center=float(fields[2]),
                rc_right=float(fields[3]),
                frames=int(fields[4]),
                replica_id=int(fields[5]),
                time_ps=parse_optional(fields[6]),
                rc_value=parse_optional(fields[7]),
                delta_to_center=parse_optional(fields[8]),
            )
        )
    if not entries:
        raise ValueError(f"No usable RC bin summary rows found in {path}")
    return metadata, entries


def load_heatmap(path: Path) -> tuple[list[str], np.ndarray]:
    lines = read_text_lines(path)
    if not lines:
        raise ValueError(f"Empty heatmap file: {path}")
    header = lines[0].strip().split("\t")
    labels: list[str] = []
    rows: list[list[float]] = []
    for raw in lines[1:]:
        line = raw.strip()
        if not line:
            continue
        fields = line.split("\t")
        labels.append(fields[0])
        rows.append([float(x) for x in fields[1:]])
    if not rows:
        return [], np.zeros((0, len(header) - 1), dtype=float)
    return labels, np.asarray(rows, dtype=float)


def select_contacts(labels: list[str], matrix: np.ndarray, top_n: int, min_peak: float) -> list[ContactTurnover]:
    selected: list[ContactTurnover] = []
    for label, values in zip(labels, matrix):
        peak = float(np.nanmax(values)) if values.size else 0.0
        if peak < min_peak:
            continue
        display, interaction_type = parse_contact_machine_label(label)
        selected.append(
            ContactTurnover(
                machine_label=label,
                display_label=display,
                interaction_type=interaction_type,
                values=values,
            )
        )
    if not selected:
        return []

    def weighted_bin_index(values: np.ndarray) -> float:
        total = float(np.sum(values))
        if total <= 0.0:
            return math.inf
        coords = np.arange(values.size, dtype=float)
        return float(np.sum(coords * values) / total)

    selected.sort(
        key=lambda item: (
            weighted_bin_index(item.values),
            -float(np.nanmax(item.values)),
            item.display_label,
        )
    )
    return selected[:top_n]


def aggregate_interaction_burden(contacts: list[ContactTurnover], n_bins: int) -> dict[str, np.ndarray]:
    burden: dict[str, np.ndarray] = defaultdict(lambda: np.zeros(n_bins, dtype=float))
    for item in contacts:
        burden[item.interaction_type] = burden[item.interaction_type] + item.values
    return dict(burden)


def contact_zero_regions(total_signal: np.ndarray) -> list[tuple[int, int]]:
    ranges: list[tuple[int, int]] = []
    start: int | None = None
    for idx, value in enumerate(total_signal):
        if value <= 0.0 and start is None:
            start = idx
        elif value > 0.0 and start is not None:
            ranges.append((start, idx - 1))
            start = None
    if start is not None:
        ranges.append((start, len(total_signal) - 1))
    return ranges


def add_contact_free_shading(ax: plt.Axes, bins: list[BinSummary], total_signal: np.ndarray) -> None:
    for start_idx, end_idx in contact_zero_regions(total_signal):
        left = bins[start_idx].rc_left
        right = bins[end_idx].rc_right
        ax.axvspan(left, right, color="#efe7d2", alpha=0.28, lw=0)


def build_turnover_summary(
    out_path: Path, bins: list[BinSummary], contacts: list[ContactTurnover]
) -> None:
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "contact",
                "interaction_type",
                "first_bin",
                "last_bin",
                "peak_bin",
                "peak_rc_nm",
                "peak_occupancy_percent",
            ]
        )
        for item in contacts:
            nz = np.flatnonzero(item.values > 0.0)
            if nz.size == 0:
                continue
            peak_idx = int(np.nanargmax(item.values))
            writer.writerow(
                [
                    item.display_label,
                    item.interaction_type,
                    int(nz[0]),
                    int(nz[-1]),
                    peak_idx,
                    f"{bins[peak_idx].rc_center:.4f}",
                    f"{float(np.nanmax(item.values)):.4f}",
                ]
            )


def plot_overview(
    out_path_base: Path,
    bins: list[BinSummary],
    all_contacts: list[ContactTurnover],
) -> None:
    centers = np.array([entry.rc_center for entry in bins], dtype=float)
    widths = np.array([entry.rc_right - entry.rc_left for entry in bins], dtype=float)
    frames = np.array([entry.frames for entry in bins], dtype=float)
    walker_ids = [entry.replica_id for entry in bins]
    total_signal = (
        np.sum(np.vstack([item.values for item in all_contacts]), axis=0)
        if all_contacts
        else np.zeros(len(bins), dtype=float)
    )
    fig, ax1 = plt.subplots(figsize=(8.8, 3.8))

    add_contact_free_shading(ax1, bins, total_signal)
    ax1.bar(
        centers,
        frames,
        width=widths * 0.92,
        color="#c9ced6",
        edgecolor="#3d4752",
        linewidth=1.0,
        zorder=2,
    )
    ax1.set_xlabel("Reaction coordinate (nm)", fontweight="bold")
    ax1.set_ylabel("Frames / bin", fontweight="bold")
    ax1.set_title("AWH RC-Binned IFP Overview", loc="left", pad=10, fontweight="bold")
    ax1.grid(axis="y", color="#d0d0d0", alpha=0.55)

    max_frames = float(np.max(frames)) if frames.size else 1.0
    scatter_y = frames + max_frames * 0.045
    scatter_colors = [WALKER_COLORS.get(replica_id, "#666666") for replica_id in walker_ids]
    ax1.scatter(
        centers,
        scatter_y,
        c=scatter_colors,
        marker="s",
        s=46,
        edgecolors="white",
        linewidths=0.7,
        zorder=3,
        clip_on=False,
    )
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="s",
            color="none",
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=0.7,
            markersize=7.8,
            label=f"Rep walker {walker_id}",
        )
        for walker_id, color in sorted(WALKER_COLORS.items())
    ]
    ax1.legend(handles=legend_handles, ncol=4, loc="upper right", handletextpad=0.3, columnspacing=0.9)
    ax1.set_xlim(bins[0].rc_left, bins[-1].rc_right)

    fig.savefig(out_path_base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_path_base.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig)


def plot_interaction_burden(
    out_path_base: Path,
    bins: list[BinSummary],
    all_contacts: list[ContactTurnover],
) -> None:
    centers = np.array([entry.rc_center for entry in bins], dtype=float)
    total_signal = (
        np.sum(np.vstack([item.values for item in all_contacts]), axis=0)
        if all_contacts
        else np.zeros(len(bins), dtype=float)
    )
    interaction_burden = aggregate_interaction_burden(all_contacts, len(bins))

    fig, ax = plt.subplots(figsize=(8.8, 4.6))

    add_contact_free_shading(ax, bins, total_signal)
    plotted_any = False
    for interaction_name in ["SaltBridge", "HBond", "Hydrophobic", "PiPi", "CationPi", "Halogen"]:
        values = interaction_burden.get(interaction_name)
        if values is None or np.nanmax(values) <= 0.0:
            continue
        plotted_any = True
        ax.plot(
            centers,
            values,
            color=INTERACTION_COLORS.get(interaction_name, "#555555"),
            lw=2.6,
            marker="o",
            ms=4.4,
            label=INTERACTION_DISPLAY.get(interaction_name, interaction_name),
        )
    ax.set_xlabel("Reaction coordinate (nm)", fontweight="bold")
    ax.set_ylabel("Summed occupancy (%)", fontweight="bold")
    ax.set_title("AWH Interaction-Class Turnover Along the Reaction Coordinate", loc="left", pad=10, fontweight="bold")
    ax.grid(axis="y", color="#d0d0d0", alpha=0.55)
    if plotted_any:
        ax.legend(loc="upper right", ncol=2, handlelength=2.4)

    fig.savefig(out_path_base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_path_base.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig)


def plot_contact_heatmap(
    out_path_base: Path,
    bins: list[BinSummary],
    selected_contacts: list[ContactTurnover],
    all_contacts: list[ContactTurnover],
) -> None:
    total_signal = (
        np.sum(np.vstack([item.values for item in all_contacts]), axis=0)
        if all_contacts
        else np.zeros(len(bins), dtype=float)
    )
    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    add_contact_free_shading(ax, bins, total_signal)

    heatmap_data = (
        np.vstack([item.values for item in selected_contacts])
        if selected_contacts
        else np.zeros((0, len(bins)), dtype=float)
    )
    if heatmap_data.size:
        vmax = float(np.nanmax(heatmap_data))
        vmax = vmax if vmax > 0.0 else 1.0
        norm = colors.PowerNorm(gamma=0.7, vmin=0.0, vmax=vmax)
        im = ax.imshow(
            heatmap_data,
            aspect="auto",
            interpolation="nearest",
            cmap="YlOrRd",
            norm=norm,
            extent=(bins[0].rc_left, bins[-1].rc_right, heatmap_data.shape[0], 0),
        )
        cbar = fig.colorbar(im, ax=ax, pad=0.015, fraction=0.05)
        cbar.set_label("Occupancy", fontweight="bold")
        yticks = np.arange(heatmap_data.shape[0]) + 0.5
        ax.set_yticks(yticks)
        ax.set_yticklabels([item.display_label for item in selected_contacts])
    else:
        ax.text(
            0.5,
            0.5,
            "No contacts passed the current plotting threshold.",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_yticks([])
    ax.set_xlabel("Reaction coordinate (nm)", fontweight="bold")
    ax.set_ylabel("Contact", fontweight="bold")
    ax.set_title("AWH Contact Heatmap Along the Reaction Coordinate", loc="left", pad=10, fontweight="bold")
    ax.set_xlim(bins[0].rc_left, bins[-1].rc_right)
    display_tick_positions: list[float] = []
    display_tick_labels: list[str] = []
    display_idx = 1
    while True:
        tick_position = bins[0].rc_left + float(display_idx - 1)
        if tick_position > bins[-1].rc_right + 1e-9:
            break
        display_tick_positions.append(tick_position)
        display_tick_labels.append(str(display_idx))
        display_idx += 1
    if display_tick_positions:
        ax.set_xticks(display_tick_positions)
        ax.set_xticklabels(display_tick_labels)

    fig.savefig(out_path_base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_path_base.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig)


def plot_contact_turnover(
    out_path_base: Path,
    bins: list[BinSummary],
    selected_contacts: list[ContactTurnover],
) -> None:
    centers = np.array([entry.rc_center for entry in bins], dtype=float)
    fig, ax = plt.subplots(figsize=(8.8, 5.0))
    ax.set_title("Top AWH Contacts Along the Reaction Coordinate", loc="left", pad=10, fontweight="bold")

    if not selected_contacts:
        ax.text(0.5, 0.5, "No contacts passed the current plotting threshold.", ha="center", va="center")
        ax.set_axis_off()
    else:
        ordered = sorted(
            selected_contacts,
            key=lambda item: (-float(np.nanmax(item.values)), item.display_label),
        )
        color_cycle = plt.get_cmap("tab10")(np.linspace(0.0, 1.0, max(3, len(ordered))))
        for color, item in zip(color_cycle, ordered):
            ax.plot(
                centers,
                item.values,
                lw=2.7,
                marker="o",
                ms=4.6,
                color=color,
                label=item.display_label,
            )
        ax.set_xlabel("Reaction coordinate (nm)", fontweight="bold")
        ax.set_ylabel("Occupancy (%)", fontweight="bold")
        ax.grid(color="#d0d0d0", alpha=0.55)
        ax.legend(loc="upper right", ncol=1, handlelength=2.0)

    fig.savefig(out_path_base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_path_base.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig)


def print_quick_summary(bins: list[BinSummary], contacts: list[ContactTurnover]) -> None:
    total_frames = sum(entry.frames for entry in bins)
    print(f"Loaded {len(bins)} RC bins spanning {bins[0].rc_left:.3f} to {bins[-1].rc_right:.3f} nm.")
    print(f"Matched frames across bins: {total_frames}.")
    if not contacts:
        print("No contacts with occupancy above the plotting threshold were found.")
        return
    ordered = sorted(contacts, key=lambda item: -float(np.nanmax(item.values)))
    print("Top contacts by peak occupancy:")
    for item in ordered[: min(6, len(ordered))]:
        peak_idx = int(np.nanargmax(item.values))
        peak_value = float(np.nanmax(item.values))
        print(
            f"  - {item.display_label}: peak {peak_value:.3f}% at bin {peak_idx} (RC {bins[peak_idx].rc_center:.3f} nm)"
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate paper-style figures for TITAN AWH RC-binned IFP outputs."
    )
    parser.add_argument(
        "run_dir",
        nargs="?",
        default=".",
        help="Directory containing awh_rc_bin_summary.dat and awh_ifp_heatmap.dat",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=12,
        help="Maximum number of contacts to show in the heatmap and turnover plot.",
    )
    parser.add_argument(
        "--min-peak",
        type=float,
        default=0.0,
        help="Minimum peak occupancy (%%) required to plot a contact.",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="Optional output directory. Defaults to <run_dir>/paper_figures",
    )
    parser.add_argument(
        "--prefix",
        default="awh_ifp_paper",
        help="Output file prefix for the generated figures.",
    )
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    outdir = Path(args.outdir).resolve() if args.outdir else run_dir / "paper_figures"
    outdir.mkdir(parents=True, exist_ok=True)

    summary_path = run_dir / "awh_rc_bin_summary.dat"
    heatmap_path = run_dir / "awh_ifp_heatmap.dat"
    if not summary_path.exists():
        raise SystemExit(f"Missing required file: {summary_path}")
    if not heatmap_path.exists():
        raise SystemExit(f"Missing required file: {heatmap_path}")

    set_publication_style()
    _, bins = load_bin_summary(summary_path)
    heatmap_labels, heatmap_matrix = load_heatmap(heatmap_path)
    all_contacts = select_contacts(
        heatmap_labels,
        heatmap_matrix,
        top_n=max(len(heatmap_labels), args.top_n),
        min_peak=args.min_peak,
    )
    selected_contacts = all_contacts[: args.top_n]

    build_turnover_summary(outdir / f"{args.prefix}_contact_summary.tsv", bins, all_contacts)
    plot_overview(outdir / f"{args.prefix}_overview", bins, all_contacts)
    plot_interaction_burden(outdir / f"{args.prefix}_interaction_classes", bins, all_contacts)
    plot_contact_heatmap(outdir / f"{args.prefix}_contact_heatmap", bins, selected_contacts, all_contacts)
    plot_contact_turnover(outdir / f"{args.prefix}_top_contacts", bins, selected_contacts)
    print_quick_summary(bins, all_contacts)
    print(f"Saved figures to {outdir}")


if __name__ == "__main__":
    main()
