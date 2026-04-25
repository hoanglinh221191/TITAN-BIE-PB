#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    from matplotlib import colors
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
    "ACE": "Ac",
    "NME": "Nm",
    "CLYS": "K",
}

NA_ONE = {
    "A": "A",
    "ADE": "A",
    "DA": "A",
    "G": "G",
    "GUA": "G",
    "DG": "G",
    "G5": "G",
    "C": "C",
    "CYT": "C",
    "DC": "C",
    "C5": "C",
    "U": "U",
    "URA": "U",
    "DT": "T",
    "T": "T",
    "THY": "T",
    "U5": "U",
}

INTERACTION_DISPLAY = {
    "HBond": "H Bond",
    "SaltBridge": "Salt Bridge",
    "PiPi": "Pi-Pi",
    "CationPi": "Cation-Pi",
    "Hydrophobic": "Hydrophobic",
    "Halogen": "Halogen",
}

PAPER_FACE = "#fbfbf8"
HEATMAP_FACE = "#fff6cc"
GRID_COLOR = "#d5d5d5"


def read_lines(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8", errors="ignore").splitlines()


def load_csv_dicts(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        return list(csv.DictReader(handle))


def load_whitespace_table(path: Path) -> list[dict[str, str]]:
    header: list[str] | None = None
    rows: list[dict[str, str]] = []
    for raw in read_lines(path):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = re.split(r"\s+", line)
        if header is None:
            header = parts
            continue
        if len(parts) < len(header):
            parts = parts + [""] * (len(header) - len(parts))
        rows.append(dict(zip(header, parts[: len(header)])))
    return rows


def numeric_rows(path: Path, comment_prefix: str = "#") -> list[list[str]]:
    rows: list[list[str]] = []
    for raw in read_lines(path):
        line = raw.strip()
        if not line or line.startswith(comment_prefix):
            continue
        rows.append(line.split())
    return rows


def biomol_kind(resname: str) -> str:
    resname = resname.upper()
    if resname in AA_ONE:
        return "protein"
    if resname in NA_ONE:
        return "nucleic"
    return "other"


def format_residue_label(chain: str, resnr: int, resname: str, compact: bool = True) -> str:
    resname_u = resname.upper()
    kind = biomol_kind(resname_u)
    if compact:
        if kind == "protein":
            code = AA_ONE.get(resname_u, resname_u[:1])
            return f"P:{code}{resnr}"
        if kind == "nucleic":
            code = NA_ONE.get(resname_u, resname_u[:1])
            return f"R:{code}{resnr}"
    return f"{chain}:{resname_u}{resnr}"


def short_residue_token(resname: str, resnr: int) -> str:
    resname_u = resname.upper()
    if resname_u in AA_ONE:
        code = AA_ONE[resname_u]
    elif resname_u in NA_ONE:
        code = NA_ONE[resname_u]
    else:
        code = resname_u
    return f"{code}{resnr}"


def pretty_interaction_name(name: str) -> str:
    if not name:
        return name
    if name in INTERACTION_DISPLAY:
        return INTERACTION_DISPLAY[name]
    spaced = re.sub(r"(?<!^)([A-Z])", r" \1", name).strip()
    return spaced.title()


def pretty_ifp_label(raw_label: str) -> str:
    label = raw_label.strip()
    match = re.match(
        r"^(?P<chain1>[^:]+):(?P<res1>[A-Za-z0-9']+?)(?P<num1>-?\d+)-"
        r"(?P<chain2>[^:]+):(?P<res2>[A-Za-z0-9']+?)(?P<num2>-?\d+):"
        r"(?P<interaction>[A-Za-z0-9]+)$",
        label,
    )
    if not match:
        return label

    left = short_residue_token(match.group("res1"), int(match.group("num1")))
    right = short_residue_token(match.group("res2"), int(match.group("num2")))
    interaction = pretty_interaction_name(match.group("interaction"))
    return f"{left}-{right}: {interaction}"


def style_contact_heatmap_axes(
    ax: plt.Axes,
    labels: list[str],
    *,
    xlabel: str,
    x_is_rc: bool,
    rc_left: float | None = None,
    rc_right: float | None = None,
) -> None:
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels([pretty_ifp_label(label) for label in labels], fontsize=15, fontweight="bold")

    if x_is_rc and rc_left is not None and rc_right is not None:
        set_heatmap_rc_axis(ax, rc_left, rc_right)
    else:
        tick_step = max(1, len(labels) // 12)
        xticks = np.arange(0, len(labels), tick_step)
        ax.set_xticks(xticks)
        ax.set_xticklabels([str(value) for value in xticks])

    ax.set_xlabel(xlabel, fontsize=18, fontweight="bold")
    ax.tick_params(axis="x", labelsize=14, width=1.2, length=6)
    ax.tick_params(axis="y", labelsize=15, width=1.2, length=5)

    for tick in ax.get_xticklabels():
        tick.set_fontweight("bold")
    for tick in ax.get_yticklabels():
        tick.set_fontweight("bold")


def parse_pdb_chain_residues(path: Path) -> dict[str, list[tuple[int, str]]]:
    chains: dict[str, list[tuple[int, str]]] = {}
    seen: set[tuple[str, int, str]] = set()
    for raw in read_lines(path):
        if not raw.startswith(("ATOM", "HETATM")):
            continue
        chain = raw[21].strip() or "?"
        try:
            resnr = int(raw[22:26])
        except ValueError:
            continue
        resname = raw[17:20].strip()
        key = (chain, resnr, resname)
        if key in seen:
            continue
        seen.add(key)
        chains.setdefault(chain, []).append((resnr, resname))
    return chains


def infer_contact_map_labels(outdir: Path, shape: tuple[int, int]) -> tuple[list[str], list[str]]:
    pdb_path = outdir / "t.pdb"
    if not pdb_path.exists():
        fallback = outdir / "awh_rep_bin_000.pdb"
        if fallback.exists():
            pdb_path = fallback
        else:
            rep_candidates = sorted(outdir.glob("awh_rep_bin_*.pdb"))
            if rep_candidates:
                pdb_path = rep_candidates[0]
    if not pdb_path.exists():
        return (
            [str(i + 1) for i in range(shape[0])],
            [str(i + 1) for i in range(shape[1])],
        )
    chains = parse_pdb_chain_residues(pdb_path)
    if len(chains) < 2:
        return (
            [str(i + 1) for i in range(shape[0])],
            [str(i + 1) for i in range(shape[1])],
        )
    chain_names = list(chains.keys())
    mol1_chain = chain_names[0]
    mol2_chain = chain_names[1]
    row_labels = [
        format_residue_label(mol1_chain, resnr, resname, compact=True)
        for resnr, resname in chains[mol1_chain][: shape[0]]
    ]
    col_labels = [
        format_residue_label(mol2_chain, resnr, resname, compact=True)
        for resnr, resname in chains[mol2_chain][: shape[1]]
    ]
    if len(row_labels) < shape[0]:
        row_labels.extend(str(i + 1) for i in range(len(row_labels), shape[0]))
    if len(col_labels) < shape[1]:
        col_labels.extend(str(i + 1) for i in range(len(col_labels), shape[1]))
    return row_labels, col_labels


def load_rmsd(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows = numeric_rows(path)
    data = np.array([[float(row[1]), float(row[2])] for row in rows], dtype=float)
    return data[:, 0], data[:, 1]


def load_rmsf(path: Path) -> tuple[list[str], np.ndarray]:
    labels: list[str] = []
    values: list[float] = []
    for raw in read_lines(path):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 6 or not parts[0].isdigit():
            continue
        labels.append(f"{parts[1]}:{parts[2]}:{parts[3]}")
        values.append(float(parts[5]))
    return labels, np.array(values, dtype=float)


def load_stat_hb(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = numeric_rows(path)
    time_ps = np.array([float(row[0]) for row in rows], dtype=float)
    mean_hb = np.array([float(row[1]) for row in rows], dtype=float)
    sd_hb = np.array([float(row[2]) for row in rows], dtype=float)
    return time_ps, mean_hb, sd_hb


def load_stat_average_hb(path: Path) -> tuple[list[str], np.ndarray]:
    labels: list[str] = []
    means: list[float] = []
    for raw in read_lines(path):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = re.split(r"\s+", line)
        if len(parts) >= 4:
            labels.append(parts[0])
            means.append(float(parts[1]))
    return labels, np.array(means, dtype=float)


def load_sc_time(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = numeric_rows(path)
    time_ps = np.array([float(row[0]) for row in rows], dtype=float)
    mean_sc = np.array([float(row[1]) for row in rows], dtype=float)
    sd_sc = np.array([float(row[2]) for row in rows], dtype=float)
    return time_ps, mean_sc, sd_sc


def load_contact_map(path: Path) -> np.ndarray:
    matrix: list[list[float]] = []
    for idx, raw in enumerate(read_lines(path)):
        line = raw.strip()
        if not line:
            continue
        parts = line.split()
        if idx == 0 and parts and parts[0].startswith("#"):
            continue
        matrix.append([float(value) for value in parts[1:]])
    return np.array(matrix, dtype=float)


def break_large_time_gaps(
    time_ps: np.ndarray, *series: np.ndarray, gap_factor: float = 50.0, min_gap_ps: float = 10000.0
) -> tuple[np.ndarray, ...]:
    if time_ps.size < 3:
        return (time_ps, *series)
    diffs = np.diff(time_ps)
    positive = diffs[diffs > 0]
    if positive.size == 0:
        return (time_ps, *series)
    typical = float(np.median(positive))
    if typical <= 0.0:
        return (time_ps, *series)

    gap_idx = np.where(diffs > max(min_gap_ps, gap_factor * typical))[0]
    if gap_idx.size == 0:
        return (time_ps, *series)

    gap_set = set(gap_idx.tolist())
    out_time: list[float] = []
    out_series: list[list[float]] = [[] for _ in series]
    for idx, value in enumerate(time_ps):
        out_time.append(float(value))
        for sidx, arr in enumerate(series):
            out_series[sidx].append(float(arr[idx]))
        if idx in gap_set:
            out_time.append(np.nan)
            for sidx in range(len(series)):
                out_series[sidx].append(np.nan)
    return (np.array(out_time, dtype=float), *[np.array(values, dtype=float) for values in out_series])


def savefig(fig: plt.Figure, path: Path, dpi: int) -> None:
    fig.patch.set_facecolor("white")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def style_axis(ax: plt.Axes, *, grid: bool = False, facecolor: str = PAPER_FACE) -> None:
    if facecolor != "none":
        ax.set_facecolor(facecolor)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.15)
    ax.spines["bottom"].set_linewidth(1.15)
    ax.tick_params(axis="both", direction="out", width=1.0, length=5)
    if grid:
        ax.grid(True, color=GRID_COLOR, linewidth=0.7, alpha=0.7)
        ax.set_axisbelow(True)


def filter_label_rows(rows: list[dict[str, str]], excluded_labels: set[str]) -> list[dict[str, str]]:
    if not excluded_labels:
        return rows
    return [row for row in rows if row.get("label", "") not in excluded_labels]


def rc_edges_from_centers(rc_axis: np.ndarray) -> np.ndarray:
    rc_axis = np.asarray(rc_axis, dtype=float)
    if rc_axis.size == 0:
        return np.array([0.0, 1.0], dtype=float)
    if rc_axis.size == 1:
        delta = 0.1
        return np.array([rc_axis[0] - delta, rc_axis[0] + delta], dtype=float)
    midpoints = 0.5 * (rc_axis[:-1] + rc_axis[1:])
    edges = np.empty(rc_axis.size + 1, dtype=float)
    edges[1:-1] = midpoints
    edges[0] = rc_axis[0] - 0.5 * (rc_axis[1] - rc_axis[0])
    edges[-1] = rc_axis[-1] + 0.5 * (rc_axis[-1] - rc_axis[-2])
    return edges


def set_heatmap_rc_axis(
    ax: plt.Axes,
    left: float,
    right: float,
    *,
    rc_min: float | None = None,
    rc_max: float | None = None,
) -> None:
    x_left = rc_min if rc_min is not None else left
    x_right = rc_max if rc_max is not None else right
    ax.set_xlim(x_left, x_right)

    tick_start = float(np.ceil(x_left))
    tick_end = float(np.floor(x_right))
    if tick_end >= tick_start:
        ticks = np.arange(tick_start, tick_end + 0.001, 1.0)
    else:
        ticks = np.linspace(x_left, x_right, num=min(6, max(2, int(np.ceil(x_right - x_left)) + 1)))
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))


def plot_rc_coverage(figdir: Path, coverage_rows: list[dict[str, str]], dpi: int) -> None:
    rc = np.array([float(row["rc_center"]) for row in coverage_rows], dtype=float)
    frames = np.array([float(row["frames"]) for row in coverage_rows], dtype=float)
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 7), sharex=True)
    style_axis(axes[0], grid=True)
    style_axis(axes[1], grid=True)
    width = np.median(np.diff(rc)) * 0.82 if rc.size > 1 else 0.12
    axes[0].bar(rc, frames, width=width, color="#d7a94b", edgecolor="black", linewidth=0.6)
    axes[0].set_ylabel("Frames per RC bin")
    if "active_walkers" in coverage_rows[0] and "dominant_fraction" in coverage_rows[0]:
        active = np.array([float(row["active_walkers"]) for row in coverage_rows], dtype=float)
        dominant = np.array([float(row["dominant_fraction"]) for row in coverage_rows], dtype=float)
        ax0r = axes[0].twinx()
        style_axis(ax0r, facecolor="none")
        ax0r.plot(rc, active, color="#0f4c81", marker="o", linewidth=1.8, label="Active walkers")
        ax0r.set_ylabel("Active walkers")
        axes[0].set_title("AWH RC coverage and walker support")

        axes[1].plot(rc, dominant, color="#8b1e3f", marker="o", linewidth=1.8)
        axes[1].fill_between(rc, 0.0, dominant, color="#8b1e3f", alpha=0.15)
        axes[1].axhline(0.5, color="0.35", linestyle="--", linewidth=1.0)
        axes[1].set_ylabel("Dominant walker fraction")
        axes[1].set_ylim(0.0, 1.0)
        axes[1].set_title("Per-bin walker dominance")
    else:
        replica = np.array([float(row.get("replica_id", 0.0)) for row in coverage_rows], dtype=float)
        delta = np.array([abs(float(row.get("delta_to_center", 0.0))) for row in coverage_rows], dtype=float)
        axes[0].set_title("AWH RC coverage")
        axes[1].plot(rc, replica, color="#0f4c81", marker="o", linewidth=1.6)
        ax1r = axes[1].twinx()
        style_axis(ax1r, facecolor="none")
        ax1r.plot(rc, delta, color="#8b1e3f", linestyle="--", linewidth=1.4)
        axes[1].set_ylabel("Representative walker")
        ax1r.set_ylabel("|Delta to bin center|")
        axes[1].set_title("Representative walker and bin-center match")
    axes[1].set_xlabel("Reaction coordinate")
    savefig(fig, figdir / "01_awh_rc_coverage.png", dpi)


def plot_walker_stability(figdir: Path, stability_rows: list[dict[str, str]], dpi: int) -> None:
    rc = np.array([float(row["rc_center"]) for row in stability_rows], dtype=float)
    walker_cols = [key for key in stability_rows[0].keys() if key.startswith("walker_") and key.endswith("_fraction")]
    matrix = np.array([[float(row[key]) for key in walker_cols] for row in stability_rows], dtype=float).T

    fig, ax = plt.subplots(figsize=(11, 6))
    style_axis(ax, facecolor=HEATMAP_FACE)
    image = ax.imshow(matrix, aspect="auto", cmap="YlOrRd", origin="lower", vmin=0.0, vmax=max(0.05, float(matrix.max())))
    fig.colorbar(image, ax=ax, label="Fraction of bin frames")
    ax.set_yticks(np.arange(len(walker_cols)))
    ax.set_yticklabels([key.replace("_fraction", "") for key in walker_cols], fontsize=8)
    xticks = np.arange(len(rc))
    tick_step = max(1, len(rc) // 12)
    ax.set_xticks(xticks[::tick_step])
    ax.set_xticklabels([f"{value:.2f}" for value in rc[::tick_step]], rotation=45, ha="right")
    ax.set_xlabel("RC bin center")
    ax.set_ylabel("Trajectory segment")
    ax.set_title("Walker stability across AWH bins")
    savefig(fig, figdir / "02_awh_walker_stability_heatmap.png", dpi)


def plot_pair_overlap(figdir: Path, pair_rows: list[dict[str, str]], dpi: int) -> None:
    n = 1 + max(max(int(row["walker_i"]), int(row["walker_j"])) for row in pair_rows)
    jaccard = np.eye(n, dtype=float)
    frame_overlap = np.eye(n, dtype=float)
    for row in pair_rows:
        i = int(row["walker_i"])
        j = int(row["walker_j"])
        jaccard[i, j] = jaccard[j, i] = float(row["bin_jaccard"])
        frame_overlap[i, j] = frame_overlap[j, i] = float(row["frame_overlap_fraction"])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    style_axis(axes[0], facecolor=HEATMAP_FACE)
    style_axis(axes[1], facecolor=HEATMAP_FACE)
    im0 = axes[0].imshow(jaccard, cmap="viridis", vmin=0.0, vmax=1.0)
    im1 = axes[1].imshow(frame_overlap, cmap="magma", vmin=0.0, vmax=max(0.1, float(frame_overlap.max())))
    fig.colorbar(im0, ax=axes[0], label="Bin Jaccard")
    fig.colorbar(im1, ax=axes[1], label="Frame overlap fraction")
    for ax, title in zip(axes, ["Walker bin overlap", "Walker frame overlap"]):
        ax.set_xlabel("Trajectory segment")
        ax.set_ylabel("Trajectory segment")
        ax.set_title(title)
    savefig(fig, figdir / "03_awh_pair_overlap.png", dpi)


def plot_structural_summary(figdir: Path, summary_rows: list[dict[str, str]], dpi: int) -> None:
    rc = np.array([float(row["rc_center"]) for row in summary_rows], dtype=float)
    hb = np.array([float(row["HBond_avg_hits_per_frame"]) for row in summary_rows], dtype=float)
    cpi = np.array([float(row["CationPi_avg_hits_per_frame"]) for row in summary_rows], dtype=float)
    hyd = np.array([float(row["Hydrophobic_avg_hits_per_frame"]) for row in summary_rows], dtype=float)
    load = np.array([float(row["contact_event_load_per_frame"]) for row in summary_rows], dtype=float)
    active = np.array([float(row["active_walkers"]) for row in summary_rows], dtype=float)
    dominant = np.array([float(row["dominant_fraction"]) for row in summary_rows], dtype=float)

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    style_axis(axes[0], grid=True)
    style_axis(axes[1], grid=True)
    axes[0].plot(rc, hb, color="#b22222", marker="o", linewidth=1.8, label="HBond")
    axes[0].plot(rc, cpi, color="#2b6cb0", marker="s", linewidth=1.8, label="Cation-Pi")
    axes[0].plot(rc, hyd, color="#c58b00", marker="^", linewidth=1.8, label="Hydrophobic")
    axes[0].plot(rc, load, color="black", linestyle="--", linewidth=1.4, label="Total contact load")
    axes[0].set_ylabel("Average hits per frame")
    axes[0].set_title("AWH structural summary by RC bin")
    axes[0].legend(loc="best", fontsize=9)

    axes[1].bar(rc, active, width=(np.median(np.diff(rc)) * 0.8 if rc.size > 1 else 0.12), color="#8fc1e3", edgecolor="black", linewidth=0.5, label="Active walkers")
    ax1r = axes[1].twinx()
    style_axis(ax1r, facecolor="none")
    ax1r.plot(rc, dominant, color="#7f1734", marker="o", linewidth=1.6, label="Dominant fraction")
    axes[1].set_xlabel("Reaction coordinate")
    axes[1].set_ylabel("Active walkers")
    ax1r.set_ylabel("Dominant walker fraction")
    axes[1].set_title("Sampling support per RC bin")
    savefig(fig, figdir / "04_awh_structural_summary.png", dpi)


def aggregate_ifp_by_bin(ifp_rows: list[dict[str, str]]) -> list[dict[str, float]]:
    bins: dict[int, dict[str, float]] = {}
    interaction_keys = [
        "HBond(%)",
        "SaltBridge(%)",
        "PiPi(%)",
        "CationPi(%)",
        "Hydrophobic(%)",
        "Halogen(%)",
    ]
    for row in ifp_rows:
        bin_id = int(row["Bin"])
        entry = bins.setdefault(
            bin_id,
            {
                "rc_center": float(row["RC_Center"]),
                "frames": float(row["Frames"]),
                "HBond": 0.0,
                "SaltBridge": 0.0,
                "PiPi": 0.0,
                "CationPi": 0.0,
                "Hydrophobic": 0.0,
                "Halogen": 0.0,
                "active_pairs": 0.0,
            },
        )
        pair_total = 0.0
        for key in interaction_keys:
            value = float(row[key])
            entry[key.replace("(%)", "")] += value
            pair_total += value
        if pair_total > 0.0:
            entry["active_pairs"] += 1.0
    return [bins[idx] for idx in sorted(bins)]


def plot_structural_summary_legacy(figdir: Path, ifp_rows: list[dict[str, str]], dpi: int) -> None:
    summary = aggregate_ifp_by_bin(ifp_rows)
    rc = np.array([row["rc_center"] for row in summary], dtype=float)
    hb = np.array([row["HBond"] for row in summary], dtype=float)
    salt = np.array([row["SaltBridge"] for row in summary], dtype=float)
    cpi = np.array([row["CationPi"] for row in summary], dtype=float)
    hyd = np.array([row["Hydrophobic"] for row in summary], dtype=float)
    total = hb + salt + cpi + hyd + np.array([row["PiPi"] + row["Halogen"] for row in summary], dtype=float)
    active_pairs = np.array([row["active_pairs"] for row in summary], dtype=float)

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    style_axis(axes[0], grid=True)
    style_axis(axes[1], grid=True)
    axes[0].plot(rc, hb, color="#b22222", marker="o", linewidth=1.8, label="HBond")
    axes[0].plot(rc, salt, color="#5b8bd0", marker="s", linewidth=1.8, label="SaltBridge")
    axes[0].plot(rc, cpi, color="#2b6cb0", marker="D", linewidth=1.8, label="Cation-Pi")
    axes[0].plot(rc, hyd, color="#c58b00", marker="^", linewidth=1.8, label="Hydrophobic")
    axes[0].plot(rc, total, color="black", linestyle="--", linewidth=1.4, label="Total IFP load")
    axes[0].set_ylabel("Summed occupancy (%)")
    axes[0].set_title("AWH interaction summary by RC bin")
    axes[0].legend(loc="best", fontsize=9)

    width = np.median(np.diff(rc)) * 0.8 if rc.size > 1 else 0.12
    axes[1].bar(rc, active_pairs, width=width, color="#8fc1e3", edgecolor="black", linewidth=0.5)
    axes[1].set_xlabel("Reaction coordinate")
    axes[1].set_ylabel("Active residue pairs")
    axes[1].set_title("Number of nonzero residue-pair contacts per RC bin")
    savefig(fig, figdir / "04_awh_structural_summary.png", dpi)


def plot_top_contacts(
    figdir: Path,
    ifp_rows: list[dict[str, str]],
    rc_axis: np.ndarray | None,
    dpi: int,
    *,
    top_n: int = 12,
    rc_min: float | None = None,
    rc_max: float | None = None,
) -> None:
    if not ifp_rows:
        return
    label_to_values: dict[str, dict[int, float]] = {}
    for row in ifp_rows:
        label = row["label"]
        bin_id = int(row["bin_id"])
        label_to_values.setdefault(label, {})[bin_id] = float(row["occupancy_pct"])

    top_labels = sorted(
        label_to_values.keys(),
        key=lambda label: max(label_to_values[label].values(), default=0.0),
        reverse=True,
    )[:top_n]
    max_bin = 1 + max(int(row["bin_id"]) for row in ifp_rows)
    matrix_width = int(rc_axis.size) if rc_axis is not None and rc_axis.size else max_bin
    matrix = np.zeros((len(top_labels), matrix_width), dtype=float)
    for row_idx, label in enumerate(top_labels):
        for bin_id, value in label_to_values[label].items():
            matrix[row_idx, bin_id] = value

    fig, ax = plt.subplots(figsize=(12.8, max(5.8, 0.62 * len(top_labels))))
    style_axis(ax, facecolor=HEATMAP_FACE)
    if rc_axis is not None and rc_axis.size >= matrix_width:
        rc_used = np.asarray(rc_axis[:matrix_width], dtype=float)
        edges = rc_edges_from_centers(rc_used)
        image = ax.imshow(
            matrix,
            aspect="auto",
            cmap="YlOrRd",
            origin="lower",
            interpolation="nearest",
            extent=[edges[0], edges[-1], -0.5, len(top_labels) - 0.5],
            vmin=0.0,
            vmax=max(1.0, float(matrix.max())),
        )
        rc_left = float(edges[0])
        rc_right = float(edges[-1])
    else:
        image = ax.imshow(
            matrix,
            aspect="auto",
            cmap="YlOrRd",
            origin="lower",
            interpolation="nearest",
            vmin=0.0,
            vmax=max(1.0, float(matrix.max())),
        )
        rc_left = None
        rc_right = None
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("Occupancy (%)", fontsize=14, fontweight="bold")
    cbar.ax.tick_params(labelsize=12, width=1.0, length=4)
    for tick in cbar.ax.get_yticklabels():
        tick.set_fontweight("bold")
    if rc_axis is not None and rc_axis.size >= matrix_width:
        style_contact_heatmap_axes(
            ax,
            top_labels,
            xlabel="Reaction coordinate (nm)",
            x_is_rc=True,
            rc_left=rc_left if rc_min is None else rc_min,
            rc_right=rc_right if rc_max is None else rc_max,
        )
    else:
        tick_step = max(1, matrix_width // 12)
        xticks = np.arange(0, matrix_width, tick_step)
        ax.set_xticks(xticks)
        ax.set_xticklabels([str(value) for value in xticks], fontsize=14, fontweight="bold")
        ax.set_yticks(np.arange(len(top_labels)))
        ax.set_yticklabels([pretty_ifp_label(label) for label in top_labels], fontsize=15, fontweight="bold")
        ax.set_xlabel("RC bin", fontsize=18, fontweight="bold")
    savefig(fig, figdir / "05_awh_top_contacts_heatmap.png", dpi)


def plot_top_contacts_from_heatmap(
    figdir: Path,
    heatmap_rows: list[dict[str, str]],
    rc_axis: np.ndarray | None,
    dpi: int,
    *,
    top_n: int = 12,
    rc_min: float | None = None,
    rc_max: float | None = None,
) -> None:
    if not heatmap_rows:
        return
    bin_cols = [key for key in heatmap_rows[0].keys() if key.startswith("bin_")]
    top_rows = sorted(
        heatmap_rows,
        key=lambda row: max(float(row[key]) for key in bin_cols),
        reverse=True,
    )[:top_n]
    matrix = np.array([[float(row[key]) for key in bin_cols] for row in top_rows], dtype=float)
    labels = [row["label"] for row in top_rows]

    fig, ax = plt.subplots(figsize=(12.8, max(5.8, 0.62 * len(labels))))
    style_axis(ax, facecolor=HEATMAP_FACE)
    if rc_axis is not None and rc_axis.size >= len(bin_cols):
        rc_used = np.asarray(rc_axis[: len(bin_cols)], dtype=float)
        edges = rc_edges_from_centers(rc_used)
        image = ax.imshow(
            matrix,
            aspect="auto",
            cmap="YlOrRd",
            origin="lower",
            interpolation="nearest",
            extent=[edges[0], edges[-1], -0.5, len(labels) - 0.5],
            vmin=0.0,
            vmax=max(0.1, float(matrix.max())),
        )
        rc_left = float(edges[0])
        rc_right = float(edges[-1])
    else:
        image = ax.imshow(
            matrix,
            aspect="auto",
            cmap="YlOrRd",
            origin="lower",
            interpolation="nearest",
            vmin=0.0,
            vmax=max(0.1, float(matrix.max())),
        )
        rc_left = None
        rc_right = None
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("Occupancy (%)", fontsize=14, fontweight="bold")
    cbar.ax.tick_params(labelsize=12, width=1.0, length=4)
    for tick in cbar.ax.get_yticklabels():
        tick.set_fontweight("bold")
    if rc_axis is not None and rc_axis.size >= len(bin_cols):
        style_contact_heatmap_axes(
            ax,
            labels,
            xlabel="Reaction coordinate (nm)",
            x_is_rc=True,
            rc_left=rc_left if rc_min is None else rc_min,
            rc_right=rc_right if rc_max is None else rc_max,
        )
    else:
        tick_step = max(1, len(bin_cols) // 12)
        xticks = np.arange(0, len(bin_cols), tick_step)
        ax.set_xticks(xticks)
        ax.set_xticklabels([str(value) for value in xticks], fontsize=14, fontweight="bold")
        ax.set_yticks(np.arange(len(labels)))
        ax.set_yticklabels([pretty_ifp_label(label) for label in labels], fontsize=15, fontweight="bold")
        ax.set_xlabel("RC bin", fontsize=18, fontweight="bold")
    savefig(fig, figdir / "05_awh_top_contacts_heatmap.png", dpi)


def plot_interface_metrics(
    figdir: Path,
    time_ps_hb: np.ndarray,
    mean_hb: np.ndarray,
    sd_hb: np.ndarray,
    hb_labels: list[str],
    hb_means: np.ndarray,
    time_ps_sc: np.ndarray | None,
    mean_sc: np.ndarray | None,
    sd_sc: np.ndarray | None,
    dpi: int,
) -> None:
    time_ps_hb, mean_hb, sd_hb = break_large_time_gaps(time_ps_hb, mean_hb, sd_hb)
    if time_ps_sc is not None and mean_sc is not None and sd_sc is not None:
        time_ps_sc, mean_sc, sd_sc = break_large_time_gaps(time_ps_sc, mean_sc, sd_sc)

    fig, axes = plt.subplots(2, 1, figsize=(11, 7))
    style_axis(axes[0], grid=True)
    style_axis(axes[1], grid=True)
    axes[0].plot(time_ps_hb / 1000.0, mean_hb, color="firebrick", linewidth=1.5, label="Mean HB")
    axes[0].fill_between(time_ps_hb / 1000.0, mean_hb - sd_hb, mean_hb + sd_hb, color="firebrick", alpha=0.2)
    axes[0].set_xlabel("Time (ns)")
    axes[0].set_ylabel("HB count")
    axes[0].set_title("Interface hydrogen bonds across trajectory segments")
    if time_ps_sc is not None and mean_sc is not None and sd_sc is not None:
        ax_sc = axes[0].twinx()
        style_axis(ax_sc, facecolor="none")
        ax_sc.plot(time_ps_sc / 1000.0, mean_sc, color="navy", linewidth=1.2, alpha=0.85, label="Mean SC contacts")
        ax_sc.fill_between(time_ps_sc / 1000.0, mean_sc - sd_sc, mean_sc + sd_sc, color="navy", alpha=0.12)
        ax_sc.set_ylabel("Sidechain contacts")

    x = np.arange(hb_means.size)
    axes[1].bar(x, hb_means, color=plt.get_cmap("tab20")(x % 20))
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(hb_labels, rotation=55, ha="right", fontsize=8)
    axes[1].set_ylabel("Mean HB")
    axes[1].set_title("Per-segment mean hydrogen bonds")
    savefig(fig, figdir / "06_interface_metrics.png", dpi)


def plot_rms(figdir: Path, rmsd_series: list[tuple[np.ndarray, np.ndarray]], rmsf_labels: list[str] | None, rmsf_stack: np.ndarray | None, walker_labels: list[str], dpi: int) -> None:
    if not rmsd_series and (rmsf_stack is None or not rmsf_stack.size):
        return
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    style_axis(axes[0], grid=True)
    style_axis(axes[1], grid=True)
    cmap = plt.get_cmap("tab20")
    for idx, (time_ps, rmsd) in enumerate(rmsd_series):
        label = walker_labels[idx] if idx < len(walker_labels) else f"seg_{idx}"
        axes[0].plot(time_ps / 1000.0, rmsd, linewidth=1.0, alpha=0.85, color=cmap(idx % 20), label=label)
    axes[0].set_xlabel("Time (ns)")
    axes[0].set_ylabel("RMSD (A)")
    axes[0].set_title("Per-segment RMSD")
    axes[0].legend(loc="best", fontsize=7, ncol=2)

    if rmsf_labels is not None and rmsf_stack is not None and rmsf_stack.size:
        mean_rmsf = rmsf_stack.mean(axis=0)
        sd_rmsf = rmsf_stack.std(axis=0)
        x = np.arange(mean_rmsf.size)
        axes[1].plot(x, mean_rmsf, color="black", linewidth=1.4)
        axes[1].fill_between(x, mean_rmsf - sd_rmsf, mean_rmsf + sd_rmsf, color="gray", alpha=0.25)
        axes[1].set_ylabel("RMSF (A)")
        axes[1].set_xlabel("Trace residue index")
        axes[1].set_title("Mean RMSF across trajectory segments")
        tick_step = max(1, mean_rmsf.size // 12)
        axes[1].set_xticks(x[::tick_step])
        axes[1].set_xticklabels(rmsf_labels[::tick_step], rotation=45, ha="right", fontsize=8)
    savefig(fig, figdir / "07_rmsd_rmsf.png", dpi)


def plot_contact_map(figdir: Path, cmap_data: np.ndarray, dpi: int, row_labels: list[str], col_labels: list[str]) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 6.2))
    style_axis(ax, facecolor=HEATMAP_FACE)
    cmap = plt.get_cmap("YlOrRd").copy()
    cmap.set_bad("white")
    vmax = float(np.nanpercentile(cmap_data, 99)) if np.isfinite(cmap_data).any() else 1.0
    vmax = max(vmax, float(np.nanmax(cmap_data)), 1e-6)
    norm = colors.PowerNorm(gamma=0.75, vmin=0.0, vmax=vmax)
    image = ax.imshow(cmap_data, cmap=cmap, norm=norm, aspect="auto", origin="lower")
    fig.colorbar(image, ax=ax, label="Average contact value")
    ax.set_xlabel("RNA residues")
    ax.set_ylabel("Protein residues")
    ax.set_title("Average contact map")
    x_step = 1 if len(col_labels) <= 20 else max(1, len(col_labels) // 14)
    y_step = 1 if len(row_labels) <= 20 else max(1, len(row_labels) // 14)
    xticks = np.arange(0, len(col_labels), x_step)
    yticks = np.arange(0, len(row_labels), y_step)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels([col_labels[i] for i in xticks], rotation=60, ha="right", fontsize=7)
    ax.set_yticklabels([row_labels[i] for i in yticks], fontsize=7)
    savefig(fig, figdir / "08_average_contact_map.png", dpi)


def plot_nucleic_profile(figdir: Path, rows: list[dict[str, str]], dpi: int) -> None:
    order = sorted(
        rows,
        key=lambda row: max(float(row["Base_Occupancy_%"]), float(row["Backbone_Occupancy_%"])),
        reverse=True,
    )[:18]
    labels = [f"{row['ResName']}{row['Protein_ResNum']}" for row in order]
    base = np.array([float(row["Base_Occupancy_%"]) for row in order], dtype=float)
    back = np.array([float(row["Backbone_Occupancy_%"]) for row in order], dtype=float)
    y = np.arange(len(order))

    fig, ax = plt.subplots(figsize=(10.5, max(4.5, 0.36 * len(order))))
    style_axis(ax, grid=True)
    ax.barh(y, back, color="#8fbcd4", label="Backbone")
    ax.barh(y, base, left=back, color="#d08b39", label="Base")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Occupancy (%)")
    ax.set_title("Protein-RNA recognition profile")
    ax.legend(loc="best")
    ax.invert_yaxis()
    savefig(fig, figdir / "09_nucleic_profile.png", dpi)


def plot_rna_annotations(figdir: Path, rows: list[dict[str, str]], dpi: int) -> None:
    order = sorted(
        rows,
        key=lambda row: max(float(row["BaseContact(%)"]), float(row["SugarContact(%)"]), float(row["PhosphateContact(%)"])),
        reverse=True,
    )
    labels = [
        f"{row['ResName_Rec']}{row['ResID_Rec']} | {row['ResName_Lig']}{row['ResID_Lig']}"
        for row in order
    ]
    base = np.array([float(row["BaseContact(%)"]) for row in order], dtype=float)
    sugar = np.array([float(row["SugarContact(%)"]) for row in order], dtype=float)
    phosphate = np.array([float(row["PhosphateContact(%)"]) for row in order], dtype=float)
    y = np.arange(len(order))

    fig, ax = plt.subplots(figsize=(10.5, max(3.8, 0.55 * len(order))))
    style_axis(ax, grid=True)
    ax.barh(y, base, color="#d08b39", label="Base")
    ax.barh(y, sugar, left=base, color="#74a57f", label="Sugar")
    ax.barh(y, phosphate, left=base + sugar, color="#5b8bd0", label="Phosphate")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Contact occupancy (%)")
    ax.set_title("RNA region annotations for top residue pairs")
    ax.legend(loc="best")
    ax.invert_yaxis()
    savefig(fig, figdir / "10_rna_region_annotations.png", dpi)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot TITAN AWH-IFP outputs.")
    parser.add_argument("--outdir", required=True, help="TITAN output directory")
    parser.add_argument("--figdir", default=None, help="Directory for generated figures")
    parser.add_argument("--dpi", type=int, default=220, help="PNG resolution")
    parser.add_argument(
        "--exclude-label",
        action="append",
        default=[],
        help="Interaction label to exclude from AWH top-contact heatmaps; may be passed multiple times.",
    )
    parser.add_argument(
        "--top-contacts-count",
        type=int,
        default=12,
        help="Number of top contacts to keep in figure 05.",
    )
    parser.add_argument("--heatmap-rc-min", type=float, default=None, help="Optional left limit for the RC heatmap axis.")
    parser.add_argument("--heatmap-rc-max", type=float, default=None, help="Optional right limit for the RC heatmap axis.")
    args = parser.parse_args()

    outdir = Path(args.outdir).expanduser().resolve()
    figdir = Path(args.figdir).expanduser().resolve() if args.figdir else (outdir / "figures_awh")
    figdir.mkdir(parents=True, exist_ok=True)

    rc_rows = load_whitespace_table(outdir / "awh_rc_bin_summary.dat")
    stability_path = outdir / "awh_walker_stability_by_bin.dat"
    pair_path = outdir / "awh_walker_pair_overlap.dat"
    structural_path = outdir / "awh_bin_structural_summary.dat"
    top_contacts_path = outdir / "awh_bin_top_contacts.csv"
    ifp_by_bin_path = outdir / "awh_ifp_by_bin.csv"
    ifp_heatmap_path = outdir / "awh_ifp_heatmap.dat"

    stability_rows = load_whitespace_table(stability_path) if stability_path.exists() else []
    pair_rows = load_whitespace_table(pair_path) if pair_path.exists() else []
    structural_rows = load_whitespace_table(structural_path) if structural_path.exists() else []
    ifp_rows = load_csv_dicts(top_contacts_path) if top_contacts_path.exists() else []
    ifp_by_bin_rows = load_csv_dicts(ifp_by_bin_path) if ifp_by_bin_path.exists() else []
    ifp_heatmap_rows = load_whitespace_table(ifp_heatmap_path) if ifp_heatmap_path.exists() else []
    excluded_labels = {label.strip() for label in args.exclude_label if label.strip()}
    ifp_rows = filter_label_rows(ifp_rows, excluded_labels)
    ifp_heatmap_rows = filter_label_rows(ifp_heatmap_rows, excluded_labels)
    rc_axis = np.array([float(row["rc_center"]) for row in (structural_rows if structural_rows else rc_rows)], dtype=float)

    plot_rc_coverage(figdir, structural_rows if structural_rows else rc_rows, args.dpi)
    if stability_rows:
        plot_walker_stability(figdir, stability_rows, args.dpi)
    if pair_rows:
        plot_pair_overlap(figdir, pair_rows, args.dpi)
    if structural_rows:
        plot_structural_summary(figdir, structural_rows, args.dpi)
    elif ifp_by_bin_rows:
        plot_structural_summary_legacy(figdir, ifp_by_bin_rows, args.dpi)
    if ifp_rows:
        plot_top_contacts(
            figdir,
            ifp_rows,
            rc_axis,
            args.dpi,
            top_n=args.top_contacts_count,
            rc_min=args.heatmap_rc_min,
            rc_max=args.heatmap_rc_max,
        )
    elif ifp_heatmap_rows:
        plot_top_contacts_from_heatmap(
            figdir,
            ifp_heatmap_rows,
            rc_axis,
            args.dpi,
            top_n=args.top_contacts_count,
            rc_min=args.heatmap_rc_min,
            rc_max=args.heatmap_rc_max,
        )

    hb_time, hb_mean, hb_sd = load_stat_hb(outdir / "stat_hydrogenbonds.dat")
    hb_labels, hb_means = load_stat_average_hb(outdir / "stat_average_HB.dat")
    sc_time, sc_mean, sc_sd = load_sc_time(outdir / "out_sc_time.dat")
    walker_list_path = outdir / "awh_walkers.txt"
    if not walker_list_path.exists():
        walker_list_path = outdir / "awh_traj_walkers.txt"
    walker_labels = [line.strip() for line in read_lines(walker_list_path) if line.strip()] if walker_list_path.exists() else []
    plot_interface_metrics(figdir, hb_time, hb_mean, hb_sd, hb_labels, hb_means, sc_time, sc_mean, sc_sd, args.dpi)

    rmsd_series: list[tuple[np.ndarray, np.ndarray]] = []
    rmsf_values: list[np.ndarray] = []
    rmsf_labels: list[str] | None = None
    rmsd_paths = sorted(outdir.glob("rmsd_*.dat"), key=lambda p: int(p.stem.split("_")[-1]))
    rmsf_paths = sorted(outdir.glob("rmsf_*.dat"), key=lambda p: int(p.stem.split("_")[-1]))
    for path in rmsd_paths:
        rmsd_series.append(load_rmsd(path))
    for path in rmsf_paths:
        labels, values = load_rmsf(path)
        if rmsf_labels is None:
            rmsf_labels = labels
        rmsf_values.append(values)
    rmsf_stack = np.vstack(rmsf_values) if rmsf_values else None
    plot_rms(figdir, rmsd_series, rmsf_labels, rmsf_stack, walker_labels, args.dpi)

    contact_map = load_contact_map(outdir / "average_contact_map.dat")
    row_labels, col_labels = infer_contact_map_labels(outdir, contact_map.shape)
    plot_contact_map(figdir, contact_map, args.dpi, row_labels, col_labels)

    nucleic_profile = load_csv_dicts(outdir / "nucleic_acid_interactions.csv")
    if nucleic_profile:
        plot_nucleic_profile(figdir, nucleic_profile, args.dpi)

    rna_annotations = load_csv_dicts(outdir / "ifp_rna_annotations.csv")
    if rna_annotations:
        plot_rna_annotations(figdir, rna_annotations, args.dpi)

    generated = sorted(figdir.glob("*.png"))
    print(f"Generated {len(generated)} figure(s) in: {figdir}")
    for path in generated:
        print(f" - {path.name}")


if __name__ == "__main__":
    main()
