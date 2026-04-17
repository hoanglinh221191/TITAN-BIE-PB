#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import string
import sys
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
    from matplotlib import colors as mcolors
except ModuleNotFoundError as exc:
    raise SystemExit(
        "Missing plotting dependency. Please install matplotlib before running this script:\n"
        "  pip install matplotlib"
    ) from exc


KB_KCAL_PER_MOL_K = 0.0019872041


@dataclass
class NearestFrameRecord:
    traj_index: int
    source: str
    frame: int
    pc1: float
    pc2: float
    distance_to_min: float
    frames_in_min_bin: int


@dataclass
class TrajectoryProjection:
    path: Path
    frames: np.ndarray
    coords: np.ndarray


def get_fes_cmap(name: str):
    if name == "nature_fes":
        return mcolors.LinearSegmentedColormap.from_list(
            "nature_fes",
            [
                "#102a43",
                "#284b63",
                "#4f6d7a",
                "#7d98a1",
                "#b8d8d8",
                "#edf2f4",
                "#fffdf8",
            ],
            N=256,
        )
    if name == "nature_fes_r":
        return get_fes_cmap("nature_fes").reversed()
    return plt.get_cmap(name)


def set_publication_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 180,
            "savefig.dpi": 400,
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
            "font.size": 10.5,
            "axes.titlesize": 10.5,
            "axes.labelsize": 11,
            "axes.labelweight": "bold",
            "axes.titleweight": "normal",
            "axes.linewidth": 0.9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "xtick.major.size": 3.8,
            "ytick.major.size": 3.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Average TITAN PCA FES from multiple pca_fel.csv files by averaging "
            "probability densities on a common PC1/PC2 grid."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help=(
            "Paths to pca_fel.csv files or directories containing pca_fel.csv. "
            "Example: run1/pca_fel.csv run2/pca_fel.csv run3/pca_fel.csv"
        ),
    )
    parser.add_argument(
        "--output",
        default="avg_pca_fes.png",
        help="Output image path. Default: avg_pca_fes.png",
    )
    parser.add_argument(
        "--pdf-output",
        default=None,
        help="Optional PDF version of the figure. Example: avg_pca_fes.pdf",
    )
    parser.add_argument(
        "--average-output",
        default=None,
        help=(
            "Optional separate image path for the mean FES only. "
            "If omitted and layout is not average-only, the script auto-creates "
            "a sibling file like *_avg_only.png."
        ),
    )
    parser.add_argument(
        "--average-pdf-output",
        default=None,
        help="Optional separate PDF path for the mean-only FES figure.",
    )
    parser.add_argument(
        "--minimum-report",
        default=None,
        help=(
            "Optional text report path for the minimum of the averaged FES. "
            "If omitted, the script auto-creates a sibling file ending with "
            "_minimum_report.txt."
        ),
    )
    parser.add_argument(
        "--layout",
        choices=("paper", "comparison", "average-only"),
        default="paper",
        help=(
            "paper: manuscript-style multi-panel figure; comparison: all runs + mean "
            "with run titles; average-only: only the mean FES panel. Default: paper"
        ),
    )
    parser.add_argument(
        "--grid-output",
        default="avg_pca_fes.dat",
        help="Output averaged FES grid path. Default: avg_pca_fes.dat",
    )
    parser.add_argument(
        "--grid-empty",
        choices=("vmax", "vmax+step", "nan"),
        default="vmax",
        help=(
            "Value written to the .dat grid for empty bins with zero probability. "
            "Default: vmax"
        ),
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=60,
        help="Number of bins per axis for the common PC grid. Default: 60",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=300.0,
        help="Temperature in Kelvin for F = -kBT ln P. Default: 300",
    )
    parser.add_argument(
        "--weighting",
        choices=("equal", "pooled"),
        default="equal",
        help=(
            "equal: each trajectory contributes the same total weight; "
            "pooled: all frames are pooled together. Default: equal"
        ),
    )
    parser.add_argument(
        "--pad",
        type=float,
        default=0.05,
        help="Fractional padding added around the global PC1/PC2 ranges. Default: 0.05",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=8.0,
        help="Maximum free-energy value shown in the color scale (kcal/mol). Default: 8.0",
    )
    parser.add_argument(
        "--contour-step",
        type=float,
        default=1.0,
        help="Contour interval in kcal/mol. Default: 1.0",
    )
    parser.add_argument(
        "--filled-levels",
        type=int,
        default=80,
        help="Number of filled contour levels. Default: 80",
    )
    parser.add_argument(
        "--cmap",
        default="nature_fes",
        help="Matplotlib colormap for FES. Default: nature_fes",
    )
    parser.add_argument(
        "--xlabel",
        default="PC1 (Å)",
        help="X-axis label. Default: PC1 (Å)",
    )
    parser.add_argument(
        "--ylabel",
        default="PC2 (Å)",
        help="Y-axis label. Default: PC2 (Å)",
    )
    parser.add_argument(
        "--no-min-marker",
        action="store_true",
        help="Disable white star marker at the free-energy minimum in each panel.",
    )
    parser.add_argument(
        "--no-equal-aspect",
        action="store_true",
        help="Use automatic aspect ratio instead of equal PC1/PC2 scaling.",
    )
    parser.add_argument(
        "--flip-x",
        nargs="*",
        type=int,
        default=[],
        help=(
            "1-based indices of input trajectories whose PC1 sign should be flipped "
            "before averaging. Example: --flip-x 2 3"
        ),
    )
    parser.add_argument(
        "--flip-y",
        nargs="*",
        type=int,
        default=[],
        help=(
            "1-based indices of input trajectories whose PC2 sign should be flipped "
            "before averaging. Example: --flip-y 2"
        ),
    )
    parser.add_argument(
        "--title",
        default="Average PCA Free Energy Surface",
        help="Figure title for the averaged panel. Default: Average PCA Free Energy Surface",
    )
    return parser.parse_args()


def resolve_input(path_str: str) -> Path:
    path = Path(path_str).expanduser().resolve()
    if path.is_dir():
        path = path / "pca_fel.csv"
    if not path.exists():
        raise FileNotFoundError(f"Input not found: {path}")
    return path


def derive_sibling_output(path: Path, suffix: str) -> Path:
    return path.with_name(f"{path.stem}{suffix}{path.suffix}")


def load_pca_fel(path: Path) -> TrajectoryProjection:
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = []
        for row in reader:
            if not row:
                continue
            try:
                frame = int(float(row["Frame"]))
                pc1 = float(row["PC1"])
                pc2 = float(row["PC2"])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(f"File {path} does not look like TITAN pca_fel.csv") from exc
            rows.append((frame, pc1, pc2))
    if not rows:
        raise ValueError(f"No PCA rows found in {path}")
    data = np.asarray(rows, dtype=float)
    return TrajectoryProjection(
        path=path,
        frames=data[:, 0].astype(int),
        coords=data[:, 1:3],
    )


def apply_requested_flips(data: np.ndarray, idx_1based: int, flip_x: set[int], flip_y: set[int]) -> np.ndarray:
    out = np.array(data, copy=True)
    if idx_1based in flip_x:
        out[:, 0] *= -1.0
    if idx_1based in flip_y:
        out[:, 1] *= -1.0
    return out


def build_common_edges(datasets: list[np.ndarray], bins: int, pad: float) -> tuple[np.ndarray, np.ndarray]:
    all_pc1 = np.concatenate([data[:, 0] for data in datasets])
    all_pc2 = np.concatenate([data[:, 1] for data in datasets])

    x_min = float(np.min(all_pc1))
    x_max = float(np.max(all_pc1))
    y_min = float(np.min(all_pc2))
    y_max = float(np.max(all_pc2))

    x_span = max(x_max - x_min, 1e-6)
    y_span = max(y_max - y_min, 1e-6)

    x_min -= pad * x_span
    x_max += pad * x_span
    y_min -= pad * y_span
    y_max += pad * y_span

    x_edges = np.linspace(x_min, x_max, bins + 1)
    y_edges = np.linspace(y_min, y_max, bins + 1)
    return x_edges, y_edges


def histogram_probability(data: np.ndarray, x_edges: np.ndarray, y_edges: np.ndarray) -> np.ndarray:
    hist, _, _ = np.histogram2d(data[:, 0], data[:, 1], bins=[x_edges, y_edges])
    total = float(hist.sum())
    if total <= 0.0:
        raise ValueError("Encountered an empty histogram while building FES")
    return hist / total


def probability_to_fes(probability: np.ndarray, temperature: float) -> np.ndarray:
    fes = np.full(probability.shape, np.nan, dtype=float)
    mask = probability > 0.0
    if not np.any(mask):
        raise ValueError("Probability grid is empty; cannot convert to FES")
    fes[mask] = -KB_KCAL_PER_MOL_K * temperature * np.log(probability[mask])
    fes[mask] -= np.nanmin(fes[mask])
    return fes


def grid_centers(x_edges: np.ndarray, y_edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return 0.5 * (x_edges[:-1] + x_edges[1:]), 0.5 * (y_edges[:-1] + y_edges[1:])


def locate_fes_minimum(
    fes: np.ndarray, x_edges: np.ndarray, y_edges: np.ndarray
) -> tuple[int, int, float, float, float]:
    if not np.any(np.isfinite(fes)):
        raise ValueError("FES grid contains no finite values")
    min_idx = np.unravel_index(np.nanargmin(fes), fes.shape)
    x_centers, y_centers = grid_centers(x_edges, y_edges)
    ix, iy = int(min_idx[0]), int(min_idx[1])
    return ix, iy, float(x_centers[ix]), float(y_centers[iy]), float(fes[ix, iy])


def count_frames_in_bin(coords: np.ndarray, x_edges: np.ndarray, y_edges: np.ndarray, ix: int, iy: int) -> int:
    x0, x1 = x_edges[ix], x_edges[ix + 1]
    y0, y1 = y_edges[iy], y_edges[iy + 1]
    in_x = (coords[:, 0] >= x0) & ((coords[:, 0] < x1) | ((ix == len(x_edges) - 2) & (coords[:, 0] <= x1)))
    in_y = (coords[:, 1] >= y0) & ((coords[:, 1] < y1) | ((iy == len(y_edges) - 2) & (coords[:, 1] <= y1)))
    return int(np.count_nonzero(in_x & in_y))


def nearest_frame_to_point(
    traj_index: int,
    dataset: TrajectoryProjection,
    target_pc1: float,
    target_pc2: float,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    ix: int,
    iy: int,
) -> NearestFrameRecord:
    deltas = dataset.coords - np.array([target_pc1, target_pc2], dtype=float)
    distances = np.sqrt(np.sum(deltas * deltas, axis=1))
    nearest_idx = int(np.argmin(distances))
    return NearestFrameRecord(
        traj_index=traj_index,
        source=(dataset.path.parent.name or dataset.path.stem),
        frame=int(dataset.frames[nearest_idx]),
        pc1=float(dataset.coords[nearest_idx, 0]),
        pc2=float(dataset.coords[nearest_idx, 1]),
        distance_to_min=float(distances[nearest_idx]),
        frames_in_min_bin=count_frames_in_bin(dataset.coords, x_edges, y_edges, ix, iy),
    )


def write_minimum_report(
    path: Path,
    datasets: list[TrajectoryProjection],
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    avg_fes: np.ndarray,
) -> None:
    ix, iy, min_pc1, min_pc2, min_f = locate_fes_minimum(avg_fes, x_edges, y_edges)
    nearest_records = [
        nearest_frame_to_point(idx, dataset, min_pc1, min_pc2, x_edges, y_edges, ix, iy)
        for idx, dataset in enumerate(datasets, start=1)
    ]
    global_best = min(nearest_records, key=lambda item: item.distance_to_min)

    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("Averaged PCA-FES minimum report\n")
        handle.write("========================================\n\n")
        handle.write(
            "Note: the star in the averaged FES marks the minimum of the averaged grid,\n"
            "not an exact raw frame. Records below identify the nearest raw frame(s)\n"
            "in PC space to that averaged-grid minimum.\n\n"
        )
        handle.write(f"Minimum grid bin index : ix={ix}, iy={iy}\n")
        handle.write(f"Minimum grid center    : PC1={min_pc1:.6f}, PC2={min_pc2:.6f}\n")
        handle.write(f"Minimum FES value      : {min_f:.6f} kcal/mol\n")
        handle.write(
            f"Bin bounds             : PC1=[{x_edges[ix]:.6f}, {x_edges[ix + 1]:.6f}], "
            f"PC2=[{y_edges[iy]:.6f}, {y_edges[iy + 1]:.6f}]\n\n"
        )

        handle.write("Nearest global frame\n")
        handle.write("----------------------------------------\n")
        handle.write(f"Trajectory index       : {global_best.traj_index}\n")
        handle.write(f"Trajectory source      : {global_best.source}\n")
        handle.write(f"Frame                  : {global_best.frame}\n")
        handle.write(f"Frame PC coordinates   : PC1={global_best.pc1:.6f}, PC2={global_best.pc2:.6f}\n")
        handle.write(f"Distance to grid min   : {global_best.distance_to_min:.6f}\n")
        handle.write(f"Frames in minimum bin  : {global_best.frames_in_min_bin}\n\n")

        handle.write("Nearest frame per trajectory\n")
        handle.write("----------------------------------------\n")
        for record in nearest_records:
            handle.write(f"Trajectory {record.traj_index}: {record.source}\n")
            handle.write(f"  Frame                : {record.frame}\n")
            handle.write(f"  Frame PC coordinates : PC1={record.pc1:.6f}, PC2={record.pc2:.6f}\n")
            handle.write(f"  Distance to grid min : {record.distance_to_min:.6f}\n")
            handle.write(f"  Frames in min bin    : {record.frames_in_min_bin}\n")
            handle.write("\n")

def write_fes_grid(
    path: Path,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    fes: np.ndarray,
    empty_mode: str,
    vmax: float,
    contour_step: float,
) -> None:
    x_centers, y_centers = grid_centers(x_edges, y_edges)
    if empty_mode == "vmax":
        empty_value = vmax
    elif empty_mode == "vmax+step":
        empty_value = vmax + max(contour_step, 0.5)
    else:
        empty_value = None

    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("# PC1_center PC2_center F_kcal_per_mol\n")
        for ix, x_val in enumerate(x_centers):
            for iy, y_val in enumerate(y_centers):
                val = fes[ix, iy]
                if np.isnan(val):
                    if empty_value is None:
                        handle.write(f"{x_val:.6f} {y_val:.6f} nan\n")
                    else:
                        handle.write(f"{x_val:.6f} {y_val:.6f} {empty_value:.6f}\n")
                else:
                    handle.write(f"{x_val:.6f} {y_val:.6f} {val:.6f}\n")


def make_panel_title(path: Path, n_frames: int, flipped_x: bool, flipped_y: bool, layout: str) -> str:
    suffix = []
    if flipped_x:
        suffix.append("flipX")
    if flipped_y:
        suffix.append("flipY")
    flip_note = ""
    if suffix:
        flip_note = f" [{', '.join(suffix)}]"
    base = (path.parent.name or path.stem).replace("_", " ")
    if layout == "paper":
        return f"{base}{flip_note}"
    return f"{base}\nframes={n_frames}{flip_note}"


def add_panel_label(ax: plt.Axes, label: str) -> None:
    if not label:
        return
    ax.text(
        0.02,
        0.98,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12.5,
        fontweight="bold",
        color="#111111",
    )


def plot_single_fes(
    ax: plt.Axes,
    fes: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    vmax: float,
    contour_step: float,
    filled_levels: int,
    cmap: str,
    xlabel: str,
    ylabel: str,
    equal_aspect: bool,
    title: str,
    panel_label: str,
    mark_minimum: bool,
):
    x_centers, y_centers = grid_centers(x_edges, y_edges)
    xx, yy = np.meshgrid(x_centers, y_centers, indexing="ij")
    masked_fes = np.ma.masked_invalid(fes)
    cmap_obj = get_fes_cmap(cmap)

    filled = ax.contourf(
        xx,
        yy,
        masked_fes,
        levels=np.linspace(0.0, vmax, max(filled_levels, 4)),
        cmap=cmap_obj,
        extend="max",
    )

    if contour_step > 0.0:
        contour_levels = np.arange(0.0, vmax + 0.5 * contour_step, contour_step)
        if len(contour_levels) >= 2:
            ax.contour(
                xx,
                yy,
                masked_fes,
                levels=contour_levels,
                colors="#1f1f1f",
                linewidths=0.38,
                alpha=0.34,
            )

    if mark_minimum and np.any(np.isfinite(fes)):
        min_idx = np.unravel_index(np.nanargmin(fes), fes.shape)
        ax.scatter(
            x_centers[min_idx[0]],
            y_centers[min_idx[1]],
            s=78,
            marker="*",
            facecolor="white",
            edgecolor="#111111",
            linewidth=0.7,
            zorder=5,
        )

    if title:
        ax.set_title(title, pad=4.0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if equal_aspect:
        ax.set_aspect("equal", adjustable="box")
    ax.tick_params(length=3.6, width=0.8)
    add_panel_label(ax, panel_label)
    return filled


def plot_fes_panels(
    input_paths: list[Path],
    datasets: list[np.ndarray],
    individual_fes: list[np.ndarray],
    avg_fes: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    output_path: Path,
    vmax: float,
    title: str,
    flip_x: set[int],
    flip_y: set[int],
    weighting: str,
    layout: str,
    contour_step: float,
    filled_levels: int,
    cmap: str,
    xlabel: str,
    ylabel: str,
    mark_minimum: bool,
    equal_aspect: bool,
) -> None:
    if layout == "average-only":
        fig, ax = plt.subplots(figsize=(4.95, 4.35), constrained_layout=False)
        mesh = plot_single_fes(
            ax=ax,
            fes=avg_fes,
            x_edges=x_edges,
            y_edges=y_edges,
            vmax=vmax,
            contour_step=contour_step,
            filled_levels=filled_levels,
            cmap=cmap,
            xlabel=xlabel,
            ylabel=ylabel,
            equal_aspect=equal_aspect,
            title="",
            panel_label="",
            mark_minimum=mark_minimum,
        )
        fig.subplots_adjust(left=0.12, right=0.83, bottom=0.12, top=0.995)
        cbar = fig.colorbar(mesh, ax=ax, shrink=0.78, pad=0.015, fraction=0.045, aspect=28)
        cbar.set_label("Free energy (kcal mol$^{-1}$)", fontsize=10.5)
        cbar.ax.tick_params(labelsize=9, length=3.0, width=0.7)
        cbar.outline.set_linewidth(0.7)
        fig.savefig(output_path, dpi=400, bbox_inches="tight")
        plt.close(fig)
        return

    n_panels = len(individual_fes) + 1
    if layout == "paper" and n_panels == 4:
        nrows, ncols = 2, 2
        fig_size = (8.6, 7.3)
    else:
        ncols = 2 if n_panels <= 4 else 3
        nrows = int(np.ceil(n_panels / ncols))
        fig_size = (4.8 * ncols, 4.2 * nrows)

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=fig_size, constrained_layout=True)
    axes_arr = np.atleast_1d(axes).ravel()

    panel_letters = iter(string.ascii_uppercase)
    mesh = None
    for i, (path, data, fes) in enumerate(zip(input_paths, datasets, individual_fes), start=1):
        ax = axes_arr[i - 1]
        mesh = plot_single_fes(
            ax=ax,
            fes=fes,
            x_edges=x_edges,
            y_edges=y_edges,
            vmax=vmax,
            contour_step=contour_step,
            filled_levels=filled_levels,
            cmap=cmap,
            xlabel=xlabel,
            ylabel=ylabel,
            equal_aspect=equal_aspect,
            title=make_panel_title(path, len(data), i in flip_x, i in flip_y, layout),
            panel_label=next(panel_letters),
            mark_minimum=mark_minimum,
        )

    avg_ax = axes_arr[len(individual_fes)]
    mesh = plot_single_fes(
        ax=avg_ax,
        fes=avg_fes,
        x_edges=x_edges,
        y_edges=y_edges,
        vmax=vmax,
        contour_step=contour_step,
        filled_levels=filled_levels,
        cmap=cmap,
        xlabel=xlabel,
        ylabel=ylabel,
        equal_aspect=equal_aspect,
        title="Mean FES",
        panel_label=next(panel_letters),
        mark_minimum=mark_minimum,
    )

    for ax in axes_arr[n_panels:]:
        ax.axis("off")

    if mesh is not None:
        cbar = fig.colorbar(mesh, ax=axes_arr[:n_panels], shrink=0.93, pad=0.015, fraction=0.03)
        cbar.set_label("Free energy (kcal mol$^{-1}$)", fontsize=10.5)
        cbar.ax.tick_params(labelsize=9, length=3.0, width=0.7)
        cbar.outline.set_linewidth(0.7)

    if layout != "paper":
        fig.suptitle(
            "TITAN PCA-FES comparison\nAverage is computed from probability densities on a shared grid",
            fontsize=12,
            fontweight="bold",
        )
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    set_publication_style()

    input_paths = [resolve_input(path_str) for path_str in args.inputs]
    flip_x = set(args.flip_x)
    flip_y = set(args.flip_y)

    print(
        "Warning: averaging separate pca_fel.csv files is only meaningful if the "
        "PC1/PC2 axes are comparable across runs.\n"
        "Best practice for three initial configurations is to run TITAN once on all "
        "three trajectories together with -i input.txt so PCA is built in one common basis."
    )

    raw_datasets = [load_pca_fel(path) for path in input_paths]
    datasets = [
        TrajectoryProjection(
            path=data.path,
            frames=data.frames,
            coords=apply_requested_flips(data.coords, idx, flip_x, flip_y),
        )
        for idx, data in enumerate(raw_datasets, start=1)
    ]
    dataset_coords = [data.coords for data in datasets]

    x_edges, y_edges = build_common_edges(dataset_coords, args.bins, args.pad)

    probabilities = [histogram_probability(data.coords, x_edges, y_edges) for data in datasets]
    individual_fes = [probability_to_fes(prob, args.temperature) for prob in probabilities]

    if args.weighting == "equal":
        avg_probability = np.mean(np.stack(probabilities, axis=0), axis=0)
    else:
        pooled_hist = np.sum(
            np.stack(probabilities, axis=0) * np.array([len(data.coords) for data in datasets])[:, None, None],
            axis=0,
        )
        avg_probability = pooled_hist / np.sum(pooled_hist)

    avg_fes = probability_to_fes(avg_probability, args.temperature)

    output_path = Path(args.output).expanduser().resolve()
    grid_output_path = Path(args.grid_output).expanduser().resolve()
    pdf_output_path = (
        None if args.pdf_output is None else Path(args.pdf_output).expanduser().resolve()
    )
    average_output_path = (
        None if args.average_output is None else Path(args.average_output).expanduser().resolve()
    )
    average_pdf_output_path = (
        None
        if args.average_pdf_output is None
        else Path(args.average_pdf_output).expanduser().resolve()
    )
    minimum_report_path = (
        None if args.minimum_report is None else Path(args.minimum_report).expanduser().resolve()
    )

    if args.layout != "average-only" and average_output_path is None:
        average_output_path = derive_sibling_output(output_path, "_avg_only")
    if pdf_output_path is not None and args.layout != "average-only" and average_pdf_output_path is None:
        average_pdf_output_path = derive_sibling_output(pdf_output_path, "_avg_only")
    if minimum_report_path is None:
        report_anchor = average_output_path if average_output_path is not None else output_path
        minimum_report_path = derive_sibling_output(report_anchor, "_minimum_report")
        minimum_report_path = minimum_report_path.with_suffix(".txt")

    write_fes_grid(
        grid_output_path,
        x_edges,
        y_edges,
        avg_fes,
        args.grid_empty,
        args.vmax,
        args.contour_step,
    )
    plot_fes_panels(
        input_paths=input_paths,
        datasets=dataset_coords,
        individual_fes=individual_fes,
        avg_fes=avg_fes,
        x_edges=x_edges,
        y_edges=y_edges,
        output_path=output_path,
        vmax=args.vmax,
        title=args.title,
        flip_x=flip_x,
        flip_y=flip_y,
        weighting=args.weighting,
        layout=args.layout,
        contour_step=args.contour_step,
        filled_levels=args.filled_levels,
        cmap=args.cmap,
        xlabel=args.xlabel,
        ylabel=args.ylabel,
        mark_minimum=not args.no_min_marker,
        equal_aspect=not args.no_equal_aspect,
    )
    if pdf_output_path is not None:
        plot_fes_panels(
            input_paths=input_paths,
            datasets=dataset_coords,
            individual_fes=individual_fes,
            avg_fes=avg_fes,
            x_edges=x_edges,
            y_edges=y_edges,
            output_path=pdf_output_path,
            vmax=args.vmax,
            title=args.title,
            flip_x=flip_x,
            flip_y=flip_y,
            weighting=args.weighting,
            layout=args.layout,
            contour_step=args.contour_step,
            filled_levels=args.filled_levels,
            cmap=args.cmap,
            xlabel=args.xlabel,
            ylabel=args.ylabel,
            mark_minimum=not args.no_min_marker,
            equal_aspect=not args.no_equal_aspect,
        )
    if average_output_path is not None:
        plot_fes_panels(
            input_paths=input_paths,
            datasets=dataset_coords,
            individual_fes=individual_fes,
            avg_fes=avg_fes,
            x_edges=x_edges,
            y_edges=y_edges,
            output_path=average_output_path,
            vmax=args.vmax,
            title=args.title,
            flip_x=flip_x,
            flip_y=flip_y,
            weighting=args.weighting,
            layout="average-only",
            contour_step=args.contour_step,
            filled_levels=args.filled_levels,
            cmap=args.cmap,
            xlabel=args.xlabel,
            ylabel=args.ylabel,
            mark_minimum=not args.no_min_marker,
            equal_aspect=not args.no_equal_aspect,
        )
    if average_pdf_output_path is not None:
        plot_fes_panels(
            input_paths=input_paths,
            datasets=dataset_coords,
            individual_fes=individual_fes,
            avg_fes=avg_fes,
            x_edges=x_edges,
            y_edges=y_edges,
            output_path=average_pdf_output_path,
            vmax=args.vmax,
            title=args.title,
            flip_x=flip_x,
            flip_y=flip_y,
            weighting=args.weighting,
            layout="average-only",
            contour_step=args.contour_step,
            filled_levels=args.filled_levels,
            cmap=args.cmap,
            xlabel=args.xlabel,
            ylabel=args.ylabel,
            mark_minimum=not args.no_min_marker,
            equal_aspect=not args.no_equal_aspect,
        )
    write_minimum_report(minimum_report_path, datasets, x_edges, y_edges, avg_fes)

    print(f"Wrote averaged FES grid to: {grid_output_path}")
    print(f"Wrote figure to: {output_path}")
    if pdf_output_path is not None:
        print(f"Wrote PDF figure to: {pdf_output_path}")
    if average_output_path is not None:
        print(f"Wrote mean-only figure to: {average_output_path}")
    if average_pdf_output_path is not None:
        print(f"Wrote mean-only PDF figure to: {average_pdf_output_path}")
    print(f"Wrote averaged-minimum report to: {minimum_report_path}")


if __name__ == "__main__":
    main()
