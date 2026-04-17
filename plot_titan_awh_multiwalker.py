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
    from matplotlib import colors
    from matplotlib.patches import FancyArrowPatch
except ModuleNotFoundError as exc:
    raise SystemExit(
        "Missing plotting dependency. Please install matplotlib before running this script:\n"
        "  pip install matplotlib"
    ) from exc


def read_lines(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8", errors="ignore").splitlines()


def numeric_rows(path: Path, comment_prefix: str = "#") -> list[list[str]]:
    rows: list[list[str]] = []
    for raw in read_lines(path):
        line = raw.strip()
        if not line or line.startswith(comment_prefix):
            continue
        rows.append(line.split())
    return rows


def load_whitespace_array(path: Path) -> np.ndarray:
    rows = numeric_rows(path)
    if not rows:
        return np.empty((0, 0))
    return np.array([[float(x) for x in row] for row in rows], dtype=float)


def load_csv_dicts(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        return list(csv.DictReader(handle))


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
    "DT": "T",
    "T": "T",
    "THY": "T",
}


def biomol_kind_from_resnames(resnames: list[str]) -> str:
    names = [name.upper() for name in resnames if name]
    if names and all(name in AA_ONE for name in names):
        return "protein"
    if names and all(name in NA_ONE for name in names):
        return "nucleic"
    return "other"


def format_residue_label(chain: str, resnr: int, resname: str, compact: bool = True, biomol_kind: str | None = None) -> str:
    resname_u = resname.upper()
    kind = biomol_kind or biomol_kind_from_resnames([resname_u])
    if compact:
        if kind == "protein":
            code = AA_ONE.get(resname_u, resname_u[:1])
            return f"P:{code}{resnr}"
        if kind == "nucleic":
            code = NA_ONE.get(resname_u, resname_u[:1])
            return f"R:{code}{resnr}"
    return f"{chain}:{resname_u}{resnr}"


def parse_dssp_residue_records(path: Path) -> list[dict[str, str | int]]:
    records: list[dict[str, str | int]] = []
    seen: set[tuple[str, int, str]] = set()
    in_residue_block = False
    for raw in read_lines(path):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("# residue_content state ="):
            in_residue_block = True
            continue
        if not in_residue_block or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 6 or not parts[0].isdigit():
            continue
        key = (parts[1], int(parts[2]), parts[3])
        if key in seen:
            continue
        seen.add(key)
        records.append({"chain": parts[1], "resnr": int(parts[2]), "resname": parts[3]})
    return records


def infer_contact_map_labels(outdir: Path, shape: tuple[int, int]) -> tuple[list[str] | None, list[str] | None]:
    rows_n, cols_n = shape
    mol1_records: list[dict[str, str | int]] = []
    mol2_records: list[dict[str, str | int]] = []

    dssp1 = outdir / "dssp_mol1_content.dat"
    dssp2 = outdir / "dssp_mol2_content.dat"
    if dssp1.exists():
        mol1_records = parse_dssp_residue_records(dssp1)
    if dssp2.exists():
        mol2_records = parse_dssp_residue_records(dssp2)

    communities = outdir / "network_communities.csv"
    if communities.exists():
        grouped: dict[str, list[dict[str, str | int]]] = {}
        for row in load_csv_dicts(communities):
            chain = row.get("Chain", "")
            if not chain:
                continue
            grouped.setdefault(chain, []).append(
                {"chain": chain, "resnr": int(row["ResNr"]), "resname": row["ResName"]}
            )
        for chain in grouped:
            grouped[chain].sort(key=lambda rec: int(rec["resnr"]))

        if not mol1_records:
            for records in grouped.values():
                if len(records) == rows_n:
                    mol1_records = records
                    break

        if not mol2_records:
            used_chain = mol1_records[0]["chain"] if mol1_records else None
            for chain, records in grouped.items():
                if chain == used_chain:
                    continue
                if len(records) == cols_n:
                    mol2_records = records
                    break

    y_labels = None
    x_labels = None
    if len(mol1_records) == rows_n:
        kind = biomol_kind_from_resnames([str(rec["resname"]) for rec in mol1_records])
        y_labels = [
            format_residue_label(str(rec["chain"]), int(rec["resnr"]), str(rec["resname"]), compact=True, biomol_kind=kind)
            for rec in mol1_records
        ]
    if len(mol2_records) == cols_n:
        kind = biomol_kind_from_resnames([str(rec["resname"]) for rec in mol2_records])
        x_labels = [
            format_residue_label(str(rec["chain"]), int(rec["resnr"]), str(rec["resname"]), compact=True, biomol_kind=kind)
            for rec in mol2_records
        ]

    return y_labels, x_labels


def load_tica_projections(path: Path) -> np.ndarray:
    data = load_whitespace_array(path)
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError(f"Unexpected TICA projection format: {path}")
    return data


def load_pca_traj(path: Path) -> np.ndarray:
    rows = load_csv_dicts(path)
    if not rows:
        return np.empty((0, 3))
    return np.array(
        [[float(row["Frame"]), float(row["PC1"]), float(row["PC2"])] for row in rows],
        dtype=float,
    )


def load_umap(path: Path) -> np.ndarray:
    rows = load_csv_dicts(path)
    if not rows:
        return np.empty((0, 3))
    return np.array(
        [[float(row["frame"]), float(row["UMAP1"]), float(row["UMAP2"])] for row in rows],
        dtype=float,
    )


def load_fes_grid(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = load_whitespace_array(path)
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError(f"Unexpected FES format: {path}")
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    x_unique = np.unique(x)
    y_unique = np.unique(y)
    if x_unique.size * y_unique.size != data.shape[0]:
        return x, y, z
    x_index = {value: idx for idx, value in enumerate(x_unique)}
    y_index = {value: idx for idx, value in enumerate(y_unique)}
    grid = np.full((y_unique.size, x_unique.size), np.nan, dtype=float)
    for xi, yi, zi in zip(x, y, z):
        grid[y_index[yi], x_index[xi]] = zi
    xx, yy = np.meshgrid(x_unique, y_unique)
    return xx, yy, grid


def parse_header_value(path: Path, prefix: str) -> str | None:
    for line in read_lines(path):
        line = line.strip()
        if line.startswith(prefix):
            return line[len(prefix) :].strip()
    return None


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
        if not parts:
            continue
        matrix.append([float(value) for value in parts[1:]])
    return np.array(matrix, dtype=float)


def load_msm_macrostates(path: Path) -> np.ndarray:
    data = load_whitespace_array(path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Unexpected macrostate format: {path}")
    return data[:, 1].astype(int)


def load_msm_stationary(path: Path) -> np.ndarray:
    return load_whitespace_array(path)[:, 1]


def load_msm_tpm(path: Path) -> tuple[np.ndarray, float | None]:
    lines = read_lines(path)
    lag_ns = None
    if lines:
        dt_match = re.search(r"dt=([0-9.eE+-]+)\s*ns", lines[0])
        lag_match = re.search(r"lag=([0-9]+)", lines[0])
        if dt_match and lag_match:
            lag_ns = float(dt_match.group(1)) * float(lag_match.group(1))
    return load_whitespace_array(path), lag_ns


def load_ifp(path: Path) -> list[dict[str, str]]:
    rows = load_csv_dicts(path)
    filtered: list[dict[str, str]] = []
    for row in rows:
        vals = []
        for key, value in row.items():
            if key.endswith("(%)"):
                try:
                    vals.append(float(value))
                except ValueError:
                    vals.append(0.0)
        if vals and max(vals) > 0.0:
            filtered.append(row)
    return filtered


def parse_dssp_content(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    states: list[str] = []
    state_names: list[str] = []
    residue_labels: list[str] = []
    by_state: dict[str, list[float]] = {}
    current_state: str | None = None

    for raw in read_lines(path):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("# residue_content state ="):
            after_eq = line.split("=", 1)[1].strip()
            current_state = after_eq.split()[0]
            by_state[current_state] = []
            continue
        if line.startswith("# state_summary columns"):
            current_state = None
            continue
        if line.startswith("#"):
            continue
        parts = line.split()
        if current_state is None and len(parts) == 4:
            states.append(parts[0])
            state_names.append(parts[1])
            continue
        if current_state is not None and len(parts) >= 6:
            if len(residue_labels) < int(parts[0]):
                residue_labels.append(f"{parts[1]}:{parts[2]}:{parts[3]}")
            by_state[current_state].append(float(parts[4]))

    matrix = np.array([by_state[state] for state in states], dtype=float)
    return state_names, residue_labels, matrix


def walker_labels_from_file(outdir: Path) -> list[str]:
    walkers = outdir / "awh_walkers.txt"
    if walkers.exists():
        labels = [line.strip() for line in read_lines(walkers) if line.strip()]
        if labels:
            return labels
    rmsd_files = sorted(outdir.glob("rmsd_*.dat"))
    return [path.stem for path in rmsd_files]


def walker_frame_counts(outdir: Path) -> list[int]:
    counts: list[int] = []
    for path in sorted(outdir.glob("rmsd_*.dat"), key=lambda p: int(p.stem.split("_")[-1])):
        header_count = parse_header_value(path, "# frames")
        counts.append(int(float(header_count)) if header_count is not None else len(numeric_rows(path)))
    if counts:
        return counts
    stat_path = outdir / "stat_average_HB.dat"
    if stat_path.exists():
        for raw in read_lines(stat_path):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            if len(parts) >= 4:
                counts.append(int(parts[3]))
    return counts


def assign_walkers(total_frames: int, counts: list[int]) -> np.ndarray:
    if not counts:
        return np.zeros(total_frames, dtype=int)
    walker = np.full(total_frames, -1, dtype=int)
    start = 0
    for idx, count in enumerate(counts):
        end = min(start + count, total_frames)
        walker[start:end] = idx
        start = end
    if start < total_frames:
        walker[start:] = len(counts) - 1
    return walker


def sample_indices(size: int, target: int) -> np.ndarray:
    if size <= target:
        return np.arange(size)
    return np.linspace(0, size - 1, target, dtype=int)


def break_large_time_gaps(time_ps: np.ndarray, *series: np.ndarray, gap_factor: float = 50.0, min_gap_ps: float = 10000.0) -> tuple[np.ndarray, ...]:
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

    out_time: list[float] = []
    out_series: list[list[float]] = [[] for _ in series]
    gap_set = set(gap_idx.tolist())
    for idx, value in enumerate(time_ps):
        out_time.append(float(value))
        for series_idx, arr in enumerate(series):
            out_series[series_idx].append(float(arr[idx]))
        if idx in gap_set:
            out_time.append(np.nan)
            for series_idx in range(len(series)):
                out_series[series_idx].append(np.nan)

    return (np.array(out_time, dtype=float), *[np.array(values, dtype=float) for values in out_series])


def aggregate_macro_model(
    micro_to_macro: np.ndarray,
    stationary: np.ndarray,
    tpm: np.ndarray,
    micro_centroids: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_macro = int(micro_to_macro.max()) + 1
    macro_pop = np.zeros(n_macro, dtype=float)
    flux = np.zeros((n_macro, n_macro), dtype=float)
    macro_xy = np.zeros((n_macro, 2), dtype=float)

    for macro in range(n_macro):
        mask = micro_to_macro == macro
        weights = stationary[mask]
        macro_pop[macro] = weights.sum()
        if np.any(mask):
            if weights.sum() > 0:
                macro_xy[macro] = np.average(micro_centroids[mask, :2], axis=0, weights=weights)
            else:
                macro_xy[macro] = micro_centroids[mask, :2].mean(axis=0)

    for i in range(tpm.shape[0]):
        a = micro_to_macro[i]
        for j in range(tpm.shape[1]):
            b = micro_to_macro[j]
            flux[a, b] += stationary[i] * tpm[i, j]

    macro_t = np.zeros_like(flux)
    for a in range(n_macro):
        if macro_pop[a] > 0:
            macro_t[a] = flux[a] / macro_pop[a]
    return macro_pop, macro_xy, macro_t


def load_timescales(path: Path) -> np.ndarray:
    return load_whitespace_array(path)


def load_its_scan(path: Path) -> tuple[np.ndarray, list[str]]:
    header = read_lines(path)[0].lstrip("#").split()
    return load_whitespace_array(path), header


def load_macrostate_traj(path: Path) -> np.ndarray:
    return load_whitespace_array(path)


def savefig(fig: plt.Figure, path: Path, dpi: int) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_tica_fes(figdir: Path, xx: np.ndarray, yy: np.ndarray, zz: np.ndarray, macro_xy: np.ndarray, macro_pop: np.ndarray, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    if xx.ndim == 2:
        mesh = ax.pcolormesh(xx, yy, zz, shading="auto", cmap="viridis_r")
    else:
        mesh = ax.tricontourf(xx, yy, zz, levels=50, cmap="viridis_r")
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("Free energy / score")
    ax.scatter(
        macro_xy[:, 0],
        macro_xy[:, 1],
        s=2500 * (macro_pop / max(macro_pop.max(), 1e-12)) + 120,
        c=np.arange(macro_xy.shape[0]),
        cmap="tab10",
        edgecolor="black",
        linewidth=1.0,
        alpha=0.9,
    )
    for idx, (xv, yv) in enumerate(macro_xy):
        ax.text(xv, yv, f"M{idx}", ha="center", va="center", fontsize=10, weight="bold")
    ax.set_xlabel("TICA1")
    ax.set_ylabel("TICA2")
    ax.set_title("TICA FES")
    savefig(fig, figdir / "01_tica_fes_macro.png", dpi)


def plot_tica_scatter(figdir: Path, tica: np.ndarray, walker_ids: np.ndarray, walker_labels: list[str], macro_traj: np.ndarray | None, dpi: int, max_scatter: int) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    idx = sample_indices(tica.shape[0], max_scatter)
    walker_cmap = plt.get_cmap("tab10")

    for walker in np.unique(walker_ids[idx]):
        mask = walker_ids[idx] == walker
        label = walker_labels[walker] if walker < len(walker_labels) else f"walker_{walker}"
        axes[0].scatter(tica[idx][mask, 1], tica[idx][mask, 2], s=6, alpha=0.6, color=walker_cmap(walker % 10), label=label)
    axes[0].set_title("TICA projections colored by walker")
    axes[0].set_xlabel("TICA1")
    axes[0].set_ylabel("TICA2")
    axes[0].legend(loc="best", fontsize=8, markerscale=2)

    if macro_traj is not None and macro_traj.shape[0] == tica.shape[0]:
        macro_ids = macro_traj[:, 2].astype(int)
        for macro in np.unique(macro_ids[idx]):
            mask = macro_ids[idx] == macro
            axes[1].scatter(tica[idx][mask, 1], tica[idx][mask, 2], s=6, alpha=0.6, color=plt.get_cmap("tab20")(macro % 20), label=f"M{macro}")
        axes[1].legend(loc="best", fontsize=8, markerscale=2, ncol=2)
        axes[1].set_title("TICA projections colored by macrostate")
    else:
        scatter = axes[1].scatter(tica[idx, 1], tica[idx, 2], c=tica[idx, 0], s=6, cmap="plasma", alpha=0.7)
        fig.colorbar(scatter, ax=axes[1], label="Frame index")
        axes[1].set_title("TICA projections colored by frame index")

    axes[1].set_xlabel("TICA1")
    axes[1].set_ylabel("TICA2")
    savefig(fig, figdir / "02_tica_scatter.png", dpi)


def plot_umap(figdir: Path, umap: np.ndarray, walker_ids: np.ndarray, walker_labels: list[str], macro_traj: np.ndarray | None, dpi: int, max_scatter: int) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    idx = sample_indices(umap.shape[0], max_scatter)
    walker_cmap = plt.get_cmap("tab10")

    for walker in np.unique(walker_ids[idx]):
        mask = walker_ids[idx] == walker
        label = walker_labels[walker] if walker < len(walker_labels) else f"walker_{walker}"
        axes[0].scatter(umap[idx][mask, 1], umap[idx][mask, 2], s=6, alpha=0.6, color=walker_cmap(walker % 10), label=label)
    axes[0].set_title("UMAP embedding colored by walker")
    axes[0].set_xlabel("UMAP1")
    axes[0].set_ylabel("UMAP2")
    axes[0].legend(loc="best", fontsize=8, markerscale=2)

    if macro_traj is not None and macro_traj.shape[0] == umap.shape[0]:
        macro_ids = macro_traj[:, 2].astype(int)
        for macro in np.unique(macro_ids[idx]):
            mask = macro_ids[idx] == macro
            axes[1].scatter(umap[idx][mask, 1], umap[idx][mask, 2], s=6, alpha=0.6, color=plt.get_cmap("tab20")(macro % 20), label=f"M{macro}")
        axes[1].legend(loc="best", fontsize=8, markerscale=2, ncol=2)
        axes[1].set_title("UMAP embedding colored by macrostate")
    else:
        scatter = axes[1].scatter(umap[idx, 1], umap[idx, 2], c=umap[idx, 0], s=6, cmap="plasma", alpha=0.7)
        fig.colorbar(scatter, ax=axes[1], label="Frame index")
        axes[1].set_title("UMAP embedding colored by frame index")

    axes[1].set_xlabel("UMAP1")
    axes[1].set_ylabel("UMAP2")
    savefig(fig, figdir / "03_umap_scatter.png", dpi)


def plot_msm_network(figdir: Path, macro_pop: np.ndarray, macro_xy: np.ndarray, macro_t: np.ndarray, dpi: int, edge_threshold: float) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    sizes = 3500 * (macro_pop / max(macro_pop.max(), 1e-12)) + 200
    edges: list[tuple[int, int, float]] = []
    for i in range(macro_t.shape[0]):
        for j in range(macro_t.shape[1]):
            if i != j and macro_t[i, j] >= edge_threshold:
                edges.append((i, j, macro_t[i, j]))
    if not edges:
        ranked: list[tuple[int, int, float]] = []
        for i in range(macro_t.shape[0]):
            for j in range(macro_t.shape[1]):
                if i != j and macro_t[i, j] > 0.0:
                    ranked.append((i, j, macro_t[i, j]))
        ranked.sort(key=lambda item: item[2], reverse=True)
        edges = ranked[:8]
    edges.sort(key=lambda item: item[2], reverse=True)

    for i, j, weight in edges[:16]:
        start = macro_xy[i]
        end = macro_xy[j]
        arrow = FancyArrowPatch(
            posA=(start[0], start[1]),
            posB=(end[0], end[1]),
            arrowstyle="-|>",
            mutation_scale=12 + 18 * weight,
            linewidth=1.0 + 6.0 * weight,
            alpha=0.25 + 0.65 * weight,
            color="0.25",
        )
        ax.add_patch(arrow)
        ax.text(0.5 * (start[0] + end[0]), 0.5 * (start[1] + end[1]), f"{weight:.2f}", fontsize=8, color="0.2")

    ax.scatter(macro_xy[:, 0], macro_xy[:, 1], s=sizes, c=np.arange(macro_xy.shape[0]), cmap="tab10", edgecolor="black", linewidth=1.2, zorder=3)
    for idx, (xv, yv) in enumerate(macro_xy):
        ax.text(xv, yv, f"M{idx}\n{macro_pop[idx]:.2f}", ha="center", va="center", fontsize=9, weight="bold", zorder=4)
    ax.set_xlabel("Macrostate centroid on TICA1")
    ax.set_ylabel("Macrostate centroid on TICA2")
    ax.set_title("MSM macrostate network")
    savefig(fig, figdir / "04_msm_macro_network.png", dpi)


def plot_msm_its(figdir: Path, its_scan: np.ndarray, its_header: list[str], timescales: np.ndarray, selected_lag_ns: float | None, dpi: int) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    lag_ns = its_scan[:, 1]
    for col in range(2, its_scan.shape[1]):
        axes[0].plot(lag_ns, its_scan[:, col], marker="o", linewidth=1.6, label=its_header[col])
    if selected_lag_ns is not None:
        axes[0].axvline(selected_lag_ns, color="black", linestyle="--", linewidth=1.0, label="selected lag")
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Lag (ns)")
    axes[0].set_ylabel("Implied timescale (ns)")
    axes[0].set_title("MSM implied timescales")
    axes[0].legend(loc="best", fontsize=8)

    axes[1].bar(timescales[:, 0].astype(int), timescales[:, 1], color="steelblue")
    axes[1].set_xlabel("Timescale index")
    axes[1].set_ylabel("Timescale (ns)")
    axes[1].set_title("Timescales at selected lag")
    savefig(fig, figdir / "05_msm_its.png", dpi)


def plot_interface_metrics(figdir: Path, time_ps_hb: np.ndarray, mean_hb: np.ndarray, sd_hb: np.ndarray, hb_labels: list[str], hb_means: np.ndarray, time_ps_sc: np.ndarray | None, mean_sc: np.ndarray | None, sd_sc: np.ndarray | None, dpi: int) -> None:
    time_ps_hb, mean_hb, sd_hb = break_large_time_gaps(time_ps_hb, mean_hb, sd_hb)
    if time_ps_sc is not None and mean_sc is not None and sd_sc is not None:
        time_ps_sc, mean_sc, sd_sc = break_large_time_gaps(time_ps_sc, mean_sc, sd_sc)

    fig, axes = plt.subplots(2, 1, figsize=(11, 7))
    axes[0].plot(time_ps_hb / 1000.0, mean_hb, color="firebrick", linewidth=1.5, label="Mean HB")
    axes[0].fill_between(time_ps_hb / 1000.0, mean_hb - sd_hb, mean_hb + sd_hb, color="firebrick", alpha=0.2, label="Mean ± SD")
    axes[0].set_xlabel("Time (ns)")
    axes[0].set_ylabel("HB count")
    axes[0].set_title("Interface hydrogen bonds across walkers")
    axes[0].legend(loc="best")
    if time_ps_sc is not None and mean_sc is not None and sd_sc is not None:
        ax_sc = axes[0].twinx()
        ax_sc.plot(time_ps_sc / 1000.0, mean_sc, color="navy", linewidth=1.3, alpha=0.8, label="Mean SC contacts")
        ax_sc.fill_between(time_ps_sc / 1000.0, mean_sc - sd_sc, mean_sc + sd_sc, color="navy", alpha=0.14)
        ax_sc.set_ylabel("Sidechain contacts")

    x = np.arange(hb_means.size)
    axes[1].bar(x, hb_means, color=plt.get_cmap("tab10")(x % 10))
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(hb_labels, rotation=20, ha="right")
    axes[1].set_ylabel("Mean HB")
    axes[1].set_title("Per-walker mean hydrogen bonds")
    savefig(fig, figdir / "06_interface_metrics.png", dpi)


def plot_rms(figdir: Path, rmsd_series: list[tuple[np.ndarray, np.ndarray]], rmsf_labels: list[str] | None, rmsf_stack: np.ndarray | None, walker_labels: list[str], dpi: int) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    cmap = plt.get_cmap("tab10")
    for idx, (time_ps, rmsd) in enumerate(rmsd_series):
        label = walker_labels[idx] if idx < len(walker_labels) else f"walker_{idx}"
        axes[0].plot(time_ps / 1000.0, rmsd, linewidth=1.0, alpha=0.85, color=cmap(idx % 10), label=label)
    axes[0].set_xlabel("Time (ns)")
    axes[0].set_ylabel("RMSD (A)")
    axes[0].set_title("Per-walker RMSD")
    axes[0].legend(loc="best", fontsize=8)

    if rmsf_labels is not None and rmsf_stack is not None and rmsf_stack.size:
        mean_rmsf = rmsf_stack.mean(axis=0)
        sd_rmsf = rmsf_stack.std(axis=0)
        x = np.arange(mean_rmsf.size)
        axes[1].plot(x, mean_rmsf, color="black", linewidth=1.5)
        axes[1].fill_between(x, mean_rmsf - sd_rmsf, mean_rmsf + sd_rmsf, color="gray", alpha=0.25)
        axes[1].set_ylabel("RMSF (A)")
        axes[1].set_xlabel("Trace residue index")
        axes[1].set_title("Mean RMSF across walkers")
        tick_step = max(1, mean_rmsf.size // 12)
        axes[1].set_xticks(x[::tick_step])
        axes[1].set_xticklabels(rmsf_labels[::tick_step], rotation=45, ha="right", fontsize=8)
    savefig(fig, figdir / "07_rmsd_rmsf.png", dpi)


def plot_contact_map(figdir: Path, cmap_data: np.ndarray, dpi: int, row_labels: list[str] | None = None, col_labels: list[str] | None = None) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    cmap = plt.get_cmap("YlOrRd").copy()
    cmap.set_bad("white")
    vmax = float(np.nanpercentile(cmap_data, 99)) if np.isfinite(cmap_data).any() else 1.0
    vmax = max(vmax, float(np.nanmax(cmap_data)), 1e-6)
    norm = colors.PowerNorm(gamma=0.75, vmin=0.0, vmax=vmax)
    image = ax.imshow(cmap_data, cmap=cmap, norm=norm, aspect="auto", origin="lower")
    fig.colorbar(image, ax=ax, label="Average contact value")
    ax.set_xlabel("Group 4 residues")
    ax.set_ylabel("Group 3 residues")
    ax.set_title("Average contact map")
    if col_labels is None:
        col_labels = [str(i + 1) for i in range(cmap_data.shape[1])]
    if row_labels is None:
        row_labels = [str(i + 1) for i in range(cmap_data.shape[0])]

    x_step = 1 if len(col_labels) <= 20 else max(1, len(col_labels) // 14)
    y_step = 1 if len(row_labels) <= 20 else max(1, len(row_labels) // 14)
    xticks = np.arange(0, len(col_labels), x_step)
    yticks = np.arange(0, len(row_labels), y_step)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels([col_labels[i] for i in xticks], rotation=60, ha="right", fontsize=7)
    ax.set_yticklabels([row_labels[i] for i in yticks], fontsize=7)
    savefig(fig, figdir / "08_average_contact_map.png", dpi)


def plot_network_communities(figdir: Path, rows: list[dict[str, str]], dpi: int) -> None:
    ordered = sorted(rows, key=lambda row: (row.get("Chain", ""), int(row.get("ResNr", "0"))))
    x = np.arange(len(ordered))
    y = np.array([int(row["CommunityId"]) for row in ordered], dtype=int)
    chains = [row.get("Chain", "?") for row in ordered]
    unique_chains = list(dict.fromkeys(chains))
    chain_to_marker = {chain: marker for chain, marker in zip(unique_chains, ["o", "s", "^", "D", "v", "P", "X"])}

    fig, ax = plt.subplots(figsize=(11, 4.8))
    for chain in unique_chains:
        idx = np.array([item == chain for item in chains])
        ax.scatter(x[idx], y[idx], c=y[idx], cmap="tab20", s=45, marker=chain_to_marker[chain], edgecolor="black", linewidth=0.3, label=f"Chain {chain}")
    ax.set_xlabel("Residue order")
    ax.set_ylabel("Community ID")
    ax.set_title("Network community assignment along sequence")
    ax.legend(loc="best", fontsize=8)
    tick_step = max(1, len(ordered) // 12)
    labels = [f"{row['Chain']}:{row['ResNr']}" for row in ordered]
    ax.set_xticks(x[::tick_step])
    ax.set_xticklabels(labels[::tick_step], rotation=45, ha="right", fontsize=8)
    savefig(fig, figdir / "09_network_communities.png", dpi)


def plot_dssp(figdir: Path, state_names: list[str], residue_labels: list[str], matrix: np.ndarray, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(12, 4.5))
    image = ax.imshow(matrix, aspect="auto", cmap="viridis", vmin=0.0, vmax=max(0.05, float(np.nanmax(matrix))))
    fig.colorbar(image, ax=ax, label="Occupancy")
    ax.set_yticks(np.arange(len(state_names)))
    ax.set_yticklabels(state_names)
    tick_step = max(1, len(residue_labels) // 14)
    ticks = np.arange(0, len(residue_labels), tick_step)
    ax.set_xticks(ticks)
    ax.set_xticklabels(residue_labels[::tick_step], rotation=45, ha="right", fontsize=8)
    ax.set_xlabel("Residue")
    ax.set_title("DSSP occupancy heatmap")
    savefig(fig, figdir / "10_dssp_mol1.png", dpi)


def plot_pca(figdir: Path, pca_fes_path: Path, pca_traj: np.ndarray, walker_ids: np.ndarray, walker_labels: list[str], dpi: int, max_scatter: int) -> None:
    xx, yy, zz = load_fes_grid(pca_fes_path)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    if xx.ndim == 2:
        mesh = axes[0].pcolormesh(xx, yy, zz, shading="auto", cmap="viridis_r")
    else:
        mesh = axes[0].tricontourf(xx, yy, zz, levels=50, cmap="viridis_r")
    fig.colorbar(mesh, ax=axes[0], label="Free energy / score")
    axes[0].set_xlabel("PC1")
    axes[0].set_ylabel("PC2")
    axes[0].set_title("PCA FES")
    idx = sample_indices(pca_traj.shape[0], max_scatter)
    for walker in np.unique(walker_ids[idx]):
        mask = walker_ids[idx] == walker
        label = walker_labels[walker] if walker < len(walker_labels) else f"walker_{walker}"
        axes[1].scatter(pca_traj[idx][mask, 1], pca_traj[idx][mask, 2], s=6, alpha=0.6, color=plt.get_cmap("tab10")(walker % 10), label=label)
    axes[1].set_xlabel("PC1")
    axes[1].set_ylabel("PC2")
    axes[1].set_title("PCA projections colored by walker")
    axes[1].legend(loc="best", fontsize=8, markerscale=2)
    savefig(fig, figdir / "11_pca_fes_scatter.png", dpi)


def plot_ifp(figdir: Path, rows: list[dict[str, str]], dpi: int) -> None:
    metric_keys = [key for key in rows[0].keys() if key.endswith("(%)")]
    row_labels = [
        f"{row['ChainID_Rec']}:{row['ResName_Rec']}{row['ResID_Rec']} | "
        f"{row['ChainID_Lig']}:{row['ResName_Lig']}{row['ResID_Lig']}"
        for row in rows
    ]
    values = np.array([[float(row[key]) for key in metric_keys] for row in rows], dtype=float)
    order = np.argsort(values.max(axis=1))[::-1][: min(25, values.shape[0])]
    values = values[order]
    row_labels = [row_labels[i] for i in order]

    fig, ax = plt.subplots(figsize=(10, max(4.5, 0.35 * len(row_labels))))
    image = ax.imshow(values, aspect="auto", cmap="YlOrRd", vmin=0.0, vmax=max(1.0, float(values.max())))
    fig.colorbar(image, ax=ax, label="Occupancy (%)")
    ax.set_xticks(np.arange(len(metric_keys)))
    ax.set_xticklabels(metric_keys, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_title("Top IFP occupancies")
    savefig(fig, figdir / "12_ifp_occupancy.png", dpi)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot key TITAN multiwalker AWH outputs.")
    parser.add_argument("--outdir", required=True, help="TITAN output directory")
    parser.add_argument("--figdir", default=None, help="Directory for generated figures")
    parser.add_argument("--dpi", type=int, default=220, help="PNG resolution")
    parser.add_argument("--max-scatter", type=int, default=30000, help="Maximum points for scatter plots")
    parser.add_argument("--edge-threshold", type=float, default=0.05, help="Minimum macrostate transition probability to draw")
    args = parser.parse_args()

    outdir = Path(args.outdir).expanduser().resolve()
    figdir = Path(args.figdir).expanduser().resolve() if args.figdir else (outdir / "figures")
    figdir.mkdir(parents=True, exist_ok=True)

    walker_labels = walker_labels_from_file(outdir)
    frame_counts = walker_frame_counts(outdir)

    tica = load_tica_projections(outdir / "tica_projections.csv")
    walker_ids = assign_walkers(tica.shape[0], frame_counts)

    macro_traj = None
    macro_traj_path = outdir / "msm_macrostate_traj.dat"
    if macro_traj_path.exists():
        macro_traj = load_macrostate_traj(macro_traj_path)

    xx, yy, zz = load_fes_grid(outdir / "tica_fes_2d.dat")
    micro_to_macro = load_msm_macrostates(outdir / "msm_macrostates.dat")
    stationary = load_msm_stationary(outdir / "msm_stationary.dat")
    tpm, lag_ns = load_msm_tpm(outdir / "msm_tpm.dat")
    micro_centroids = load_whitespace_array(outdir / "msm_centroids.dat")
    macro_pop, macro_xy, macro_t = aggregate_macro_model(micro_to_macro, stationary, tpm, micro_centroids)

    plot_tica_fes(figdir, xx, yy, zz, macro_xy, macro_pop, args.dpi)
    plot_tica_scatter(figdir, tica, walker_ids, walker_labels, macro_traj, args.dpi, args.max_scatter)

    umap_path = outdir / "titan_umap_embedding.csv"
    if umap_path.exists():
        umap = load_umap(umap_path)
        plot_umap(figdir, umap, walker_ids[: umap.shape[0]], walker_labels, macro_traj, args.dpi, args.max_scatter)

    plot_msm_network(figdir, macro_pop, macro_xy, macro_t, args.dpi, args.edge_threshold)
    its_scan, its_header = load_its_scan(outdir / "msm_its_scan.dat")
    timescales = load_timescales(outdir / "msm_timescales.dat")
    plot_msm_its(figdir, its_scan, its_header, timescales, lag_ns, args.dpi)

    hb_time, hb_mean, hb_sd = load_stat_hb(outdir / "stat_hydrogenbonds.dat")
    hb_labels, hb_means = load_stat_average_hb(outdir / "stat_average_HB.dat")
    sc_time = sc_mean = sc_sd = None
    sc_path = outdir / "out_sc_time.dat"
    if sc_path.exists():
        sc_time, sc_mean, sc_sd = load_sc_time(sc_path)
    plot_interface_metrics(figdir, hb_time, hb_mean, hb_sd, hb_labels, hb_means, sc_time, sc_mean, sc_sd, args.dpi)

    rmsd_series: list[tuple[np.ndarray, np.ndarray]] = []
    rmsf_values: list[np.ndarray] = []
    rmsf_labels: list[str] | None = None
    for path in sorted(outdir.glob("rmsd_*.dat"), key=lambda p: int(p.stem.split("_")[-1])):
        rmsd_series.append(load_rmsd(path))
    for path in sorted(outdir.glob("rmsf_*.dat"), key=lambda p: int(p.stem.split("_")[-1])):
        labels, values = load_rmsf(path)
        if rmsf_labels is None:
            rmsf_labels = labels
        rmsf_values.append(values)
    rmsf_stack = np.vstack(rmsf_values) if rmsf_values else None
    plot_rms(figdir, rmsd_series, rmsf_labels, rmsf_stack, walker_labels, args.dpi)

    contact_map = load_contact_map(outdir / "average_contact_map.dat")
    contact_rows, contact_cols = infer_contact_map_labels(outdir, contact_map.shape)
    plot_contact_map(figdir, contact_map, args.dpi, row_labels=contact_rows, col_labels=contact_cols)

    communities_path = outdir / "network_communities.csv"
    if communities_path.exists():
        plot_network_communities(figdir, load_csv_dicts(communities_path), args.dpi)

    dssp_path = outdir / "dssp_mol1_content.dat"
    if dssp_path.exists():
        state_names, residue_labels, dssp_matrix = parse_dssp_content(dssp_path)
        if dssp_matrix.size:
            plot_dssp(figdir, state_names, residue_labels, dssp_matrix, args.dpi)

    pca_fes_path = outdir / "pca_fes_2d.dat"
    pca_traj_path = outdir / "pca_fel.csv"
    if pca_fes_path.exists() and pca_traj_path.exists():
        pca_traj = load_pca_traj(pca_traj_path)
        if pca_traj.size:
            plot_pca(figdir, pca_fes_path, pca_traj, walker_ids[: pca_traj.shape[0]], walker_labels, args.dpi, args.max_scatter)

    ifp_path = outdir / "ifp_occupancy.csv"
    if ifp_path.exists():
        ifp_rows = load_ifp(ifp_path)
        if ifp_rows:
            plot_ifp(figdir, ifp_rows, args.dpi)

    generated = sorted(figdir.glob("*.png"))
    print(f"Generated {len(generated)} figure(s) in: {figdir}")
    for path in generated:
        print(f" - {path.name}")


if __name__ == "__main__":
    main()
