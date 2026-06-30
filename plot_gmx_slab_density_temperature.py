#!/usr/bin/env python3
"""Plot dense- and dilute-phase protein slab density versus temperature."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_TEMPERATURES = (290, 300, 310, 320, 330, 340, 350, 360, 370)


def read_xvg(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows: list[tuple[float, float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line[0] in "#@":
                continue
            fields = line.split()
            if len(fields) >= 2:
                rows.append((float(fields[0]), float(fields[1])))
    if len(rows) < 3:
        raise ValueError(f"Profile has fewer than three data rows: {path}")
    data = np.asarray(rows, dtype=float)
    return data[:, 0], data[:, 1]


def plateau_densities(
    z_nm: np.ndarray,
    density: np.ndarray,
    dense_halfwidth_nm: float,
    dilute_edge_width_nm: float,
) -> tuple[float, float]:
    center = 0.5 * (float(z_nm[0]) + float(z_nm[-1]))
    dense_mask = np.abs(z_nm - center) <= dense_halfwidth_nm
    dilute_mask = (z_nm <= z_nm[0] + dilute_edge_width_nm) | (
        z_nm >= z_nm[-1] - dilute_edge_width_nm
    )
    if not np.any(dense_mask) or not np.any(dilute_mask):
        raise ValueError("Dense or dilute plateau selection contains no bins")
    return float(np.mean(density[dense_mask])), float(np.mean(density[dilute_mask]))


def mean_sd(values: list[float]) -> tuple[float, float]:
    array = np.asarray(values, dtype=float)
    return float(np.mean(array)), float(np.std(array, ddof=1)) if array.size > 1 else 0.0


def save_figure(fig: plt.Figure, stem: Path, dpi: int) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(
        stem.with_suffix(".tiff"),
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--temperatures", type=int, nargs="+", default=DEFAULT_TEMPERATURES)
    parser.add_argument("--dense-halfwidth-nm", type=float, default=10.0)
    parser.add_argument("--dilute-edge-width-nm", type=float, default=15.0)
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, float | int]] = []

    for temperature in args.temperatures:
        directory = args.root / str(temperature) / "gmx_density_1000ns"
        full_path = directory / "density_mass_protein.xvg"
        block_paths = sorted((directory / "blocks").glob("block*_density_mass_protein.xvg"))
        if not full_path.is_file():
            print(f"Skipping {temperature} K: {full_path} is missing")
            continue
        if not block_paths:
            raise FileNotFoundError(f"No block protein-density profiles found in {directory / 'blocks'}")

        z_nm, density = read_xvg(full_path)
        full_dense, full_dilute = plateau_densities(
            z_nm, density, args.dense_halfwidth_nm, args.dilute_edge_width_nm
        )
        block_dense: list[float] = []
        block_dilute: list[float] = []
        for block_path in block_paths:
            block_z, block_density = read_xvg(block_path)
            dense, dilute = plateau_densities(
                block_z,
                block_density,
                args.dense_halfwidth_nm,
                args.dilute_edge_width_nm,
            )
            block_dense.append(dense)
            block_dilute.append(dilute)

        dense_mean, dense_sd = mean_sd(block_dense)
        dilute_mean, dilute_sd = mean_sd(block_dilute)
        records.append(
            {
                "Temperature_K": temperature,
                "Blocks": len(block_paths),
                "Dense_full_kg_m3": full_dense,
                "Dense_block_mean_kg_m3": dense_mean,
                "Dense_block_sd_kg_m3": dense_sd,
                "Dilute_full_kg_m3": full_dilute,
                "Dilute_block_mean_kg_m3": dilute_mean,
                "Dilute_block_sd_kg_m3": dilute_sd,
            }
        )

    if not records:
        raise RuntimeError("No complete GROMACS protein-density datasets were found")

    csv_path = args.output_dir / "gmx_protein_dense_dilute_vs_temperature.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    temperature = np.asarray([row["Temperature_K"] for row in records], dtype=float)
    dense = np.asarray([row["Dense_block_mean_kg_m3"] for row in records], dtype=float)
    dense_sd = np.asarray([row["Dense_block_sd_kg_m3"] for row in records], dtype=float)
    dilute = np.asarray([row["Dilute_block_mean_kg_m3"] for row in records], dtype=float)
    dilute_sd = np.asarray([row["Dilute_block_sd_kg_m3"] for row in records], dtype=float)

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 8,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.8,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.0))
    dense_color = "#b84a32"
    dilute_color = "#177e89"
    axes[0].errorbar(
        temperature,
        dense,
        yerr=dense_sd,
        marker="o",
        markersize=4.2,
        linewidth=1.6,
        capsize=2.5,
        color=dense_color,
    )
    axes[1].errorbar(
        temperature,
        dilute,
        yerr=dilute_sd,
        marker="o",
        markersize=4.2,
        linewidth=1.6,
        capsize=2.5,
        color=dilute_color,
    )
    axes[0].set_title("Dense phase")
    axes[1].set_title("Dilute phase")
    if np.allclose(dilute, 0.0) and np.allclose(dilute_sd, 0.0):
        axes[1].set_ylim(-0.002, 0.05)
        axes[1].set_yticks(np.linspace(0.0, 0.05, 6))
        axes[1].text(
            0.5,
            0.88,
            "No protein sampled\nin edge plateau",
            transform=axes[1].transAxes,
            ha="center",
            va="top",
            fontsize=8,
            color="#4b4b4b",
        )
    for label, axis in zip("AB", axes):
        axis.set_xlabel("Temperature (K)")
        axis.set_ylabel(r"Protein mass density (kg m$^{-3}$)")
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.text(-0.16, 1.04, label, transform=axis.transAxes, fontsize=12, fontweight="bold")
    fig.suptitle("Phosphorylated N-protein/RNA slab density", fontsize=10, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93), w_pad=2.4)
    save_figure(fig, args.output_dir / "fig_gmx_protein_dense_dilute_vs_temperature", args.dpi)

    print(csv_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
