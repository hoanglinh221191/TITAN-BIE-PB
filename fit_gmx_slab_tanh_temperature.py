#!/usr/bin/env python3
"""Fit symmetric tanh slab profiles and plot coexistence densities versus temperature."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


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
    if len(rows) < 5:
        raise ValueError(f"Profile has fewer than five data rows: {path}")
    data = np.asarray(rows, dtype=float)
    return data[:, 0], data[:, 1]


def symmetric_slab_tanh(
    z_nm: np.ndarray,
    rho_dense: float,
    rho_dilute: float,
    halfwidth_nm: float,
    interface_width_nm: float,
    center_nm: float,
) -> np.ndarray:
    radial_z = np.abs(z_nm - center_nm)
    return rho_dilute + 0.5 * (rho_dense - rho_dilute) * (
        1.0 - np.tanh((radial_z - halfwidth_nm) / interface_width_nm)
    )


def fit_profile(z_nm: np.ndarray, density: np.ndarray) -> tuple[np.ndarray, float, np.ndarray]:
    center_nm = 0.5 * (float(z_nm[0]) + float(z_nm[-1]))
    radial_z = np.abs(z_nm - center_nm)
    edge_mask = radial_z >= np.quantile(radial_z, 0.8)
    rho_dilute0 = max(0.0, float(np.mean(density[edge_mask])))
    rho_dense0 = max(float(np.percentile(density, 90)), rho_dilute0 + 1.0e-6)
    half_level = rho_dilute0 + 0.5 * (rho_dense0 - rho_dilute0)
    occupied = radial_z[density >= half_level]
    halfwidth0 = float(np.max(occupied)) if occupied.size else 0.1 * float(np.ptp(z_nm))
    maximum = max(float(np.max(density)), 1.0)

    def model(z: np.ndarray, dense: float, dilute: float, halfwidth: float, interface: float):
        return symmetric_slab_tanh(z, dense, dilute, halfwidth, interface, center_nm)

    params, _ = curve_fit(
        model,
        z_nm,
        density,
        p0=(rho_dense0, rho_dilute0, max(halfwidth0, 0.5), 1.0),
        bounds=(
            (0.0, 0.0, 0.05, 0.05),
            (3.0 * maximum, maximum, 0.5 * float(np.ptp(z_nm)), 20.0),
        ),
        maxfev=100000,
    )
    fitted = model(z_nm, *params)
    residual = float(np.sum((density - fitted) ** 2))
    total = float(np.sum((density - np.mean(density)) ** 2))
    r_squared = 1.0 - residual / total if total > 0.0 else float("nan")
    return params, r_squared, fitted


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


def configure_style() -> None:
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--temperatures", type=int, nargs="+", default=DEFAULT_TEMPERATURES)
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, float | int]] = []
    profiles: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for temperature in args.temperatures:
        directory = args.root / str(temperature) / "gmx_density_1000ns"
        full_path = directory / "density_mass_protein.xvg"
        block_paths = sorted((directory / "blocks").glob("block*_density_mass_protein.xvg"))
        if not full_path.is_file():
            print(f"Skipping {temperature} K: {full_path} is missing")
            continue
        if not block_paths:
            raise FileNotFoundError(f"No block profiles found in {directory / 'blocks'}")

        z_nm, density = read_xvg(full_path)
        full_params, full_r2, fitted = fit_profile(z_nm, density)
        profiles[temperature] = (z_nm, density, fitted)
        block_params: list[np.ndarray] = []
        block_r2: list[float] = []
        for block_path in block_paths:
            block_z, block_density = read_xvg(block_path)
            params, r_squared, _ = fit_profile(block_z, block_density)
            block_params.append(params)
            block_r2.append(r_squared)

        dense_mean, dense_sd = mean_sd([float(p[0]) for p in block_params])
        dilute_mean, dilute_sd = mean_sd([float(p[1]) for p in block_params])
        halfwidth_mean, halfwidth_sd = mean_sd([float(p[2]) for p in block_params])
        interface_mean, interface_sd = mean_sd([float(p[3]) for p in block_params])
        records.append(
            {
                "Temperature_K": temperature,
                "Blocks": len(block_params),
                "Dense_full_kg_m3": float(full_params[0]),
                "Dense_block_mean_kg_m3": dense_mean,
                "Dense_block_sd_kg_m3": dense_sd,
                "Dilute_full_kg_m3": float(full_params[1]),
                "Dilute_block_mean_kg_m3": dilute_mean,
                "Dilute_block_sd_kg_m3": dilute_sd,
                "Halfwidth_full_nm": float(full_params[2]),
                "Halfwidth_block_mean_nm": halfwidth_mean,
                "Halfwidth_block_sd_nm": halfwidth_sd,
                "Interface_full_nm": float(full_params[3]),
                "Interface_block_mean_nm": interface_mean,
                "Interface_block_sd_nm": interface_sd,
                "Full_fit_R2": full_r2,
                "Block_fit_R2_mean": float(np.mean(block_r2)),
                "Block_fit_R2_min": float(np.min(block_r2)),
            }
        )

    if not records:
        raise RuntimeError("No complete GROMACS protein-density datasets were found")

    csv_path = args.output_dir / "gmx_protein_tanh_fit_vs_temperature.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    configure_style()
    temperature = np.asarray([row["Temperature_K"] for row in records], dtype=float)
    dense = np.asarray([row["Dense_block_mean_kg_m3"] for row in records], dtype=float)
    dense_sd = np.asarray([row["Dense_block_sd_kg_m3"] for row in records], dtype=float)
    dilute = np.asarray([row["Dilute_block_mean_kg_m3"] for row in records], dtype=float)
    dilute_sd = np.asarray([row["Dilute_block_sd_kg_m3"] for row in records], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.0))
    axes[0].errorbar(temperature, dense, yerr=dense_sd, marker="o", markersize=4.2,
                     linewidth=1.6, capsize=2.5, color="#b84a32")
    axes[1].errorbar(temperature, dilute, yerr=dilute_sd, marker="o", markersize=4.2,
                     linewidth=1.6, capsize=2.5, color="#177e89")
    axes[0].set_title(r"Dense plateau, $\rho_{dense}$")
    axes[1].set_title(r"Dilute plateau, $\rho_{dilute}$")
    for label, axis in zip("AB", axes):
        axis.set_xlabel("Temperature (K)")
        axis.set_ylabel(r"Protein mass density (kg m$^{-3}$)")
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.text(-0.16, 1.04, label, transform=axis.transAxes, fontsize=12, fontweight="bold")
    fig.suptitle("Symmetric tanh fit of phosphorylated N-protein/RNA slab", fontsize=10,
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93), w_pad=2.4)
    save_figure(fig, args.output_dir / "fig_gmx_tanh_dense_dilute_vs_temperature", args.dpi)

    colors = plt.cm.inferno(np.linspace(0.12, 0.88, len(records)))
    fig, axes = plt.subplots(2, 3, figsize=(8.2, 5.0), sharex=False, sharey=True)
    for axis, color, row in zip(axes.flat, colors, records):
        temp = int(row["Temperature_K"])
        z_nm, density, fitted = profiles[temp]
        center = 0.5 * (float(z_nm[0]) + float(z_nm[-1]))
        axis.plot(z_nm - center, density, color=color, linewidth=0.8, alpha=0.45, label="GROMACS")
        axis.plot(z_nm - center, fitted, color=color, linewidth=1.7, label="Tanh fit")
        axis.set_title(f"{temp} K, $R^2$={float(row['Full_fit_R2']):.3f}")
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    for axis in axes[-1, :]:
        axis.set_xlabel("Distance from slab center (nm)")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"Protein density (kg m$^{-3}$)")
    axes[0, 0].legend(frameon=False, fontsize=7)
    fig.suptitle("Protein slab profiles and symmetric tanh fits", fontsize=10, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95), h_pad=1.8, w_pad=1.5)
    save_figure(fig, args.output_dir / "fig_gmx_tanh_profile_diagnostics", args.dpi)

    print(csv_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
