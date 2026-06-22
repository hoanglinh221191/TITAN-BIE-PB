#!/usr/bin/env python3
"""Fit TITAN slab density profiles and draw a preliminary binodal figure."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


TEMPERATURES = (290, 300, 310, 320, 330, 340)
KJMOL_TO_KCALMOL = 1.0 / 4.184
MOLECULES_PER_NM3_TO_MM = 1660.539067
PROTEIN_BB_PER_CHAIN = 419.0


def tanh_slab(distance, rho_dense, rho_dilute, halfwidth, interface_width):
    width = np.maximum(interface_width, 1.0e-6)
    return rho_dilute + 0.5 * (rho_dense - rho_dilute) * (
        1.0 - np.tanh((distance - halfwidth) / width)
    )


def read_profile(path: Path):
    rows = list(csv.DictReader(path.open("r", encoding="utf-8-sig")))
    distance = np.asarray([float(row["BinCenter_nm"]) for row in rows])
    frames = np.asarray([float(row["Frames"]) for row in rows])
    protein = np.asarray([float(row["ProteinDensity"]) for row in rows])
    rna = np.asarray([float(row["RNADensity"]) for row in rows])
    mask = frames > 0.0
    return distance[mask], protein[mask], rna[mask]


def read_key_values(path: Path):
    values = {}
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        fields = raw.split()
        if len(fields) >= 2:
            values[fields[0]] = fields[1]
    return values


def fit_density(distance: np.ndarray, density: np.ndarray):
    maximum = float(np.max(density))
    central = density[distance <= min(3.0, float(np.max(distance)) * 0.2)]
    rho_dense0 = float(np.mean(central)) if central.size else maximum
    tail = density[distance >= 0.7 * float(np.max(distance))]
    rho_dilute0 = max(0.0, float(np.mean(tail)) if tail.size else 0.0)
    half_level = rho_dilute0 + 0.5 * (rho_dense0 - rho_dilute0)
    halfwidth0 = float(distance[np.argmin(np.abs(density - half_level))])
    p0 = [max(rho_dense0, maximum * 0.5), rho_dilute0, halfwidth0, 1.0]
    upper_density = max(1.0, 3.0 * maximum)
    try:
        params, covariance = curve_fit(
            tanh_slab,
            distance,
            density,
            p0=p0,
            bounds=(
                [0.0, 0.0, 0.0, 0.05],
                [upper_density, upper_density, float(np.max(distance)), 20.0],
            ),
            maxfev=50000,
        )
        errors = np.sqrt(np.maximum(0.0, np.diag(covariance)))
    except Exception:
        params = np.asarray([np.nan, np.nan, np.nan, np.nan])
        errors = np.asarray([np.nan, np.nan, np.nan, np.nan])
    fitted = tanh_slab(distance, *params) if np.all(np.isfinite(params)) else np.full_like(density, np.nan)
    residual = float(np.sum((density - fitted) ** 2))
    total = float(np.sum((density - np.mean(density)) ** 2))
    r2 = 1.0 - residual / total if total > 0.0 else np.nan
    return params, errors, r2, fitted


def save_figure(fig, stem: Path, dpi: int):
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
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    results = []
    profiles = {}
    for temperature in TEMPERATURES:
        directory = args.root / str(temperature) / "titan_binodal_1000ns"
        distance, protein, rna = read_profile(directory / "llps_slab_profiles.csv")
        convergence = read_key_values(directory / "llps_convergence_report.dat")
        p_par, p_err, p_r2, p_fit = fit_density(distance, protein)
        r_par, r_err, r_r2, r_fit = fit_density(distance, rna)
        profiles[temperature] = (distance, protein, p_fit, rna, r_fit)
        results.append(
            {
                "Temperature_K": temperature,
                "Frames": int(float(convergence.get("frames_used", 0))),
                "TimeSpan_ns": float(convergence.get("time_span_ns", "nan")),
                "StableDensePhase": convergence.get("stable_dense_phase", "unknown"),
                "ProteinDense_nodes_nm3": p_par[0],
                "ProteinDense_SEfit": p_err[0],
                "ProteinDilute_nodes_nm3": p_par[1],
                "ProteinDilute_SEfit": p_err[1],
                "ProteinDense_mM_monomer": p_par[0] / PROTEIN_BB_PER_CHAIN * MOLECULES_PER_NM3_TO_MM,
                "ProteinDilute_mM_monomer": p_par[1] / PROTEIN_BB_PER_CHAIN * MOLECULES_PER_NM3_TO_MM,
                "ProteinOrderParameter_nodes_nm3": p_par[0] - p_par[1],
                "ProteinHalfWidth_nm": p_par[2],
                "ProteinInterfaceWidth_nm": p_par[3],
                "ProteinFit_R2": p_r2,
                "RNADense_nodes_nm3": r_par[0],
                "RNADilute_nodes_nm3": r_par[1],
                "RNAHalfWidth_nm": r_par[2],
                "RNAInterfaceWidth_nm": r_par[3],
                "RNAFit_R2": r_r2,
            }
        )

    csv_path = args.output_dir / "slab_binodal_fit_summary.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)

    colors = plt.cm.inferno(np.linspace(0.12, 0.88, len(TEMPERATURES)))
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.35))
    for color, temperature in zip(colors, TEMPERATURES):
        distance, protein, fitted, _, _ = profiles[temperature]
        axes[0].plot(distance, protein, color=color, alpha=0.35, linewidth=0.8)
        axes[0].plot(distance, fitted, color=color, linewidth=1.7, label=f"{temperature} K")
    axes[0].set_xlabel("Distance from slab center (nm)")
    axes[0].set_ylabel(r"Protein BB density (nodes nm$^{-3}$)")
    axes[0].set_title("Slab density profiles")
    axes[0].legend(frameon=False, fontsize=7, ncol=2)

    temperature = np.asarray([row["Temperature_K"] for row in results], dtype=float)
    dense = np.asarray([row["ProteinDense_mM_monomer"] for row in results], dtype=float)
    dilute = np.asarray([row["ProteinDilute_mM_monomer"] for row in results], dtype=float)
    axes[1].plot(dilute, temperature, "o-", color="#287271", label="Dilute branch")
    axes[1].plot(dense, temperature, "o-", color="#b24c32", label="Dense branch")
    axes[1].fill_betweenx(temperature, dilute, dense, color="#d8c7a1", alpha=0.25)
    axes[1].set_xlabel("N-protein monomer-equivalent (mM)")
    axes[1].set_ylabel("Temperature (K)")
    axes[1].set_title("Preliminary coexistence curve")
    axes[1].legend(frameon=False, fontsize=7)

    halfwidth = np.asarray([row["ProteinHalfWidth_nm"] for row in results], dtype=float)
    interface = np.asarray([row["ProteinInterfaceWidth_nm"] for row in results], dtype=float)
    axes[2].plot(temperature, halfwidth, "o-", color="#315b7d", label="Slab half-width")
    axes[2].plot(temperature, interface, "s-", color="#d08b32", label="Interface width")
    axes[2].set_xlabel("Temperature (K)")
    axes[2].set_ylabel("Length (nm)")
    axes[2].set_title("Slab morphology")
    axes[2].legend(frameon=False, fontsize=7)

    for label, axis in zip("ABC", axes):
        axis.text(-0.17, 1.06, label, transform=axis.transAxes, fontsize=12, fontweight="bold")
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    fig.suptitle("Phosphorylated N-protein/RNA slab, 1000 ns to trajectory end", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save_figure(fig, args.output_dir / "fig_slab_preliminary_binodal", args.dpi)
    print(csv_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
