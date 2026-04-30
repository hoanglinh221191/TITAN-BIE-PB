#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gc
import json
import math
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Sequence

SCRIPT_DIR = Path(__file__).resolve().parent
LOCAL_PYDEPS = SCRIPT_DIR / ".pydeps"
if LOCAL_PYDEPS.exists():
    sys.path.insert(0, str(LOCAL_PYDEPS))

try:
    import numpy as np
except ModuleNotFoundError as exc:
    raise SystemExit(
        "Missing numerical dependency. Please install numpy before running this script:\n"
        "  pip install numpy"
    ) from exc

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import colors as mcolors
    from matplotlib import gridspec
    from matplotlib.lines import Line2D
except ModuleNotFoundError as exc:
    raise SystemExit(
        "Missing plotting dependency. Please install matplotlib before running this script:\n"
        "  pip install matplotlib"
    ) from exc


SYSTEM_ORDER = ("phos", "unphos")
RADIAL_SYSTEM_ORDER = ("unphos", "phos")
SYSTEM_DISPLAY = {
    "phos": "Phosphorylated",
    "unphos": "Unphosphorylated (WT)",
}
SYSTEM_COLORS = {
    "phos": "#b54d32",
    "unphos": "#1e6f73",
}
SYSTEM_LIGHT = {
    "phos": "#e9c6bc",
    "unphos": "#b8d9da",
}
SIDECHAIN_COLORS = {
    "protein": "#4c78a8",
    "rna": "#f58518",
}
DENSE_RG_COLORS = {
    ("phos", "protein"): "#b54d32",
    ("phos", "rna"): "#e59d2f",
    ("unphos", "protein"): "#1e6f73",
    ("unphos", "rna"): "#4178b8",
}
PANEL_LABELS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
NATURAL_SORT_RE = re.compile(r"(\d+)")
DISCOVERY_RE = re.compile(r"^10dimer(?P<system>unphos|phos)_conf(?P<conf>\d+)_800$", re.IGNORECASE)
NODE_LABEL_RE = re.compile(r"(?P<resname>[A-Za-z]+)(?P<resid>-?\d+)(?:_(?P<chain>[^_]+))?$")
LIQUIDITY_REPORT_KEYS = (
    "tau_ns",
    "mean_num_condensates",
    "cv_num_condensates",
    "mean_giant_cluster_prochains",
    "cv_giant_cluster_prochains",
    "mean_giant_cluster_rnachains",
    "cv_giant_cluster_rnachains",
    "mean_rg_pro_dense",
    "cv_rg_pro_dense",
    "mean_rg_rna_dense",
    "cv_rg_rna_dense",
    "mean_sticker_contact_lifetime_ns",
    "mean_protein_rna_contact_lifetime_ns",
    "protein_neighbor_persistence",
    "rna_neighbor_persistence",
    "global_neighbor_persistence",
    "protein_rewiring_rate_per_ns",
    "rna_rewiring_rate_per_ns",
    "global_rewiring_rate_per_ns",
    "new_contacts_per_ns",
    "lost_contacts_per_ns",
    "net_change_per_ns",
    "protein_diffusion_nm2_per_ns",
    "rna_diffusion_nm2_per_ns",
)
RNA_COMPACTION_BINNED_COLUMNS = [
    "BoundProteinChainsBin",
    "FrameCount",
    "RgMean_nm",
    "RgStd_nm",
    "EndToEndMean_nm",
    "EndToEndStd_nm",
    "AsphericityMean",
    "AsphericityStd",
    "LongRangeContactFractionMean",
    "LongRangeContactFractionStd",
]
RNA_CHAIN_COMPACTION_COLUMNS = [
    "RNAChainID",
    "FrameCount",
    "InsideDenseFraction",
    "MeanBoundProteinChains",
    "MeanProteinRNABeadContacts",
    "MeanRg_nm",
    "StdRg_nm",
    "MeanEndToEnd_nm",
    "StdEndToEnd_nm",
    "MeanRelativeShapeAnisotropy",
    "MeanLongRangeContactFraction",
]
RNA_INTRACHAIN_COLUMNS = [
    "RNAChainID",
    "FrameCount",
    "MeanIntrachainContacts",
    "StdIntrachainContacts",
    "MeanLongRangeIntrachainContacts",
    "StdLongRangeIntrachainContacts",
    "MeanLongRangeContactFraction",
    "StdLongRangeContactFraction",
]
RNA_COMPACTION_TRACE_COLUMNS = [
    "Time_ns",
    "Frame",
    "RNAChainID",
    "InsideDense",
    "ClusterID",
    "BoundProteinChains",
    "ProteinRNABeadContacts",
    "Rg_nm",
    "EndToEnd_nm",
    "Asphericity_nm2",
    "NormalizedAsphericity",
    "Acylindricity_nm2",
    "RelativeShapeAnisotropy",
    "AspectRatio",
    "IntrachainContacts",
    "LongRangeIntrachainContacts",
    "LongRangeContactFraction",
]
RNA_ENERGY_BINNED_COLUMNS = [
    "BoundProteinChainsBin",
    "FrameCount",
    "EelecMean_kJmol",
    "EelecStd_kJmol",
    "EvdwMean_kJmol",
    "EvdwStd_kJmol",
    "EtotalMean_kJmol",
    "EtotalStd_kJmol",
    "AttractionStrengthMean_kJmol",
    "AttractionStrengthStd_kJmol",
    "RgMean_nm",
    "RgStd_nm",
    "EndToEndMean_nm",
    "EndToEndStd_nm",
    "RelativeShapeAnisotropyMean",
    "RelativeShapeAnisotropyStd",
    "LongRangeContactFractionMean",
    "LongRangeContactFractionStd",
]
RNA_BINDING_ENERGY_TRACE_COLUMNS = [
    "Time_ns",
    "Frame",
    "RNAChainID",
    "InsideDense",
    "ClusterID",
    "BoundProteinChains",
    "ProteinRNABeadContacts",
    "Rg_nm",
    "EndToEnd_nm",
    "RelativeShapeAnisotropy",
    "LongRangeContactFraction",
    "Eelec_ProteinRNA_kJmol",
    "Evdw_ProteinRNA_kJmol",
    "Etotal_ProteinRNA_kJmol",
]
RNA_COMPACTION_STATE_HISTOGRAM_COLUMNS = [
    "Scope",
    "StateID",
    "FrameCount",
    "Probability",
    "RgBin",
    "RgMin_nm",
    "RgMax_nm",
    "LongRangeBin",
    "LongRangeMin",
    "LongRangeMax",
]
RNA_TOPOLOGY_TRACE_COLUMNS = [
    "Time_ns",
    "Frame",
    "RNAChainID",
    "InsideDense",
    "ClusterID",
    "BoundProteinChains",
    "TopologyState",
    "AnyBridge",
    "NumBridgeProteins",
    "AnyClamp",
    "NumClampPairs",
    "BridgeResidueSpanMax",
    "BridgeGeomSpanMax_nm",
    "LongRangeIntrachainContacts",
    "ClampSupportedLongRangeContacts",
    "BridgeTotal_kJmol",
    "NonBridgeTotal_kJmol",
    "BridgeAttractionStrength_kJmol",
    "NonBridgeAttractionStrength_kJmol",
    "PRNAAttractionStrength_kJmol",
    "LocalPPStickerContacts",
    "LocalPPElec_kJmol",
    "LocalPPVdw_kJmol",
    "LocalPPTotal_kJmol",
    "LocalPPAttractionStrength_kJmol",
    "LocalPPtoPRNAContactRatio",
    "LocalPPtoPRNAAttractionRatio",
    "LocalPPtoBridgeAttractionRatio",
]
RNA_TOPOLOGY_BINNED_COLUMNS = [
    "Scope",
    "BoundProteinChainsBin",
    "FrameCount",
    "AnyBridgeFraction",
    "AnyClampFraction",
    "LongRangeGivenBound",
    "LongRangeGivenBridge",
    "LongRangeGivenClamp",
    "MeanNumBridgeProteins",
    "MeanNumClampPairs",
    "MeanBridgeResidueSpanMax",
    "MeanBridgeGeomSpanMax_nm",
    "MeanClampSupportedLongRangeContacts",
    "MeanBridgeTotal_kJmol",
    "MeanNonBridgeTotal_kJmol",
    "MeanBridgeAttractionStrength_kJmol",
    "MeanNonBridgeAttractionStrength_kJmol",
    "MeanPRNAAttractionStrength_kJmol",
    "MeanLocalPPElec_kJmol",
    "MeanLocalPPVdw_kJmol",
    "MeanLocalPPTotal_kJmol",
    "MeanLocalPPAttractionStrength_kJmol",
    "MeanLocalPPStickerContacts",
    "MeanLocalPPtoPRNAContactRatio",
    "MeanLocalPPtoPRNAAttractionRatio",
    "MeanLocalPPtoBridgeAttractionRatio",
    "MeanLongRangeContactFractionWhenBound",
]


@dataclass
class Table:
    columns: dict[str, np.ndarray]
    order: list[str]
    x_name: str | None = None
    x_label: str | None = None

    def row_count(self) -> int:
        if not self.order:
            return 0
        return int(self.columns[self.order[0]].size)

    def has_columns(self, names: Sequence[str]) -> bool:
        return all(name in self.columns for name in names)

    def as_rows(self) -> list[list[object]]:
        n_rows = self.row_count()
        return [[self.columns[name][idx] for name in self.order] for idx in range(n_rows)]


@dataclass
class DXGrid:
    counts: tuple[int, int, int]
    origin: tuple[float, float, float]
    deltas: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]
    label: str
    data: np.ndarray | None = None
    projection_xy_data: np.ndarray | None = None
    radial_profile_data: Table | None = None
    shape_metrics: dict[str, float] = field(default_factory=dict)

    def compatible_with(self, other: "DXGrid") -> bool:
        return (
            self.counts == other.counts
            and np.allclose(self.origin, other.origin, atol=1e-6)
            and np.allclose(self.deltas, other.deltas, atol=1e-6)
        )

    def projection_xy(self) -> np.ndarray:
        if self.projection_xy_data is not None:
            return self.projection_xy_data
        if self.data is None:
            return np.zeros((self.counts[0], self.counts[1]), dtype=np.float32)
        self.projection_xy_data = np.sum(self.data, axis=2, dtype=np.float64).astype(np.float32, copy=False)
        return self.projection_xy_data

    @property
    def voxel_size(self) -> float:
        return float(np.linalg.norm(np.asarray(self.deltas[0], dtype=float)))

    @property
    def spacing_x(self) -> float:
        return float(np.linalg.norm(np.asarray(self.deltas[0], dtype=float)))

    @property
    def spacing_y(self) -> float:
        return float(np.linalg.norm(np.asarray(self.deltas[1], dtype=float)))

    @property
    def spacing_z(self) -> float:
        return float(np.linalg.norm(np.asarray(self.deltas[2], dtype=float)))

    @property
    def extent_xy_nm(self) -> tuple[float, float, float, float]:
        x0 = float(self.origin[0])
        y0 = float(self.origin[1])
        x1 = x0 + float(self.counts[0]) * self.spacing_x
        y1 = y0 + float(self.counts[1]) * self.spacing_y
        return (x0, x1, y0, y1)

    @property
    def has_full_data(self) -> bool:
        return self.data is not None


@dataclass
class NucleicInteraction:
    resnum: int
    resname: str
    base_frames: float
    backbone_frames: float
    base_occupancy: float
    backbone_occupancy: float
    recognition_type: str

    @property
    def total_occupancy(self) -> float:
        return self.base_occupancy + self.backbone_occupancy


@dataclass
class SidechainContact:
    chain: str
    resnum: int
    resname: str
    score: float

    @property
    def label(self) -> str:
        return f"{self.chain}:{self.resname}{self.resnum}"


@dataclass
class NetworkBundle:
    nodes: dict[str, dict[str, str]] = field(default_factory=dict)
    corr_edges: dict[tuple[str, str], float] = field(default_factory=dict)
    gc_edges: dict[tuple[str, str], float] = field(default_factory=dict)
    te_edges: dict[tuple[str, str], float] = field(default_factory=dict)
    mi_edges: dict[tuple[str, str], float] = field(default_factory=dict)

    def available(self) -> bool:
        return bool(self.nodes or self.corr_edges or self.gc_edges or self.te_edges or self.mi_edges)


@dataclass
class NetworkAnchorBundle:
    labels: dict[str, str] = field(default_factory=dict)
    gc_anchor_scores: dict[str, float] = field(default_factory=dict)
    te_source_scores: dict[str, float] = field(default_factory=dict)
    te_sink_scores: dict[str, float] = field(default_factory=dict)

    def available(self) -> bool:
        return bool(self.gc_anchor_scores or self.te_source_scores or self.te_sink_scores)


@dataclass
class DynamicsBundle:
    summary: dict[str, float] = field(default_factory=dict)
    table: Table | None = None

    def available(self) -> bool:
        return bool(self.summary) or (self.table is not None and self.table.row_count() > 0)


@dataclass
class ReplicateData:
    system: str
    replicate_label: str
    source_dir: Path
    frame_weight: int
    llps_timeseries: Table | None = None
    llps_autocorr: Table | None = None
    llps_degree: Table | None = None
    llps_summary: dict[str, object] = field(default_factory=dict)
    dbscan_timeseries: Table | None = None
    polymer_scaling: Table | None = None
    density_protein: DXGrid | None = None
    density_rna: DXGrid | None = None
    nucleic_interactions: list[NucleicInteraction] = field(default_factory=list)
    sidechain_protein: list[SidechainContact] = field(default_factory=list)
    sidechain_rna: list[SidechainContact] = field(default_factory=list)
    network_bundle: NetworkBundle | None = None
    network_anchor_bundle: NetworkAnchorBundle | None = None
    tda_points: np.ndarray | None = None
    density_shape_metrics: dict[str, dict[str, float]] = field(default_factory=dict)
    sticker_contact_lifetimes: DynamicsBundle | None = None
    protein_rna_contact_lifetimes: DynamicsBundle | None = None
    neighbor_rewiring: DynamicsBundle | None = None
    interface_turnover: DynamicsBundle | None = None
    dense_msd_protein: DynamicsBundle | None = None
    dense_msd_rna: DynamicsBundle | None = None
    liquidity_summary: dict[str, object] = field(default_factory=dict)
    rna_compaction_summary: dict[str, object] = field(default_factory=dict)
    rna_compaction_binned: Table | None = None
    rna_energy_compaction_summary: dict[str, object] = field(default_factory=dict)
    rna_energy_compaction_binned: Table | None = None
    rna_information_summary: dict[str, object] = field(default_factory=dict)
    rna_compaction_state_histogram: Table | None = None
    rna_topology_summary: dict[str, object] = field(default_factory=dict)
    rna_topology_binned: Table | None = None
    rna_topology_trace: Table | None = None
    rna_chain_compaction_summary: Table | None = None
    rna_intrachain_contacts_summary: Table | None = None
    rna_compaction_trace: Table | None = None
    rna_binding_energy_trace: Table | None = None
    warnings: list[str] = field(default_factory=list)


@dataclass
class AggregatedSystem:
    system: str
    replicates: list[ReplicateData]
    warnings: list[str] = field(default_factory=list)
    llps_timeseries: Table | None = None
    llps_autocorr: Table | None = None
    llps_degree: Table | None = None
    llps_summary: dict[str, object] = field(default_factory=dict)
    dbscan_timeseries: Table | None = None
    polymer_scaling: Table | None = None
    density_protein: DXGrid | None = None
    density_rna: DXGrid | None = None
    nucleic_interactions: list[NucleicInteraction] = field(default_factory=list)
    sidechain_protein: list[SidechainContact] = field(default_factory=list)
    sidechain_rna: list[SidechainContact] = field(default_factory=list)
    network_bundle: NetworkBundle | None = None
    network_anchor_bundle: NetworkAnchorBundle | None = None
    tda_points: np.ndarray | None = None
    density_shape_metrics: dict[str, dict[str, float]] = field(default_factory=dict)
    sticker_contact_lifetimes: DynamicsBundle | None = None
    protein_rna_contact_lifetimes: DynamicsBundle | None = None
    neighbor_rewiring: DynamicsBundle | None = None
    interface_turnover: DynamicsBundle | None = None
    dense_msd_protein: DynamicsBundle | None = None
    dense_msd_rna: DynamicsBundle | None = None
    liquidity_summary: dict[str, object] = field(default_factory=dict)
    rna_compaction_summary: dict[str, object] = field(default_factory=dict)
    rna_compaction_binned: Table | None = None
    rna_energy_compaction_summary: dict[str, object] = field(default_factory=dict)
    rna_energy_compaction_binned: Table | None = None
    rna_information_summary: dict[str, object] = field(default_factory=dict)
    rna_compaction_state_histogram: Table | None = None
    rna_topology_summary: dict[str, object] = field(default_factory=dict)
    rna_topology_binned: Table | None = None
    rna_topology_trace: Table | None = None
    rna_chain_compaction_summary: Table | None = None
    rna_intrachain_contacts_summary: Table | None = None
    rna_compaction_trace: Table | None = None
    rna_binding_energy_trace: Table | None = None
    rg_rna_dense_values: np.ndarray | None = None
    steady_start_ns: float = 0.0
    interface_mode: str = "none"
    liquidity_metrics: dict[str, float] = field(default_factory=dict)

    @property
    def total_weight(self) -> int:
        return sum(max(1, rep.frame_weight) for rep in self.replicates)


@dataclass(frozen=True)
class RuntimeOptions:
    include_density: bool
    include_radial_density: bool
    include_shape_tensor: bool
    include_network: bool
    include_network_anchors: bool
    include_tda: bool
    export_density_dx: bool
    max_tda_points: int
    max_time_ns: float | None
    steady_state_subdir: str | None
    density_subdir: str | None
    radial_density_bin_nm: float
    radial_density_smooth_window: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate TITAN LLPS outputs for two systems (phos vs unphos), average three "
            "replicates per system on a common time grid, and export paper-ready figures."
        )
    )
    parser.add_argument(
        "--root-dir",
        default=".",
        help=(
            "Root directory containing system replicate folders such as "
            "10dimerphos_conf0_800 and 10dimerunphos_conf2_800."
        ),
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory where averaged tables, figures, and summary text will be written.",
    )
    parser.add_argument(
        "--ensemble-manifest",
        default=None,
        help=(
            "Read llps_ensemble_manifest.tsv for the LLPS paper bundle aggregator. "
            "Each row should include condition_id, replica_id, output_dir, temperature, "
            "salt, protein_conc, rna_conc, rna_length, and label."
        ),
    )
    parser.add_argument(
        "--fig7-from-averaged-dir",
        "--contact-dynamics-from-averaged-dir",
        dest="contact_dynamics_from_averaged_dir",
        default=None,
        help=(
            "Read existing phos_avg/unphos_avg dynamics CSV files from this directory and "
            "export contact-dynamics figures. Useful when raw replicate folders are not available."
        ),
    )
    parser.add_argument(
        "--group",
        action="append",
        default=[],
        metavar="SYSTEM=DIR1,DIR2,DIR3",
        help=(
            "Manual override for replicate grouping. Example: "
            "--group phos=run_a,run_b,run_c --group unphos=run_d,run_e,run_f"
        ),
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
        help="Raster export DPI for PNG/TIFF outputs. Default: 600",
    )
    parser.add_argument(
        "--steady-start-ns",
        default="auto",
        help=(
            "Steady-state analysis window start in ns, or 'auto' to infer it from percolation "
            "time / half-trajectory duration. Default: auto"
        ),
    )
    parser.add_argument(
        "--topn",
        type=int,
        default=20,
        help="Number of top residues / edges to show in interface and supplement figures. Default: 20",
    )
    parser.add_argument(
        "--protein-residue-min",
        type=int,
        default=1,
        help="Minimum protein residue index to show in Figure 3 sequence profile. Default: 1",
    )
    parser.add_argument(
        "--protein-residue-max",
        type=int,
        default=419,
        help="Maximum protein residue index to show in Figure 3 sequence profile. Default: 419",
    )
    parser.add_argument(
        "--protein-profile-smooth-window",
        type=int,
        default=9,
        help=(
            "Moving-average window used to smooth the Figure 3 protein residue profile. "
            "Use 1 to disable smoothing. Default: 9"
        ),
    )
    parser.add_argument(
        "--max-time-ns",
        type=float,
        default=2500.0,
        help=(
            "Maximum simulation time kept in averaged llps/dbscan time-series. "
            "Replicates shorter than this are used only up to their last frame; later points "
            "are averaged from the remaining replicate(s). Use 0 or a negative value to disable. "
            "Default: 2500"
        ),
    )
    parser.add_argument(
        "--skip-density",
        action="store_true",
        help="Skip loading density DX files. Figure 2 will be exported as a density placeholder.",
    )
    parser.add_argument(
        "--steady-state-subdir",
        default="1000ns_",
        help=(
            "Subdirectory inside each replicate folder for steady-state analyses "
            "(density, radial density, shape tensor, interface, polymer scaling, contact dynamics). "
            "Default: 1000ns_. Use 'none' to read all files from the replicate root."
        ),
    )
    parser.add_argument(
        "--density-subdir",
        default=None,
        help=(
            "Subdirectory inside each replicate folder where steady-state density DX files are preferred. "
            "Default: same as --steady-state-subdir. Use 'none' to read density files only from the replicate root."
        ),
    )
    parser.add_argument(
        "--no-radial-density",
        action="store_true",
        help="Disable radial density profile calculation from llps_density_protein/rna.dx.",
    )
    parser.add_argument(
        "--no-shape-tensor",
        action="store_true",
        help="Disable condensate shape tensor/asphericity calculation from density DX arrays.",
    )
    parser.add_argument(
        "--radial-density-bin-nm",
        type=float,
        default=1.0,
        help="Radial bin width in nm for protein/RNA density profiles. Default: 1.0",
    )
    parser.add_argument(
        "--radial-density-smooth-window",
        type=int,
        default=3,
        help="Moving-average window for radial density curves. Use 1 to disable. Default: 3",
    )
    parser.add_argument(
        "--radial-density-from-averaged-dir",
        default=None,
        help=(
            "Read existing phos_avg/unphos_avg radial density CSV files from this directory "
            "and export only radial density figures."
        ),
    )
    parser.add_argument(
        "--export-density-dx",
        action="store_true",
        help="Retain full 3D density grids and write averaged .dx files. This uses substantially more RAM.",
    )
    parser.add_argument(
        "--with-network",
        action="store_true",
        help="Load network/allostery files and build Supplementary Figure S2. Disabled by default to save RAM.",
    )
    parser.add_argument(
        "--with-network-anchors",
        action="store_true",
        help=(
            "Low-memory network mode: stream network_gc_edges.csv/network_te_edges.csv to compute "
            "dynamic anchor rankings without retaining full edge maps. Prefer this over --with-network "
            "on 16 GB RAM machines."
        ),
    )
    parser.add_argument(
        "--with-tda",
        action="store_true",
        help="Load TDA point clouds and build Supplementary Figure S3. Disabled by default to save RAM.",
    )
    parser.add_argument(
        "--max-tda-points",
        type=int,
        default=200000,
        help=(
            "Maximum pooled TDA points to retain per replicate when --with-tda is enabled. "
            "Use 0 to keep all points. Default: 200000"
        ),
    )
    return parser.parse_args()


def build_runtime_options(args: argparse.Namespace) -> RuntimeOptions:
    if args.skip_density and args.export_density_dx:
        raise SystemExit("--export-density-dx cannot be used together with --skip-density.")
    if args.max_tda_points < 0:
        raise SystemExit("--max-tda-points must be >= 0.")
    if args.radial_density_bin_nm <= 0.0:
        raise SystemExit("--radial-density-bin-nm must be > 0.")
    if args.radial_density_smooth_window < 1:
        raise SystemExit("--radial-density-smooth-window must be >= 1.")
    max_time_ns = float(args.max_time_ns) if float(args.max_time_ns) > 0.0 else None
    steady_state_subdir_raw = str(args.steady_state_subdir).strip()
    steady_state_subdir = None
    if steady_state_subdir_raw and steady_state_subdir_raw.lower() not in {"none", "root", "."}:
        steady_state_subdir = steady_state_subdir_raw.strip("/\\")
    density_subdir_raw = str(args.density_subdir).strip()
    density_subdir = None
    if not density_subdir_raw:
        density_subdir = steady_state_subdir
    elif density_subdir_raw.lower() not in {"none", "root", "."}:
        density_subdir = density_subdir_raw.strip("/\\")
    return RuntimeOptions(
        include_density=not args.skip_density,
        include_radial_density=(not args.skip_density and not args.no_radial_density),
        include_shape_tensor=(not args.skip_density and not args.no_shape_tensor),
        include_network=bool(args.with_network),
        include_network_anchors=bool(args.with_network_anchors or args.with_network),
        include_tda=bool(args.with_tda),
        export_density_dx=bool(args.export_density_dx),
        max_tda_points=int(args.max_tda_points),
        max_time_ns=max_time_ns,
        steady_state_subdir=steady_state_subdir,
        density_subdir=density_subdir,
        radial_density_bin_nm=float(args.radial_density_bin_nm),
        radial_density_smooth_window=int(args.radial_density_smooth_window),
    )


def natural_key(text: str) -> list[object]:
    parts = NATURAL_SORT_RE.split(text)
    key: list[object] = []
    for part in parts:
        if part.isdigit():
            key.append(int(part))
        else:
            key.append(part.lower())
    return key


def parse_ns_from_subdir_name(name: str | None) -> float | None:
    if not name:
        return None
    match = re.search(r"(\d+(?:\.\d+)?)\s*ns", name, flags=re.IGNORECASE)
    if not match:
        return None
    return float(match.group(1))


def read_text_lines(path: Path) -> list[str]:
    lines = path.read_text(encoding="utf-8-sig", errors="ignore").splitlines()
    if lines:
        lines[0] = lines[0].lstrip("\ufeff")
    return lines


def safe_float(text: str) -> float:
    stripped = text.strip()
    if not stripped or stripped.lower() in {"na", "nan", "null", "none"}:
        return float("nan")
    return float(stripped)


def safe_int(text: str) -> int:
    return int(float(text.strip()))


def unique_sorted(values: Iterable[float]) -> np.ndarray:
    arr = np.asarray(list(values), dtype=float)
    if arr.size == 0:
        return np.asarray([], dtype=float)
    return np.unique(np.round(arr, 12))


def set_publication_style(dpi: int) -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 180,
            "savefig.dpi": dpi,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
            "font.size": 8.0,
            "axes.titlesize": 9.0,
            "axes.labelsize": 9.0,
            "axes.labelweight": "normal",
            "axes.titleweight": "bold",
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "xtick.major.size": 3.2,
            "ytick.major.size": 3.2,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "legend.frameon": False,
            "legend.fontsize": 8.0,
            "grid.color": "#d9d9d9",
            "grid.linewidth": 0.6,
            "grid.alpha": 0.7,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.14,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
        ha="left",
    )


def add_panel_label_in_margin(fig: plt.Figure, ax: plt.Axes, label: str, *, fontsize: float = 12.0) -> None:
    bbox = ax.get_position()
    fig.text(
        bbox.x0 - 0.055,
        bbox.y1 + 0.006,
        label,
        fontsize=fontsize,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def sanitize_lag_ns(values: np.ndarray) -> np.ndarray:
    clean = np.asarray(values, dtype=float)
    if clean.size == 0:
        return clean
    if np.nanmedian(clean) < 0.0 or np.nanmax(clean) <= 0.0:
        clean = np.abs(clean)
    return clean


def simplify_node_label(raw: str, fallback_label: str | None = None) -> str:
    label = fallback_label or raw
    match = NODE_LABEL_RE.fullmatch(raw.strip())
    if not match:
        return label
    resname = match.group("resname")
    resid = match.group("resid")
    chain = match.group("chain") or ""
    compact = f"{resname}{resid}"
    if chain:
        compact += f"_{chain}"
    return compact


def pair_label(pair: tuple[str, str], nodes: dict[str, dict[str, str]]) -> str:
    left = simplify_node_label(pair[0], nodes.get(pair[0], {}).get("Label"))
    right = simplify_node_label(pair[1], nodes.get(pair[1], {}).get("Label"))
    return f"{left} vs {right}"


def load_whitespace_table(path: Path, columns: Sequence[str], x_name: str | None = None, x_label: str | None = None) -> Table | None:
    rows: list[list[float]] = []
    for raw in read_text_lines(path):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < len(columns):
            continue
        try:
            rows.append([safe_float(fields[idx]) for idx in range(len(columns))])
        except ValueError:
            continue
    if not rows:
        return None
    matrix = np.asarray(rows, dtype=float)
    data = {name: matrix[:, idx] for idx, name in enumerate(columns)}
    return Table(columns=data, order=list(columns), x_name=x_name, x_label=x_label)


def load_csv_numeric_table(path: Path, columns: Sequence[str], x_name: str | None = None, x_label: str | None = None) -> Table | None:
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            try:
                rows.append([safe_float(row[name]) for name in columns])
            except (ValueError, KeyError, TypeError):
                continue
    if not rows:
        return None
    matrix = np.asarray(rows, dtype=float)
    data = {name: matrix[:, idx] for idx, name in enumerate(columns)}
    return Table(columns=data, order=list(columns), x_name=x_name, x_label=x_label)


def load_csv_mixed_table(
    path: Path,
    columns: Sequence[str],
    *,
    numeric_columns: Sequence[str],
    x_name: str | None = None,
    x_label: str | None = None,
) -> Table | None:
    numeric_set = set(numeric_columns)
    collected: dict[str, list[object]] = {name: [] for name in columns}
    with path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            parsed_row: dict[str, object] = {}
            valid = True
            for name in columns:
                raw_value = row.get(name)
                if raw_value is None:
                    valid = False
                    break
                raw_value = raw_value.strip()
                if name in numeric_set:
                    try:
                        parsed_row[name] = safe_float(raw_value)
                    except ValueError:
                        valid = False
                        break
                else:
                    parsed_row[name] = raw_value
            if not valid:
                continue
            for name in columns:
                collected[name].append(parsed_row[name])
    if not collected[columns[0]]:
        return None
    data: dict[str, np.ndarray] = {}
    for name in columns:
        if name in numeric_set:
            data[name] = np.asarray(collected[name], dtype=float)
        else:
            data[name] = np.asarray(collected[name], dtype=object)
    return Table(columns=data, order=list(columns), x_name=x_name or columns[0], x_label=x_label)


def load_summary(path: Path) -> dict[str, object]:
    summary: dict[str, object] = {}
    for raw in read_text_lines(path):
        line = raw.strip()
        if not line or ":" not in line:
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip()
        if key == "percolated":
            summary[key] = value.lower() == "yes"
            continue
        try:
            summary[key] = safe_float(value)
        except ValueError:
            summary[key] = value
    return summary


def load_lifetime_distribution(path: Path) -> DynamicsBundle | None:
    summary: dict[str, float] = {}
    rows: list[tuple[float, float]] = []
    for raw in read_text_lines(path):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            fields = line[1:].strip().split()
            if len(fields) >= 2:
                key = fields[0].strip().lower()
                try:
                    summary[key] = safe_float(fields[1])
                except ValueError:
                    pass
            continue
        fields = line.split()
        if len(fields) < 2:
            continue
        try:
            rows.append((safe_float(fields[0]), safe_float(fields[1])))
        except ValueError:
            continue
    table = None
    if rows:
        matrix = np.asarray(rows, dtype=float)
        order = np.argsort(matrix[:, 0])
        matrix = matrix[order]
        table = Table(
            columns={"lifetime_ns": matrix[:, 0], "count": matrix[:, 1]},
            order=["lifetime_ns", "count"],
            x_name="lifetime_ns",
            x_label="Lifetime (ns)",
        )
    if not summary and table is None:
        return None
    return DynamicsBundle(summary=summary, table=table)


def load_keyvalue_timeseries(path: Path, columns: Sequence[str], x_name: str = "time_ns") -> DynamicsBundle | None:
    summary: dict[str, float] = {}
    rows: list[list[float]] = []
    for raw in read_text_lines(path):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            continue
        if ":" in line:
            key, value = line.split(":", 1)
            try:
                summary[key.strip()] = safe_float(value)
            except ValueError:
                pass
            continue
        fields = line.split()
        if len(fields) < len(columns):
            continue
        try:
            rows.append([safe_float(fields[idx]) for idx in range(len(columns))])
        except ValueError:
            continue
    table = None
    if rows:
        matrix = np.asarray(rows, dtype=float)
        order = np.argsort(matrix[:, 0])
        matrix = matrix[order]
        table = Table(
            columns={name: matrix[:, idx] for idx, name in enumerate(columns)},
            order=list(columns),
            x_name=x_name,
            x_label="Time (ns)",
        )
    if not summary and table is None:
        return None
    return DynamicsBundle(summary=summary, table=table)


def load_dense_msd(path: Path) -> DynamicsBundle | None:
    summary: dict[str, float] = {}
    rows: list[tuple[float, float, float]] = []
    for raw in read_text_lines(path):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            fields = line[1:].strip().split()
            if not fields:
                continue
            key = fields[0].strip()
            try:
                if key == "diffusion_nm2_per_ns" and len(fields) >= 2:
                    summary[key] = safe_float(fields[1])
                elif key == "fit_window_ns" and len(fields) >= 3:
                    summary["fit_start_ns"] = safe_float(fields[1])
                    summary["fit_end_ns"] = safe_float(fields[2])
            except ValueError:
                pass
            continue
        fields = line.split()
        if len(fields) < 3:
            continue
        try:
            rows.append((safe_float(fields[0]), safe_float(fields[1]), safe_float(fields[2])))
        except ValueError:
            continue
    table = None
    if rows:
        matrix = np.asarray(rows, dtype=float)
        order = np.argsort(matrix[:, 0])
        matrix = matrix[order]
        table = Table(
            columns={"lag_ns": matrix[:, 0], "msd_nm2": matrix[:, 1], "count": matrix[:, 2]},
            order=["lag_ns", "msd_nm2", "count"],
            x_name="lag_ns",
            x_label="Lag time (ns)",
        )
    if not summary and table is None:
        return None
    return DynamicsBundle(summary=summary, table=table)


def load_saved_dynamics_csv(path: Path) -> DynamicsBundle | None:
    if not path.exists():
        return None
    summary: dict[str, float] = {}
    header: list[str] | None = None
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        reader = csv.reader(handle)
        for raw_row in reader:
            row = [cell.strip() for cell in raw_row]
            if not row or not any(row):
                continue
            first = row[0]
            if first.startswith("#"):
                key = first[1:].strip()
                if len(row) >= 2:
                    try:
                        summary[key] = safe_float(row[1])
                    except ValueError:
                        pass
                continue
            if header is None:
                header = row
                continue
            if len(row) < len(header):
                continue
            try:
                rows.append([safe_float(row[idx]) for idx in range(len(header))])
            except ValueError:
                continue

    table = None
    if header is not None and rows:
        matrix = np.asarray(rows, dtype=float)
        x_name = "time_ns" if "time_ns" in header else ("lifetime_ns" if "lifetime_ns" in header else header[0])
        x_idx = header.index(x_name) if x_name in header else 0
        order_idx = np.argsort(matrix[:, x_idx])
        matrix = matrix[order_idx]
        table = Table(
            columns={name: matrix[:, idx] for idx, name in enumerate(header)},
            order=list(header),
            x_name=x_name,
            x_label="Time (ns)" if x_name == "time_ns" else "Lifetime (ns)",
        )
    if not summary and table is None:
        return None
    return DynamicsBundle(summary=summary, table=table)


def normalize_optional_subdir(raw_value: str | None) -> str | None:
    if raw_value is None:
        return None
    value = str(raw_value).strip()
    if not value or value.lower() in {"none", "root", "."}:
        return None
    return value.strip("/\\")


def load_averaged_interface_systems(avg_dir: Path) -> dict[str, AggregatedSystem]:
    systems: dict[str, AggregatedSystem] = {}
    for system in SYSTEM_ORDER:
        prefix = system
        aggregated = AggregatedSystem(system=system, replicates=[])
        aggregated.sticker_contact_lifetimes = load_saved_dynamics_csv(
            avg_dir / f"{prefix}_avg_llps_sticker_contact_lifetimes.csv"
        )
        aggregated.protein_rna_contact_lifetimes = load_saved_dynamics_csv(
            avg_dir / f"{prefix}_avg_llps_protein_rna_contact_lifetimes.csv"
        )
        aggregated.neighbor_rewiring = load_saved_dynamics_csv(avg_dir / f"{prefix}_avg_llps_neighbor_rewiring.csv")
        aggregated.interface_turnover = load_saved_dynamics_csv(avg_dir / f"{prefix}_avg_llps_interface_turnover.csv")
        systems[system] = aggregated
    return systems


def load_nucleic_interactions(path: Path) -> list[NucleicInteraction]:
    entries: list[NucleicInteraction] = []
    with path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            try:
                entries.append(
                    NucleicInteraction(
                        resnum=safe_int(row["Protein_ResNum"]),
                        resname=row["ResName"].strip(),
                        base_frames=safe_float(row["Base_Contact_Frames"]),
                        backbone_frames=safe_float(row["Backbone_Contact_Frames"]),
                        base_occupancy=safe_float(row["Base_Occupancy_%"]),
                        backbone_occupancy=safe_float(row["Backbone_Occupancy_%"]),
                        recognition_type=row["Recognition_Type"].strip(),
                    )
                )
            except (KeyError, ValueError):
                continue
    return entries


def load_sidechain_contacts(path: Path) -> list[SidechainContact]:
    contacts: list[SidechainContact] = []
    for raw in read_text_lines(path):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 4:
            continue
        try:
            contacts.append(
                SidechainContact(
                    chain=fields[0],
                    resnum=safe_int(fields[1]),
                    resname=fields[2],
                    score=safe_float(fields[3]),
                )
            )
        except ValueError:
            continue
    return contacts


def load_pair_value_csv(path: Path, value_column: str) -> dict[tuple[str, str], float]:
    pairs: dict[tuple[str, str], float] = {}
    with path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            try:
                source = row["Source"].strip()
                target = row["Target"].strip()
                value = safe_float(row[value_column])
            except (KeyError, ValueError):
                continue
            pairs[(source, target)] = value
    return pairs


def load_network_nodes(path: Path) -> dict[str, dict[str, str]]:
    nodes: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            node_id = row.get("Id", "").strip()
            if node_id:
                nodes[node_id] = row
    return nodes


def load_network_bundle(directory: Path) -> NetworkBundle | None:
    bundle = NetworkBundle()

    node_path = directory / "network_nodes.csv"
    corr_path = directory / "network_edges.csv"
    gc_path = directory / "network_gc_edges.csv"
    te_path = directory / "network_te_edges.csv"
    mi_path = directory / "raw_mi_edges.csv"
    corr_knn_gc = directory / "corr_knn_gc_pairs.csv"
    corr_knn_te = directory / "corr_knn_te_pairs.csv"
    corr_knn_mi = directory / "corr_knn_mi_pairs.csv"

    if node_path.exists():
        bundle.nodes = load_network_nodes(node_path)
    if corr_path.exists():
        bundle.corr_edges = load_pair_value_csv(corr_path, "Correlation")
    if gc_path.exists():
        bundle.gc_edges = load_pair_value_csv(gc_path, "GC")
    elif corr_knn_gc.exists():
        bundle.gc_edges = load_pair_value_csv(corr_knn_gc, "GC")
    if te_path.exists():
        bundle.te_edges = load_pair_value_csv(te_path, "Net_TE")
    elif corr_knn_te.exists():
        bundle.te_edges = load_pair_value_csv(corr_knn_te, "Net_TE")
    if mi_path.exists():
        bundle.mi_edges = load_pair_value_csv(mi_path, "MI")
    elif corr_knn_mi.exists():
        bundle.mi_edges = load_pair_value_csv(corr_knn_mi, "MI")

    return bundle if bundle.available() else None


def load_network_anchor_labels(path: Path) -> dict[str, str]:
    labels: dict[str, str] = {}
    with path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            node_id = row.get("Id", "").strip()
            if not node_id:
                continue
            label = row.get("Label", "").strip() or None
            labels[node_id] = simplify_node_label(node_id, label)
    return labels


def stream_gc_anchor_scores(path: Path, value_column: str) -> dict[str, float]:
    scores: dict[str, float] = defaultdict(float)
    with path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            source = row.get("Source", "").strip()
            target = row.get("Target", "").strip()
            if not source or not target:
                continue
            try:
                value = abs(safe_float(row.get(value_column) or ""))
            except ValueError:
                continue
            if not np.isfinite(value) or value <= 0.0:
                continue
            scores[source] += value
            scores[target] += value
    return dict(scores)


def stream_te_anchor_scores(path: Path, value_column: str) -> tuple[dict[str, float], dict[str, float]]:
    outgoing: dict[str, float] = defaultdict(float)
    incoming: dict[str, float] = defaultdict(float)
    with path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            source = row.get("Source", "").strip()
            target = row.get("Target", "").strip()
            if not source or not target:
                continue
            try:
                net = safe_float(row.get(value_column) or "")
            except ValueError:
                continue
            if not np.isfinite(net) or net == 0.0:
                continue
            if net > 0.0:
                outgoing[source] += net
                incoming[target] += net
            else:
                outgoing[target] += -net
                incoming[source] += -net
    nodes = set(outgoing) | set(incoming)
    source_scores = {
        node: score
        for node in nodes
        if (score := max(0.0, outgoing.get(node, 0.0) - incoming.get(node, 0.0))) > 0.0
    }
    sink_scores = {
        node: score
        for node in nodes
        if (score := max(0.0, incoming.get(node, 0.0) - outgoing.get(node, 0.0))) > 0.0
    }
    return source_scores, sink_scores


def load_network_anchor_bundle(directory: Path | None) -> NetworkAnchorBundle | None:
    if directory is None:
        return None
    bundle = NetworkAnchorBundle()

    node_path = directory / "network_nodes.csv"
    gc_path = directory / "network_gc_edges.csv"
    te_path = directory / "network_te_edges.csv"
    corr_knn_gc = directory / "corr_knn_gc_pairs.csv"
    corr_knn_te = directory / "corr_knn_te_pairs.csv"

    if node_path.exists():
        bundle.labels = load_network_anchor_labels(node_path)
    if gc_path.exists():
        bundle.gc_anchor_scores = stream_gc_anchor_scores(gc_path, "GC")
    elif corr_knn_gc.exists():
        bundle.gc_anchor_scores = stream_gc_anchor_scores(corr_knn_gc, "GC")
    if te_path.exists():
        bundle.te_source_scores, bundle.te_sink_scores = stream_te_anchor_scores(te_path, "Net_TE")
    elif corr_knn_te.exists():
        bundle.te_source_scores, bundle.te_sink_scores = stream_te_anchor_scores(corr_knn_te, "Net_TE")

    return bundle if bundle.available() or bundle.labels else None


def load_tda_points(path: Path, max_points: int = 0) -> np.ndarray | None:
    rows: list[list[float]] = []
    rng = np.random.default_rng(0)
    seen = 0
    with path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            try:
                point = [safe_float(row["X"]), safe_float(row["Y"]), safe_float(row["Z"])]
            except (KeyError, ValueError):
                continue
            seen += 1
            if max_points <= 0 or len(rows) < max_points:
                rows.append(point)
                continue
            swap_idx = int(rng.integers(0, seen))
            if swap_idx < max_points:
                rows[swap_idx] = point
    if not rows:
        return None
    return np.asarray(rows, dtype=np.float32)


def load_dx_grid(path: Path, keep_full_data: bool = False) -> DXGrid | None:
    counts: tuple[int, int, int] | None = None
    origin: tuple[float, float, float] | None = None
    deltas: list[tuple[float, float, float]] = []
    expected_items: int | None = None
    label = path.stem
    reading_values = False
    projection: np.ndarray | None = None
    flat_values: np.ndarray | None = None
    value_idx = 0

    with path.open("r", encoding="utf-8-sig", errors="ignore") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("object 1 class gridpositions counts"):
                parts = line.split()
                counts = (int(parts[-3]), int(parts[-2]), int(parts[-1]))
                continue
            if line.startswith("origin "):
                parts = line.split()
                origin = (float(parts[1]), float(parts[2]), float(parts[3]))
                continue
            if line.startswith("delta "):
                parts = line.split()
                deltas.append((float(parts[1]), float(parts[2]), float(parts[3])))
                continue
            if "data follows" in line and line.startswith("object 3 class array"):
                items_match = re.search(r"items\s+(\d+)", line)
                if items_match:
                    expected_items = int(items_match.group(1))
                    reading_values = True
                    if counts is None:
                        return None
                    projection = np.zeros((counts[0], counts[1]), dtype=np.float64)
                    if keep_full_data:
                        flat_values = np.empty(expected_items, dtype=np.float32)
                continue
            if reading_values and line.startswith("object "):
                label_match = re.search(r'"([^"]+)"', line)
                if label_match:
                    label = label_match.group(1)
                break
            if not reading_values:
                continue
            if counts is None or expected_items is None or projection is None:
                return None
            _, ny, nz = counts
            for token in line.split():
                if value_idx >= expected_items:
                    break
                try:
                    value = float(token)
                except ValueError:
                    continue
                if flat_values is not None:
                    flat_values[value_idx] = value
                x_idx = value_idx // (ny * nz)
                y_idx = (value_idx // nz) % ny
                projection[x_idx, y_idx] += value
                value_idx += 1

    if counts is None or origin is None or len(deltas) != 3 or expected_items is None or projection is None:
        return None

    if value_idx < expected_items:
        return None

    data = flat_values.reshape(counts) if flat_values is not None else None
    return DXGrid(
        counts=counts,
        origin=origin,
        deltas=(deltas[0], deltas[1], deltas[2]),
        label=label,
        data=data,
        projection_xy_data=projection.astype(np.float32),
    )


def radial_density_profile_from_dx(grid: DXGrid, bin_width_nm: float) -> Table | None:
    if grid.data is None or grid.data.size == 0 or bin_width_nm <= 0.0:
        return None
    data = np.asarray(grid.data, dtype=np.float64)
    if not np.any(np.isfinite(data)):
        return None

    nx, ny, nz = grid.counts
    x = float(grid.origin[0]) + (np.arange(nx, dtype=float) + 0.5) * grid.spacing_x
    y = float(grid.origin[1]) + (np.arange(ny, dtype=float) + 0.5) * grid.spacing_y
    z = float(grid.origin[2]) + (np.arange(nz, dtype=float) + 0.5) * grid.spacing_z

    positive = np.where(np.isfinite(data) & (data > 0.0), data, 0.0)
    total_density = float(np.sum(positive, dtype=np.float64))
    if total_density > 0.0:
        x_weight = np.sum(positive, axis=(1, 2), dtype=np.float64)
        y_weight = np.sum(positive, axis=(0, 2), dtype=np.float64)
        z_weight = np.sum(positive, axis=(0, 1), dtype=np.float64)
        center = np.asarray(
            [
                float(np.dot(x, x_weight) / total_density),
                float(np.dot(y, y_weight) / total_density),
                float(np.dot(z, z_weight) / total_density),
            ],
            dtype=float,
        )
    else:
        center = np.asarray([(x[0] + x[-1]) * 0.5, (y[0] + y[-1]) * 0.5, (z[0] + z[-1]) * 0.5], dtype=float)

    y2z2 = (y[:, None] - center[1]) ** 2 + (z[None, :] - center[2]) ** 2
    corner_distances: list[float] = []
    for xx in (x[0], x[-1]):
        for yy in (y[0], y[-1]):
            for zz in (z[0], z[-1]):
                corner_distances.append(float(np.linalg.norm(np.asarray([xx, yy, zz], dtype=float) - center)))
    n_bins = max(1, int(math.ceil(max(corner_distances) / bin_width_nm)))
    density_sum = np.zeros(n_bins, dtype=np.float64)
    voxel_count = np.zeros(n_bins, dtype=np.float64)

    for ix, xx in enumerate(x):
        radii = np.sqrt((xx - center[0]) ** 2 + y2z2)
        bin_index = np.floor(radii / bin_width_nm).astype(np.int64, copy=False)
        mask = (bin_index >= 0) & (bin_index < n_bins)
        if not np.any(mask):
            continue
        flat_bins = bin_index[mask].ravel()
        flat_values = data[ix][mask].ravel()
        density_sum += np.bincount(flat_bins, weights=flat_values, minlength=n_bins)
        voxel_count += np.bincount(flat_bins, minlength=n_bins)

    valid = voxel_count > 0
    if not np.any(valid):
        return None
    r_nm = (np.arange(n_bins, dtype=float) + 0.5) * bin_width_nm
    density = np.divide(density_sum, voxel_count, out=np.zeros_like(density_sum), where=valid)
    return Table(
        columns={
            "r_nm": r_nm[valid],
            "density_nm3": density[valid],
            "voxel_count": voxel_count[valid],
        },
        order=["r_nm", "density_nm3", "voxel_count"],
        x_name="r_nm",
        x_label="Distance from condensate center, r (nm)",
    )


def shape_tensor_metrics_from_dx(grid: DXGrid) -> dict[str, float]:
    if grid.data is None or grid.data.size == 0:
        return {}
    data = np.asarray(grid.data, dtype=np.float64)
    density = np.where(np.isfinite(data) & (data > 0.0), data, 0.0)
    total_density = float(np.sum(density, dtype=np.float64))
    if total_density <= 0.0:
        return {}

    nx, ny, nz = grid.counts
    x = float(grid.origin[0]) + (np.arange(nx, dtype=float) + 0.5) * grid.spacing_x
    y = float(grid.origin[1]) + (np.arange(ny, dtype=float) + 0.5) * grid.spacing_y
    z = float(grid.origin[2]) + (np.arange(nz, dtype=float) + 0.5) * grid.spacing_z
    x_weight = np.sum(density, axis=(1, 2), dtype=np.float64)
    y_weight = np.sum(density, axis=(0, 2), dtype=np.float64)
    z_weight = np.sum(density, axis=(0, 1), dtype=np.float64)
    center = np.asarray(
        [
            float(np.dot(x, x_weight) / total_density),
            float(np.dot(y, y_weight) / total_density),
            float(np.dot(z, z_weight) / total_density),
        ],
        dtype=float,
    )

    dy = y - center[1]
    dz = z - center[2]
    dy_grid = dy[:, None]
    dz_grid = dz[None, :]
    tensor = np.zeros((3, 3), dtype=np.float64)
    for ix, xx in enumerate(x):
        plane = density[ix]
        dx = float(xx - center[0])
        plane_sum = float(np.sum(plane, dtype=np.float64))
        tensor[0, 0] += dx * dx * plane_sum
        tensor[0, 1] += dx * float(np.sum(plane * dy_grid, dtype=np.float64))
        tensor[0, 2] += dx * float(np.sum(plane * dz_grid, dtype=np.float64))
        tensor[1, 1] += float(np.sum(plane * dy_grid * dy_grid, dtype=np.float64))
        tensor[1, 2] += float(np.sum(plane * dy_grid * dz_grid, dtype=np.float64))
        tensor[2, 2] += float(np.sum(plane * dz_grid * dz_grid, dtype=np.float64))
    tensor[1, 0] = tensor[0, 1]
    tensor[2, 0] = tensor[0, 2]
    tensor[2, 1] = tensor[1, 2]
    tensor /= total_density

    eigenvalues = np.linalg.eigvalsh(tensor)
    eigenvalues = np.sort(np.maximum(eigenvalues, 0.0))[::-1]
    lam1, lam2, lam3 = (float(eigenvalues[0]), float(eigenvalues[1]), float(eigenvalues[2]))
    rg2 = lam1 + lam2 + lam3
    if rg2 <= 0.0:
        relative_shape_anisotropy = 0.0
        normalized_asphericity = 0.0
        aspect_ratio = 0.0
    else:
        mean_lam = rg2 / 3.0
        relative_shape_anisotropy = float(1.5 * np.sum((eigenvalues - mean_lam) ** 2) / (rg2 * rg2))
        normalized_asphericity = (lam1 - 0.5 * (lam2 + lam3)) / rg2
        aspect_ratio = math.sqrt(lam1 / lam3) if lam3 > 0.0 else float("inf")
    return {
        "center_x_nm": float(center[0]),
        "center_y_nm": float(center[1]),
        "center_z_nm": float(center[2]),
        "lambda1_nm2": lam1,
        "lambda2_nm2": lam2,
        "lambda3_nm2": lam3,
        "principal_radius1_nm": math.sqrt(lam1),
        "principal_radius2_nm": math.sqrt(lam2),
        "principal_radius3_nm": math.sqrt(lam3),
        "rg_nm": math.sqrt(rg2),
        "asphericity_nm2": lam1 - 0.5 * (lam2 + lam3),
        "normalized_asphericity": normalized_asphericity,
        "acylindricity_nm2": lam2 - lam3,
        "relative_shape_anisotropy": relative_shape_anisotropy,
        "aspect_ratio": aspect_ratio,
        "density_mass": total_density,
    }


def maybe_attach_radial_profile(grid: DXGrid | None, options: RuntimeOptions) -> DXGrid | None:
    if grid is None:
        return grid
    if options.include_radial_density:
        grid.radial_profile_data = radial_density_profile_from_dx(grid, options.radial_density_bin_nm)
    if options.include_shape_tensor:
        grid.shape_metrics = shape_tensor_metrics_from_dx(grid)
    if not options.export_density_dx:
        grid.data = None
    return grid


def load_saved_table_csv(path: Path, x_name: str | None = None) -> Table | None:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            return None
        header = [item.strip() for item in header]
        rows: list[list[float]] = []
        for row in reader:
            if len(row) < len(header):
                continue
            try:
                rows.append([safe_float(row[idx]) for idx in range(len(header))])
            except ValueError:
                continue
    if not header or not rows:
        return None
    matrix = np.asarray(rows, dtype=float)
    return Table(
        columns={name: matrix[:, idx] for idx, name in enumerate(header)},
        order=list(header),
        x_name=x_name or header[0],
        x_label="Distance from condensate center, r (nm)" if (x_name or header[0]) == "r_nm" else None,
    )


def load_averaged_radial_density_systems(avg_dir: Path) -> dict[str, AggregatedSystem]:
    systems: dict[str, AggregatedSystem] = {}
    for system in SYSTEM_ORDER:
        aggregated = AggregatedSystem(system=system, replicates=[])
        protein_profile = load_saved_table_csv(avg_dir / f"{system}_avg_radial_density_protein.csv", x_name="r_nm")
        rna_profile = load_saved_table_csv(avg_dir / f"{system}_avg_radial_density_rna.csv", x_name="r_nm")
        if protein_profile is not None:
            aggregated.density_protein = DXGrid((0, 0, 0), (0.0, 0.0, 0.0), ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)), "protein")
            aggregated.density_protein.radial_profile_data = protein_profile
        if rna_profile is not None:
            aggregated.density_rna = DXGrid((0, 0, 0), (0.0, 0.0, 0.0), ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)), "rna")
            aggregated.density_rna.radial_profile_data = rna_profile
        systems[system] = aggregated
    return systems


def resolve_density_file(directory: Path, filename: str, options: RuntimeOptions) -> Path | None:
    candidates: list[Path] = []
    if options.density_subdir:
        candidates.append(directory / options.density_subdir / filename)
    candidates.append(directory / filename)
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def resolve_steady_file(directory: Path, filename: str, options: RuntimeOptions) -> Path | None:
    candidates: list[Path] = []
    if options.steady_state_subdir:
        candidates.append(directory / options.steady_state_subdir / filename)
    candidates.append(directory / filename)
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def resolve_steady_dir_for_any(directory: Path, filenames: Sequence[str], options: RuntimeOptions) -> Path:
    if options.steady_state_subdir:
        candidate_dir = directory / options.steady_state_subdir
        if any((candidate_dir / filename).exists() for filename in filenames):
            return candidate_dir
    return directory


def discover_groups(root_dir: Path) -> dict[str, list[tuple[int, Path]]]:
    grouped: dict[str, list[tuple[int, Path]]] = {system: [] for system in SYSTEM_ORDER}
    for candidate in root_dir.rglob("*"):
        if not candidate.is_dir():
            continue
        match = DISCOVERY_RE.fullmatch(candidate.name)
        if not match:
            continue
        system = match.group("system").lower()
        conf = int(match.group("conf"))
        grouped[system].append((conf, candidate))
    return grouped


def parse_group_overrides(raw_groups: Sequence[str], base_dir: Path) -> dict[str, list[tuple[int, Path]]]:
    grouped: dict[str, list[tuple[int, Path]]] = {}
    for entry in raw_groups:
        if "=" not in entry:
            raise SystemExit(f"Invalid --group value '{entry}'. Expected format SYSTEM=DIR1,DIR2,DIR3")
        key, value = entry.split("=", 1)
        system = key.strip().lower()
        if system not in SYSTEM_ORDER:
            raise SystemExit(f"Unsupported system '{system}' in --group. Expected phos or unphos.")
        paths = [item.strip() for item in value.split(",") if item.strip()]
        if len(paths) != 3:
            raise SystemExit(f"--group {system}=... must contain exactly 3 directories.")
        resolved: list[tuple[int, Path]] = []
        for idx, raw_path in enumerate(paths):
            directory = Path(raw_path)
            if not directory.is_absolute():
                directory = (base_dir / directory).resolve()
            else:
                directory = directory.resolve()
            resolved.append((idx, directory))
        grouped[system] = resolved
    return grouped


def resolve_replicate_dirs(root_dir: Path, group_args: Sequence[str]) -> dict[str, list[tuple[int, Path]]]:
    grouped = parse_group_overrides(group_args, root_dir) if group_args else discover_groups(root_dir)
    for system in SYSTEM_ORDER:
        dirs = grouped.get(system, [])
        if len(dirs) != 3:
            found = ", ".join(path.name for _, path in sorted(dirs, key=lambda item: item[0])) or "none"
            raise SystemExit(
                f"Expected exactly 3 replicate directories for system '{system}', found {len(dirs)} ({found})."
            )
        seen_ids: set[int] = set()
        for conf, path in dirs:
            if conf in seen_ids:
                raise SystemExit(f"Duplicate replicate id conf{conf} detected for system '{system}'.")
            if not path.exists() or not path.is_dir():
                raise SystemExit(f"Replicate directory does not exist: {path}")
            seen_ids.add(conf)
        grouped[system] = sorted(dirs, key=lambda item: item[0])
    return grouped


def load_replicate(system: str, conf_id: int, directory: Path, options: RuntimeOptions) -> ReplicateData:
    llps_autocorr_path = resolve_steady_file(directory, "llps_autocorr.dat", options)
    llps_degree_path = resolve_steady_file(directory, "llps_degree_dist.dat", options)
    polymer_scaling_path = resolve_steady_file(directory, "polymer_scaling.dat", options)
    summary_path = resolve_steady_file(directory, "llps_summary.dat", options)
    nucleic_interactions_path = resolve_steady_file(directory, "nucleic_acid_interactions.csv", options)
    sidechain_protein_path = resolve_steady_file(directory, "out_list_res_sidechain_contact_Protein.dat", options)
    sidechain_rna_path = resolve_steady_file(directory, "out_list_res_sidechain_contact_RNA.dat", options)
    sticker_contact_lifetimes_path = resolve_steady_file(directory, "llps_sticker_contact_lifetimes.dat", options)
    protein_rna_contact_lifetimes_path = resolve_steady_file(directory, "llps_protein_rna_contact_lifetimes.dat", options)
    neighbor_rewiring_path = resolve_steady_file(directory, "llps_neighbor_rewiring.dat", options)
    interface_turnover_path = resolve_steady_file(directory, "llps_interface_turnover.dat", options)
    dense_msd_protein_path = resolve_steady_file(directory, "llps_dense_msd_protein.dat", options)
    dense_msd_rna_path = resolve_steady_file(directory, "llps_dense_msd_rna.dat", options)
    liquidity_summary_path = resolve_steady_file(directory, "llps_liquidity_summary.dat", options)
    rna_compaction_summary_path = resolve_steady_file(
        directory, "llps_rna_compaction_binding_summary.dat", options
    )
    rna_compaction_binned_path = resolve_steady_file(
        directory, "llps_rna_compaction_binding_binned.csv", options
    )
    rna_energy_compaction_summary_path = resolve_steady_file(
        directory, "llps_rna_energy_compaction_summary.dat", options
    )
    rna_energy_compaction_binned_path = resolve_steady_file(
        directory, "llps_rna_energy_compaction_binned.csv", options
    )
    rna_information_summary_path = resolve_steady_file(
        directory, "llps_rna_information_summary.dat", options
    )
    rna_compaction_state_histogram_path = resolve_steady_file(
        directory, "llps_rna_compaction_state_histogram.csv", options
    )
    rna_topology_summary_path = resolve_steady_file(
        directory, "llps_rna_topology_summary.dat", options
    )
    rna_topology_binned_path = resolve_steady_file(
        directory, "llps_rna_topology_binding_binned.csv", options
    )
    rna_topology_trace_path = resolve_steady_file(
        directory, "llps_rna_topology_trace.csv", options
    )
    rna_chain_compaction_summary_path = resolve_steady_file(
        directory, "llps_rna_chain_compaction_summary.csv", options
    )
    rna_intrachain_contacts_summary_path = resolve_steady_file(
        directory, "llps_rna_intrachain_contacts_summary.csv", options
    )
    rna_compaction_trace_path = resolve_steady_file(
        directory, "llps_rna_compaction_trace.csv", options
    )
    rna_binding_energy_trace_path = resolve_steady_file(
        directory, "llps_rna_binding_energy_trace.csv", options
    )
    network_dir = resolve_steady_dir_for_any(
        directory,
        ["network_edges.csv", "network_gc_edges.csv", "network_te_edges.csv", "raw_mi_edges.csv", "network_nodes.csv"],
        options,
    )
    tda_path = resolve_steady_file(directory, "tda_global_pointcloud.csv", options)

    llps_timeseries = load_whitespace_table(
        directory / "llps_timeseries.dat",
        ["frame", "time_ns", "gcc_fraction", "gcc_size", "gcc_rg_nm", "mean_degree", "rho_dense", "rho_dilute"],
        x_name="time_ns",
        x_label="Time (ns)",
    ) if (directory / "llps_timeseries.dat").exists() else None

    llps_autocorr = load_whitespace_table(
        llps_autocorr_path,
        ["lag_ns", "ct"],
        x_name="lag_ns",
        x_label="Lag (ns)",
    ) if llps_autocorr_path is not None else None
    if llps_autocorr is not None:
        llps_autocorr.columns["lag_ns"] = sanitize_lag_ns(llps_autocorr.columns["lag_ns"])
        order = np.argsort(llps_autocorr.columns["lag_ns"])
        for name in llps_autocorr.order:
            llps_autocorr.columns[name] = llps_autocorr.columns[name][order]

    llps_degree = load_whitespace_table(
        llps_degree_path,
        ["k", "pk_mean", "pk_std"],
        x_name="k",
        x_label="k",
    ) if llps_degree_path is not None else None

    dbscan_timeseries = load_csv_numeric_table(
        directory / "llps_dbscan_timeseries.csv",
        [
            "Time_ps",
            "TotalBeads",
            "NumCondensates",
            "GiantCluster_ProChains",
            "GiantCluster_RnaChains",
            "NoiseBeads",
            "Rg_Pro_Dilute",
            "Rg_Pro_Dense",
            "Rg_RNA_Dilute",
            "Rg_RNA_Dense",
        ],
        x_name="Time_ps",
        x_label="Time (ps)",
    ) if (directory / "llps_dbscan_timeseries.csv").exists() else None

    polymer_scaling = load_whitespace_table(
        polymer_scaling_path,
        ["s", "mean_R2_nm2", "R_rms_nm", "count"],
        x_name="s",
        x_label="Sequence separation s",
    ) if polymer_scaling_path is not None else None

    summary = load_summary(summary_path) if summary_path is not None else {}
    density_protein_path = resolve_density_file(directory, "llps_density_protein.dx", options) if options.include_density else None
    density_protein = (
        load_dx_grid(
            density_protein_path,
            keep_full_data=(options.export_density_dx or options.include_radial_density or options.include_shape_tensor),
        )
        if density_protein_path is not None
        else None
    )
    density_protein = maybe_attach_radial_profile(density_protein, options)
    density_rna_path = resolve_density_file(directory, "llps_density_rna.dx", options) if options.include_density else None
    density_rna = (
        load_dx_grid(
            density_rna_path,
            keep_full_data=(options.export_density_dx or options.include_radial_density or options.include_shape_tensor),
        )
        if density_rna_path is not None
        else None
    )
    density_rna = maybe_attach_radial_profile(density_rna, options)
    nucleic_interactions = load_nucleic_interactions(nucleic_interactions_path) if nucleic_interactions_path is not None else []
    sidechain_protein = load_sidechain_contacts(sidechain_protein_path) if sidechain_protein_path is not None else []
    sidechain_rna = load_sidechain_contacts(sidechain_rna_path) if sidechain_rna_path is not None else []
    network_bundle = (
        load_network_bundle(network_dir)
        if options.include_network and network_dir is not None
        else None
    )
    network_anchor_bundle = (
        load_network_anchor_bundle(network_dir)
        if options.include_network_anchors and network_dir is not None
        else None
    )
    tda_points = (
        load_tda_points(tda_path, max_points=options.max_tda_points)
        if options.include_tda and tda_path is not None
        else None
    )
    sticker_contact_lifetimes = (
        load_lifetime_distribution(sticker_contact_lifetimes_path)
        if sticker_contact_lifetimes_path is not None
        else None
    )
    protein_rna_contact_lifetimes = (
        load_lifetime_distribution(protein_rna_contact_lifetimes_path)
        if protein_rna_contact_lifetimes_path is not None
        else None
    )
    neighbor_rewiring = (
        load_keyvalue_timeseries(
            neighbor_rewiring_path,
            [
                "time_ns",
                "dt_ns",
                "protein_persistence",
                "rna_persistence",
                "global_persistence",
                "protein_changes",
                "rna_changes",
                "global_changes",
            ],
        )
        if neighbor_rewiring_path is not None
        else None
    )
    interface_turnover = (
        load_keyvalue_timeseries(
            interface_turnover_path,
            [
                "time_ns",
                "dt_ns",
                "new_contacts",
                "lost_contacts",
                "new_per_ns",
                "lost_per_ns",
                "net_per_ns",
            ],
        )
        if interface_turnover_path is not None
        else None
    )
    dense_msd_protein = load_dense_msd(dense_msd_protein_path) if dense_msd_protein_path is not None else None
    dense_msd_rna = load_dense_msd(dense_msd_rna_path) if dense_msd_rna_path is not None else None
    liquidity_summary = load_summary(liquidity_summary_path) if liquidity_summary_path is not None else {}
    rna_compaction_summary = load_summary(rna_compaction_summary_path) if rna_compaction_summary_path is not None else {}
    rna_compaction_binned = (
        load_csv_numeric_table(
            rna_compaction_binned_path,
            RNA_COMPACTION_BINNED_COLUMNS,
            x_name="BoundProteinChainsBin",
            x_label="Bound protein chains",
        )
        if rna_compaction_binned_path is not None
        else None
    )
    rna_energy_compaction_summary = (
        load_summary(rna_energy_compaction_summary_path)
        if rna_energy_compaction_summary_path is not None
        else {}
    )
    rna_energy_compaction_binned = (
        load_csv_numeric_table(
            rna_energy_compaction_binned_path,
            RNA_ENERGY_BINNED_COLUMNS,
            x_name="BoundProteinChainsBin",
            x_label="Bound protein chains",
        )
        if rna_energy_compaction_binned_path is not None
        else None
    )
    rna_information_summary = (
        load_summary(rna_information_summary_path) if rna_information_summary_path is not None else {}
    )
    rna_compaction_state_histogram = (
        load_csv_mixed_table(
            rna_compaction_state_histogram_path,
            RNA_COMPACTION_STATE_HISTOGRAM_COLUMNS,
            numeric_columns=[
                "StateID",
                "FrameCount",
                "Probability",
                "RgBin",
                "RgMin_nm",
                "RgMax_nm",
                "LongRangeBin",
                "LongRangeMin",
                "LongRangeMax",
            ],
            x_name="StateID",
            x_label="Compaction state ID",
        )
        if rna_compaction_state_histogram_path is not None
        else None
    )
    rna_topology_summary = (
        load_summary(rna_topology_summary_path) if rna_topology_summary_path is not None else {}
    )
    rna_topology_binned = (
        load_csv_mixed_table(
            rna_topology_binned_path,
            RNA_TOPOLOGY_BINNED_COLUMNS,
            numeric_columns=[name for name in RNA_TOPOLOGY_BINNED_COLUMNS if name != "Scope"],
            x_name="BoundProteinChainsBin",
            x_label="Bound protein chains",
        )
        if rna_topology_binned_path is not None
        else None
    )
    rna_topology_trace = (
        load_csv_numeric_table(
            rna_topology_trace_path,
            RNA_TOPOLOGY_TRACE_COLUMNS,
            x_name="Time_ns",
            x_label="Time (ns)",
        )
        if rna_topology_trace_path is not None
        else None
    )
    rna_chain_compaction_summary = (
        load_csv_numeric_table(
            rna_chain_compaction_summary_path,
            RNA_CHAIN_COMPACTION_COLUMNS,
            x_name="RNAChainID",
            x_label="RNA chain ID",
        )
        if rna_chain_compaction_summary_path is not None
        else None
    )
    rna_intrachain_contacts_summary = (
        load_csv_numeric_table(
            rna_intrachain_contacts_summary_path,
            RNA_INTRACHAIN_COLUMNS,
            x_name="RNAChainID",
            x_label="RNA chain ID",
        )
        if rna_intrachain_contacts_summary_path is not None
        else None
    )
    rna_compaction_trace = (
        load_csv_numeric_table(
            rna_compaction_trace_path,
            RNA_COMPACTION_TRACE_COLUMNS,
            x_name="Time_ns",
            x_label="Time (ns)",
        )
        if rna_compaction_trace_path is not None
        else None
    )
    rna_binding_energy_trace = (
        load_csv_numeric_table(
            rna_binding_energy_trace_path,
            RNA_BINDING_ENERGY_TRACE_COLUMNS,
            x_name="Time_ns",
            x_label="Time (ns)",
        )
        if rna_binding_energy_trace_path is not None
        else None
    )

    if llps_timeseries is not None:
        frame_weight = max(1, llps_timeseries.row_count())
    elif dbscan_timeseries is not None:
        frame_weight = max(1, dbscan_timeseries.row_count())
    else:
        frame_weight = 1

    replicate = ReplicateData(
        system=system,
        replicate_label=f"conf{conf_id}",
        source_dir=directory,
        frame_weight=frame_weight,
        llps_timeseries=llps_timeseries,
        llps_autocorr=llps_autocorr,
        llps_degree=llps_degree,
        llps_summary=summary,
        dbscan_timeseries=dbscan_timeseries,
        polymer_scaling=polymer_scaling,
        density_protein=density_protein,
        density_rna=density_rna,
        nucleic_interactions=nucleic_interactions,
        sidechain_protein=sidechain_protein,
        sidechain_rna=sidechain_rna,
        network_bundle=network_bundle,
        network_anchor_bundle=network_anchor_bundle,
        tda_points=tda_points,
        sticker_contact_lifetimes=sticker_contact_lifetimes,
        protein_rna_contact_lifetimes=protein_rna_contact_lifetimes,
        neighbor_rewiring=neighbor_rewiring,
        interface_turnover=interface_turnover,
        dense_msd_protein=dense_msd_protein,
        dense_msd_rna=dense_msd_rna,
        liquidity_summary=liquidity_summary,
        rna_compaction_summary=rna_compaction_summary,
        rna_compaction_binned=rna_compaction_binned,
        rna_energy_compaction_summary=rna_energy_compaction_summary,
        rna_energy_compaction_binned=rna_energy_compaction_binned,
        rna_information_summary=rna_information_summary,
        rna_compaction_state_histogram=rna_compaction_state_histogram,
        rna_topology_summary=rna_topology_summary,
        rna_topology_binned=rna_topology_binned,
        rna_topology_trace=rna_topology_trace,
        rna_chain_compaction_summary=rna_chain_compaction_summary,
        rna_intrachain_contacts_summary=rna_intrachain_contacts_summary,
        rna_compaction_trace=rna_compaction_trace,
        rna_binding_energy_trace=rna_binding_energy_trace,
    )

    if llps_timeseries is None:
        replicate.warnings.append("Missing llps_timeseries.dat")
    if dbscan_timeseries is None:
        replicate.warnings.append("Missing llps_dbscan_timeseries.csv")
    if not summary:
        replicate.warnings.append("Missing llps_summary.dat")
    if options.include_network and network_bundle is None:
        replicate.warnings.append("Missing network/allostery files for full --with-network mode")
    if options.include_network_anchors and network_anchor_bundle is None:
        replicate.warnings.append("Missing network GC/TE anchor files for --with-network-anchors")
    if sticker_contact_lifetimes is None:
        replicate.warnings.append("Missing llps_sticker_contact_lifetimes.dat")
    if protein_rna_contact_lifetimes is None:
        replicate.warnings.append("Missing llps_protein_rna_contact_lifetimes.dat")
    if neighbor_rewiring is None:
        replicate.warnings.append("Missing llps_neighbor_rewiring.dat")
    if interface_turnover is None:
        replicate.warnings.append("Missing llps_interface_turnover.dat")
    if dense_msd_protein is None:
        replicate.warnings.append("Missing llps_dense_msd_protein.dat")
    if dense_msd_rna is None:
        replicate.warnings.append("Missing llps_dense_msd_rna.dat")
    if not rna_compaction_summary:
        replicate.warnings.append("Missing llps_rna_compaction_binding_summary.dat")
    if rna_compaction_binned is None:
        replicate.warnings.append("Missing llps_rna_compaction_binding_binned.csv")
    if not rna_energy_compaction_summary:
        replicate.warnings.append("Missing llps_rna_energy_compaction_summary.dat")
    if rna_energy_compaction_binned is None:
        replicate.warnings.append("Missing llps_rna_energy_compaction_binned.csv")
    if not rna_information_summary:
        replicate.warnings.append("Missing llps_rna_information_summary.dat")
    if rna_compaction_state_histogram is None:
        replicate.warnings.append("Missing llps_rna_compaction_state_histogram.csv")
    if not rna_topology_summary:
        replicate.warnings.append("Missing llps_rna_topology_summary.dat")
    if rna_topology_binned is None:
        replicate.warnings.append("Missing llps_rna_topology_binding_binned.csv")
    if rna_topology_trace is None:
        replicate.warnings.append("Missing llps_rna_topology_trace.csv")
    if rna_chain_compaction_summary is None:
        replicate.warnings.append("Missing llps_rna_chain_compaction_summary.csv")
    if rna_intrachain_contacts_summary is None:
        replicate.warnings.append("Missing llps_rna_intrachain_contacts_summary.csv")
    if rna_binding_energy_trace is None:
        replicate.warnings.append("Missing llps_rna_binding_energy_trace.csv")
    return replicate


def weighted_mean(values: Sequence[float], weights: Sequence[float]) -> float:
    value_arr = np.asarray(values, dtype=float)
    weight_arr = np.asarray(weights, dtype=float)
    if value_arr.size == 0 or np.sum(weight_arr) <= 0.0:
        return float("nan")
    return float(np.sum(value_arr * weight_arr) / np.sum(weight_arr))


def combine_weighted_mean_std(
    means: Sequence[float],
    stds: Sequence[float],
    weights: Sequence[float],
) -> tuple[float, float]:
    mean_arr = np.asarray(means, dtype=float)
    std_arr = np.asarray(stds, dtype=float)
    weight_arr = np.asarray(weights, dtype=float)
    mask = np.isfinite(mean_arr) & np.isfinite(weight_arr) & (weight_arr > 0.0)
    if not np.any(mask):
        return float("nan"), float("nan")
    mean_arr = mean_arr[mask]
    std_arr = np.where(np.isfinite(std_arr[mask]), std_arr[mask], 0.0)
    weight_arr = weight_arr[mask]
    total_weight = float(np.sum(weight_arr))
    mean_value = float(np.sum(mean_arr * weight_arr) / total_weight)
    variance = float(np.sum(weight_arr * (std_arr * std_arr + (mean_arr - mean_value) ** 2)) / total_weight)
    return mean_value, math.sqrt(max(0.0, variance))


def aggregate_time_aligned(
    replicates: Sequence[ReplicateData],
    attr_name: str,
    required_columns: Sequence[str],
    max_x: float | None = None,
) -> tuple[Table | None, list[str], tuple[float, float] | None]:
    warnings: list[str] = []
    usable: list[tuple[Table, float, str]] = []
    for rep in replicates:
        table = getattr(rep, attr_name)
        if table is None:
            warnings.append(f"{rep.system}/{rep.replicate_label}: missing {attr_name}")
            continue
        if not table.has_columns(required_columns):
            warnings.append(f"{rep.system}/{rep.replicate_label}: incomplete {attr_name}")
            continue
        x = np.asarray(table.columns[required_columns[0]], dtype=float)
        if x.size == 0:
            warnings.append(f"{rep.system}/{rep.replicate_label}: empty {attr_name}")
            continue
        usable.append((table, float(rep.frame_weight), rep.replicate_label))

    if not usable:
        return None, warnings, None

    x_name = required_columns[0]
    lo = min(float(np.min(table.columns[x_name])) for table, _, _ in usable)
    hi = max(float(np.max(table.columns[x_name])) for table, _, _ in usable)
    if max_x is not None:
        hi = min(hi, float(max_x))
    if hi < lo:
        warnings.append(f"{attr_name}: invalid aggregation range after applying time cap")
        return None, warnings, None

    grid_parts: list[np.ndarray] = []
    for table, _, _ in usable:
        x = np.asarray(table.columns[x_name], dtype=float)
        mask = (x >= lo - 1e-9) & (x <= hi + 1e-9)
        if np.any(mask):
            grid_parts.append(x[mask])
    if grid_parts:
        grid = unique_sorted(np.concatenate(grid_parts))
    else:
        grid = np.asarray([], dtype=float)
    if grid.size == 0:
        if math.isclose(lo, hi, rel_tol=0.0, abs_tol=1e-9):
            grid = np.asarray([lo], dtype=float)
        else:
            grid = np.linspace(lo, hi, num=200, dtype=float)
            warnings.append(f"{attr_name}: fell back to evenly spaced interpolation grid")

    aggregated: dict[str, np.ndarray] = {x_name: grid}
    coverage_weight = np.zeros_like(grid, dtype=float)
    coverage_counts = np.zeros_like(grid, dtype=int)
    for table, weight, _ in usable:
        source_x = np.asarray(table.columns[x_name], dtype=float)
        order = np.argsort(source_x)
        unique_x = np.unique(source_x[order])
        if unique_x.size == 0:
            continue
        in_range = (grid >= unique_x[0] - 1e-9) & (grid <= unique_x[-1] + 1e-9)
        coverage_weight[in_range] += weight
        coverage_counts[in_range] += 1
    valid_mask = coverage_weight > 0.0
    if not np.all(valid_mask):
        warnings.append(
            f"{attr_name}: dropped {int(np.count_nonzero(~valid_mask))} grid point(s) without replicate coverage"
        )
        grid = grid[valid_mask]
        coverage_weight = coverage_weight[valid_mask]
        coverage_counts = coverage_counts[valid_mask]
        aggregated[x_name] = grid
    if coverage_counts.size and int(np.min(coverage_counts)) < len(usable):
        warnings.append(
            f"{attr_name}: replicate coverage decreases at late time points; later averages use the remaining trajectories only"
        )
    for column in required_columns[1:]:
        accum = np.zeros_like(grid, dtype=float)
        for table, weight, _ in usable:
            source_x = np.asarray(table.columns[x_name], dtype=float)
            source_y = np.asarray(table.columns[column], dtype=float)
            order = np.argsort(source_x)
            source_x = source_x[order]
            source_y = source_y[order]
            unique_x, unique_idx = np.unique(source_x, return_index=True)
            source_y = source_y[unique_idx]
            if unique_x.size == 0:
                continue
            in_range = (grid >= unique_x[0] - 1e-9) & (grid <= unique_x[-1] + 1e-9)
            if not np.any(in_range):
                continue
            if unique_x.size == 1:
                interp_y = np.full(int(np.count_nonzero(in_range)), source_y[0], dtype=float)
            else:
                interp_y = np.interp(grid[in_range], unique_x, source_y)
            accum[in_range] += interp_y * weight
        aggregated[column] = accum / np.maximum(coverage_weight, 1.0)

    x_label = usable[0][0].x_label
    return Table(columns=aggregated, order=list(required_columns), x_name=x_name, x_label=x_label), warnings, (lo, hi)


def aggregate_degree_distribution(replicates: Sequence[ReplicateData]) -> tuple[Table | None, list[str]]:
    warnings: list[str] = []
    usable: list[tuple[Table, float]] = []
    for rep in replicates:
        table = rep.llps_degree
        if table is None:
            warnings.append(f"{rep.system}/{rep.replicate_label}: missing llps_degree")
            continue
        usable.append((table, float(rep.frame_weight)))

    if not usable:
        return None, warnings

    k_values = unique_sorted(np.concatenate([table.columns["k"] for table, _ in usable]))
    mean_values = np.zeros_like(k_values, dtype=float)
    std_values = np.zeros_like(k_values, dtype=float)
    total_weights = np.zeros_like(k_values, dtype=float)

    for idx, kval in enumerate(k_values):
        means: list[float] = []
        stds: list[float] = []
        weights: list[float] = []
        for table, weight in usable:
            k_table = np.asarray(table.columns["k"], dtype=float)
            match = np.where(np.isclose(k_table, kval))[0]
            if match.size == 0:
                continue
            means.append(float(table.columns["pk_mean"][match[0]]))
            stds.append(float(table.columns["pk_std"][match[0]]))
            weights.append(weight)
        if not weights:
            continue
        mean_val = weighted_mean(means, weights)
        variance = 0.0
        total_w = float(np.sum(weights))
        for mean_i, std_i, weight_i in zip(means, stds, weights):
            variance += weight_i * (std_i * std_i + (mean_i - mean_val) ** 2)
        mean_values[idx] = mean_val
        std_values[idx] = math.sqrt(max(variance / max(total_w, 1.0), 0.0))
        total_weights[idx] = total_w

    table = Table(
        columns={"k": k_values, "pk_mean": mean_values, "pk_std": std_values},
        order=["k", "pk_mean", "pk_std"],
        x_name="k",
        x_label="k",
    )
    return table, warnings


def aggregate_summary(replicates: Sequence[ReplicateData]) -> tuple[dict[str, object], list[str]]:
    warnings: list[str] = []
    if not replicates:
        return {}, warnings

    weights = [float(rep.frame_weight) for rep in replicates]
    summary: dict[str, object] = {}

    percolated_votes = [bool(rep.llps_summary.get("percolated", False)) for rep in replicates]
    summary["percolated"] = sum(1 for vote in percolated_votes if vote) >= math.ceil(len(percolated_votes) / 2)

    for key in ("n_macronodes", "dimer_size"):
        values = [rep.llps_summary.get(key) for rep in replicates if key in rep.llps_summary]
        if values:
            numeric = [int(float(value)) for value in values if not math.isnan(float(value))]
            if numeric:
                summary[key] = max(set(numeric), key=numeric.count)

    for key in ("mean_gcc_rg_nm", "mean_gcc_fraction"):
        values: list[float] = []
        value_weights: list[float] = []
        for rep, weight in zip(replicates, weights):
            value = rep.llps_summary.get(key)
            if value is None:
                continue
            val = float(value)
            if math.isnan(val):
                continue
            values.append(val)
            value_weights.append(weight)
        if values:
            summary[key] = weighted_mean(values, value_weights)

    positive_tau_values: list[float] = []
    positive_tau_weights: list[float] = []
    tau_missing = 0
    for rep, weight in zip(replicates, weights):
        value = rep.llps_summary.get("tau_ns")
        if value is None:
            tau_missing += 1
            continue
        tau = float(value)
        if tau > 0:
            positive_tau_values.append(tau)
            positive_tau_weights.append(weight)
        else:
            tau_missing += 1
    if positive_tau_values:
        summary["tau_ns"] = weighted_mean(positive_tau_values, positive_tau_weights)
    else:
        summary["tau_ns"] = 0.0
        warnings.append("No positive tau_ns values found across replicates.")
    if tau_missing:
        warnings.append(f"tau_ns unavailable or non-positive in {tau_missing} replicate(s).")

    perc_values: list[float] = []
    perc_weights: list[float] = []
    for rep, weight in zip(replicates, weights):
        if not bool(rep.llps_summary.get("percolated", False)):
            continue
        value = rep.llps_summary.get("percolation_time_ns")
        if value is None:
            continue
        val = float(value)
        if math.isnan(val) or val < 0:
            continue
        perc_values.append(val)
        perc_weights.append(weight)
    if perc_values:
        summary["percolation_time_ns"] = weighted_mean(perc_values, perc_weights)
    else:
        summary["percolation_time_ns"] = -1.0
        warnings.append("No valid percolation_time_ns available across percolated replicates.")

    return summary, warnings


def aggregate_polymer_scaling(replicates: Sequence[ReplicateData]) -> tuple[Table | None, list[str]]:
    warnings: list[str] = []
    usable: list[Table] = []
    for rep in replicates:
        if rep.polymer_scaling is None:
            warnings.append(f"{rep.system}/{rep.replicate_label}: missing polymer_scaling")
            continue
        usable.append(rep.polymer_scaling)
    if not usable:
        return None, warnings

    s_values = unique_sorted(np.concatenate([table.columns["s"] for table in usable]))
    mean_r2 = np.zeros_like(s_values, dtype=float)
    rms_nm = np.zeros_like(s_values, dtype=float)
    counts = np.zeros_like(s_values, dtype=float)

    for idx, sval in enumerate(s_values):
        weighted_r2 = 0.0
        total_count = 0.0
        for table in usable:
            s_table = np.asarray(table.columns["s"], dtype=float)
            match = np.where(np.isclose(s_table, sval))[0]
            if match.size == 0:
                continue
            count = float(table.columns["count"][match[0]])
            weighted_r2 += float(table.columns["mean_R2_nm2"][match[0]]) * count
            total_count += count
        if total_count <= 0.0:
            continue
        mean_r2[idx] = weighted_r2 / total_count
        rms_nm[idx] = math.sqrt(max(mean_r2[idx], 0.0))
        counts[idx] = total_count

    table = Table(
        columns={"s": s_values, "mean_R2_nm2": mean_r2, "R_rms_nm": rms_nm, "count": counts},
        order=["s", "mean_R2_nm2", "R_rms_nm", "count"],
        x_name="s",
        x_label="Sequence separation s",
    )
    return table, warnings


def aggregate_radial_profile_tables(usable: Sequence[tuple[Table | None, float]]) -> Table | None:
    tables = [(table, weight) for table, weight in usable if table is not None and table.has_columns(["r_nm", "density_nm3"])]
    if not tables:
        return None
    r_values = unique_sorted(float(r) for table, _ in tables for r in table.columns["r_nm"])
    if r_values.size == 0:
        return None
    density_mean = np.zeros_like(r_values, dtype=float)
    voxel_counts = np.zeros_like(r_values, dtype=float)
    replicate_counts = np.zeros_like(r_values, dtype=float)
    indexed: list[tuple[dict[float, int], Table, float]] = []
    for table, weight in tables:
        lookup = {round(float(x), 12): idx for idx, x in enumerate(table.columns["r_nm"])}
        indexed.append((lookup, table, weight))
    for idx, rval in enumerate(r_values):
        key = round(float(rval), 12)
        values: list[float] = []
        weights: list[float] = []
        voxels: list[float] = []
        for lookup, table, weight in indexed:
            row_idx = lookup.get(key)
            if row_idx is None:
                continue
            density = float(table.columns["density_nm3"][row_idx])
            if not np.isfinite(density):
                continue
            values.append(density)
            weights.append(weight)
            if "voxel_count" in table.columns:
                voxels.append(float(table.columns["voxel_count"][row_idx]))
        density_mean[idx] = weighted_mean(values, weights) if values else float("nan")
        voxel_counts[idx] = float(np.sum(voxels)) if voxels else 0.0
        replicate_counts[idx] = float(len(values))
    return Table(
        columns={
            "r_nm": r_values,
            "density_nm3": density_mean,
            "voxel_count": voxel_counts,
            "n_replicates": replicate_counts,
        },
        order=["r_nm", "density_nm3", "voxel_count", "n_replicates"],
        x_name="r_nm",
        x_label="Distance from condensate center, r (nm)",
    )


def aggregate_shape_metrics(usable: Sequence[tuple[dict[str, float], float]]) -> dict[str, float]:
    keys = sorted({key for metrics, _ in usable for key in metrics})
    aggregated: dict[str, float] = {}
    for key in keys:
        values: list[float] = []
        weights: list[float] = []
        for metrics, weight in usable:
            value = metrics.get(key)
            if value is None or not np.isfinite(value):
                continue
            values.append(float(value))
            weights.append(float(weight))
        if values:
            aggregated[key] = weighted_mean(values, weights)
    if aggregated:
        aggregated["replicates_used"] = float(
            sum(1 for metrics, _ in usable if metrics and any(np.isfinite(value) for value in metrics.values()))
        )
    return aggregated


def projection_density_nm2(grid: DXGrid) -> np.ndarray:
    return np.asarray(grid.projection_xy(), dtype=np.float64) * grid.spacing_z


def projection_center_xy_nm(grid: DXGrid, projection_nm2: np.ndarray) -> tuple[float, float]:
    nx, ny = projection_nm2.shape
    x = float(grid.origin[0]) + (np.arange(nx, dtype=float) + 0.5) * grid.spacing_x
    y = float(grid.origin[1]) + (np.arange(ny, dtype=float) + 0.5) * grid.spacing_y
    positive = np.where(np.isfinite(projection_nm2) & (projection_nm2 > 0.0), projection_nm2, 0.0)
    total = float(np.sum(positive, dtype=np.float64))
    if total <= 0.0:
        return float(np.mean(x)), float(np.mean(y))
    x_weight = np.sum(positive, axis=1, dtype=np.float64)
    y_weight = np.sum(positive, axis=0, dtype=np.float64)
    return float(np.dot(x, x_weight) / total), float(np.dot(y, y_weight) / total)


def resample_projection_2d(
    projection_nm2: np.ndarray,
    x_src: np.ndarray,
    y_src: np.ndarray,
    x_target: np.ndarray,
    y_target: np.ndarray,
) -> np.ndarray:
    if projection_nm2.size == 0 or x_src.size < 2 or y_src.size < 2 or x_target.size == 0 or y_target.size == 0:
        return np.full((x_target.size, y_target.size), np.nan, dtype=np.float64)
    tmp = np.empty((x_target.size, y_src.size), dtype=np.float64)
    for j in range(y_src.size):
        tmp[:, j] = np.interp(x_target, x_src, projection_nm2[:, j], left=np.nan, right=np.nan)
    out = np.empty((x_target.size, y_target.size), dtype=np.float64)
    for i in range(x_target.size):
        row = tmp[i, :]
        if np.count_nonzero(np.isfinite(row)) < 2:
            out[i, :] = np.nan
            continue
        out[i, :] = np.interp(y_target, y_src, row, left=np.nan, right=np.nan)
    return out


def aggregate_centered_projection_grids(
    usable: Sequence[tuple[DXGrid, float, str]],
    attr_name: str,
    warnings: list[str],
) -> DXGrid | None:
    prepared: list[tuple[np.ndarray, np.ndarray, np.ndarray, float]] = []
    spacing_x: list[float] = []
    spacing_y: list[float] = []
    for grid, weight, _ in usable:
        projection = projection_density_nm2(grid)
        if projection.size == 0:
            continue
        center_x, center_y = projection_center_xy_nm(grid, projection)
        nx, ny = projection.shape
        x_src = float(grid.origin[0]) + (np.arange(nx, dtype=float) + 0.5) * grid.spacing_x - center_x
        y_src = float(grid.origin[1]) + (np.arange(ny, dtype=float) + 0.5) * grid.spacing_y - center_y
        if x_src.size < 2 or y_src.size < 2:
            continue
        prepared.append((projection, x_src, y_src, float(weight)))
        spacing_x.append(grid.spacing_x)
        spacing_y.append(grid.spacing_y)

    if not prepared:
        return None

    dx = float(np.median(spacing_x)) if spacing_x else 1.0
    dy = float(np.median(spacing_y)) if spacing_y else 1.0
    if dx <= 0.0 or dy <= 0.0:
        warnings.append(f"{attr_name}: invalid DX spacing; density projection could not be center-aligned")
        return None

    x_min = min(float(x_src[0]) for _, x_src, _, _ in prepared)
    x_max = max(float(x_src[-1]) for _, x_src, _, _ in prepared)
    y_min = min(float(y_src[0]) for _, _, y_src, _ in prepared)
    y_max = max(float(y_src[-1]) for _, _, y_src, _ in prepared)
    if x_max <= x_min or y_max <= y_min:
        warnings.append(f"{attr_name}: invalid centered projection extent")
        return None

    x_target = np.arange(x_min, x_max + 0.5 * dx, dx, dtype=float)
    y_target = np.arange(y_min, y_max + 0.5 * dy, dy, dtype=float)
    accum = np.zeros((x_target.size, y_target.size), dtype=np.float64)
    weight_accum = np.zeros_like(accum)
    coverage_counts = np.zeros_like(accum, dtype=np.int16)

    for projection, x_src, y_src, weight in prepared:
        resampled = resample_projection_2d(projection, x_src, y_src, x_target, y_target)
        valid = np.isfinite(resampled)
        if not np.any(valid):
            continue
        accum[valid] += resampled[valid] * weight
        weight_accum[valid] += weight
        coverage_counts[valid] += 1

    valid = weight_accum > 0.0
    if not np.any(valid):
        warnings.append(f"{attr_name}: centered projection has no overlapping finite pixels")
        return None
    if int(np.min(coverage_counts[valid])) < len(prepared):
        warnings.append(
            f"{attr_name}: center-aligned density projection edge pixels use partial replicate coverage"
        )

    averaged = np.zeros_like(accum, dtype=np.float64)
    averaged[valid] = accum[valid] / weight_accum[valid]
    return DXGrid(
        counts=(int(x_target.size), int(y_target.size), 1),
        origin=(float(x_target[0] - 0.5 * dx), float(y_target[0] - 0.5 * dy), 0.0),
        deltas=((dx, 0.0, 0.0), (0.0, dy, 0.0), (0.0, 0.0, 1.0)),
        label=f"{attr_name}_centered_projection",
        data=None,
        projection_xy_data=averaged.astype(np.float32),
    )


def aggregate_density_grids(replicates: Sequence[ReplicateData], attr_name: str) -> tuple[DXGrid | None, list[str]]:
    warnings: list[str] = []
    shape_usable: list[tuple[DXGrid, float, str]] = []
    reference: DXGrid | None = None
    for rep in replicates:
        grid = getattr(rep, attr_name)
        if grid is None:
            warnings.append(f"{rep.system}/{rep.replicate_label}: missing {attr_name}")
            continue
        weight = float(rep.frame_weight)
        shape_usable.append((grid, weight, rep.replicate_label))
        if reference is None:
            reference = grid
        elif not reference.compatible_with(grid):
            warnings.append(
                f"{rep.system}/{rep.replicate_label}: incompatible {attr_name} grid; 2D projection is center-aligned before averaging, and full 3D DX averaging skips this replicate"
            )

    if not shape_usable:
        return None, warnings

    radial_profile = aggregate_radial_profile_tables([(grid.radial_profile_data, weight) for grid, weight, _ in shape_usable])
    shape_metrics = aggregate_shape_metrics([(grid.shape_metrics, weight) for grid, weight, _ in shape_usable])
    aggregated = aggregate_centered_projection_grids(shape_usable, attr_name, warnings)
    if aggregated is None:
        ref = shape_usable[0][0]
        warnings.append(f"{attr_name}: falling back to the first available density projection")
        aggregated = DXGrid(
            counts=ref.counts,
            origin=ref.origin,
            deltas=ref.deltas,
            label=ref.label,
            data=None,
            projection_xy_data=ref.projection_xy().astype(np.float32, copy=False),
        )
    aggregated.radial_profile_data = radial_profile
    aggregated.shape_metrics = shape_metrics
    return aggregated, warnings


def classify_recognition_type(base_occupancy: float, backbone_occupancy: float) -> str:
    if base_occupancy > backbone_occupancy * 1.5:
        return "Highly_Specific"
    if backbone_occupancy > base_occupancy * 1.5:
        return "Non_Specific_Backbone"
    return "Mixed_Binding"


def aggregate_nucleic_interactions(replicates: Sequence[ReplicateData]) -> tuple[list[NucleicInteraction], list[str]]:
    warnings: list[str] = []
    merged: dict[tuple[int, str], dict[str, float]] = defaultdict(lambda: {"base_frames": 0.0, "backbone_frames": 0.0})
    total_frames = float(sum(max(1, rep.frame_weight) for rep in replicates))
    any_rows = False
    for rep in replicates:
        if not rep.nucleic_interactions:
            warnings.append(f"{rep.system}/{rep.replicate_label}: empty or missing nucleic_acid_interactions")
            continue
        any_rows = True
        for entry in rep.nucleic_interactions:
            key = (entry.resnum, entry.resname)
            merged[key]["base_frames"] += entry.base_frames
            merged[key]["backbone_frames"] += entry.backbone_frames

    if not any_rows:
        return [], warnings

    aggregated: list[NucleicInteraction] = []
    for (resnum, resname), values in sorted(merged.items(), key=lambda item: item[0][0]):
        base_occ = (values["base_frames"] / max(total_frames, 1.0)) * 100.0
        back_occ = (values["backbone_frames"] / max(total_frames, 1.0)) * 100.0
        if base_occ <= 5.0 and back_occ <= 5.0:
            continue
        aggregated.append(
            NucleicInteraction(
                resnum=resnum,
                resname=resname,
                base_frames=values["base_frames"],
                backbone_frames=values["backbone_frames"],
                base_occupancy=base_occ,
                backbone_occupancy=back_occ,
                recognition_type=classify_recognition_type(base_occ, back_occ),
            )
        )
    aggregated.sort(key=lambda item: item.total_occupancy, reverse=True)
    return aggregated, warnings


def aggregate_sidechain_contacts(
    replicates: Sequence[ReplicateData], attr_name: str
) -> tuple[list[SidechainContact], list[str]]:
    warnings: list[str] = []
    merged: dict[tuple[str, int, str], float] = defaultdict(float)
    for rep in replicates:
        contacts = getattr(rep, attr_name)
        if not contacts:
            warnings.append(f"{rep.system}/{rep.replicate_label}: empty or missing {attr_name}")
            continue
        for entry in contacts:
            merged[(entry.chain, entry.resnum, entry.resname)] += entry.score
    aggregated = [
        SidechainContact(chain=chain, resnum=resnum, resname=resname, score=score)
        for (chain, resnum, resname), score in merged.items()
    ]
    aggregated.sort(key=lambda item: item.score, reverse=True)
    return aggregated, warnings


def aggregate_network_bundle(replicates: Sequence[ReplicateData]) -> tuple[NetworkBundle | None, list[str]]:
    warnings: list[str] = []
    bundles = [rep.network_bundle for rep in replicates if rep.network_bundle is not None]
    if not bundles:
        warnings.append("No network/allostery files found across replicates.")
        return None, warnings

    aggregated = NetworkBundle()
    node_union: dict[str, dict[str, str]] = {}
    weights = {rep.replicate_label: float(rep.frame_weight) for rep in replicates}

    for rep in replicates:
        bundle = rep.network_bundle
        if bundle is None:
            warnings.append(f"{rep.system}/{rep.replicate_label}: missing network bundle")
            continue
        for node_id, row in bundle.nodes.items():
            node_union.setdefault(node_id, row)

    aggregated.nodes = node_union

    def merge_edge_dict(selector: str) -> dict[tuple[str, str], float]:
        merged: dict[tuple[str, str], float] = {}
        key_union: set[tuple[str, str]] = set()
        for rep in replicates:
            bundle = rep.network_bundle
            if bundle is None:
                continue
            key_union.update(getattr(bundle, selector).keys())
        for key in key_union:
            total = 0.0
            total_weight = 0.0
            for rep in replicates:
                bundle = rep.network_bundle
                if bundle is None:
                    continue
                value = getattr(bundle, selector).get(key, 0.0)
                weight = weights[rep.replicate_label]
                total += value * weight
                total_weight += weight
            merged[key] = total / max(total_weight, 1.0)
        return merged

    aggregated.corr_edges = merge_edge_dict("corr_edges")
    aggregated.gc_edges = merge_edge_dict("gc_edges")
    aggregated.te_edges = merge_edge_dict("te_edges")
    aggregated.mi_edges = merge_edge_dict("mi_edges")
    return aggregated, warnings


def aggregate_network_anchor_bundle(
    replicates: Sequence[ReplicateData],
) -> tuple[NetworkAnchorBundle | None, list[str]]:
    warnings: list[str] = []
    bundles = [rep.network_anchor_bundle for rep in replicates if rep.network_anchor_bundle is not None]
    if not bundles:
        warnings.append("No low-memory network anchor files found across replicates.")
        return None, warnings

    aggregated = NetworkAnchorBundle()
    for rep in replicates:
        bundle = rep.network_anchor_bundle
        if bundle is None:
            warnings.append(f"{rep.system}/{rep.replicate_label}: missing network anchor bundle")
            continue
        for node_id, label in bundle.labels.items():
            aggregated.labels.setdefault(node_id, label)

    weights = {rep.replicate_label: float(max(1, rep.frame_weight)) for rep in replicates}

    def merge_score_dict(selector: str) -> dict[str, float]:
        key_union: set[str] = set()
        for rep in replicates:
            bundle = rep.network_anchor_bundle
            if bundle is not None:
                key_union.update(getattr(bundle, selector).keys())
        merged: dict[str, float] = {}
        for key in key_union:
            total = 0.0
            total_weight = 0.0
            for rep in replicates:
                bundle = rep.network_anchor_bundle
                if bundle is None:
                    continue
                total += float(getattr(bundle, selector).get(key, 0.0)) * weights[rep.replicate_label]
                total_weight += weights[rep.replicate_label]
            value = total / max(total_weight, 1.0)
            if np.isfinite(value) and value > 0.0:
                merged[key] = value
        return merged

    aggregated.gc_anchor_scores = merge_score_dict("gc_anchor_scores")
    aggregated.te_source_scores = merge_score_dict("te_source_scores")
    aggregated.te_sink_scores = merge_score_dict("te_sink_scores")
    return aggregated if aggregated.available() or aggregated.labels else None, warnings


def aggregate_tda_points(replicates: Sequence[ReplicateData]) -> tuple[np.ndarray | None, list[str]]:
    warnings: list[str] = []
    pooled: list[np.ndarray] = []
    for rep in replicates:
        if rep.tda_points is None:
            warnings.append(f"{rep.system}/{rep.replicate_label}: missing tda_global_pointcloud.csv")
            continue
        pooled.append(rep.tda_points)
    if not pooled:
        return None, warnings
    return np.concatenate(pooled, axis=0).astype(np.float32, copy=False), warnings


def weighted_quantile_from_hist(values: np.ndarray, counts: np.ndarray, quantile: float) -> float:
    if values.size == 0 or counts.size == 0 or float(np.sum(counts)) <= 0.0:
        return float("nan")
    order = np.argsort(values)
    sorted_values = values[order]
    sorted_counts = counts[order]
    cutoff = float(np.sum(sorted_counts)) * float(np.clip(quantile, 0.0, 1.0))
    cumulative = np.cumsum(sorted_counts)
    idx = int(np.searchsorted(cumulative, cutoff, side="left"))
    idx = min(max(idx, 0), sorted_values.size - 1)
    return float(sorted_values[idx])


def aggregate_lifetime_bundle(
    replicates: Sequence[ReplicateData], attr_name: str
) -> tuple[DynamicsBundle | None, list[str]]:
    warnings: list[str] = []
    bins: dict[float, float] = defaultdict(float)
    usable = 0
    for rep in replicates:
        bundle = getattr(rep, attr_name)
        if bundle is None or bundle.table is None:
            warnings.append(f"{rep.system}/{rep.replicate_label}: missing {attr_name}")
            continue
        usable += 1
        lifetimes = bundle.table.columns.get("lifetime_ns", np.asarray([], dtype=float))
        counts = bundle.table.columns.get("count", np.asarray([], dtype=float))
        for lifetime, count in zip(lifetimes, counts):
            if np.isfinite(lifetime) and np.isfinite(count) and count > 0:
                bins[round(float(lifetime), 6)] += float(count)
    if not bins:
        return None, warnings
    values = np.asarray(sorted(bins), dtype=float)
    counts = np.asarray([bins[round(float(value), 6)] for value in values], dtype=float)
    total = float(np.sum(counts))
    mean = float(np.sum(values * counts) / total) if total > 0.0 else float("nan")
    summary = {
        "total_events": total,
        "mean_ns": mean,
        "median_ns": weighted_quantile_from_hist(values, counts, 0.5),
        "p90_ns": weighted_quantile_from_hist(values, counts, 0.9),
        "max_ns": float(np.max(values)) if values.size else float("nan"),
        "replicates_used": float(usable),
    }
    table = Table(
        columns={"lifetime_ns": values, "count": counts},
        order=["lifetime_ns", "count"],
        x_name="lifetime_ns",
        x_label="Lifetime (ns)",
    )
    return DynamicsBundle(summary=summary, table=table), warnings


def aggregate_dynamics_bundle(
    replicates: Sequence[ReplicateData],
    attr_name: str,
    required_columns: Sequence[str],
    max_x: float | None = None,
) -> tuple[DynamicsBundle | None, list[str], tuple[float, float] | None]:
    warnings: list[str] = []
    usable: list[tuple[DynamicsBundle, float, str]] = []
    for rep in replicates:
        bundle = getattr(rep, attr_name)
        if bundle is None:
            warnings.append(f"{rep.system}/{rep.replicate_label}: missing {attr_name}")
            continue
        if bundle.table is not None and bundle.table.has_columns(required_columns):
            usable.append((bundle, float(max(1, rep.frame_weight)), rep.replicate_label))
        elif bundle.summary:
            usable.append((bundle, float(max(1, rep.frame_weight)), rep.replicate_label))
        else:
            warnings.append(f"{rep.system}/{rep.replicate_label}: empty {attr_name}")
    if not usable:
        return None, warnings, None

    summary: dict[str, float] = {}
    for key in sorted({key for bundle, _, _ in usable for key in bundle.summary}):
        values: list[float] = []
        weights: list[float] = []
        for bundle, weight, _ in usable:
            value = bundle.summary.get(key)
            if value is not None and np.isfinite(value):
                values.append(float(value))
                weights.append(weight)
        if values:
            summary[key] = weighted_mean(values, weights)

    table_bundles = [(bundle, weight, label) for bundle, weight, label in usable if bundle.table is not None]
    if not table_bundles:
        return DynamicsBundle(summary=summary, table=None), warnings, None

    x_col_name = required_columns[0]
    x_values = unique_sorted(
        x
        for bundle, _, _ in table_bundles
        for x in bundle.table.columns[x_col_name]
        if max_x is None or x <= max_x + 1e-9
    )
    if x_values.size == 0:
        return DynamicsBundle(summary=summary, table=None), warnings, None

    indexed: list[tuple[dict[float, int], DynamicsBundle, float]] = []
    for bundle, weight, _ in table_bundles:
        lookup = {round(float(x), 12): idx for idx, x in enumerate(bundle.table.columns[x_col_name])}
        indexed.append((lookup, bundle, weight))

    out_columns: dict[str, list[float]] = {name: [] for name in required_columns}
    replicate_counts: list[float] = []
    for x in x_values:
        key = round(float(x), 12)
        out_columns[x_col_name].append(float(x))
        active_count = 0
        for name in required_columns[1:]:
            values: list[float] = []
            weights: list[float] = []
            for lookup, bundle, weight in indexed:
                idx = lookup.get(key)
                if idx is None:
                    continue
                value = float(bundle.table.columns[name][idx])
                if np.isfinite(value):
                    values.append(value)
                    weights.append(weight)
            if values:
                active_count = max(active_count, len(values))
            out_columns[name].append(weighted_mean(values, weights))
        replicate_counts.append(float(active_count))
    out_columns["n_replicates"] = replicate_counts
    order = list(required_columns) + ["n_replicates"]
    table = Table(
        columns={name: np.asarray(values, dtype=float) for name, values in out_columns.items()},
        order=order,
        x_name=x_col_name,
        x_label="Time (ns)",
    )
    return DynamicsBundle(summary=summary, table=table), warnings, (float(np.min(x_values)), float(np.max(x_values)))


def fit_msd_diffusion(table: Table | None, fit_start_ns: float | None = None, fit_end_ns: float | None = None) -> dict[str, float]:
    if table is None or not table.has_columns(["lag_ns", "msd_nm2", "count"]):
        return {}
    lag = np.asarray(table.columns["lag_ns"], dtype=float)
    msd = np.asarray(table.columns["msd_nm2"], dtype=float)
    count = np.asarray(table.columns["count"], dtype=float)
    mask = np.isfinite(lag) & np.isfinite(msd) & np.isfinite(count) & (lag > 0.0) & (count > 0.0)
    if np.count_nonzero(mask) < 2:
        return {}
    lag = lag[mask]
    msd = msd[mask]
    count = count[mask]
    if fit_start_ns is None or fit_end_ns is None or fit_end_ns <= fit_start_ns:
        unique_lags = np.unique(np.round(lag, 12))
        positive_steps = np.diff(unique_lags)
        positive_steps = positive_steps[positive_steps > 0.0]
        dt = float(np.min(positive_steps)) if positive_steps.size else float(np.min(lag))
        fit_start_ns = max(float(np.min(lag)), 5.0 * dt)
        fit_end_ns = 0.25 * float(np.max(lag))
        if fit_end_ns <= fit_start_ns:
            fit_start_ns = float(np.min(lag))
            fit_end_ns = float(np.max(lag))
    fit_mask = (lag >= fit_start_ns - 1e-9) & (lag <= fit_end_ns + 1e-9)
    if np.count_nonzero(fit_mask) < 2:
        fit_mask = np.ones_like(lag, dtype=bool)
    try:
        coeffs = np.polyfit(lag[fit_mask], msd[fit_mask], deg=1, w=np.sqrt(np.maximum(count[fit_mask], 1.0)))
    except (ValueError, np.linalg.LinAlgError):
        return {}
    slope = float(coeffs[0])
    intercept = float(coeffs[1])
    return {
        "diffusion_nm2_per_ns": slope / 6.0,
        "fit_start_ns": float(fit_start_ns),
        "fit_end_ns": float(fit_end_ns),
        "fit_slope_nm2_per_ns": slope,
        "fit_intercept_nm2": intercept,
    }


def aggregate_dense_msd_bundle(
    replicates: Sequence[ReplicateData], attr_name: str, max_x: float | None = None
) -> tuple[DynamicsBundle | None, list[str], tuple[float, float] | None]:
    warnings: list[str] = []
    usable: list[tuple[DynamicsBundle, str]] = []
    for rep in replicates:
        bundle = getattr(rep, attr_name)
        if bundle is None or bundle.table is None or not bundle.table.has_columns(["lag_ns", "msd_nm2", "count"]):
            warnings.append(f"{rep.system}/{rep.replicate_label}: missing {attr_name}")
            continue
        usable.append((bundle, rep.replicate_label))
    if not usable:
        return None, warnings, None

    x_values = unique_sorted(
        lag
        for bundle, _ in usable
        for lag in bundle.table.columns["lag_ns"]
        if max_x is None or lag <= max_x + 1e-9
    )
    if x_values.size == 0:
        return None, warnings, None

    indexed: list[tuple[dict[float, int], DynamicsBundle]] = []
    for bundle, _ in usable:
        lookup = {round(float(x), 12): idx for idx, x in enumerate(bundle.table.columns["lag_ns"])}
        indexed.append((lookup, bundle))

    lag_out: list[float] = []
    msd_out: list[float] = []
    count_out: list[float] = []
    nrep_out: list[float] = []
    for lag_value in x_values:
        key = round(float(lag_value), 12)
        weighted_sum = 0.0
        total_count = 0.0
        active = 0
        for lookup, bundle in indexed:
            idx = lookup.get(key)
            if idx is None:
                continue
            msd_value = float(bundle.table.columns["msd_nm2"][idx])
            count_value = float(bundle.table.columns["count"][idx])
            if not (np.isfinite(msd_value) and np.isfinite(count_value) and count_value > 0.0):
                continue
            weighted_sum += msd_value * count_value
            total_count += count_value
            active += 1
        if total_count <= 0.0:
            continue
        lag_out.append(float(lag_value))
        msd_out.append(weighted_sum / total_count)
        count_out.append(total_count)
        nrep_out.append(float(active))
    if not lag_out:
        return None, warnings, None

    table = Table(
        columns={
            "lag_ns": np.asarray(lag_out, dtype=float),
            "msd_nm2": np.asarray(msd_out, dtype=float),
            "count": np.asarray(count_out, dtype=float),
            "n_replicates": np.asarray(nrep_out, dtype=float),
        },
        order=["lag_ns", "msd_nm2", "count", "n_replicates"],
        x_name="lag_ns",
        x_label="Lag time (ns)",
    )
    summary = fit_msd_diffusion(table)
    if not summary:
        diffusion_values: list[float] = []
        weights: list[float] = []
        for bundle, _ in usable:
            diffusion = bundle.summary.get("diffusion_nm2_per_ns")
            if diffusion is not None and np.isfinite(diffusion):
                total_count = 1.0
                if bundle.table is not None and "count" in bundle.table.columns:
                    total_count = max(1.0, float(np.nansum(bundle.table.columns["count"])))
                diffusion_values.append(float(diffusion))
                weights.append(total_count)
        if diffusion_values:
            summary["diffusion_nm2_per_ns"] = weighted_mean(diffusion_values, weights)
    summary["replicates_used"] = float(len(usable))
    return DynamicsBundle(summary=summary, table=table), warnings, (float(min(lag_out)), float(max(lag_out)))


def aggregate_rna_compaction_binned(replicates: Sequence[ReplicateData]) -> tuple[Table | None, list[str]]:
    warnings: list[str] = []
    usable: list[tuple[Table, str]] = []
    for rep in replicates:
        table = rep.rna_compaction_binned
        if table is None or not table.has_columns(RNA_COMPACTION_BINNED_COLUMNS):
            warnings.append(f"{rep.system}/{rep.replicate_label}: missing rna_compaction_binned")
            continue
        usable.append((table, rep.replicate_label))
    if not usable:
        return None, warnings

    bins = unique_sorted(
        bin_value
        for table, _ in usable
        for bin_value in table.columns["BoundProteinChainsBin"]
    )
    columns: dict[str, np.ndarray] = {
        "BoundProteinChainsBin": bins.astype(float),
        "FrameCount": np.zeros_like(bins, dtype=float),
        "NReplicates": np.zeros_like(bins, dtype=float),
    }
    for name in RNA_COMPACTION_BINNED_COLUMNS:
        if name not in columns:
            columns[name] = np.full_like(bins, np.nan, dtype=float)

    metric_pairs = [
        ("RgMean_nm", "RgStd_nm"),
        ("EndToEndMean_nm", "EndToEndStd_nm"),
        ("AsphericityMean", "AsphericityStd"),
        ("LongRangeContactFractionMean", "LongRangeContactFractionStd"),
    ]
    for idx, bin_value in enumerate(bins):
        weights: list[float] = []
        rows: list[tuple[Table, int]] = []
        for table, _ in usable:
            source_bins = np.asarray(table.columns["BoundProteinChainsBin"], dtype=float)
            match = np.where(np.isclose(source_bins, bin_value))[0]
            if match.size == 0:
                continue
            row_idx = int(match[0])
            frame_count = float(table.columns["FrameCount"][row_idx])
            if np.isfinite(frame_count) and frame_count > 0.0:
                weights.append(frame_count)
                rows.append((table, row_idx))
        if not rows:
            continue
        columns["FrameCount"][idx] = float(np.sum(weights))
        columns["NReplicates"][idx] = float(len(rows))
        for mean_name, std_name in metric_pairs:
            means = [float(table.columns[mean_name][row_idx]) for table, row_idx in rows]
            stds = [float(table.columns[std_name][row_idx]) for table, row_idx in rows]
            mean_value, std_value = combine_weighted_mean_std(means, stds, weights)
            columns[mean_name][idx] = mean_value
            columns[std_name][idx] = std_value

    order = RNA_COMPACTION_BINNED_COLUMNS + ["NReplicates"]
    return Table(columns=columns, order=order, x_name="BoundProteinChainsBin", x_label="Bound protein chains"), warnings


def aggregate_rna_energy_compaction_binned(
    replicates: Sequence[ReplicateData],
) -> tuple[Table | None, list[str]]:
    warnings: list[str] = []
    usable: list[tuple[Table, str]] = []
    for rep in replicates:
        table = rep.rna_energy_compaction_binned
        if table is None or not table.has_columns(RNA_ENERGY_BINNED_COLUMNS):
            warnings.append(f"{rep.system}/{rep.replicate_label}: missing rna_energy_compaction_binned")
            continue
        usable.append((table, rep.replicate_label))
    if not usable:
        return None, warnings

    bins = unique_sorted(
        bin_value
        for table, _ in usable
        for bin_value in table.columns["BoundProteinChainsBin"]
    )
    columns: dict[str, np.ndarray] = {
        "BoundProteinChainsBin": bins.astype(float),
        "FrameCount": np.zeros_like(bins, dtype=float),
        "NReplicates": np.zeros_like(bins, dtype=float),
    }
    for name in RNA_ENERGY_BINNED_COLUMNS:
        if name not in columns:
            columns[name] = np.full_like(bins, np.nan, dtype=float)

    metric_pairs = [
        ("EelecMean_kJmol", "EelecStd_kJmol"),
        ("EvdwMean_kJmol", "EvdwStd_kJmol"),
        ("EtotalMean_kJmol", "EtotalStd_kJmol"),
        ("AttractionStrengthMean_kJmol", "AttractionStrengthStd_kJmol"),
        ("RgMean_nm", "RgStd_nm"),
        ("EndToEndMean_nm", "EndToEndStd_nm"),
        ("RelativeShapeAnisotropyMean", "RelativeShapeAnisotropyStd"),
        ("LongRangeContactFractionMean", "LongRangeContactFractionStd"),
    ]
    for idx, bin_value in enumerate(bins):
        weights: list[float] = []
        rows: list[tuple[Table, int]] = []
        for table, _ in usable:
            source_bins = np.asarray(table.columns["BoundProteinChainsBin"], dtype=float)
            match = np.where(np.isclose(source_bins, bin_value))[0]
            if match.size == 0:
                continue
            row_idx = int(match[0])
            frame_count = float(table.columns["FrameCount"][row_idx])
            if np.isfinite(frame_count) and frame_count > 0.0:
                weights.append(frame_count)
                rows.append((table, row_idx))
        if not rows:
            continue
        columns["FrameCount"][idx] = float(np.sum(weights))
        columns["NReplicates"][idx] = float(len(rows))
        for mean_name, std_name in metric_pairs:
            means = [float(table.columns[mean_name][row_idx]) for table, row_idx in rows]
            stds = [float(table.columns[std_name][row_idx]) for table, row_idx in rows]
            mean_value, std_value = combine_weighted_mean_std(means, stds, weights)
            columns[mean_name][idx] = mean_value
            columns[std_name][idx] = std_value

    order = RNA_ENERGY_BINNED_COLUMNS + ["NReplicates"]
    return Table(columns=columns, order=order, x_name="BoundProteinChainsBin", x_label="Bound protein chains"), warnings


def aggregate_rna_chain_compaction_table(
    replicates: Sequence[ReplicateData],
) -> tuple[Table | None, list[str]]:
    warnings: list[str] = []
    usable: list[tuple[Table, str]] = []
    for rep in replicates:
        table = rep.rna_chain_compaction_summary
        if table is None or not table.has_columns(RNA_CHAIN_COMPACTION_COLUMNS):
            warnings.append(f"{rep.system}/{rep.replicate_label}: missing rna_chain_compaction_summary")
            continue
        usable.append((table, rep.replicate_label))
    if not usable:
        return None, warnings

    chain_ids = unique_sorted(chain for table, _ in usable for chain in table.columns["RNAChainID"])
    columns: dict[str, np.ndarray] = {
        "RNAChainID": chain_ids.astype(float),
        "FrameCount": np.zeros_like(chain_ids, dtype=float),
        "NReplicates": np.zeros_like(chain_ids, dtype=float),
    }
    for name in RNA_CHAIN_COMPACTION_COLUMNS:
        if name not in columns:
            columns[name] = np.full_like(chain_ids, np.nan, dtype=float)

    mean_columns = [
        "InsideDenseFraction",
        "MeanBoundProteinChains",
        "MeanProteinRNABeadContacts",
        "MeanRelativeShapeAnisotropy",
        "MeanLongRangeContactFraction",
    ]
    metric_pairs = [("MeanRg_nm", "StdRg_nm"), ("MeanEndToEnd_nm", "StdEndToEnd_nm")]

    for idx, chain_id in enumerate(chain_ids):
        weights: list[float] = []
        rows: list[tuple[Table, int]] = []
        for table, _ in usable:
            source_ids = np.asarray(table.columns["RNAChainID"], dtype=float)
            match = np.where(np.isclose(source_ids, chain_id))[0]
            if match.size == 0:
                continue
            row_idx = int(match[0])
            frame_count = float(table.columns["FrameCount"][row_idx])
            if np.isfinite(frame_count) and frame_count > 0.0:
                weights.append(frame_count)
                rows.append((table, row_idx))
        if not rows:
            continue
        columns["FrameCount"][idx] = float(np.sum(weights))
        columns["NReplicates"][idx] = float(len(rows))
        for name in mean_columns:
            values = [float(table.columns[name][row_idx]) for table, row_idx in rows]
            columns[name][idx] = weighted_mean(values, weights)
        for mean_name, std_name in metric_pairs:
            means = [float(table.columns[mean_name][row_idx]) for table, row_idx in rows]
            stds = [float(table.columns[std_name][row_idx]) for table, row_idx in rows]
            mean_value, std_value = combine_weighted_mean_std(means, stds, weights)
            columns[mean_name][idx] = mean_value
            columns[std_name][idx] = std_value

    order = RNA_CHAIN_COMPACTION_COLUMNS + ["NReplicates"]
    return Table(columns=columns, order=order, x_name="RNAChainID", x_label="RNA chain ID"), warnings


def aggregate_rna_intrachain_contacts_table(
    replicates: Sequence[ReplicateData],
) -> tuple[Table | None, list[str]]:
    warnings: list[str] = []
    usable: list[tuple[Table, str]] = []
    for rep in replicates:
        table = rep.rna_intrachain_contacts_summary
        if table is None or not table.has_columns(RNA_INTRACHAIN_COLUMNS):
            warnings.append(f"{rep.system}/{rep.replicate_label}: missing rna_intrachain_contacts_summary")
            continue
        usable.append((table, rep.replicate_label))
    if not usable:
        return None, warnings

    chain_ids = unique_sorted(chain for table, _ in usable for chain in table.columns["RNAChainID"])
    columns: dict[str, np.ndarray] = {
        "RNAChainID": chain_ids.astype(float),
        "FrameCount": np.zeros_like(chain_ids, dtype=float),
        "NReplicates": np.zeros_like(chain_ids, dtype=float),
    }
    for name in RNA_INTRACHAIN_COLUMNS:
        if name not in columns:
            columns[name] = np.full_like(chain_ids, np.nan, dtype=float)

    metric_pairs = [
        ("MeanIntrachainContacts", "StdIntrachainContacts"),
        ("MeanLongRangeIntrachainContacts", "StdLongRangeIntrachainContacts"),
        ("MeanLongRangeContactFraction", "StdLongRangeContactFraction"),
    ]
    for idx, chain_id in enumerate(chain_ids):
        weights: list[float] = []
        rows: list[tuple[Table, int]] = []
        for table, _ in usable:
            source_ids = np.asarray(table.columns["RNAChainID"], dtype=float)
            match = np.where(np.isclose(source_ids, chain_id))[0]
            if match.size == 0:
                continue
            row_idx = int(match[0])
            frame_count = float(table.columns["FrameCount"][row_idx])
            if np.isfinite(frame_count) and frame_count > 0.0:
                weights.append(frame_count)
                rows.append((table, row_idx))
        if not rows:
            continue
        columns["FrameCount"][idx] = float(np.sum(weights))
        columns["NReplicates"][idx] = float(len(rows))
        for mean_name, std_name in metric_pairs:
            means = [float(table.columns[mean_name][row_idx]) for table, row_idx in rows]
            stds = [float(table.columns[std_name][row_idx]) for table, row_idx in rows]
            mean_value, std_value = combine_weighted_mean_std(means, stds, weights)
            columns[mean_name][idx] = mean_value
            columns[std_name][idx] = std_value

    order = RNA_INTRACHAIN_COLUMNS + ["NReplicates"]
    return Table(columns=columns, order=order, x_name="RNAChainID", x_label="RNA chain ID"), warnings


def aggregate_rna_compaction_state_histogram(
    replicates: Sequence[ReplicateData],
) -> tuple[Table | None, list[str]]:
    warnings: list[str] = []
    merged: dict[tuple[str, int], dict[str, object]] = {}
    scope_totals: dict[str, float] = defaultdict(float)
    usable = 0
    for rep in replicates:
        table = rep.rna_compaction_state_histogram
        if table is None or not table.has_columns(RNA_COMPACTION_STATE_HISTOGRAM_COLUMNS):
            warnings.append(f"{rep.system}/{rep.replicate_label}: missing rna_compaction_state_histogram")
            continue
        usable += 1
        scopes = np.asarray(table.columns["Scope"], dtype=object)
        state_ids = np.asarray(table.columns["StateID"], dtype=float)
        frame_counts = np.asarray(table.columns["FrameCount"], dtype=float)
        rg_bins = np.asarray(table.columns["RgBin"], dtype=float)
        rg_min = np.asarray(table.columns["RgMin_nm"], dtype=float)
        rg_max = np.asarray(table.columns["RgMax_nm"], dtype=float)
        long_bins = np.asarray(table.columns["LongRangeBin"], dtype=float)
        long_min = np.asarray(table.columns["LongRangeMin"], dtype=float)
        long_max = np.asarray(table.columns["LongRangeMax"], dtype=float)
        for idx in range(table.row_count()):
            scope = str(scopes[idx])
            state_id = int(round(float(state_ids[idx])))
            frame_count = float(frame_counts[idx])
            if not np.isfinite(frame_count) or frame_count <= 0.0:
                continue
            key = (scope, state_id)
            entry = merged.setdefault(
                key,
                {
                    "Scope": scope,
                    "StateID": float(state_id),
                    "FrameCount": 0.0,
                    "RgBin": float(rg_bins[idx]),
                    "RgMin_nm": float(rg_min[idx]),
                    "RgMax_nm": float(rg_max[idx]),
                    "LongRangeBin": float(long_bins[idx]),
                    "LongRangeMin": float(long_min[idx]),
                    "LongRangeMax": float(long_max[idx]),
                },
            )
            entry["FrameCount"] = float(entry["FrameCount"]) + frame_count
            scope_totals[scope] += frame_count
    if not merged:
        return None, warnings
    rows = sorted(merged.values(), key=lambda row: (str(row["Scope"]), int(round(float(row["StateID"])))))
    columns: dict[str, np.ndarray] = {
        "Scope": np.asarray([str(row["Scope"]) for row in rows], dtype=object),
        "StateID": np.asarray([float(row["StateID"]) for row in rows], dtype=float),
        "FrameCount": np.asarray([float(row["FrameCount"]) for row in rows], dtype=float),
        "Probability": np.asarray(
            [
                float(row["FrameCount"]) / max(scope_totals.get(str(row["Scope"]), 0.0), 1.0)
                for row in rows
            ],
            dtype=float,
        ),
        "RgBin": np.asarray([float(row["RgBin"]) for row in rows], dtype=float),
        "RgMin_nm": np.asarray([float(row["RgMin_nm"]) for row in rows], dtype=float),
        "RgMax_nm": np.asarray([float(row["RgMax_nm"]) for row in rows], dtype=float),
        "LongRangeBin": np.asarray([float(row["LongRangeBin"]) for row in rows], dtype=float),
        "LongRangeMin": np.asarray([float(row["LongRangeMin"]) for row in rows], dtype=float),
        "LongRangeMax": np.asarray([float(row["LongRangeMax"]) for row in rows], dtype=float),
    }
    if usable < len(replicates):
        warnings.append(
            "rna_compaction_state_histogram: aggregated using the replicate(s) that produced information-proxy histograms"
        )
    return Table(columns=columns, order=list(RNA_COMPACTION_STATE_HISTOGRAM_COLUMNS), x_name="StateID", x_label="Compaction state ID"), warnings


def pool_rna_compaction_trace(
    replicates: Sequence[ReplicateData], max_time_ns: float | None = None
) -> tuple[Table | None, list[str]]:
    warnings: list[str] = []
    columns: dict[str, list[np.ndarray]] = {name: [] for name in RNA_COMPACTION_TRACE_COLUMNS}
    total_rows = 0
    for rep in replicates:
        table = rep.rna_compaction_trace
        if table is None:
            continue
        if not table.has_columns(RNA_COMPACTION_TRACE_COLUMNS):
            warnings.append(f"{rep.system}/{rep.replicate_label}: incomplete rna_compaction_trace")
            continue
        times = np.asarray(table.columns["Time_ns"], dtype=float)
        mask = np.isfinite(times)
        if max_time_ns is not None:
            mask &= times <= float(max_time_ns) + 1e-9
        if not np.any(mask):
            continue
        for name in RNA_COMPACTION_TRACE_COLUMNS:
            columns[name].append(np.asarray(table.columns[name], dtype=float)[mask])
        total_rows += int(np.count_nonzero(mask))
    if total_rows == 0:
        return None, warnings
    merged = {name: np.concatenate(parts) for name, parts in columns.items() if parts}
    order_idx = np.argsort(merged["Time_ns"])
    merged = {name: values[order_idx] for name, values in merged.items()}
    return Table(columns=merged, order=list(RNA_COMPACTION_TRACE_COLUMNS), x_name="Time_ns", x_label="Time (ns)"), warnings


def pool_rna_binding_energy_trace(
    replicates: Sequence[ReplicateData], max_time_ns: float | None = None
) -> tuple[Table | None, list[str]]:
    warnings: list[str] = []
    columns: dict[str, list[np.ndarray]] = {name: [] for name in RNA_BINDING_ENERGY_TRACE_COLUMNS}
    total_rows = 0
    for rep in replicates:
        table = rep.rna_binding_energy_trace
        if table is None:
            continue
        if not table.has_columns(RNA_BINDING_ENERGY_TRACE_COLUMNS):
            warnings.append(f"{rep.system}/{rep.replicate_label}: incomplete rna_binding_energy_trace")
            continue
        times = np.asarray(table.columns["Time_ns"], dtype=float)
        mask = np.isfinite(times)
        if max_time_ns is not None:
            mask &= times <= float(max_time_ns) + 1e-9
        if not np.any(mask):
            continue
        for name in RNA_BINDING_ENERGY_TRACE_COLUMNS:
            columns[name].append(np.asarray(table.columns[name], dtype=float)[mask])
        total_rows += int(np.count_nonzero(mask))
    if total_rows == 0:
        return None, warnings
    merged = {name: np.concatenate(parts) for name, parts in columns.items() if parts}
    order_idx = np.argsort(merged["Time_ns"])
    merged = {name: values[order_idx] for name, values in merged.items()}
    return Table(columns=merged, order=list(RNA_BINDING_ENERGY_TRACE_COLUMNS), x_name="Time_ns", x_label="Time (ns)"), warnings


def aggregate_rna_topology_binned(
    replicates: Sequence[ReplicateData],
) -> tuple[Table | None, list[str]]:
    warnings: list[str] = []
    usable: list[tuple[Table, str]] = []
    for rep in replicates:
        table = rep.rna_topology_binned
        if table is None or not table.has_columns(RNA_TOPOLOGY_BINNED_COLUMNS):
            warnings.append(f"{rep.system}/{rep.replicate_label}: missing rna_topology_binned")
            continue
        usable.append((table, rep.replicate_label))
    if not usable:
        return None, warnings

    keys: list[tuple[str, float]] = []
    seen: set[tuple[str, float]] = set()
    for table, _ in usable:
        scopes = np.asarray(table.columns["Scope"], dtype=object)
        bins = np.asarray(table.columns["BoundProteinChainsBin"], dtype=float)
        for scope_value, bin_value in zip(scopes, bins):
            if not np.isfinite(bin_value):
                continue
            key = (str(scope_value), float(round(float(bin_value), 12)))
            if key not in seen:
                seen.add(key)
                keys.append(key)
    keys.sort(key=lambda item: (item[0], item[1]))
    if not keys:
        return None, warnings

    columns: dict[str, np.ndarray] = {
        "Scope": np.asarray([scope for scope, _ in keys], dtype=object),
        "BoundProteinChainsBin": np.asarray([bin_value for _, bin_value in keys], dtype=float),
        "FrameCount": np.zeros(len(keys), dtype=float),
        "NReplicates": np.zeros(len(keys), dtype=float),
    }
    for name in RNA_TOPOLOGY_BINNED_COLUMNS:
        if name not in columns:
            columns[name] = np.full(len(keys), np.nan, dtype=float)

    metric_columns = [
        name
        for name in RNA_TOPOLOGY_BINNED_COLUMNS
        if name not in {"Scope", "BoundProteinChainsBin", "FrameCount"}
    ]
    for idx, (scope, bin_value) in enumerate(keys):
        rows: list[tuple[Table, int]] = []
        weights: list[float] = []
        for table, _ in usable:
            scopes = np.asarray(table.columns["Scope"], dtype=object)
            bins = np.asarray(table.columns["BoundProteinChainsBin"], dtype=float)
            mask = np.asarray([str(value) == scope for value in scopes], dtype=bool)
            mask &= np.isclose(bins, bin_value)
            matches = np.where(mask)[0]
            if matches.size == 0:
                continue
            row_idx = int(matches[0])
            frame_count = float(table.columns["FrameCount"][row_idx])
            if np.isfinite(frame_count) and frame_count > 0.0:
                rows.append((table, row_idx))
                weights.append(frame_count)
        if not rows:
            continue
        columns["FrameCount"][idx] = float(np.sum(weights))
        columns["NReplicates"][idx] = float(len(rows))
        for name in metric_columns:
            values = np.asarray(
                [float(table.columns[name][row_idx]) for table, row_idx in rows],
                dtype=float,
            )
            weight_arr = np.asarray(weights, dtype=float)
            mask = np.isfinite(values) & np.isfinite(weight_arr) & (weight_arr > 0.0)
            if np.any(mask):
                columns[name][idx] = float(np.sum(values[mask] * weight_arr[mask]) / np.sum(weight_arr[mask]))

    order = RNA_TOPOLOGY_BINNED_COLUMNS + ["NReplicates"]
    return Table(columns=columns, order=order, x_name="BoundProteinChainsBin", x_label="Bound protein chains"), warnings


def pool_rna_topology_trace(
    replicates: Sequence[ReplicateData], max_time_ns: float | None = None
) -> tuple[Table | None, list[str]]:
    warnings: list[str] = []
    columns: dict[str, list[np.ndarray]] = {name: [] for name in RNA_TOPOLOGY_TRACE_COLUMNS}
    total_rows = 0
    for rep in replicates:
        table = rep.rna_topology_trace
        if table is None:
            continue
        if not table.has_columns(RNA_TOPOLOGY_TRACE_COLUMNS):
            warnings.append(f"{rep.system}/{rep.replicate_label}: incomplete rna_topology_trace")
            continue
        times = np.asarray(table.columns["Time_ns"], dtype=float)
        mask = np.isfinite(times)
        if max_time_ns is not None:
            mask &= times <= float(max_time_ns) + 1e-9
        if not np.any(mask):
            continue
        for name in RNA_TOPOLOGY_TRACE_COLUMNS:
            columns[name].append(np.asarray(table.columns[name], dtype=float)[mask])
        total_rows += int(np.count_nonzero(mask))
    if total_rows == 0:
        return None, warnings
    merged = {name: np.concatenate(parts) for name, parts in columns.items() if parts}
    order_idx = np.argsort(merged["Time_ns"])
    merged = {name: values[order_idx] for name, values in merged.items()}
    return Table(columns=merged, order=list(RNA_TOPOLOGY_TRACE_COLUMNS), x_name="Time_ns", x_label="Time (ns)"), warnings


def aggregate_numeric_summary(
    replicates: Sequence[ReplicateData], attr_name: str
) -> tuple[dict[str, object], list[str]]:
    warnings: list[str] = []
    usable: list[tuple[dict[str, object], float]] = []
    for rep in replicates:
        summary = getattr(rep, attr_name)
        if not summary:
            continue
        usable.append((summary, float(max(1, rep.frame_weight))))
    if not usable:
        return {}, warnings
    out: dict[str, object] = {}
    for key in sorted({key for summary, _ in usable for key in summary}):
        values: list[float] = []
        weights: list[float] = []
        for summary, weight in usable:
            value = summary.get(key)
            if isinstance(value, bool):
                continue
            try:
                numeric = float(value)
            except (TypeError, ValueError):
                continue
            if np.isfinite(numeric):
                values.append(numeric)
                weights.append(weight)
        if values:
            out[key] = weighted_mean(values, weights)
    return out, warnings


def infer_steady_start(system: AggregatedSystem, arg_value: str) -> float:
    if arg_value.lower() != "auto":
        try:
            return max(0.0, float(arg_value))
        except ValueError as exc:
            raise SystemExit(f"Invalid --steady-start-ns value: {arg_value}") from exc

    perc_time = float(system.llps_summary.get("percolation_time_ns", -1.0) or -1.0)
    time_candidates: list[float] = []
    if perc_time > 0:
        time_candidates.append(perc_time)

    max_duration = 0.0
    if system.llps_timeseries is not None and system.llps_timeseries.row_count():
        max_duration = max(max_duration, float(np.max(system.llps_timeseries.columns["time_ns"])))
    if system.dbscan_timeseries is not None and system.dbscan_timeseries.row_count():
        max_duration = max(max_duration, float(np.max(system.dbscan_timeseries.columns["Time_ps"])) / 1000.0)
    if max_duration > 0:
        time_candidates.append(max_duration * 0.5)

    return max(time_candidates) if time_candidates else 0.0


def slim_replicate_metadata(replicates: Sequence[ReplicateData]) -> list[ReplicateData]:
    slimmed: list[ReplicateData] = []
    for rep in replicates:
        slimmed.append(
            ReplicateData(
                system=rep.system,
                replicate_label=rep.replicate_label,
                source_dir=rep.source_dir,
                frame_weight=rep.frame_weight,
                warnings=list(rep.warnings),
            )
        )
    return slimmed


def aggregate_system(
    system: str,
    replicates: Sequence[ReplicateData],
    steady_start_arg: str,
    options: RuntimeOptions,
) -> AggregatedSystem:
    aggregated = AggregatedSystem(system=system, replicates=slim_replicate_metadata(replicates))

    aggregated.llps_timeseries, warnings, _ = aggregate_time_aligned(
        replicates,
        "llps_timeseries",
        ["time_ns", "gcc_fraction", "gcc_size", "gcc_rg_nm", "mean_degree", "rho_dense", "rho_dilute"],
        max_x=options.max_time_ns,
    )
    aggregated.warnings.extend(warnings)

    aggregated.llps_autocorr, warnings, _ = aggregate_time_aligned(
        replicates,
        "llps_autocorr",
        ["lag_ns", "ct"],
    )
    aggregated.warnings.extend(warnings)

    aggregated.dbscan_timeseries, warnings, _ = aggregate_time_aligned(
        replicates,
        "dbscan_timeseries",
        [
            "Time_ps",
            "TotalBeads",
            "NumCondensates",
            "GiantCluster_ProChains",
            "GiantCluster_RnaChains",
            "NoiseBeads",
            "Rg_Pro_Dilute",
            "Rg_Pro_Dense",
            "Rg_RNA_Dilute",
            "Rg_RNA_Dense",
        ],
        max_x=(options.max_time_ns * 1000.0) if options.max_time_ns is not None else None,
    )
    aggregated.warnings.extend(warnings)

    aggregated.llps_degree, warnings = aggregate_degree_distribution(replicates)
    aggregated.warnings.extend(warnings)

    aggregated.llps_summary, warnings = aggregate_summary(replicates)
    aggregated.warnings.extend(warnings)

    aggregated.polymer_scaling, warnings = aggregate_polymer_scaling(replicates)
    aggregated.warnings.extend(warnings)

    if options.include_density:
        aggregated.density_protein, warnings = aggregate_density_grids(replicates, "density_protein")
        aggregated.warnings.extend(warnings)
        aggregated.density_rna, warnings = aggregate_density_grids(replicates, "density_rna")
        aggregated.warnings.extend(warnings)
        if aggregated.density_protein is not None and aggregated.density_protein.shape_metrics:
            aggregated.density_shape_metrics["protein"] = dict(aggregated.density_protein.shape_metrics)
        if aggregated.density_rna is not None and aggregated.density_rna.shape_metrics:
            aggregated.density_shape_metrics["rna"] = dict(aggregated.density_rna.shape_metrics)
    else:
        aggregated.warnings.append("Density loading skipped by flag (--skip-density).")

    aggregated.nucleic_interactions, warnings = aggregate_nucleic_interactions(replicates)
    aggregated.warnings.extend(warnings)
    aggregated.sidechain_protein, warnings = aggregate_sidechain_contacts(replicates, "sidechain_protein")
    aggregated.warnings.extend(warnings)
    aggregated.sidechain_rna, warnings = aggregate_sidechain_contacts(replicates, "sidechain_rna")
    aggregated.warnings.extend(warnings)

    if aggregated.nucleic_interactions:
        aggregated.interface_mode = "nucleic"
    elif aggregated.sidechain_protein or aggregated.sidechain_rna:
        aggregated.interface_mode = "sidechain"
    else:
        aggregated.interface_mode = "none"

    if options.include_network:
        aggregated.network_bundle, warnings = aggregate_network_bundle(replicates)
        aggregated.warnings.extend(warnings)
    else:
        aggregated.warnings.append("Network/allostery loading skipped by default. Use --with-network to enable it.")
    if options.include_network_anchors:
        aggregated.network_anchor_bundle, warnings = aggregate_network_anchor_bundle(replicates)
        aggregated.warnings.extend(warnings)
    else:
        aggregated.warnings.append(
            "Low-memory network-anchor loading skipped by default. Use --with-network-anchors to enable Fig25 without full edge maps."
        )
    if options.include_tda:
        aggregated.tda_points, warnings = aggregate_tda_points(replicates)
        aggregated.warnings.extend(warnings)
    else:
        aggregated.warnings.append("TDA loading skipped by default. Use --with-tda to enable it.")

    aggregated.sticker_contact_lifetimes, warnings = aggregate_lifetime_bundle(
        replicates, "sticker_contact_lifetimes"
    )
    aggregated.warnings.extend(warnings)
    aggregated.protein_rna_contact_lifetimes, warnings = aggregate_lifetime_bundle(
        replicates, "protein_rna_contact_lifetimes"
    )
    aggregated.warnings.extend(warnings)
    aggregated.neighbor_rewiring, warnings, _ = aggregate_dynamics_bundle(
        replicates,
        "neighbor_rewiring",
        [
            "time_ns",
            "dt_ns",
            "protein_persistence",
            "rna_persistence",
            "global_persistence",
            "protein_changes",
            "rna_changes",
            "global_changes",
        ],
        max_x=options.max_time_ns,
    )
    aggregated.warnings.extend(warnings)
    aggregated.interface_turnover, warnings, _ = aggregate_dynamics_bundle(
        replicates,
        "interface_turnover",
        ["time_ns", "dt_ns", "new_contacts", "lost_contacts", "new_per_ns", "lost_per_ns", "net_per_ns"],
        max_x=options.max_time_ns,
    )
    aggregated.warnings.extend(warnings)
    aggregated.dense_msd_protein, warnings, _ = aggregate_dense_msd_bundle(
        replicates,
        "dense_msd_protein",
        max_x=options.max_time_ns,
    )
    aggregated.warnings.extend(warnings)
    aggregated.dense_msd_rna, warnings, _ = aggregate_dense_msd_bundle(
        replicates,
        "dense_msd_rna",
        max_x=options.max_time_ns,
    )
    aggregated.warnings.extend(warnings)
    aggregated.liquidity_summary, warnings = aggregate_numeric_summary(replicates, "liquidity_summary")
    aggregated.warnings.extend(warnings)
    aggregated.rna_compaction_summary, warnings = aggregate_numeric_summary(replicates, "rna_compaction_summary")
    aggregated.warnings.extend(warnings)
    aggregated.rna_compaction_binned, warnings = aggregate_rna_compaction_binned(replicates)
    aggregated.warnings.extend(warnings)
    aggregated.rna_energy_compaction_summary, warnings = aggregate_numeric_summary(
        replicates,
        "rna_energy_compaction_summary",
    )
    aggregated.warnings.extend(warnings)
    aggregated.rna_energy_compaction_binned, warnings = aggregate_rna_energy_compaction_binned(replicates)
    aggregated.warnings.extend(warnings)
    aggregated.rna_information_summary, warnings = aggregate_numeric_summary(replicates, "rna_information_summary")
    aggregated.warnings.extend(warnings)
    aggregated.rna_compaction_state_histogram, warnings = aggregate_rna_compaction_state_histogram(replicates)
    aggregated.warnings.extend(warnings)
    aggregated.rna_topology_summary, warnings = aggregate_numeric_summary(replicates, "rna_topology_summary")
    aggregated.warnings.extend(warnings)
    aggregated.rna_topology_binned, warnings = aggregate_rna_topology_binned(replicates)
    aggregated.warnings.extend(warnings)
    aggregated.rna_chain_compaction_summary, warnings = aggregate_rna_chain_compaction_table(replicates)
    aggregated.warnings.extend(warnings)
    aggregated.rna_intrachain_contacts_summary, warnings = aggregate_rna_intrachain_contacts_table(replicates)
    aggregated.warnings.extend(warnings)
    aggregated.rna_compaction_trace, warnings = pool_rna_compaction_trace(
        replicates,
        max_time_ns=options.max_time_ns,
    )
    aggregated.warnings.extend(warnings)
    aggregated.rna_binding_energy_trace, warnings = pool_rna_binding_energy_trace(
        replicates,
        max_time_ns=options.max_time_ns,
    )
    aggregated.warnings.extend(warnings)
    aggregated.rna_topology_trace, warnings = pool_rna_topology_trace(
        replicates,
        max_time_ns=options.max_time_ns,
    )
    aggregated.warnings.extend(warnings)

    aggregated.steady_start_ns = infer_steady_start(aggregated, steady_start_arg)
    aggregated.rg_rna_dense_values = collect_dense_rg_values(
        replicates,
        "Rg_RNA_Dense",
        aggregated.steady_start_ns,
        options.max_time_ns,
    )
    aggregated.liquidity_metrics = compute_liquidity_metrics(
        replicates,
        aggregated.steady_start_ns,
        options.max_time_ns,
        float(aggregated.llps_summary.get("tau_ns", 0.0) or 0.0),
    )
    if aggregated.sticker_contact_lifetimes is not None:
        aggregated.liquidity_metrics["mean_sticker_contact_lifetime_ns"] = float(
            aggregated.sticker_contact_lifetimes.summary.get("mean_ns", float("nan"))
        )
    if aggregated.protein_rna_contact_lifetimes is not None:
        aggregated.liquidity_metrics["mean_protein_rna_contact_lifetime_ns"] = float(
            aggregated.protein_rna_contact_lifetimes.summary.get("mean_ns", float("nan"))
        )
    if aggregated.neighbor_rewiring is not None:
        for key, value in aggregated.neighbor_rewiring.summary.items():
            aggregated.liquidity_metrics[key] = float(value)
    if aggregated.interface_turnover is not None:
        for key, value in aggregated.interface_turnover.summary.items():
            aggregated.liquidity_metrics[key] = float(value)
    for key, value in aggregated.liquidity_summary.items():
        try:
            numeric_value = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(numeric_value):
            aggregated.liquidity_metrics[key] = numeric_value
    if aggregated.dense_msd_protein is not None:
        diffusion = aggregated.dense_msd_protein.summary.get("diffusion_nm2_per_ns")
        if diffusion is not None and np.isfinite(diffusion):
            aggregated.liquidity_metrics["protein_diffusion_nm2_per_ns"] = float(diffusion)
    if aggregated.dense_msd_rna is not None:
        diffusion = aggregated.dense_msd_rna.summary.get("diffusion_nm2_per_ns")
        if diffusion is not None and np.isfinite(diffusion):
            aggregated.liquidity_metrics["rna_diffusion_nm2_per_ns"] = float(diffusion)
    return aggregated


def filtered_table_after_time(table: Table | None, threshold: float, time_column: str) -> Table | None:
    if table is None or table.row_count() == 0 or time_column not in table.columns:
        return None
    mask = np.asarray(table.columns[time_column], dtype=float) >= threshold
    if not np.any(mask):
        return None
    columns = {name: values[mask] for name, values in table.columns.items()}
    return Table(columns=columns, order=list(table.order), x_name=table.x_name, x_label=table.x_label)


def compute_plateau_mean(table: Table | None, time_threshold: float, time_column: str, value_column: str) -> float | None:
    filtered = filtered_table_after_time(table, time_threshold, time_column)
    if filtered is None or value_column not in filtered.columns:
        return None
    return float(np.mean(filtered.columns[value_column]))


def values_after_time(
    table: Table | None,
    threshold: float,
    time_column: str,
    value_column: str,
    max_time: float | None = None,
) -> np.ndarray:
    if table is None or table.row_count() == 0 or time_column not in table.columns or value_column not in table.columns:
        return np.asarray([], dtype=float)
    x = np.asarray(table.columns[time_column], dtype=float)
    y = np.asarray(table.columns[value_column], dtype=float)
    mask = x >= threshold
    if max_time is not None:
        mask &= x <= max_time
    if not np.any(mask):
        return np.asarray([], dtype=float)
    return y[mask]


def compute_mean_sd_cv(values: np.ndarray) -> tuple[float, float, float] | None:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return None
    mean = float(np.mean(arr))
    sd = float(np.std(arr))
    cv = float(sd / mean) if mean != 0.0 else float("nan")
    return mean, sd, cv


def compute_liquidity_metrics(
    replicates: Sequence[ReplicateData],
    steady_start_ns: float,
    max_time_ns: float | None,
    tau_ns: float,
) -> dict[str, float]:
    metrics: dict[str, float] = {}
    metrics["tau_ns"] = float(tau_ns)
    max_time_ps = (max_time_ns * 1000.0) if max_time_ns is not None else None
    threshold_ps = steady_start_ns * 1000.0
    definitions = {
        "num_condensates": "NumCondensates",
        "giant_cluster_prochains": "GiantCluster_ProChains",
        "giant_cluster_rnachains": "GiantCluster_RnaChains",
        "rg_pro_dense": "Rg_Pro_Dense",
        "rg_rna_dense": "Rg_RNA_Dense",
    }
    for prefix, column in definitions.items():
        pooled: list[np.ndarray] = []
        for rep in replicates:
            values = values_after_time(rep.dbscan_timeseries, threshold_ps, "Time_ps", column, max_time_ps)
            if values.size:
                pooled.append(values)
        if not pooled:
            continue
        combined = np.concatenate(pooled)
        stats = compute_mean_sd_cv(combined)
        if stats is None:
            continue
        mean, sd, cv = stats
        metrics[f"mean_{prefix}"] = mean
        metrics[f"sd_{prefix}"] = sd
        metrics[f"cv_{prefix}"] = cv
        metrics[f"n_{prefix}"] = float(combined.size)
    return metrics


def collect_dense_rg_values(
    replicates: Sequence[ReplicateData],
    column: str,
    steady_start_ns: float,
    max_time_ns: float | None,
) -> np.ndarray | None:
    max_time_ps = (max_time_ns * 1000.0) if max_time_ns is not None else None
    threshold_ps = steady_start_ns * 1000.0
    pooled: list[np.ndarray] = []
    for rep in replicates:
        values = values_after_time(rep.dbscan_timeseries, threshold_ps, "Time_ps", column, max_time_ps)
        if values.size == 0:
            continue
        clean = values[np.isfinite(values) & (values > 0.0)]
        if clean.size:
            pooled.append(clean)
    if not pooled:
        return None
    return np.concatenate(pooled)


def save_csv_table(table: Table, path: Path) -> None:
    def format_value(value: object) -> object:
        if isinstance(value, (str, np.str_)):
            return str(value)
        if isinstance(value, (bytes, bytearray)):
            return value.decode("utf-8", errors="ignore")
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            return str(value)
        if np.isfinite(numeric) and abs(numeric - round(numeric)) < 1e-9:
            return int(round(numeric))
        return f"{numeric:.8f}"

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(table.order)
        for row in table.as_rows():
            writer.writerow([format_value(value) for value in row])


def save_summary_text(system: AggregatedSystem, path: Path) -> None:
    plateau_time = system.steady_start_ns
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"system: {system.system}\n")
        handle.write(f"display_name: {SYSTEM_DISPLAY[system.system]}\n")
        handle.write(f"replicate_count: {len(system.replicates)}\n")
        handle.write(f"total_weighted_frames: {system.total_weight}\n")
        handle.write(f"steady_start_ns: {plateau_time:.6f}\n")
        for rep in system.replicates:
            handle.write(f"replicate: {rep.replicate_label} dir={rep.source_dir} frame_weight={rep.frame_weight}\n")
        for key in ("percolated", "percolation_time_ns", "tau_ns", "n_macronodes", "dimer_size", "mean_gcc_rg_nm", "mean_gcc_fraction"):
            if key in system.llps_summary:
                handle.write(f"{key}: {system.llps_summary[key]}\n")
        plateau_metrics = {
            "mean_num_condensates_after_steady": compute_plateau_mean(system.dbscan_timeseries, plateau_time * 1000.0, "Time_ps", "NumCondensates"),
            "mean_giant_prochains_after_steady": compute_plateau_mean(system.dbscan_timeseries, plateau_time * 1000.0, "Time_ps", "GiantCluster_ProChains"),
            "mean_giant_rnachains_after_steady": compute_plateau_mean(system.dbscan_timeseries, plateau_time * 1000.0, "Time_ps", "GiantCluster_RnaChains"),
            "mean_rg_pro_dense_after_steady": compute_plateau_mean(system.dbscan_timeseries, plateau_time * 1000.0, "Time_ps", "Rg_Pro_Dense"),
            "mean_rg_rna_dense_after_steady": compute_plateau_mean(system.dbscan_timeseries, plateau_time * 1000.0, "Time_ps", "Rg_RNA_Dense"),
        }
        for key, value in plateau_metrics.items():
            if value is not None:
                handle.write(f"{key}: {value:.6f}\n")
        for key in LIQUIDITY_REPORT_KEYS:
            if key in system.liquidity_metrics:
                handle.write(f"{key}: {system.liquidity_metrics[key]:.6f}\n")
        if system.rna_compaction_summary:
            handle.write("rna_compaction_binding_summary:\n")
            for key in sorted(system.rna_compaction_summary):
                handle.write(f"  {key}: {system.rna_compaction_summary[key]}\n")
        if system.rna_energy_compaction_summary:
            handle.write("rna_energy_compaction_summary:\n")
            for key in sorted(system.rna_energy_compaction_summary):
                handle.write(f"  {key}: {system.rna_energy_compaction_summary[key]}\n")
        if system.rna_information_summary:
            handle.write("rna_information_summary:\n")
            for key in sorted(system.rna_information_summary):
                handle.write(f"  {key}: {system.rna_information_summary[key]}\n")
        if system.rna_topology_summary:
            handle.write("rna_topology_summary:\n")
            for key in sorted(system.rna_topology_summary):
                handle.write(f"  {key}: {system.rna_topology_summary[key]}\n")
        if system.warnings:
            handle.write("warnings:\n")
            for warning in sorted(set(system.warnings)):
                handle.write(f"  - {warning}\n")


def save_density_dx(grid: DXGrid, path: Path) -> None:
    if grid.data is None:
        raise ValueError(f"Full 3D density data was not retained for {path.name}.")
    nx, ny, nz = grid.counts
    ntotal = nx * ny * nz
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"object 1 class gridpositions counts {nx} {ny} {nz}\n")
        handle.write(f"origin {grid.origin[0]} {grid.origin[1]} {grid.origin[2]}\n")
        for delta in grid.deltas:
            handle.write(f"delta {delta[0]} {delta[1]} {delta[2]}\n")
        handle.write(f"object 2 class gridconnections counts {nx} {ny} {nz}\n")
        handle.write(f"object 3 class array type double rank 0 items {ntotal} data follows\n")
        flat = grid.data.reshape(-1)
        for idx in range(0, flat.size, 3):
            chunk = flat[idx : idx + 3]
            handle.write(" ".join(f"{value:.6f}" for value in chunk) + "\n")
        handle.write(f'object "{grid.label}" class field\n')


def save_nucleic_interactions(entries: Sequence[NucleicInteraction], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "Protein_ResNum",
                "ResName",
                "Base_Contact_Frames",
                "Backbone_Contact_Frames",
                "Base_Occupancy_%",
                "Backbone_Occupancy_%",
                "Recognition_Type",
            ]
        )
        for entry in entries:
            writer.writerow(
                [
                    entry.resnum,
                    entry.resname,
                    f"{entry.base_frames:.6f}",
                    f"{entry.backbone_frames:.6f}",
                    f"{entry.base_occupancy:.6f}",
                    f"{entry.backbone_occupancy:.6f}",
                    entry.recognition_type,
                ]
            )


def save_sidechain_contacts(entries: Sequence[SidechainContact], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Chain", "ResNum", "ResName", "Score"])
        for entry in entries:
            writer.writerow([entry.chain, entry.resnum, entry.resname, f"{entry.score:.6f}"])


def save_network_bundle(bundle: NetworkBundle, prefix: str, outdir: Path) -> None:
    node_path = outdir / f"{prefix}_avg_network_nodes.csv"
    with node_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = ["Id", "Label", "Type", "AtomName"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for node_id, row in sorted(bundle.nodes.items(), key=lambda item: natural_key(item[0])):
            writer.writerow(
                {
                    "Id": node_id,
                    "Label": row.get("Label", ""),
                    "Type": row.get("Type", ""),
                    "AtomName": row.get("AtomName", ""),
                }
            )

    def write_edge_csv(name: str, header: str, mapping: dict[tuple[str, str], float]) -> None:
        path = outdir / f"{prefix}_avg_{name}.csv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["Source", "Target", header])
            for (source, target), value in sorted(mapping.items(), key=lambda item: (natural_key(item[0][0]), natural_key(item[0][1]))):
                writer.writerow([source, target, f"{value:.8f}"])

    write_edge_csv("network_edges", "Correlation", bundle.corr_edges)
    write_edge_csv("network_gc_edges", "GC", bundle.gc_edges)
    write_edge_csv("network_te_edges", "Net_TE", bundle.te_edges)
    write_edge_csv("raw_mi_edges", "MI", bundle.mi_edges)


def save_tda_points(points: np.ndarray, path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["X", "Y", "Z"])
        for x, y, z in points:
            writer.writerow([f"{x:.8f}", f"{y:.8f}", f"{z:.8f}"])


def save_dynamics_bundle(bundle: DynamicsBundle, path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        for key, value in sorted(bundle.summary.items()):
            writer.writerow([f"# {key}", f"{value:.8f}"])
        if bundle.table is not None:
            writer.writerow(bundle.table.order)
            for row in bundle.table.as_rows():
                writer.writerow([f"{value:.8f}" for value in row])


def save_key_value_summary(summary: dict[str, object], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for key in sorted(summary):
            value = summary[key]
            if isinstance(value, float):
                if np.isfinite(value):
                    handle.write(f"{key}: {value:.8f}\n")
                else:
                    handle.write(f"{key}: NA\n")
            else:
                handle.write(f"{key}: {value}\n")


SHAPE_METRIC_EXPORT_ORDER = [
    "center_x_nm",
    "center_y_nm",
    "center_z_nm",
    "rg_nm",
    "lambda1_nm2",
    "lambda2_nm2",
    "lambda3_nm2",
    "principal_radius1_nm",
    "principal_radius2_nm",
    "principal_radius3_nm",
    "asphericity_nm2",
    "normalized_asphericity",
    "acylindricity_nm2",
    "relative_shape_anisotropy",
    "aspect_ratio",
    "density_mass",
    "replicates_used",
]


def save_density_shape_metrics(system: AggregatedSystem, path: Path) -> None:
    keys = [key for key in SHAPE_METRIC_EXPORT_ORDER if any(key in metrics for metrics in system.density_shape_metrics.values())]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["species"] + keys)
        for species in ("protein", "rna"):
            metrics = system.density_shape_metrics.get(species, {})
            if not metrics:
                continue
            writer.writerow(
                [species]
                + [
                    f"{float(metrics[key]):.8f}" if key in metrics and np.isfinite(metrics[key]) else ""
                    for key in keys
                ]
            )


def save_density_shape_compare(systems: dict[str, AggregatedSystem], path: Path) -> None:
    keys = [key for key in SHAPE_METRIC_EXPORT_ORDER if any(key in system.density_shape_metrics.get(species, {}) for system in systems.values() for species in ("protein", "rna"))]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["system", "display_name", "species"] + keys)
        for system_key in SYSTEM_ORDER:
            system = systems[system_key]
            for species in ("protein", "rna"):
                metrics = system.density_shape_metrics.get(species, {})
                if not metrics:
                    continue
                writer.writerow(
                    [system_key, SYSTEM_DISPLAY[system_key], species]
                    + [
                        f"{float(metrics[key]):.8f}" if key in metrics and np.isfinite(metrics[key]) else ""
                        for key in keys
                    ]
                )


def save_liquidity_metrics_compare(systems: dict[str, AggregatedSystem], path: Path) -> None:
    metric_rows = [
        ("tau_ns", "Contact relaxation time tau (ns)"),
        ("cv_num_condensates", "CV of condensate count"),
        ("cv_giant_cluster_prochains", "CV of giant-cluster protein chains"),
        ("cv_giant_cluster_rnachains", "CV of giant-cluster RNA chains"),
        ("cv_rg_pro_dense", "CV of dense-phase protein Rg"),
        ("cv_rg_rna_dense", "CV of dense-phase RNA Rg"),
        ("mean_sticker_contact_lifetime_ns", "Mean sticker-contact lifetime (ns)"),
        ("mean_protein_rna_contact_lifetime_ns", "Mean protein-RNA contact lifetime (ns)"),
        ("global_neighbor_persistence", "Mean neighbor persistence"),
        ("global_rewiring_rate_per_ns", "Global rewiring rate (events/ns)"),
        ("new_contacts_per_ns", "New protein-RNA contacts per ns"),
        ("lost_contacts_per_ns", "Lost protein-RNA contacts per ns"),
        ("protein_diffusion_nm2_per_ns", r"Protein dense-phase diffusion (nm2/ns)"),
        ("rna_diffusion_nm2_per_ns", r"RNA dense-phase diffusion (nm2/ns)"),
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "label", "phos", "unphos"])
        for key, label in metric_rows:
            writer.writerow(
                [
                    key,
                    label,
                    f"{systems['phos'].liquidity_metrics.get(key, float('nan')):.8f}",
                    f"{systems['unphos'].liquidity_metrics.get(key, float('nan')):.8f}",
                ]
            )


def save_diffusion_metrics_compare(systems: dict[str, AggregatedSystem], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "species",
                "phos_diffusion_nm2_per_ns",
                "unphos_diffusion_nm2_per_ns",
                "phos_over_unphos",
                "phos_fit_start_ns",
                "phos_fit_end_ns",
                "unphos_fit_start_ns",
                "unphos_fit_end_ns",
            ]
        )
        for species, attr_name in (("protein", "dense_msd_protein"), ("rna", "dense_msd_rna")):
            phos_bundle = getattr(systems["phos"], attr_name)
            unphos_bundle = getattr(systems["unphos"], attr_name)
            phos_d = float(phos_bundle.summary.get("diffusion_nm2_per_ns", float("nan"))) if phos_bundle else float("nan")
            unphos_d = float(unphos_bundle.summary.get("diffusion_nm2_per_ns", float("nan"))) if unphos_bundle else float("nan")
            ratio = phos_d / unphos_d if np.isfinite(phos_d) and np.isfinite(unphos_d) and unphos_d != 0.0 else float("nan")
            row = [
                species,
                phos_d,
                unphos_d,
                ratio,
                float(phos_bundle.summary.get("fit_start_ns", float("nan"))) if phos_bundle else float("nan"),
                float(phos_bundle.summary.get("fit_end_ns", float("nan"))) if phos_bundle else float("nan"),
                float(unphos_bundle.summary.get("fit_start_ns", float("nan"))) if unphos_bundle else float("nan"),
                float(unphos_bundle.summary.get("fit_end_ns", float("nan"))) if unphos_bundle else float("nan"),
            ]
            writer.writerow([row[0]] + [f"{value:.8f}" if np.isfinite(value) else "" for value in row[1:]])


def rg_histogram_table(systems: dict[str, AggregatedSystem], bins: int = 32) -> Table | None:
    values_by_system: dict[str, np.ndarray] = {}
    for system in SYSTEM_ORDER:
        values = systems[system].rg_rna_dense_values
        if values is None:
            continue
        clean = np.asarray(values, dtype=float)
        clean = clean[np.isfinite(clean) & (clean > 0.0)]
        if clean.size:
            values_by_system[system] = clean
    if not values_by_system:
        return None
    combined = np.concatenate(list(values_by_system.values()))
    lo = float(np.min(combined))
    hi = float(np.max(combined))
    if hi <= lo:
        hi = lo + 1.0
    edges = np.linspace(lo, hi, max(4, int(bins)) + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    columns: dict[str, np.ndarray] = {
        "bin_left_nm": edges[:-1],
        "bin_right_nm": edges[1:],
        "bin_center_nm": centers,
    }
    order = ["bin_left_nm", "bin_right_nm", "bin_center_nm"]
    for system in SYSTEM_ORDER:
        values = values_by_system.get(system, np.asarray([], dtype=float))
        if values.size:
            counts, _ = np.histogram(values, bins=edges)
            density, _ = np.histogram(values, bins=edges, density=True)
        else:
            counts = np.zeros(centers.size, dtype=int)
            density = np.zeros(centers.size, dtype=float)
        columns[f"{system}_count"] = counts.astype(float)
        columns[f"{system}_probability_density"] = density.astype(float)
        order.extend([f"{system}_count", f"{system}_probability_density"])
    return Table(columns=columns, order=order, x_name="bin_center_nm", x_label="Dense-phase RNA Rg (nm)")


def save_rg_rna_dense_histogram_compare(systems: dict[str, AggregatedSystem], path: Path) -> None:
    table = rg_histogram_table(systems)
    if table is not None:
        save_csv_table(table, path)


def save_rna_compaction_binned_compare(systems: dict[str, AggregatedSystem], path: Path) -> None:
    metric_columns = [
        "FrameCount",
        "RgMean_nm",
        "RgStd_nm",
        "EndToEndMean_nm",
        "EndToEndStd_nm",
        "AsphericityMean",
        "AsphericityStd",
        "LongRangeContactFractionMean",
        "LongRangeContactFractionStd",
        "NReplicates",
    ]
    bins = unique_sorted(
        bin_value
        for system in SYSTEM_ORDER
        if systems[system].rna_compaction_binned is not None
        for bin_value in systems[system].rna_compaction_binned.columns["BoundProteinChainsBin"]
    )
    if bins.size == 0:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["BoundProteinChainsBin"] + [f"{system}_{name}" for system in SYSTEM_ORDER for name in metric_columns])
        for bin_value in bins:
            row: list[str | int | float] = [int(round(float(bin_value)))]
            for system in SYSTEM_ORDER:
                table = systems[system].rna_compaction_binned
                if table is None:
                    row.extend([""] * len(metric_columns))
                    continue
                source_bins = np.asarray(table.columns["BoundProteinChainsBin"], dtype=float)
                match = np.where(np.isclose(source_bins, bin_value))[0]
                if match.size == 0:
                    row.extend([""] * len(metric_columns))
                    continue
                idx = int(match[0])
                row.extend(
                    f"{float(table.columns[name][idx]):.8f}" if name in table.columns and np.isfinite(table.columns[name][idx]) else ""
                    for name in metric_columns
                )
            writer.writerow(row)


def save_rna_energy_compaction_binned_compare(systems: dict[str, AggregatedSystem], path: Path) -> None:
    metric_columns = [
        "FrameCount",
        "EelecMean_kJmol",
        "EelecStd_kJmol",
        "EvdwMean_kJmol",
        "EvdwStd_kJmol",
        "EtotalMean_kJmol",
        "EtotalStd_kJmol",
        "AttractionStrengthMean_kJmol",
        "AttractionStrengthStd_kJmol",
        "RgMean_nm",
        "RgStd_nm",
        "EndToEndMean_nm",
        "EndToEndStd_nm",
        "RelativeShapeAnisotropyMean",
        "RelativeShapeAnisotropyStd",
        "LongRangeContactFractionMean",
        "LongRangeContactFractionStd",
        "NReplicates",
    ]
    bins = unique_sorted(
        bin_value
        for system in SYSTEM_ORDER
        if systems[system].rna_energy_compaction_binned is not None
        for bin_value in systems[system].rna_energy_compaction_binned.columns["BoundProteinChainsBin"]
    )
    if bins.size == 0:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["BoundProteinChainsBin"] + [f"{system}_{name}" for system in SYSTEM_ORDER for name in metric_columns])
        for bin_value in bins:
            row: list[str | int | float] = [int(round(float(bin_value)))]
            for system in SYSTEM_ORDER:
                table = systems[system].rna_energy_compaction_binned
                if table is None:
                    row.extend([""] * len(metric_columns))
                    continue
                source_bins = np.asarray(table.columns["BoundProteinChainsBin"], dtype=float)
                match = np.where(np.isclose(source_bins, bin_value))[0]
                if match.size == 0:
                    row.extend([""] * len(metric_columns))
                    continue
                idx = int(match[0])
                row.extend(
                    f"{float(table.columns[name][idx]):.8f}" if name in table.columns and np.isfinite(table.columns[name][idx]) else ""
                    for name in metric_columns
                )
            writer.writerow(row)


def save_rna_information_summary_compare(systems: dict[str, AggregatedSystem], path: Path) -> None:
    preferred = [
        "inside_dense_only_H_compaction_bits",
        "inside_dense_only_H_compaction_norm",
        "inside_dense_only_H_partner_bits",
        "inside_dense_only_H_partner_norm",
        "inside_dense_only_H_partner_given_binding_bits",
        "inside_dense_only_MI_binding_compaction_bits",
        "inside_dense_only_MI_binding_compaction_norm",
        "inside_dense_only_H_compaction_given_binding_bits",
        "inside_dense_only_transition_entropy_rate_bits_per_step",
        "inside_dense_only_transition_entropy_rate_bits_per_ns",
        "all_frames_H_compaction_bits",
        "all_frames_MI_binding_compaction_bits",
        "contacting_only_H_partner_bits",
        "contacting_only_H_partner_given_binding_bits",
    ]
    key_union = {
        key
        for system in SYSTEM_ORDER
        for key in systems[system].rna_information_summary
        if key not in {"primary_scope"}
    }
    ordered_keys = [key for key in preferred if key in key_union]
    ordered_keys.extend(sorted(key for key in key_union if key not in ordered_keys))
    if not ordered_keys:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "phos", "unphos"])
        for key in ordered_keys:
            row = [key]
            for system in SYSTEM_ORDER:
                value = systems[system].rna_information_summary.get(key)
                try:
                    numeric = float(value) if value is not None else float("nan")
                    row.append(f"{numeric:.8f}" if np.isfinite(numeric) else "")
                except (TypeError, ValueError):
                    row.append("" if value is None else str(value))
            writer.writerow(row)


def save_rna_topology_summary_compare(systems: dict[str, AggregatedSystem], path: Path) -> None:
    preferred = [
        "inside_dense_only_P_any_bridge_given_bound",
        "inside_dense_only_P_any_clamp_given_bound",
        "inside_dense_only_P_longrange_given_bound",
        "inside_dense_only_P_longrange_given_bridge",
        "inside_dense_only_P_longrange_given_clamp",
        "inside_dense_only_bridge_enrichment_over_bound",
        "inside_dense_only_clamp_enrichment_over_bound",
        "inside_dense_only_bridge_productivity_index",
        "inside_dense_only_clamp_productivity_index",
        "inside_dense_only_mean_bridge_attraction_strength_kjmol",
        "inside_dense_only_mean_nonbridge_attraction_strength_kjmol",
        "inside_dense_only_mean_prna_attraction_strength_kjmol",
        "inside_dense_only_mean_local_pp_attraction_strength_kjmol",
        "inside_dense_only_mean_local_pp_to_prna_attraction_ratio",
        "inside_dense_only_mean_local_pp_to_bridge_attraction_ratio",
    ]
    key_union = {
        key
        for system in SYSTEM_ORDER
        for key in systems[system].rna_topology_summary
    }
    ordered_keys = [key for key in preferred if key in key_union]
    ordered_keys.extend(sorted(key for key in key_union if key not in ordered_keys))
    if not ordered_keys:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "phos", "unphos"])
        for key in ordered_keys:
            row = [key]
            for system in SYSTEM_ORDER:
                value = systems[system].rna_topology_summary.get(key)
                try:
                    numeric = float(value) if value is not None else float("nan")
                    row.append(f"{numeric:.8f}" if np.isfinite(numeric) else "")
                except (TypeError, ValueError):
                    row.append("" if value is None else str(value))
            writer.writerow(row)


def save_network_anchor_compare(systems: dict[str, AggregatedSystem], path: Path, topn: int) -> None:
    rows: list[list[object]] = []
    metrics = ["gc_anchor", "te_source", "te_sink"]
    for system in SYSTEM_ORDER:
        for metric in metrics:
            entries = network_anchor_entries(systems[system], metric, topn)
            for rank, (label, value) in enumerate(entries, start=1):
                rows.append([system, SYSTEM_DISPLAY[system], metric, rank, label, f"{value:.8f}"])
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["system", "display_name", "metric", "rank", "node_label", "score"])
        writer.writerows(rows)


def save_rna_chain_compaction_compare(systems: dict[str, AggregatedSystem], path: Path) -> None:
    metric_columns = [
        "FrameCount",
        "InsideDenseFraction",
        "MeanBoundProteinChains",
        "MeanProteinRNABeadContacts",
        "MeanRg_nm",
        "StdRg_nm",
        "MeanEndToEnd_nm",
        "StdEndToEnd_nm",
        "MeanRelativeShapeAnisotropy",
        "MeanLongRangeContactFraction",
        "NReplicates",
    ]
    chain_ids = unique_sorted(
        chain_id
        for system in SYSTEM_ORDER
        if systems[system].rna_chain_compaction_summary is not None
        for chain_id in systems[system].rna_chain_compaction_summary.columns["RNAChainID"]
    )
    if chain_ids.size == 0:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["RNAChainID"] + [f"{system}_{name}" for system in SYSTEM_ORDER for name in metric_columns])
        for chain_id in chain_ids:
            row: list[str | int | float] = [int(round(float(chain_id)))]
            for system in SYSTEM_ORDER:
                table = systems[system].rna_chain_compaction_summary
                if table is None:
                    row.extend([""] * len(metric_columns))
                    continue
                source_ids = np.asarray(table.columns["RNAChainID"], dtype=float)
                match = np.where(np.isclose(source_ids, chain_id))[0]
                if match.size == 0:
                    row.extend([""] * len(metric_columns))
                    continue
                idx = int(match[0])
                row.extend(
                    f"{float(table.columns[name][idx]):.8f}" if name in table.columns and np.isfinite(table.columns[name][idx]) else ""
                    for name in metric_columns
                )
            writer.writerow(row)


def save_figure_all_formats(fig: plt.Figure, base_path: Path, dpi: int) -> list[str]:
    notes: list[str] = []
    png_path = base_path.with_suffix(".png")
    pdf_path = base_path.with_suffix(".pdf")
    tiff_path = base_path.with_suffix(".tiff")
    for axis in fig.axes:
        axis.grid(False, which="both", axis="both")
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    try:
        fig.savefig(
            tiff_path,
            dpi=dpi,
            bbox_inches="tight",
            facecolor="white",
            pil_kwargs={"compression": "tiff_lzw"},
        )
    except Exception as exc:  # pragma: no cover - depends on Pillow backend
        fig.savefig(tiff_path, dpi=dpi, bbox_inches="tight", facecolor="white")
        notes.append(f"TIFF LZW compression unavailable for {tiff_path.name}: {exc}")
    finally:
        plt.close(fig)
    return notes


def system_color(system: str) -> str:
    return SYSTEM_COLORS[system]


def shaded_color(system: str) -> str:
    return SYSTEM_LIGHT[system]


def dense_rg_color(system: str, molecule: str) -> str:
    return DENSE_RG_COLORS[(system, molecule)]


def residue_axis(residue_min: int, residue_max: int) -> np.ndarray:
    return np.arange(residue_min, residue_max + 1, dtype=int)


def draw_missing_axis(ax: plt.Axes, message: str) -> None:
    ax.axis("off")
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=9, color="#555555")


def scope_mask_from_table(table: Table, preferred_scopes: Sequence[str]) -> np.ndarray:
    if "Scope" not in table.columns:
        return np.ones(table.row_count(), dtype=bool)
    scopes = np.asarray(table.columns["Scope"], dtype=object)
    for scope in preferred_scopes:
        mask = np.asarray([str(value) == scope for value in scopes], dtype=bool)
        if np.any(mask):
            return mask
    return np.zeros(table.row_count(), dtype=bool)


def add_steady_window(ax: plt.Axes, start_ns: float, x_unit: str) -> None:
    if start_ns <= 0:
        return
    if x_unit == "ps":
        start_value = start_ns * 1000.0
    else:
        start_value = start_ns
    xmin, xmax = ax.get_xlim()
    if xmax <= start_value:
        return
    ax.axvspan(start_value, xmax, color="#f0f0f0", alpha=0.35, zorder=0)


def plot_metric_compare(ax: plt.Axes, systems: dict[str, AggregatedSystem], table_attr: str, x_col: str, y_col: str, title: str, ylabel: str, x_scale: float = 1.0) -> None:
    plotted = False
    for system in SYSTEM_ORDER:
        table = getattr(systems[system], table_attr)
        if table is None or y_col not in table.columns or x_col not in table.columns:
            continue
        x = np.asarray(table.columns[x_col], dtype=float) / x_scale
        y = np.asarray(table.columns[y_col], dtype=float)
        if x.size == 0 or y.size == 0:
            continue
        ax.plot(x, y, color=system_color(system), linewidth=1.6, label=SYSTEM_DISPLAY[system])
        plotted = True
    if not plotted:
        draw_missing_axis(ax, f"{title}\nNot available")
        return
    ax.set_title(title)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(ylabel)
    ax.grid(True, axis="y")
    ax.legend(loc="best")
    steady = max(systems[system].steady_start_ns for system in SYSTEM_ORDER)
    add_steady_window(ax, steady, "ns")


def build_density_projection(grid: DXGrid) -> np.ndarray:
    projection = grid.projection_xy()
    if projection.size == 0:
        return projection
    return (projection * grid.spacing_z).T


def plot_density_panel(ax: plt.Axes, grid: DXGrid | None, title: str, cmap: mcolors.Colormap, vmax: float | None) -> None:
    if grid is None:
        draw_missing_axis(ax, f"{title}\nNot available")
        return
    projection = build_density_projection(grid)
    image = ax.imshow(
        projection,
        origin="lower",
        cmap=cmap,
        aspect="equal",
        vmax=vmax,
        extent=grid.extent_xy_nm,
    )
    ax.set_title(title, pad=7)
    ax.set_xlabel("Centered X (nm)")
    ax.set_ylabel("Centered Y (nm)")
    ax.set_box_aspect(1)
    return image


def fit_polymer_scaling(table: Table | None) -> tuple[float, float] | None:
    if table is None or table.row_count() < 3:
        return None
    s = np.asarray(table.columns["s"], dtype=float)
    r_rms = np.asarray(table.columns["R_rms_nm"], dtype=float)
    mask = (s > 0) & (r_rms > 0)
    if np.count_nonzero(mask) < 3:
        return None
    coeffs = np.polyfit(np.log10(s[mask]), np.log10(r_rms[mask]), deg=1)
    return float(coeffs[0]), float(coeffs[1])


def plot_interface_nucleic(ax: plt.Axes, system: AggregatedSystem, topn: int) -> None:
    entries = system.nucleic_interactions[:topn]
    if not entries:
        draw_missing_axis(ax, f"{SYSTEM_DISPLAY[system.system]}\nNo nucleic-acid recognition data")
        return
    labels = [f"{entry.resname}{entry.resnum}" for entry in entries][::-1]
    base = [entry.base_occupancy for entry in entries][::-1]
    backbone = [entry.backbone_occupancy for entry in entries][::-1]
    ypos = np.arange(len(entries))
    ax.barh(ypos, backbone, color="#86b6f6", label="Backbone")
    ax.barh(ypos, base, left=backbone, color=system_color(system.system), label="Base")
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Occupancy (%)")
    ax.set_title(SYSTEM_DISPLAY[system.system])
    ax.grid(True, axis="x")
    handles, labels_ = ax.get_legend_handles_labels()
    unique = dict(zip(labels_, handles))
    ax.legend(unique.values(), unique.keys(), loc="lower right")


def plot_interface_sidechain(
    ax: plt.Axes,
    entries: Sequence[SidechainContact],
    title: str,
    color: str,
    topn: int,
    show_xlabel: bool = True,
) -> None:
    if not entries:
        draw_missing_axis(ax, f"{title}\nNo sidechain-contact data")
        return
    subset = list(entries[:topn])[::-1]
    ypos = np.arange(len(subset))
    ax.barh(ypos, [entry.score for entry in subset], color=color)
    ax.set_yticks(ypos)
    ax.set_yticklabels([entry.label for entry in subset])
    if show_xlabel:
        ax.set_xlabel("Contact score (a.u.)", labelpad=6)
    else:
        ax.set_xlabel("")
    ax.set_title(title, pad=8)
    ax.grid(True, axis="x")
    ax.tick_params(axis="y", labelsize=7.5)


def protein_profile_from_nucleic(
    entries: Sequence[NucleicInteraction],
    residue_min: int,
    residue_max: int,
) -> np.ndarray:
    profile = np.zeros(residue_max - residue_min + 1, dtype=float)
    for entry in entries:
        if residue_min <= entry.resnum <= residue_max:
            profile[entry.resnum - residue_min] += entry.total_occupancy
    return profile


def protein_profile_from_sidechain(
    entries: Sequence[SidechainContact],
    residue_min: int,
    residue_max: int,
) -> np.ndarray:
    profile = np.zeros(residue_max - residue_min + 1, dtype=float)
    for entry in entries:
        if residue_min <= entry.resnum <= residue_max:
            profile[entry.resnum - residue_min] += entry.score
    return profile


def smooth_profile(values: np.ndarray, window: int) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0 or window <= 1:
        return arr.copy()
    window = min(int(window), int(arr.size))
    if window <= 1:
        return arr.copy()
    if window % 2 == 0:
        window = max(1, window - 1)
    if window <= 1:
        return arr.copy()
    kernel = np.ones(window, dtype=float) / float(window)
    pad_left = window // 2
    pad_right = window - 1 - pad_left
    padded = np.pad(arr, (pad_left, pad_right), mode="edge")
    smoothed = np.convolve(padded, kernel, mode="valid")
    return smoothed.astype(float, copy=False)


def normalize_profile_percent(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return arr.copy()
    vmax = float(np.max(arr))
    if vmax <= 0.0:
        return np.zeros_like(arr, dtype=float)
    return (arr / vmax) * 100.0


def plot_interface_profile(
    ax_main: plt.Axes,
    ax_diff: plt.Axes,
    systems: dict[str, AggregatedSystem],
    residue_min: int,
    residue_max: int,
    smooth_window: int,
) -> None:
    residues = residue_axis(residue_min, residue_max)
    use_nucleic = all(system.interface_mode == "nucleic" for system in systems.values())
    normalized_profiles: dict[str, np.ndarray] = {}
    for system_key in SYSTEM_ORDER:
        system = systems[system_key]
        if use_nucleic:
            y_raw = protein_profile_from_nucleic(system.nucleic_interactions, residue_min, residue_max)
        else:
            y_raw = protein_profile_from_sidechain(system.sidechain_protein, residue_min, residue_max)
        y_smooth = smooth_profile(y_raw, smooth_window)
        y_norm = normalize_profile_percent(y_smooth)
        normalized_profiles[system_key] = y_norm

    plotted = False
    for system_key in SYSTEM_ORDER:
        y = normalized_profiles.get(system_key)
        if y is None or not np.any(y):
            continue
        ax_main.plot(
            residues,
            y,
            color=system_color(system_key),
            linewidth=1.8,
            label=SYSTEM_DISPLAY[system_key],
        )
        plotted = True
    if not plotted:
        draw_missing_axis(ax_main, "Protein residue-wise RNA interaction profile\nNot available")
        draw_missing_axis(ax_diff, "phos - unphos profile\nNot available")
        return

    ax_main.set_xlim(residue_min, residue_max)
    ax_main.set_ylim(0.0, 100.0)
    if use_nucleic:
        ax_main.set_ylabel("Normalized RNA interaction (%)")
        ax_main.set_title("Protein residue-wise RNA recognition profile", pad=8)
        note = f"Smoothed ({max(1, smooth_window)}-res moving average); each system normalized to 100%"
    else:
        ax_main.set_ylabel("Normalized RNA interaction (%)")
        ax_main.set_title("Protein residue-wise RNA interface profile", pad=8)
        note = f"Sidechain-contact fallback; smoothed ({max(1, smooth_window)}-res window) and normalized to 100%"
    major_step = 50 if (residue_max - residue_min) >= 200 else 25
    ax_main.set_xticks(np.arange(residue_min, residue_max + 1, major_step))
    ax_main.set_yticks(np.arange(0, 101, 20))
    ax_main.grid(True, axis="y")
    ax_main.margins(x=0.0)
    ax_main.legend(loc="upper right", frameon=True, framealpha=0.95)
    ax_main.text(
        0.015,
        0.97,
        note,
        transform=ax_main.transAxes,
        ha="left",
        va="top",
        fontsize=7.5,
        color="#555555",
        bbox={"facecolor": "white", "edgecolor": "#dddddd", "alpha": 0.9, "boxstyle": "round,pad=0.2"},
    )

    phos = normalized_profiles.get("phos", np.zeros_like(residues, dtype=float))
    unphos = normalized_profiles.get("unphos", np.zeros_like(residues, dtype=float))
    delta = phos - unphos
    ax_diff.plot(residues, delta, color="#333333", linewidth=1.6)
    ax_diff.fill_between(residues, 0.0, delta, where=(delta >= 0.0), color=SYSTEM_LIGHT["phos"], alpha=0.75)
    ax_diff.fill_between(residues, 0.0, delta, where=(delta < 0.0), color=SYSTEM_LIGHT["unphos"], alpha=0.75)
    ax_diff.axhline(0.0, color="#666666", linewidth=0.9)
    ax_diff.set_xlim(residue_min, residue_max)
    ax_diff.set_ylim(-100.0, 100.0)
    ax_diff.set_xticks(np.arange(residue_min, residue_max + 1, major_step))
    ax_diff.set_yticks(np.arange(-100, 101, 50))
    ax_diff.set_xlabel("Protein residue index")
    ax_diff.set_ylabel(r"$\Delta$ normalized interaction (%)")
    ax_diff.set_title("Difference profile: phos - unphos", pad=8)
    ax_diff.grid(True, axis="y")
    ax_diff.margins(x=0.0)


def compute_tau_bar_values(systems: dict[str, AggregatedSystem]) -> tuple[list[str], list[float], list[str]]:
    labels: list[str] = []
    values: list[float] = []
    colors: list[str] = []
    for system in SYSTEM_ORDER:
        labels.append(SYSTEM_DISPLAY[system])
        values.append(float(systems[system].llps_summary.get("tau_ns", 0.0) or 0.0))
        colors.append(system_color(system))
    return labels, values, colors


def liquidity_metric_values(systems: dict[str, AggregatedSystem], metric_key: str) -> tuple[list[str], list[float], list[str]]:
    labels: list[str] = []
    values: list[float] = []
    colors: list[str] = []
    for system in SYSTEM_ORDER:
        labels.append(SYSTEM_DISPLAY[system])
        values.append(float(systems[system].liquidity_metrics.get(metric_key, 0.0) or 0.0))
        colors.append(system_color(system))
    return labels, values, colors


def plot_liquidity_bar(
    ax: plt.Axes,
    systems: dict[str, AggregatedSystem],
    metric_key: str,
    title: str,
    ylabel: str,
) -> None:
    labels, values, colors = liquidity_metric_values(systems, metric_key)
    xpos = np.arange(len(labels))
    ax.bar(xpos, values, color=colors, width=0.6)
    ax.set_xticks(xpos)
    ax.set_xticklabels(labels, rotation=0)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.grid(True, axis="y")
    upper = max(values + [1.0])
    if metric_key.startswith("cv_"):
        upper = max(upper * 1.25, 0.05)
    else:
        upper = max(upper * 1.20, 1.0)
    ax.set_ylim(0.0, upper)
    for idx, value in enumerate(values):
        ax.text(idx, value + upper * 0.03, f"{value:.3f}", ha="center", va="bottom", fontsize=8)


def plot_lifetime_survival(
    ax: plt.Axes,
    systems: dict[str, AggregatedSystem],
    attr_name: str,
    title: str,
    *,
    log_x: bool = False,
    show_legend: bool = True,
) -> None:
    plotted = False
    for system in SYSTEM_ORDER:
        bundle = getattr(systems[system], attr_name)
        if bundle is None or bundle.table is None:
            continue
        x = np.asarray(bundle.table.columns["lifetime_ns"], dtype=float)
        counts = np.asarray(bundle.table.columns["count"], dtype=float)
        mask = np.isfinite(x) & np.isfinite(counts) & (counts > 0)
        if log_x:
            mask &= x > 0.0
        if not np.any(mask):
            continue
        x = x[mask]
        counts = counts[mask]
        order = np.argsort(x)
        x = x[order]
        counts = counts[order]
        survival = np.cumsum(counts[::-1])[::-1] / max(float(np.sum(counts)), 1.0)
        ax.step(x, survival, where="post", color=system_color(system), linewidth=1.7, label=SYSTEM_DISPLAY[system])
        mean = float(bundle.summary.get("mean_ns", float("nan")))
        if np.isfinite(mean) and (not log_x or mean > 0.0):
            ax.axvline(mean, color=system_color(system), linewidth=1.0, alpha=0.55)
        plotted = True
    if plotted:
        ax.set_title(title)
        if log_x:
            ax.set_xscale("log")
            ax.set_xlabel("Lifetime (ns, log scale)")
        else:
            ax.set_xlabel("Lifetime (ns)")
        ax.set_ylabel("Survival fraction")
        ax.set_ylim(0.0, 1.05)
        ax.grid(True, axis="both")
        if show_legend:
            ax.legend(loc="best", frameon=False)
    else:
        draw_missing_axis(ax, f"{title}\nNot available")


def plot_grouped_summary_bars(
    ax: plt.Axes,
    systems: dict[str, AggregatedSystem],
    source_attr: str,
    metrics: Sequence[tuple[str, str]],
    title: str,
    ylabel: str,
    *,
    show_legend: bool = True,
) -> None:
    x = np.arange(len(metrics), dtype=float)
    width = 0.34
    plotted = False
    for idx, system in enumerate(SYSTEM_ORDER):
        bundle = getattr(systems[system], source_attr)
        values = [
            float(bundle.summary.get(key, float("nan"))) if bundle is not None else float("nan")
            for key, _ in metrics
        ]
        if any(np.isfinite(values)):
            offset = (idx - 0.5) * width
            ax.bar(x + offset, values, width=width, color=system_color(system), label=SYSTEM_DISPLAY[system])
            plotted = True
    if plotted:
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xticks(x)
        ax.set_xticklabels([label for _, label in metrics], rotation=25, ha="right")
        ax.grid(True, axis="y")
        if show_legend:
            ax.legend(loc="best", frameon=False)
    else:
        draw_missing_axis(ax, f"{title}\nNot available")


def add_system_figure_legend(fig: plt.Figure, *, y: float = 0.94) -> None:
    handles = [
        Line2D([0], [0], color=system_color(system), linewidth=1.8, label=SYSTEM_DISPLAY[system])
        for system in SYSTEM_ORDER
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, y), ncol=2, frameon=False)


def plot_dynamics_line(
    ax: plt.Axes,
    systems: dict[str, AggregatedSystem],
    source_attr: str,
    y_col: str,
    title: str,
    ylabel: str,
) -> None:
    plotted = False
    for system in SYSTEM_ORDER:
        bundle = getattr(systems[system], source_attr)
        if bundle is None or bundle.table is None or y_col not in bundle.table.columns:
            continue
        x = np.asarray(bundle.table.columns["time_ns"], dtype=float)
        y = np.asarray(bundle.table.columns[y_col], dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        if not np.any(mask):
            continue
        ax.plot(x[mask], y[mask], color=system_color(system), linewidth=1.5, label=SYSTEM_DISPLAY[system])
        plotted = True
    if plotted:
        ax.set_title(title)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel(ylabel)
        ax.grid(True, axis="both")
        ax.legend(loc="best")
    else:
        draw_missing_axis(ax, f"{title}\nNot available")


def top_network_strength(bundle: NetworkBundle | None, topn: int) -> list[tuple[str, float]]:
    if bundle is None or not bundle.corr_edges:
        return []
    strength: dict[str, float] = defaultdict(float)
    for (source, target), corr in bundle.corr_edges.items():
        strength[source] += abs(corr)
        strength[target] += abs(corr)
    ranking = sorted(strength.items(), key=lambda item: item[1], reverse=True)[:topn]
    return [(simplify_node_label(node, bundle.nodes.get(node, {}).get("Label")), value) for node, value in ranking]


def network_node_anchor_scores(
    bundle: NetworkBundle | None,
    selector: str,
    mode: str,
    topn: int,
) -> list[tuple[str, float]]:
    if bundle is None:
        return []
    mapping = getattr(bundle, selector)
    if not mapping:
        return []
    scores: dict[str, float] = defaultdict(float)
    if mode == "gc_strength":
        for (source, target), value in mapping.items():
            score = abs(float(value))
            scores[source] += score
            scores[target] += score
    elif mode == "te_source":
        outgoing: dict[str, float] = defaultdict(float)
        incoming: dict[str, float] = defaultdict(float)
        for (source, target), value in mapping.items():
            net = float(value)
            if net >= 0.0:
                outgoing[source] += net
                incoming[target] += net
            else:
                outgoing[target] += -net
                incoming[source] += -net
        for node in set(outgoing) | set(incoming):
            scores[node] = max(0.0, outgoing.get(node, 0.0) - incoming.get(node, 0.0))
    elif mode == "te_sink":
        outgoing = defaultdict(float)
        incoming = defaultdict(float)
        for (source, target), value in mapping.items():
            net = float(value)
            if net >= 0.0:
                outgoing[source] += net
                incoming[target] += net
            else:
                outgoing[target] += -net
                incoming[source] += -net
        for node in set(outgoing) | set(incoming):
            scores[node] = max(0.0, incoming.get(node, 0.0) - outgoing.get(node, 0.0))
    else:
        return []
    ranking = sorted(
        ((node, value) for node, value in scores.items() if np.isfinite(value) and value > 0.0),
        key=lambda item: item[1],
        reverse=True,
    )[:topn]
    return [
        (simplify_node_label(node, bundle.nodes.get(node, {}).get("Label")), value)
        for node, value in ranking
    ]


def network_anchor_entries(system: AggregatedSystem, metric: str, topn: int) -> list[tuple[str, float]]:
    bundle = system.network_anchor_bundle
    if bundle is not None:
        if metric == "gc_anchor":
            mapping = bundle.gc_anchor_scores
        elif metric == "te_source":
            mapping = bundle.te_source_scores
        elif metric == "te_sink":
            mapping = bundle.te_sink_scores
        else:
            mapping = {}
        if mapping:
            ranking = sorted(
                ((node, value) for node, value in mapping.items() if np.isfinite(value) and value > 0.0),
                key=lambda item: item[1],
                reverse=True,
            )[:topn]
            return [
                (bundle.labels.get(node, simplify_node_label(node, None)), float(value))
                for node, value in ranking
            ]

    fallback_specs = {
        "gc_anchor": ("gc_edges", "gc_strength"),
        "te_source": ("te_edges", "te_source"),
        "te_sink": ("te_edges", "te_sink"),
    }
    selector_mode = fallback_specs.get(metric)
    if selector_mode is None:
        return []
    selector, mode = selector_mode
    return network_node_anchor_scores(system.network_bundle, selector, mode, topn)


ANCHOR_LABEL_RE = re.compile(r"^(?P<name>[A-Za-z]{1,4})(?P<idx>\d+)(?:[_:\-].*)?$")
RNA_ANCHOR_NAMES = {"A", "C", "G", "U", "DA", "DC", "DG", "DT", "RA", "RC", "RG", "RU"}


def parse_anchor_label(label: str) -> tuple[str, int | None, str]:
    compact = label.strip()
    match = ANCHOR_LABEL_RE.match(compact)
    if not match:
        return "other", None, compact
    name = match.group("name").upper()
    try:
        index = int(match.group("idx"))
    except ValueError:
        index = None
    molecule = "rna" if name in RNA_ANCHOR_NAMES else "protein"
    return molecule, index, compact


def collect_network_anchor_points(
    systems: dict[str, AggregatedSystem],
    metric: str,
    molecule: str,
    topn: int,
) -> dict[str, list[dict[str, float | int | str]]]:
    collected: dict[str, list[dict[str, float | int | str]]] = {system: [] for system in SYSTEM_ORDER}
    for system in SYSTEM_ORDER:
        entries = network_anchor_entries(systems[system], metric, topn)
        max_score = max((score for _, score in entries), default=0.0)
        if max_score <= 0.0:
            continue
        for rank, (label, score) in enumerate(entries, start=1):
            parsed_molecule, index, clean_label = parse_anchor_label(label)
            if parsed_molecule != molecule or index is None:
                continue
            collected[system].append(
                {
                    "label": clean_label,
                    "index": int(index),
                    "score": float(score),
                    "relative_score": 100.0 * float(score) / max_score,
                    "rank": int(rank),
                }
            )
    return collected


def top_edges(bundle: NetworkBundle | None, selector: str, topn: int) -> list[tuple[str, float]]:
    if bundle is None:
        return []
    mapping = getattr(bundle, selector)
    if not mapping:
        return []
    ranking = sorted(mapping.items(), key=lambda item: item[1], reverse=True)[:topn]
    return [(pair_label(pair, bundle.nodes), value) for pair, value in ranking]


def plot_ranked_horizontal(ax: plt.Axes, entries: Sequence[tuple[str, float]], title: str, color: str) -> None:
    if not entries:
        draw_missing_axis(ax, f"{title}\nNot available")
        return
    labels = [label for label, _ in entries][::-1]
    values = [value for _, value in entries][::-1]
    ypos = np.arange(len(entries))
    ax.barh(ypos, values, color=color)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels)
    ax.set_title(title)
    ax.grid(True, axis="x")


def make_fig1_timeseries(systems: dict[str, AggregatedSystem], outdir: Path, dpi: int) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4))
    fig.subplots_adjust(wspace=0.25, hspace=0.34)
    panel_axes = axes.ravel()
    plot_metric_compare(panel_axes[0], systems, "dbscan_timeseries", "Time_ps", "NumCondensates", "Condensates vs time", "Count", x_scale=1000.0)
    plot_metric_compare(panel_axes[1], systems, "dbscan_timeseries", "Time_ps", "GiantCluster_ProChains", "Protein chains in giant condensate", "Chains", x_scale=1000.0)
    plot_metric_compare(panel_axes[2], systems, "dbscan_timeseries", "Time_ps", "GiantCluster_RnaChains", "RNA chains in giant condensate", "Chains", x_scale=1000.0)

    ax = panel_axes[3]
    plotted = False
    for system in SYSTEM_ORDER:
        table = systems[system].dbscan_timeseries
        if table is None:
            continue
        x = np.asarray(table.columns["Time_ps"], dtype=float) / 1000.0
        if "Rg_Pro_Dense" in table.columns:
            ax.plot(
                x,
                table.columns["Rg_Pro_Dense"],
                color=dense_rg_color(system, "protein"),
                linewidth=1.7,
                label=f"{SYSTEM_DISPLAY[system]} protein",
            )
            plotted = True
        if "Rg_RNA_Dense" in table.columns:
            ax.plot(
                x,
                table.columns["Rg_RNA_Dense"],
                color=dense_rg_color(system, "rna"),
                linewidth=1.7,
                label=f"{SYSTEM_DISPLAY[system]} RNA",
            )
            plotted = True
    if plotted:
        ax.set_title("Dense-phase radius of gyration")
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Rg (nm)")
        ax.grid(True, axis="y")
        legend = ax.legend(
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            ncol=1,
            frameon=True,
            framealpha=0.95,
            borderpad=0.35,
            labelspacing=0.35,
            handlelength=2.0,
        )
        legend.get_frame().set_facecolor("white")
        legend.get_frame().set_edgecolor("#cccccc")
        steady = max(systems[system].steady_start_ns for system in SYSTEM_ORDER)
        add_steady_window(ax, steady, "ns")
    else:
        draw_missing_axis(ax, "Dense-phase radius of gyration\nNot available")

    for idx, axis in enumerate(panel_axes):
        add_panel_label(axis, PANEL_LABELS[idx])
    fig.suptitle("Figure 1. LLPS kinetics: phos vs unphos", fontsize=10, fontweight="bold", y=1.02)
    return save_figure_all_formats(fig, outdir / "fig1_llps_timeseries_compare", dpi)


def make_fig2_density_polymer(systems: dict[str, AggregatedSystem], outdir: Path, dpi: int) -> list[str]:
    fig = plt.figure(figsize=(7.2, 5.9))
    gs = fig.add_gridspec(
        2,
        3,
        width_ratios=[1.0, 1.0, 0.055],
        height_ratios=[1.0, 1.0],
        left=0.10,
        right=0.90,
        bottom=0.09,
        top=0.88,
        hspace=0.42,
        wspace=0.16,
    )
    axes = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 0]),
        fig.add_subplot(gs[1, 1]),
    ]
    cax_protein = fig.add_subplot(gs[0, 2])
    cax_rna = fig.add_subplot(gs[1, 2])

    protein_vmax = max(
        [float(np.max(build_density_projection(system.density_protein))) for system in systems.values() if system.density_protein is not None]
        or [1.0]
    )
    rna_vmax = max(
        [float(np.max(build_density_projection(system.density_rna))) for system in systems.values() if system.density_rna is not None]
        or [1.0]
    )

    protein_images = []
    rna_images = []
    protein_images.append(
        plot_density_panel(
            axes[0],
            systems["phos"].density_protein,
            "Protein, phosphorylated",
            plt.get_cmap("magma"),
            protein_vmax,
        )
    )
    protein_images.append(
        plot_density_panel(
            axes[1],
            systems["unphos"].density_protein,
            "Protein, unphosphorylated (WT)",
            plt.get_cmap("magma"),
            protein_vmax,
        )
    )
    rna_images.append(
        plot_density_panel(
            axes[2],
            systems["phos"].density_rna,
            "RNA, phosphorylated",
            plt.get_cmap("viridis"),
            rna_vmax,
        )
    )
    rna_images.append(
        plot_density_panel(
            axes[3],
            systems["unphos"].density_rna,
            "RNA, unphosphorylated (WT)",
            plt.get_cmap("viridis"),
            rna_vmax,
        )
    )

    valid_protein_image = next((img for img in protein_images if img is not None), None)
    valid_rna_image = next((img for img in rna_images if img is not None), None)
    if valid_protein_image is not None:
        fig.colorbar(
            valid_protein_image,
            cax=cax_protein,
            label=r"Projected bead density (nm$^{-2}$)",
        )
    else:
        cax_protein.axis("off")
    if valid_rna_image is not None:
        fig.colorbar(
            valid_rna_image,
            cax=cax_rna,
            label=r"Projected bead density (nm$^{-2}$)",
        )
    else:
        cax_rna.axis("off")

    extents = [
        grid.extent_xy_nm
        for system in systems.values()
        for grid in (system.density_protein, system.density_rna)
        if grid is not None
    ]
    if extents:
        max_x = max(max(abs(extent[0]), abs(extent[1])) for extent in extents)
        max_y = max(max(abs(extent[2]), abs(extent[3])) for extent in extents)
        for ax in axes:
            ax.set_xlim(-max_x, max_x)
            ax.set_ylim(-max_y, max_y)

    for idx, axis in enumerate(axes):
        add_panel_label_in_margin(fig, axis, PANEL_LABELS[idx], fontsize=13.0)
    fig.suptitle(
        "Figure 2. XY projections of protein/RNA bead density integrated along z",
        fontsize=10,
        fontweight="bold",
        y=0.965,
    )
    return save_figure_all_formats(fig, outdir / "fig2_density_polymer_compare", dpi)


def radial_density_table(system: AggregatedSystem, species: str) -> Table | None:
    grid = system.density_protein if species == "protein" else system.density_rna
    if grid is None:
        return None
    table = grid.radial_profile_data
    if table is None or not table.has_columns(["r_nm", "density_nm3"]):
        return None
    return table


def make_radial_density_figure(
    systems: dict[str, AggregatedSystem],
    outdir: Path,
    dpi: int,
    species: str,
    figure_number: int,
    smooth_window: int,
) -> list[str]:
    species_label = "Protein" if species == "protein" else "RNA"
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    fig.subplots_adjust(top=0.78, bottom=0.18, left=0.14, right=0.96)
    plotted = False
    max_y = 0.0
    for system in RADIAL_SYSTEM_ORDER:
        table = radial_density_table(systems[system], species)
        if table is None:
            continue
        r_nm = np.asarray(table.columns["r_nm"], dtype=float)
        density = np.asarray(table.columns["density_nm3"], dtype=float)
        mask = np.isfinite(r_nm) & np.isfinite(density)
        if not np.any(mask):
            continue
        r_nm = r_nm[mask]
        density = density[mask]
        order = np.argsort(r_nm)
        r_nm = r_nm[order]
        density = density[order]
        density_plot = smooth_profile(density, smooth_window)
        ax.plot(r_nm, density_plot, color=system_color(system), linewidth=1.8, label=SYSTEM_DISPLAY[system])
        max_y = max(max_y, float(np.max(density_plot)) if density_plot.size else 0.0)
        plotted = True

    if plotted:
        ax.axhline(0.0, color="#999999", linewidth=0.7, alpha=0.6)
        ax.set_title(f"{species_label} radial density profile")
        ax.set_xlabel(r"Distance from condensate center, $r$ (nm)")
        ax.set_ylabel(rf"{species_label} density (nm$^{{-3}}$)")
        ax.set_xlim(left=0.0)
        ax.set_ylim(bottom=0.0, top=max(max_y * 1.15, 1e-6))
        ax.grid(True, axis="y")
        ax.legend(loc="upper right", frameon=False)
    else:
        draw_missing_axis(ax, f"{species_label} radial density profile\nDensity DX not available")

    fig.suptitle(f"Figure {figure_number}. {species_label} radial density: phosphorylated vs WT", fontsize=10, fontweight="bold", y=0.995)
    suffix = "protein" if species == "protein" else "rna"
    return save_figure_all_formats(fig, outdir / f"fig{figure_number}_radial_density_{suffix}_compare", dpi)


def plot_shape_metric_bar(ax: plt.Axes, systems: dict[str, AggregatedSystem], metric_key: str, title: str, ylabel: str) -> None:
    species_order = ("protein", "rna")
    x = np.arange(len(species_order), dtype=float)
    width = 0.34
    plotted = False
    for idx, system_key in enumerate(SYSTEM_ORDER):
        values: list[float] = []
        for species in species_order:
            value = systems[system_key].density_shape_metrics.get(species, {}).get(metric_key, float("nan"))
            values.append(float(value) if np.isfinite(value) else float("nan"))
        if any(np.isfinite(values)):
            ax.bar(x + (idx - 0.5) * width, values, width=width, color=system_color(system_key), label=SYSTEM_DISPLAY[system_key])
            plotted = True
    if plotted:
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xticks(x)
        ax.set_xticklabels(["Protein", "RNA"])
        ax.grid(True, axis="y")
    else:
        draw_missing_axis(ax, f"{title}\nDensity shape metrics not available")


def make_fig10_shape_asymmetry(systems: dict[str, AggregatedSystem], outdir: Path, dpi: int) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.5))
    fig.subplots_adjust(wspace=0.34, hspace=0.42, top=0.84)
    specs = [
        ("normalized_asphericity", "Normalized asphericity", "Asphericity / Rg$^2$"),
        ("relative_shape_anisotropy", "Relative shape anisotropy", r"$\kappa^2$"),
        ("acylindricity_nm2", "Acylindricity", r"nm$^2$"),
        ("aspect_ratio", "Principal-axis aspect ratio", r"$\sqrt{\lambda_1/\lambda_3}$"),
    ]
    for idx, (metric_key, title, ylabel) in enumerate(specs):
        axis = axes.ravel()[idx]
        plot_shape_metric_bar(axis, systems, metric_key, title, ylabel)
        add_panel_label(axis, PANEL_LABELS[idx])
    add_system_figure_legend(fig, y=0.94)
    fig.suptitle("Figure 10. Condensate shape asymmetry from density tensors", fontsize=10, fontweight="bold", y=0.995)
    return save_figure_all_formats(fig, outdir / "fig10_density_shape_asymmetry_compare", dpi)


def plot_dense_msd_panel(ax: plt.Axes, systems: dict[str, AggregatedSystem], species: str) -> None:
    attr_name = "dense_msd_protein" if species == "protein" else "dense_msd_rna"
    species_label = "Protein" if species == "protein" else "RNA"
    plotted = False
    diffusion_notes: list[str] = []
    for system in SYSTEM_ORDER:
        bundle = getattr(systems[system], attr_name)
        if bundle is None or bundle.table is None or not bundle.table.has_columns(["lag_ns", "msd_nm2", "count"]):
            continue
        lag = np.asarray(bundle.table.columns["lag_ns"], dtype=float)
        msd = np.asarray(bundle.table.columns["msd_nm2"], dtype=float)
        count = np.asarray(bundle.table.columns["count"], dtype=float)
        mask = np.isfinite(lag) & np.isfinite(msd) & np.isfinite(count) & (lag > 0.0) & (count > 0.0)
        if np.count_nonzero(mask) < 2:
            continue
        diffusion = float(bundle.summary.get("diffusion_nm2_per_ns", float("nan")))
        label = SYSTEM_DISPLAY[system]
        if np.isfinite(diffusion):
            diffusion_notes.append(f"{SYSTEM_DISPLAY[system]}: D={diffusion:.3g} nm$^2$/ns")
        ax.plot(lag[mask], msd[mask], color=system_color(system), linewidth=1.7, label=label)
        slope = float(bundle.summary.get("fit_slope_nm2_per_ns", float("nan")))
        intercept = float(bundle.summary.get("fit_intercept_nm2", float("nan")))
        fit_start = float(bundle.summary.get("fit_start_ns", float("nan")))
        fit_end = float(bundle.summary.get("fit_end_ns", float("nan")))
        if np.isfinite(slope) and np.isfinite(intercept) and np.isfinite(fit_start) and np.isfinite(fit_end) and fit_end > fit_start:
            fit_x = np.asarray([fit_start, fit_end], dtype=float)
            fit_y = slope * fit_x + intercept
            ax.plot(fit_x, fit_y, color=system_color(system), linewidth=1.0, linestyle="--", alpha=0.85)
        plotted = True
    if plotted:
        ax.set_title(f"{species_label} dense-phase MSD")
        ax.set_xlabel("Lag time (ns)")
        ax.set_ylabel(r"MSD (nm$^2$)")
        ax.set_xlim(left=0.0)
        ax.set_ylim(bottom=0.0)
        if diffusion_notes:
            ax.text(
                0.03,
                0.97,
                "\n".join(diffusion_notes),
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=7.2,
                bbox={
                    "boxstyle": "round,pad=0.25",
                    "facecolor": "white",
                    "edgecolor": "#dddddd",
                    "alpha": 0.92,
                },
            )
    else:
        draw_missing_axis(ax, f"{species_label} dense-phase MSD\nNot available")


def make_fig11_dense_msd(systems: dict[str, AggregatedSystem], outdir: Path, dpi: int) -> list[str]:
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))
    fig.subplots_adjust(wspace=0.32, top=0.74)
    plot_dense_msd_panel(axes[0], systems, "protein")
    plot_dense_msd_panel(axes[1], systems, "rna")
    for idx, axis in enumerate(axes):
        add_panel_label(axis, PANEL_LABELS[idx])
    add_system_figure_legend(fig, y=0.89)
    fig.suptitle("Figure 11. Dense-phase MSD and diffusion", fontsize=10, fontweight="bold", y=0.995)
    return save_figure_all_formats(fig, outdir / "fig11_dense_msd_diffusion_compare", dpi)


def make_fig12_rg_rna_dense_histogram(systems: dict[str, AggregatedSystem], outdir: Path, dpi: int) -> list[str]:
    fig, ax = plt.subplots(figsize=(7.2, 4.1))
    fig.subplots_adjust(top=0.78, bottom=0.18, left=0.14, right=0.96)
    table = rg_histogram_table(systems)
    if table is None:
        draw_missing_axis(ax, "Dense-phase RNA Rg distribution\nNot available")
    else:
        centers = np.asarray(table.columns["bin_center_nm"], dtype=float)
        plotted = False
        for system in SYSTEM_ORDER:
            density = np.asarray(table.columns[f"{system}_probability_density"], dtype=float)
            values = systems[system].rg_rna_dense_values
            clean = np.asarray([], dtype=float) if values is None else np.asarray(values, dtype=float)
            clean = clean[np.isfinite(clean) & (clean > 0.0)]
            if clean.size == 0 or not np.any(np.isfinite(density)):
                continue
            mean = float(np.mean(clean))
            median = float(np.median(clean))
            label = f"{SYSTEM_DISPLAY[system]} median={median:.2f} nm"
            ax.step(centers, density, where="mid", color=system_color(system), linewidth=1.8, label=label)
            ax.fill_between(centers, density, step="mid", color=system_color(system), alpha=0.12)
            ax.axvline(mean, color=system_color(system), linestyle="--", linewidth=1.0, alpha=0.8)
            plotted = True
        if plotted:
            ax.set_title("Distribution of dense-phase RNA radius of gyration")
            ax.set_xlabel("Dense-phase RNA Rg (nm)")
            ax.set_ylabel("Probability density")
            ax.grid(True, axis="y")
            legend = ax.legend(loc="best", frameon=True, framealpha=0.94)
            legend.get_frame().set_facecolor("white")
            legend.get_frame().set_edgecolor("#cccccc")
        else:
            draw_missing_axis(ax, "Dense-phase RNA Rg distribution\nNot available")
    add_panel_label(ax, "A")
    fig.suptitle("Figure 12. Dense-phase RNA compaction distribution", fontsize=10, fontweight="bold", y=0.985)
    return save_figure_all_formats(fig, outdir / "fig12_rg_rna_dense_histogram_compare", dpi)


def plot_rna_binned_metric(
    ax: plt.Axes,
    systems: dict[str, AggregatedSystem],
    mean_column: str,
    std_column: str | None,
    title: str,
    ylabel: str,
) -> None:
    plotted = False
    for system in SYSTEM_ORDER:
        table = systems[system].rna_compaction_binned
        if table is None or mean_column not in table.columns:
            continue
        x = np.asarray(table.columns["BoundProteinChainsBin"], dtype=float)
        y = np.asarray(table.columns[mean_column], dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        if not np.any(mask):
            continue
        order = np.argsort(x[mask])
        x_plot = x[mask][order]
        y_plot = y[mask][order]
        ax.plot(
            x_plot,
            y_plot,
            marker="o",
            markersize=3.8,
            linewidth=1.7,
            color=system_color(system),
            label=SYSTEM_DISPLAY[system],
        )
        if std_column and std_column in table.columns:
            yerr = np.asarray(table.columns[std_column], dtype=float)[mask][order]
            valid_err = np.isfinite(yerr)
            if np.any(valid_err):
                lower = y_plot - np.where(valid_err, yerr, 0.0)
                upper = y_plot + np.where(valid_err, yerr, 0.0)
                ax.fill_between(x_plot, lower, upper, color=system_color(system), alpha=0.13, linewidth=0)
        plotted = True
    if not plotted:
        draw_missing_axis(ax, f"{title}\nNot available")
        return
    ax.set_title(title)
    ax.set_xlabel("Protein chains bound to one RNA chain")
    ax.set_ylabel(ylabel)
    ax.set_xlim(left=-0.25)
    ax.grid(True, axis="both")


def make_fig13_rna_binding_compaction_binned(
    systems: dict[str, AggregatedSystem], outdir: Path, dpi: int
) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.8))
    fig.subplots_adjust(hspace=0.42, wspace=0.34, top=0.86)
    panels = [
        ("RgMean_nm", "RgStd_nm", "RNA radius of gyration", "Rg (nm)"),
        ("EndToEndMean_nm", "EndToEndStd_nm", "RNA end-to-end distance", "End-to-end (nm)"),
        ("AsphericityMean", "AsphericityStd", "RNA asphericity", r"Asphericity (nm$^2$)"),
        (
            "LongRangeContactFractionMean",
            "LongRangeContactFractionStd",
            "Long-range intrachain contacts",
            "Long-range contact fraction",
        ),
    ]
    for idx, (mean_col, std_col, title, ylabel) in enumerate(panels):
        axis = axes.ravel()[idx]
        plot_rna_binned_metric(axis, systems, mean_col, std_col, title, ylabel)
        add_panel_label(axis, PANEL_LABELS[idx])
    add_system_figure_legend(fig, y=0.94)
    fig.suptitle(
        "Figure 13. RNA compaction binned by protein binding",
        fontsize=10,
        fontweight="bold",
        y=0.995,
    )
    return save_figure_all_formats(fig, outdir / "fig13_rna_binding_compaction_binned_compare", dpi)


def plot_rna_chain_metric(
    ax: plt.Axes,
    systems: dict[str, AggregatedSystem],
    column: str,
    title: str,
    ylabel: str,
    source_attr: str = "rna_chain_compaction_summary",
) -> None:
    plotted = False
    for system in SYSTEM_ORDER:
        table = getattr(systems[system], source_attr)
        if table is None or column not in table.columns:
            continue
        x = np.asarray(table.columns["RNAChainID"], dtype=float)
        y = np.asarray(table.columns[column], dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        if not np.any(mask):
            continue
        order = np.argsort(x[mask])
        ax.plot(
            x[mask][order],
            y[mask][order],
            marker="o",
            markersize=3.2,
            linewidth=1.35,
            color=system_color(system),
            label=SYSTEM_DISPLAY[system],
        )
        plotted = True
    if not plotted:
        draw_missing_axis(ax, f"{title}\nNot available")
        return
    ax.set_title(title)
    ax.set_xlabel("RNA chain ID")
    ax.set_ylabel(ylabel)
    ax.grid(True, axis="y")


def make_fig14_rna_chain_compaction(
    systems: dict[str, AggregatedSystem], outdir: Path, dpi: int
) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.8), sharex=True)
    fig.subplots_adjust(hspace=0.34, wspace=0.34, top=0.86)
    panels = [
        ("MeanRg_nm", "Mean RNA Rg per chain", "Mean Rg (nm)", "rna_chain_compaction_summary"),
        ("MeanBoundProteinChains", "Protein chains bound per RNA", "Mean bound protein chains", "rna_chain_compaction_summary"),
        ("MeanRelativeShapeAnisotropy", "RNA chain shape anisotropy", "Relative shape anisotropy", "rna_chain_compaction_summary"),
        (
            "MeanLongRangeContactFraction",
            "RNA long-range intrachain contacts",
            "Long-range contact fraction",
            "rna_intrachain_contacts_summary",
        ),
    ]
    for idx, (column, title, ylabel, attr_name) in enumerate(panels):
        axis = axes.ravel()[idx]
        plot_rna_chain_metric(axis, systems, column, title, ylabel, attr_name)
        add_panel_label(axis, PANEL_LABELS[idx])
    add_system_figure_legend(fig, y=0.94)
    fig.suptitle(
        "Figure 14. Per-RNA-chain compaction and protein binding",
        fontsize=10,
        fontweight="bold",
        y=0.995,
    )
    return save_figure_all_formats(fig, outdir / "fig14_rna_chain_compaction_compare", dpi)


def sampled_trace_xy(
    table: Table,
    x_column: str,
    y_column: str,
    max_points: int = 60000,
) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(table.columns[x_column], dtype=float)
    y = np.asarray(table.columns[y_column], dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if x.size > max_points:
        idx = np.linspace(0, x.size - 1, max_points, dtype=int)
        x = x[idx]
        y = y[idx]
    return x, y


def plot_rna_trace_scatter(
    ax: plt.Axes,
    systems: dict[str, AggregatedSystem],
    x_column: str,
    y_column: str,
    title: str,
    xlabel: str,
    ylabel: str,
) -> None:
    plotted = False
    for system in SYSTEM_ORDER:
        table = systems[system].rna_compaction_trace
        if table is None or x_column not in table.columns or y_column not in table.columns:
            continue
        x, y = sampled_trace_xy(table, x_column, y_column)
        if x.size == 0:
            continue
        ax.scatter(
            x,
            y,
            s=5,
            alpha=0.16,
            color=system_color(system),
            edgecolors="none",
            rasterized=True,
            label=SYSTEM_DISPLAY[system],
        )
        plotted = True
    if not plotted:
        draw_missing_axis(ax, f"{title}\nDetailed trace not available")
        return
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, axis="both")


def make_fig15_rna_compaction_trace(
    systems: dict[str, AggregatedSystem], outdir: Path, dpi: int
) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.8))
    fig.subplots_adjust(hspace=0.42, wspace=0.34, top=0.86)
    panels = [
        ("BoundProteinChains", "Rg_nm", "Binding vs RNA Rg", "Bound protein chains", "RNA Rg (nm)"),
        (
            "ProteinRNABeadContacts",
            "Rg_nm",
            "Bead contacts vs RNA Rg",
            "Protein-RNA bead contacts",
            "RNA Rg (nm)",
        ),
        (
            "BoundProteinChains",
            "RelativeShapeAnisotropy",
            "Binding vs RNA anisotropy",
            "Bound protein chains",
            "Relative shape anisotropy",
        ),
        (
            "BoundProteinChains",
            "LongRangeContactFraction",
            "Binding vs RNA intrachain contacts",
            "Bound protein chains",
            "Long-range contact fraction",
        ),
    ]
    for idx, (x_col, y_col, title, xlabel, ylabel) in enumerate(panels):
        axis = axes.ravel()[idx]
        plot_rna_trace_scatter(axis, systems, x_col, y_col, title, xlabel, ylabel)
        add_panel_label(axis, PANEL_LABELS[idx])
    add_system_figure_legend(fig, y=0.94)
    fig.suptitle(
        "Figure 15. Frame-level RNA compaction vs protein binding",
        fontsize=10,
        fontweight="bold",
        y=0.995,
    )
    return save_figure_all_formats(fig, outdir / "fig15_rna_compaction_binding_trace_compare", dpi)


def plot_rna_energy_binned_metric(
    ax: plt.Axes,
    systems: dict[str, AggregatedSystem],
    mean_column: str,
    std_column: str | None,
    title: str,
    ylabel: str,
) -> None:
    plotted = False
    for system in SYSTEM_ORDER:
        table = systems[system].rna_energy_compaction_binned
        if table is None or mean_column not in table.columns:
            continue
        x = np.asarray(table.columns["BoundProteinChainsBin"], dtype=float)
        y = np.asarray(table.columns[mean_column], dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        if not np.any(mask):
            continue
        order = np.argsort(x[mask])
        x_plot = x[mask][order]
        y_plot = y[mask][order]
        ax.plot(
            x_plot,
            y_plot,
            marker="o",
            markersize=3.8,
            linewidth=1.7,
            color=system_color(system),
        )
        if std_column and std_column in table.columns:
            yerr = np.asarray(table.columns[std_column], dtype=float)[mask][order]
            valid_err = np.isfinite(yerr)
            if np.any(valid_err):
                lower = y_plot - np.where(valid_err, yerr, 0.0)
                upper = y_plot + np.where(valid_err, yerr, 0.0)
                ax.fill_between(x_plot, lower, upper, color=system_color(system), alpha=0.13, linewidth=0)
        plotted = True
    if not plotted:
        draw_missing_axis(ax, f"{title}\nNot available")
        return
    ax.set_title(title)
    ax.set_xlabel("Protein chains bound to one RNA chain")
    ax.set_ylabel(ylabel)
    ax.set_xlim(left=-0.25)


def trace_metric_from_table(table: Table, metric_name: str) -> np.ndarray:
    if metric_name == "AttractionStrength_kJmol":
        values = -np.asarray(table.columns["Etotal_ProteinRNA_kJmol"], dtype=float)
    elif metric_name == "AttractionPerBoundChain_kJmol":
        total = -np.asarray(table.columns["Etotal_ProteinRNA_kJmol"], dtype=float)
        denom = np.asarray(table.columns["BoundProteinChains"], dtype=float)
        values = np.full_like(total, np.nan, dtype=float)
        mask = np.isfinite(total) & np.isfinite(denom) & (denom > 0.0)
        values[mask] = total[mask] / denom[mask]
    elif metric_name == "AttractionPerBeadContact_kJmol":
        total = -np.asarray(table.columns["Etotal_ProteinRNA_kJmol"], dtype=float)
        denom = np.asarray(table.columns["ProteinRNABeadContacts"], dtype=float)
        values = np.full_like(total, np.nan, dtype=float)
        mask = np.isfinite(total) & np.isfinite(denom) & (denom > 0.0)
        values[mask] = total[mask] / denom[mask]
    else:
        values = np.asarray(table.columns[metric_name], dtype=float)
    return values.astype(float, copy=False)


def binding_energy_trace_xy(
    table: Table,
    y_column: str,
    *,
    x_metric: str = "AttractionStrength_kJmol",
    inside_dense_only: bool = True,
    contacting_only: bool = True,
    max_points: int = 45000,
) -> tuple[np.ndarray, np.ndarray]:
    x = trace_metric_from_table(table, x_metric)
    y = trace_metric_from_table(table, y_column)
    mask = np.isfinite(x) & np.isfinite(y)
    if inside_dense_only and "InsideDense" in table.columns:
        mask &= np.asarray(table.columns["InsideDense"], dtype=float) > 0.5
    if contacting_only and "BoundProteinChains" in table.columns:
        mask &= np.asarray(table.columns["BoundProteinChains"], dtype=float) > 0.0
    x = x[mask]
    y = y[mask]
    if x.size > max_points:
        idx = np.linspace(0, x.size - 1, max_points, dtype=int)
        x = x[idx]
        y = y[idx]
    return x, y


def binned_trend_line(x: np.ndarray, y: np.ndarray, bins: int = 18) -> tuple[np.ndarray, np.ndarray]:
    if x.size < 8 or y.size < 8:
        return np.asarray([], dtype=float), np.asarray([], dtype=float)
    xmin = float(np.min(x))
    xmax = float(np.max(x))
    if not np.isfinite(xmin) or not np.isfinite(xmax) or xmax <= xmin:
        return np.asarray([], dtype=float), np.asarray([], dtype=float)
    edges = np.linspace(xmin, xmax, max(4, int(bins)) + 1)
    centers: list[float] = []
    medians: list[float] = []
    for idx in range(edges.size - 1):
        lo = edges[idx]
        hi = edges[idx + 1]
        if idx == edges.size - 2:
            mask = (x >= lo) & (x <= hi)
        else:
            mask = (x >= lo) & (x < hi)
        if np.count_nonzero(mask) < 5:
            continue
        centers.append(0.5 * (lo + hi))
        medians.append(float(np.median(y[mask])))
    return np.asarray(centers, dtype=float), np.asarray(medians, dtype=float)


def plot_energy_trace_relationship(
    ax: plt.Axes,
    systems: dict[str, AggregatedSystem],
    y_column: str,
    title: str,
    ylabel: str,
    *,
    x_metric: str = "AttractionStrength_kJmol",
    xlabel: str = r"Protein-RNA attraction strength, $-E_{\mathrm{total}}$ (kJ/mol)",
    inside_dense_only: bool = True,
    contacting_only: bool = True,
) -> None:
    plotted = False
    derived_x_metrics = {"AttractionStrength_kJmol", "AttractionPerBoundChain_kJmol", "AttractionPerBeadContact_kJmol"}
    for system in SYSTEM_ORDER:
        table = systems[system].rna_binding_energy_trace
        if table is None:
            continue
        needed = ["Etotal_ProteinRNA_kJmol", y_column]
        if x_metric not in derived_x_metrics:
            needed.append(x_metric)
        if not table.has_columns(needed):
            continue
        x, y = binding_energy_trace_xy(
            table,
            y_column,
            x_metric=x_metric,
            inside_dense_only=inside_dense_only,
            contacting_only=contacting_only,
        )
        if x.size == 0:
            continue
        ax.scatter(
            x,
            y,
            s=5,
            alpha=0.12,
            color=system_color(system),
            edgecolors="none",
            rasterized=True,
        )
        trend_x, trend_y = binned_trend_line(x, y)
        if trend_x.size:
            ax.plot(trend_x, trend_y, color=system_color(system), linewidth=2.0)
        plotted = True
    if not plotted:
        draw_missing_axis(ax, f"{title}\nEnergy-trace data not available")
        return
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(left=0.0 if x_metric == "AttractionStrength_kJmol" else None)


def make_fig16_rna_energy_compaction(
    systems: dict[str, AggregatedSystem], outdir: Path, dpi: int
) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.9))
    fig.subplots_adjust(hspace=0.44, wspace=0.36, top=0.85, bottom=0.12)
    panels = [
        ("EelecMean_kJmol", "EelecStd_kJmol", "Protein-RNA electrostatics vs binding", r"$E_{\mathrm{elec}}$ (kJ/mol)"),
        (
            "AttractionStrengthMean_kJmol",
            "AttractionStrengthStd_kJmol",
            "Protein-RNA attraction vs binding",
            r"Attraction strength, $-E_{\mathrm{total}}$ (kJ/mol)",
        ),
    ]
    for idx, (mean_col, std_col, title, ylabel) in enumerate(panels):
        axis = axes.ravel()[idx]
        plot_rna_energy_binned_metric(axis, systems, mean_col, std_col, title, ylabel)
        add_panel_label(axis, PANEL_LABELS[idx])

    plot_energy_trace_relationship(
        axes[1, 0],
        systems,
        "Rg_nm",
        "RNA size vs protein-RNA attraction",
        "RNA Rg (nm)",
    )
    add_panel_label(axes[1, 0], PANEL_LABELS[2])
    plot_energy_trace_relationship(
        axes[1, 1],
        systems,
        "LongRangeContactFraction",
        "RNA tertiary compaction vs protein-RNA attraction",
        "Long-range contact fraction",
    )
    add_panel_label(axes[1, 1], PANEL_LABELS[3])
    add_system_figure_legend(fig, y=0.935)
    fig.text(
        0.5,
        0.04,
        "Panels C-D use contacting RNA-chain frames inside the giant condensate; solid lines show binned medians.",
        ha="center",
        va="center",
        fontsize=7.3,
        color="#555555",
    )
    fig.suptitle(
        "Figure 16. Protein-RNA interaction energy versus RNA compaction",
        fontsize=10,
        fontweight="bold",
        y=0.992,
    )
    return save_figure_all_formats(fig, outdir / "fig16_rna_energy_compaction_compare", dpi)


def mean_trace_metric(
    system: AggregatedSystem,
    metric_name: str,
    *,
    inside_dense_only: bool = True,
    contacting_only: bool = True,
) -> float:
    table = system.rna_binding_energy_trace
    if table is None or "Etotal_ProteinRNA_kJmol" not in table.columns:
        return float("nan")
    values = trace_metric_from_table(table, metric_name)
    mask = np.isfinite(values)
    if inside_dense_only and "InsideDense" in table.columns:
        mask &= np.asarray(table.columns["InsideDense"], dtype=float) > 0.5
    if contacting_only and "BoundProteinChains" in table.columns:
        mask &= np.asarray(table.columns["BoundProteinChains"], dtype=float) > 0.0
    if not np.any(mask):
        return float("nan")
    return float(np.mean(values[mask]))


def mean_dense_rna_rg(system: AggregatedSystem) -> float:
    values = system.rg_rna_dense_values
    if values is None:
        return float("nan")
    clean = np.asarray(values, dtype=float)
    clean = clean[np.isfinite(clean) & (clean > 0.0)]
    if clean.size == 0:
        return float("nan")
    return float(np.mean(clean))


def plot_scalar_comparison_bar(
    ax: plt.Axes,
    values: dict[str, float],
    title: str,
    ylabel: str,
) -> None:
    labels = [SYSTEM_DISPLAY[system] for system in SYSTEM_ORDER]
    arr = np.asarray([float(values.get(system, float("nan"))) for system in SYSTEM_ORDER], dtype=float)
    mask = np.isfinite(arr)
    if not np.any(mask):
        draw_missing_axis(ax, f"{title}\nNot available")
        return
    xpos = np.arange(len(SYSTEM_ORDER), dtype=float)
    ax.bar(xpos, np.where(mask, arr, 0.0), color=[system_color(system) for system in SYSTEM_ORDER], width=0.6)
    ax.set_xticks(xpos)
    ax.set_xticklabels(labels)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    finite = arr[mask]
    ymin = min(0.0, float(np.min(finite)))
    ymax = float(np.max(finite))
    span = ymax - ymin
    if span <= 1e-12:
        span = max(abs(ymax), 1.0)
    ax.set_ylim(ymin - 0.10 * span, ymax + 0.18 * span)
    for idx, value in enumerate(arr):
        if not np.isfinite(value):
            continue
        offset = 0.04 * span if value >= 0.0 else -0.06 * span
        va = "bottom" if value >= 0.0 else "top"
        ax.text(idx, value + offset, f"{value:.2f}", ha="center", va=va, fontsize=8)


def summary_metric_values(
    systems: dict[str, AggregatedSystem],
    source_attr: str,
    metric_key: str,
    *,
    transform: Callable[[float], float] | None = None,
) -> dict[str, float]:
    out: dict[str, float] = {}
    for system in SYSTEM_ORDER:
        source = getattr(systems[system], source_attr)
        if not source:
            out[system] = float("nan")
            continue
        try:
            value = float(source.get(metric_key, float("nan")))
        except AttributeError:
            value = float("nan")
        if transform is not None and np.isfinite(value):
            value = float(transform(value))
        out[system] = value
    return out


def binding_explained_compaction_fraction(system: AggregatedSystem, scope: str = "inside_dense_only") -> float:
    summary = system.rna_information_summary
    if not summary:
        return float("nan")
    try:
        h_compaction = float(summary.get(f"{scope}_H_compaction_bits", float("nan")))
        h_conditional = float(summary.get(f"{scope}_H_compaction_given_binding_bits", float("nan")))
    except (TypeError, ValueError):
        return float("nan")
    if not np.isfinite(h_compaction) or not np.isfinite(h_conditional) or h_compaction <= 0.0:
        return float("nan")
    value = 1.0 - (h_conditional / h_compaction)
    return float(np.clip(value, 0.0, 1.0))


def plot_summary_metric_bar(
    ax: plt.Axes,
    systems: dict[str, AggregatedSystem],
    source_attr: str,
    metric_key: str,
    title: str,
    ylabel: str,
    *,
    transform: Callable[[float], float] | None = None,
) -> None:
    values = summary_metric_values(systems, source_attr, metric_key, transform=transform)
    plot_scalar_comparison_bar(ax, values, title, ylabel)


def plot_summary_metric_grouped_bars(
    ax: plt.Axes,
    systems: dict[str, AggregatedSystem],
    source_attr: str,
    metrics: Sequence[tuple[str, str, Callable[[float], float] | None]],
    title: str,
    ylabel: str,
) -> None:
    x = np.arange(len(metrics), dtype=float)
    width = 0.34
    plotted = False
    for idx, system in enumerate(SYSTEM_ORDER):
        source = getattr(systems[system], source_attr)
        values: list[float] = []
        for key, _, transform in metrics:
            try:
                value = float(source.get(key, float("nan"))) if source else float("nan")
            except AttributeError:
                value = float("nan")
            if transform is not None and np.isfinite(value):
                value = float(transform(value))
            values.append(value)
        if any(np.isfinite(values)):
            offset = (idx - 0.5) * width
            ax.bar(x + offset, values, width=width, color=system_color(system), label=SYSTEM_DISPLAY[system])
            plotted = True
    if not plotted:
        draw_missing_axis(ax, f"{title}\nNot available")
        return
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xticks(x)
    ax.set_xticklabels([label for _, label, _ in metrics], rotation=18, ha="right")
    ymin, ymax = ax.get_ylim()
    lower = min(0.0, ymin)
    upper = max(ymax, 0.05)
    ax.set_ylim(lower, upper * 1.10 if upper > 0 else 1.0)


def make_fig17_mechanistic_summary(
    systems: dict[str, AggregatedSystem], outdir: Path, dpi: int
) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.5))
    fig.subplots_adjust(hspace=0.40, wspace=0.34, top=0.86, bottom=0.10)
    protein_diffusion = {
        system: float(systems[system].dense_msd_protein.summary.get("diffusion_nm2_per_ns", float("nan")))
        if systems[system].dense_msd_protein is not None
        else float("nan")
        for system in SYSTEM_ORDER
    }
    rna_diffusion = {
        system: float(systems[system].dense_msd_rna.summary.get("diffusion_nm2_per_ns", float("nan")))
        if systems[system].dense_msd_rna is not None
        else float("nan")
        for system in SYSTEM_ORDER
    }
    dense_rna_rg = {system: mean_dense_rna_rg(systems[system]) for system in SYSTEM_ORDER}
    attraction = {
        system: mean_trace_metric(
            systems[system],
            "AttractionStrength_kJmol",
            inside_dense_only=True,
            contacting_only=True,
        )
        for system in SYSTEM_ORDER
    }
    panels = [
        (protein_diffusion, "Protein dense-phase diffusion", r"Diffusion (nm$^2$/ns)"),
        (rna_diffusion, "RNA dense-phase diffusion", r"Diffusion (nm$^2$/ns)"),
        (dense_rna_rg, "Mean dense-phase RNA size", "Mean RNA Rg (nm)"),
        (
            attraction,
            "Protein-RNA attraction in dense phase",
            r"Mean $-E_{\mathrm{total}}$ (kJ/mol)",
        ),
    ]
    for idx, (values, title, ylabel) in enumerate(panels):
        axis = axes.ravel()[idx]
        plot_scalar_comparison_bar(axis, values, title, ylabel)
        add_panel_label(axis, PANEL_LABELS[idx])
    fig.text(
        0.5,
        0.045,
        "Panel D averages contacting RNA-chain frames inside the giant condensate only.",
        ha="center",
        va="center",
        fontsize=7.3,
        color="#555555",
    )
    fig.suptitle(
        "Figure 17. Mechanistic summary: weaker attraction, larger RNA, faster condensates",
        fontsize=10,
        fontweight="bold",
        y=0.992,
    )
    return save_figure_all_formats(fig, outdir / "fig17_mechanistic_summary_compare", dpi)


def get_trace_table(system: AggregatedSystem, preferred_attr: str, fallback_attr: str | None = None) -> Table | None:
    table = getattr(system, preferred_attr)
    if table is not None:
        return table
    if fallback_attr is not None:
        return getattr(system, fallback_attr)
    return None


def trace_mask(
    table: Table,
    *,
    inside_dense_only: bool = False,
    contacting_only: bool = False,
) -> np.ndarray:
    first_column = table.order[0]
    mask = np.isfinite(np.asarray(table.columns[first_column], dtype=float))
    if inside_dense_only and "InsideDense" in table.columns:
        mask &= np.asarray(table.columns["InsideDense"], dtype=float) > 0.5
    if contacting_only and "BoundProteinChains" in table.columns:
        mask &= np.asarray(table.columns["BoundProteinChains"], dtype=float) > 0.0
    return mask


def compaction_state_landscape(
    system: AggregatedSystem,
    preferred_scopes: Sequence[str] = ("inside_dense_only", "contacting_only", "all_frames"),
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    table = system.rna_compaction_state_histogram
    if table is None or not table.has_columns(RNA_COMPACTION_STATE_HISTOGRAM_COLUMNS):
        return None
    mask = scope_mask_from_table(table, preferred_scopes)
    if not np.any(mask):
        return None
    rg_bins = np.asarray(table.columns["RgBin"], dtype=float)[mask].astype(int)
    rg_min = np.asarray(table.columns["RgMin_nm"], dtype=float)[mask]
    rg_max = np.asarray(table.columns["RgMax_nm"], dtype=float)[mask]
    long_bins = np.asarray(table.columns["LongRangeBin"], dtype=float)[mask].astype(int)
    long_min = np.asarray(table.columns["LongRangeMin"], dtype=float)[mask]
    long_max = np.asarray(table.columns["LongRangeMax"], dtype=float)[mask]
    probabilities = np.asarray(table.columns["Probability"], dtype=float)[mask]
    if rg_bins.size == 0 or long_bins.size == 0:
        return None
    unique_rg = np.unique(rg_bins)
    unique_long = np.unique(long_bins)
    x_edges = np.asarray(
        [np.min(rg_min[rg_bins == bin_id]) for bin_id in unique_rg]
        + [np.max(rg_max[rg_bins == unique_rg[-1]])],
        dtype=float,
    )
    y_edges = np.asarray(
        [np.min(long_min[long_bins == bin_id]) for bin_id in unique_long]
        + [np.max(long_max[long_bins == unique_long[-1]])],
        dtype=float,
    )
    matrix = np.zeros((unique_long.size, unique_rg.size), dtype=float)
    for rg_bin, long_bin, probability in zip(rg_bins, long_bins, probabilities):
        if not np.isfinite(probability):
            continue
        x_idx = int(np.where(unique_rg == rg_bin)[0][0])
        y_idx = int(np.where(unique_long == long_bin)[0][0])
        matrix[y_idx, x_idx] += float(probability)
    return x_edges, y_edges, matrix


def plot_compaction_state_landscape(
    ax: plt.Axes,
    system: AggregatedSystem,
    *,
    preferred_scopes: Sequence[str] = ("inside_dense_only", "contacting_only", "all_frames"),
    vmax: float | None = None,
    x_limits: tuple[float, float] | None = None,
    y_limits: tuple[float, float] | None = None,
) -> object | None:
    landscape = compaction_state_landscape(system, preferred_scopes=preferred_scopes)
    if landscape is None:
        draw_missing_axis(ax, f"{SYSTEM_DISPLAY[system.system]}\nNo RNA-state histogram")
        return None
    x_edges, y_edges, matrix = landscape
    image = ax.pcolormesh(
        x_edges,
        y_edges,
        matrix,
        shading="auto",
        cmap="YlOrBr",
        vmin=0.0,
        vmax=vmax,
    )
    ax.set_title(SYSTEM_DISPLAY[system.system], pad=7)
    ax.set_xlabel("RNA Rg (nm)")
    ax.set_ylabel("Long-range contact fraction")
    if x_limits is not None:
        ax.set_xlim(*x_limits)
    if y_limits is not None:
        ax.set_ylim(*y_limits)
    return image


def plot_integer_trace_distribution(
    ax: plt.Axes,
    systems: dict[str, AggregatedSystem],
    column: str,
    title: str,
    xlabel: str,
    *,
    preferred_attr: str = "rna_binding_energy_trace",
    fallback_attr: str | None = "rna_compaction_trace",
    inside_dense_only: bool = False,
    contacting_only: bool = False,
) -> None:
    plotted = False
    max_value = 0
    distributions: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for system_key in SYSTEM_ORDER:
        table = get_trace_table(systems[system_key], preferred_attr, fallback_attr)
        if table is None or column not in table.columns:
            continue
        values = np.asarray(table.columns[column], dtype=float)
        mask = trace_mask(table, inside_dense_only=inside_dense_only, contacting_only=contacting_only)
        mask &= np.isfinite(values)
        if not np.any(mask):
            continue
        ints = np.rint(values[mask]).astype(int)
        ints = ints[ints >= 0]
        if ints.size == 0:
            continue
        counts = np.bincount(ints)
        probs = counts.astype(float) / max(float(np.sum(counts)), 1.0)
        x = np.arange(counts.size, dtype=float)
        distributions[system_key] = (x, probs)
        max_value = max(max_value, int(np.max(ints)))
        plotted = True
    if not plotted:
        draw_missing_axis(ax, f"{title}\nNot available")
        return
    for system_key in SYSTEM_ORDER:
        entry = distributions.get(system_key)
        if entry is None:
            continue
        x, probs = entry
        ax.step(x, probs, where="mid", color=system_color(system_key), linewidth=1.9, label=SYSTEM_DISPLAY[system_key])
        ax.fill_between(x, probs, step="mid", color=system_color(system_key), alpha=0.10)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Probability")
    ax.set_xlim(-0.25, max_value + 0.25)


def make_fig18_compaction_efficiency(
    systems: dict[str, AggregatedSystem], outdir: Path, dpi: int
) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 6.0))
    fig.subplots_adjust(hspace=0.50, wspace=0.42, top=0.83, bottom=0.15)

    plot_integer_trace_distribution(
        axes[0, 0],
        systems,
        "BoundProteinChains",
        "RNA binding-state occupancy",
        "Protein chains bound to one RNA chain",
        inside_dense_only=True,
        contacting_only=False,
    )
    add_panel_label_in_margin(fig, axes[0, 0], PANEL_LABELS[0], fontsize=13.0)

    plot_rna_binned_metric(
        axes[0, 1],
        systems,
        "RgMean_nm",
        "RgStd_nm",
        "RNA size at matched valency",
        "RNA Rg (nm)",
    )
    add_panel_label_in_margin(fig, axes[0, 1], PANEL_LABELS[1], fontsize=13.0)

    plot_rna_binned_metric(
        axes[1, 0],
        systems,
        "LongRangeContactFractionMean",
        "LongRangeContactFractionStd",
        "RNA tertiary compaction at matched valency",
        "Long-range contact fraction",
    )
    add_panel_label_in_margin(fig, axes[1, 0], PANEL_LABELS[2], fontsize=13.0)

    plot_energy_trace_relationship(
        axes[1, 1],
        systems,
        "LongRangeContactFraction",
        "Compaction efficiency per bound-chain attraction",
        "Long-range contact fraction",
        x_metric="AttractionPerBoundChain_kJmol",
        xlabel="Attraction per bound chain,\n" + r"$-E_{\mathrm{total}}/N_{\mathrm{bound}}$ (kJ/mol)",
        inside_dense_only=True,
        contacting_only=True,
    )
    add_panel_label_in_margin(fig, axes[1, 1], PANEL_LABELS[3], fontsize=13.0)

    add_system_figure_legend(fig, y=0.92)
    fig.text(
        0.5,
        0.04,
        "Panels B-C compare matched protein valency; panel D tests how efficiently interaction energy is converted into RNA tertiary compaction.",
        ha="center",
        va="center",
        fontsize=7.2,
        color="#555555",
    )
    fig.suptitle(
        "Figure 18. Phosphorylation reduces RNA compaction efficiency",
        fontsize=10,
        fontweight="bold",
        y=0.985,
    )
    return save_figure_all_formats(fig, outdir / "fig18_compaction_efficiency_compare", dpi)


def make_fig19_dynamic_binding_efficiency(
    systems: dict[str, AggregatedSystem], outdir: Path, dpi: int
) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 6.0))
    fig.subplots_adjust(hspace=0.58, wspace=0.46, top=0.83, bottom=0.17)

    plot_summary_metric_bar(
        axes[0, 0],
        systems,
        "liquidity_metrics",
        "mean_protein_rna_contact_lifetime_ns",
        "Protein-RNA contact lifetime",
        "Mean lifetime (ns)",
    )
    add_panel_label_in_margin(fig, axes[0, 0], PANEL_LABELS[0], fontsize=13.0)

    plot_grouped_summary_bars(
        axes[0, 1],
        systems,
        "interface_turnover",
        [
            ("new_contacts_per_ns", "New"),
            ("lost_contacts_per_ns", "Lost"),
            ("net_change_per_ns", "Net"),
        ],
        "Protein-RNA interface turnover",
        "Contacts/ns",
        show_legend=False,
    )
    add_panel_label_in_margin(fig, axes[0, 1], PANEL_LABELS[1], fontsize=13.0)

    plot_summary_metric_grouped_bars(
        axes[1, 0],
        systems,
        "rna_energy_compaction_summary",
        [
            (
                "inside_dense_only_spearman_AttractionStrength_kJmol_vs_LongRangeContactFraction",
                "Attraction -> long-range",
                None,
            ),
            (
                "inside_dense_only_spearman_AttractionStrength_kJmol_vs_Rg_nm",
                "Attraction -> collapse",
                lambda value: -value,
            ),
        ],
        "Attraction-compaction coupling",
        "Coupling strength (Spearman)",
    )
    add_panel_label_in_margin(fig, axes[1, 0], PANEL_LABELS[2], fontsize=13.0)

    plot_summary_metric_grouped_bars(
        axes[1, 1],
        systems,
        "rna_compaction_summary",
        [
            (
                "inside_dense_only_spearman_BoundProteinChains_vs_LongRangeContactFraction",
                "Valency -> long-range",
                None,
            ),
            (
                "inside_dense_only_spearman_BoundProteinChains_vs_Rg_nm",
                "Valency -> collapse",
                lambda value: -value,
            ),
        ],
        "Valency-compaction coupling",
        "Coupling strength (Spearman)",
    )
    add_panel_label_in_margin(fig, axes[1, 1], PANEL_LABELS[3], fontsize=13.0)

    add_system_figure_legend(fig, y=0.92)
    fig.text(
        0.5,
        0.05,
        "This figure uses dynamic/ensemble proxies rather than a direct entropy calculation: WT combines faster contact turnover with stronger coupling between binding and RNA compaction.",
        ha="center",
        va="center",
        fontsize=7.1,
        color="#555555",
    )
    fig.suptitle(
        "Figure 19. Dynamic binding can be more compacting than persistent binding",
        fontsize=10,
        fontweight="bold",
        y=0.985,
    )
    return save_figure_all_formats(fig, outdir / "fig19_dynamic_binding_efficiency_compare", dpi)


def make_fig20_main_mechanism(
    systems: dict[str, AggregatedSystem], outdir: Path, dpi: int
) -> list[str]:
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.7))
    fig.subplots_adjust(wspace=0.42, top=0.77, bottom=0.20, left=0.08, right=0.98)

    plot_rna_binned_metric(
        axes[0],
        systems,
        "RgMean_nm",
        "RgStd_nm",
        "Larger RNA at matched valency",
        "RNA Rg (nm)",
    )
    add_panel_label_in_margin(fig, axes[0], "A", fontsize=13.0)

    plot_rna_binned_metric(
        axes[1],
        systems,
        "LongRangeContactFractionMean",
        "LongRangeContactFractionStd",
        "Lower tertiary compaction at matched valency",
        "Long-range contact fraction",
    )
    add_panel_label_in_margin(fig, axes[1], "B", fontsize=13.0)

    plot_summary_metric_bar(
        axes[2],
        systems,
        "liquidity_metrics",
        "mean_protein_rna_contact_lifetime_ns",
        "Longer protein-RNA lifetime",
        "Mean lifetime (ns)",
    )
    add_panel_label_in_margin(fig, axes[2], "C", fontsize=13.0)

    add_system_figure_legend(fig, y=0.89)
    fig.text(
        0.5,
        0.055,
        "At matched protein valency, phosphorylated RNA remains larger and less compact even though protein-RNA contacts persist longer.",
        ha="center",
        va="center",
        fontsize=7.2,
        color="#555555",
    )
    fig.suptitle(
        "Figure 20. Persistent but less productive RNA binding in the phosphorylated system",
        fontsize=10,
        fontweight="bold",
        y=0.985,
    )
    return save_figure_all_formats(fig, outdir / "fig20_main_mechanism_compare", dpi)


def make_fig21_rna_information_entropy(
    systems: dict[str, AggregatedSystem], outdir: Path, dpi: int
) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.9))
    fig.subplots_adjust(hspace=0.52, wspace=0.42, top=0.83, bottom=0.17)

    plot_summary_metric_bar(
        axes[0, 0],
        systems,
        "rna_information_summary",
        "inside_dense_only_H_compaction_bits",
        "RNA conformational entropy",
        "Entropy (bits)",
    )
    add_panel_label_in_margin(fig, axes[0, 0], PANEL_LABELS[0], fontsize=13.0)

    plot_summary_metric_bar(
        axes[0, 1],
        systems,
        "rna_information_summary",
        "inside_dense_only_H_partner_given_binding_bits",
        "Partner entropy at matched valency",
        "Entropy (bits)",
    )
    add_panel_label_in_margin(fig, axes[0, 1], PANEL_LABELS[1], fontsize=13.0)

    plot_summary_metric_bar(
        axes[1, 0],
        systems,
        "rna_information_summary",
        "inside_dense_only_MI_binding_compaction_bits",
        "Binding-compaction mutual information",
        "Mutual information (bits)",
    )
    add_panel_label_in_margin(fig, axes[1, 0], PANEL_LABELS[2], fontsize=13.0)

    plot_summary_metric_bar(
        axes[1, 1],
        systems,
        "rna_information_summary",
        "inside_dense_only_H_compaction_given_binding_bits",
        "Residual RNA entropy after conditioning on binding",
        "Conditional entropy (bits)",
    )
    add_panel_label_in_margin(fig, axes[1, 1], PANEL_LABELS[3], fontsize=13.0)

    fig.text(
        0.5,
        0.05,
        "All panels use the inside-dense-only scope. Lower mutual information together with higher residual compaction entropy is consistent with binding that is less effective at specifying RNA compaction.",
        ha="center",
        va="center",
        fontsize=7.1,
        color="#555555",
    )
    fig.suptitle(
        "Figure 21. RNA information-theoretic entropy proxies in the condensate",
        fontsize=10,
        fontweight="bold",
        y=0.985,
    )
    return save_figure_all_formats(fig, outdir / "fig21_rna_information_entropy_compare", dpi)


def make_fig22_rna_information_landscape(
    systems: dict[str, AggregatedSystem], outdir: Path, dpi: int
) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(7.8, 6.2))
    fig.subplots_adjust(hspace=0.48, wspace=0.40, top=0.84, bottom=0.15)

    landscapes = {
        system_key: compaction_state_landscape(system, preferred_scopes=("inside_dense_only", "contacting_only", "all_frames"))
        for system_key, system in systems.items()
    }
    valid_landscapes = [entry for entry in landscapes.values() if entry is not None]
    vmax = None
    x_limits = None
    y_limits = None
    if valid_landscapes:
        vmax = max(float(np.nanmax(matrix)) for _, _, matrix in valid_landscapes if matrix.size)
        x_limits = (
            min(float(edges[0]) for edges, _, _ in valid_landscapes),
            max(float(edges[-1]) for edges, _, _ in valid_landscapes),
        )
        y_limits = (
            min(float(edges[0]) for _, edges, _ in valid_landscapes),
            max(float(edges[-1]) for _, edges, _ in valid_landscapes),
        )

    artist = plot_compaction_state_landscape(
        axes[0, 0],
        systems["phos"],
        vmax=vmax,
        x_limits=x_limits,
        y_limits=y_limits,
    )
    add_panel_label_in_margin(fig, axes[0, 0], PANEL_LABELS[0], fontsize=13.0)

    second_artist = plot_compaction_state_landscape(
        axes[0, 1],
        systems["unphos"],
        vmax=vmax,
        x_limits=x_limits,
        y_limits=y_limits,
    )
    add_panel_label_in_margin(fig, axes[0, 1], PANEL_LABELS[1], fontsize=13.0)
    if artist is not None or second_artist is not None:
        reference = artist if artist is not None else second_artist
        fig.colorbar(reference, ax=[axes[0, 0], axes[0, 1]], fraction=0.045, pad=0.04, label="State probability")

    plot_summary_metric_bar(
        axes[1, 0],
        systems,
        "rna_information_summary",
        "inside_dense_only_transition_entropy_rate_bits_per_ns",
        "RNA-state transition entropy rate",
        "Entropy rate (bits/ns)",
    )
    add_panel_label_in_margin(fig, axes[1, 0], PANEL_LABELS[2], fontsize=13.0)

    plot_scalar_comparison_bar(
        axes[1, 1],
        {
            system: binding_explained_compaction_fraction(systems[system], scope="inside_dense_only")
            for system in SYSTEM_ORDER
        },
        "Compaction entropy explained by binding",
        "Explained fraction",
    )
    add_panel_label_in_margin(fig, axes[1, 1], PANEL_LABELS[3], fontsize=13.0)

    fig.text(
        0.5,
        0.045,
        "Top panels show the inside-dense RNA compaction-state landscape over Rg and long-range intrachain contacts; bottom panels summarize dynamic state exploration and the fraction of RNA compaction entropy explained by protein binding.",
        ha="center",
        va="center",
        fontsize=7.1,
        color="#555555",
    )
    fig.suptitle(
        "Figure 22. RNA state landscape and entropy-rate proxies",
        fontsize=10,
        fontweight="bold",
        y=0.985,
    )
    return save_figure_all_formats(fig, outdir / "fig22_rna_information_landscape_compare", dpi)


def plot_topology_binned_metric(
    ax: plt.Axes,
    systems: dict[str, AggregatedSystem],
    y_col: str,
    title: str,
    ylabel: str,
    *,
    preferred_scopes: Sequence[str] = ("inside_dense_only", "contacting_only", "all_frames"),
    min_bound: float = 1.0,
) -> None:
    plotted = False
    xmax = 0.0
    for system in SYSTEM_ORDER:
        table = systems[system].rna_topology_binned
        if table is None or y_col not in table.columns:
            continue
        scope_mask = scope_mask_from_table(table, preferred_scopes)
        x = np.asarray(table.columns["BoundProteinChainsBin"], dtype=float)
        y = np.asarray(table.columns[y_col], dtype=float)
        frames = np.asarray(table.columns["FrameCount"], dtype=float)
        mask = scope_mask & np.isfinite(x) & np.isfinite(y) & np.isfinite(frames)
        mask &= frames > 0.0
        mask &= x >= min_bound
        if not np.any(mask):
            continue
        order = np.argsort(x[mask])
        x_plot = x[mask][order]
        y_plot = y[mask][order]
        ax.plot(
            x_plot,
            y_plot,
            color=system_color(system),
            linewidth=1.8,
            marker="o",
            markersize=3.6,
            label=SYSTEM_DISPLAY[system],
        )
        xmax = max(xmax, float(np.nanmax(x_plot)))
        plotted = True
    if not plotted:
        draw_missing_axis(ax, f"{title}\nNot available")
        return
    ax.set_title(title)
    ax.set_xlabel("Bound protein chains")
    ax.set_ylabel(ylabel)
    ax.set_xlim(0.75, max(1.25, xmax + 0.25))
    ax.legend(loc="best", frameon=False)


def make_fig23_rna_bridge_clamp_topology(
    systems: dict[str, AggregatedSystem], outdir: Path, dpi: int
) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 6.0))
    fig.subplots_adjust(hspace=0.56, wspace=0.44, top=0.83, bottom=0.15)

    plot_summary_metric_grouped_bars(
        axes[0, 0],
        systems,
        "rna_topology_summary",
        [
            ("inside_dense_only_P_any_bridge_given_bound", "Bridge | bound", None),
            ("inside_dense_only_P_any_clamp_given_bound", "Clamp | bound", None),
        ],
        "Productive topology probability",
        "Probability",
    )
    add_panel_label_in_margin(fig, axes[0, 0], PANEL_LABELS[0], fontsize=13.0)

    plot_summary_metric_grouped_bars(
        axes[0, 1],
        systems,
        "rna_topology_summary",
        [
            ("inside_dense_only_P_longrange_given_bound", "Bound", None),
            ("inside_dense_only_P_longrange_given_bridge", "Bridge", None),
            ("inside_dense_only_P_longrange_given_clamp", "Clamp", None),
        ],
        "Loop closure conditional probability",
        "P(long-range RNA contact)",
    )
    add_panel_label_in_margin(fig, axes[0, 1], PANEL_LABELS[1], fontsize=13.0)

    plot_summary_metric_grouped_bars(
        axes[1, 0],
        systems,
        "rna_topology_summary",
        [
            ("inside_dense_only_bridge_productivity_index", "Bridge", None),
            ("inside_dense_only_clamp_productivity_index", "Clamp", None),
        ],
        "Bridge/clamp productivity",
        "Long-range contact fraction",
    )
    add_panel_label_in_margin(fig, axes[1, 0], PANEL_LABELS[2], fontsize=13.0)

    plot_summary_metric_grouped_bars(
        axes[1, 1],
        systems,
        "rna_topology_summary",
        [
            ("inside_dense_only_bridge_enrichment_over_bound", "Bridge", None),
            ("inside_dense_only_clamp_enrichment_over_bound", "Clamp", None),
        ],
        "Loop-closure enrichment",
        "Enrichment over bound",
    )
    axes[1, 1].axhline(1.0, color="#777777", linewidth=0.9, linestyle=":")
    add_panel_label_in_margin(fig, axes[1, 1], PANEL_LABELS[3], fontsize=13.0)

    add_system_figure_legend(fig, y=0.92)
    fig.text(
        0.5,
        0.045,
        "Inside-dense topology metrics distinguish local decorating contacts from bridge/clamp states that support RNA loop closure.",
        ha="center",
        va="center",
        fontsize=7.1,
        color="#555555",
    )
    fig.suptitle(
        "Figure 23. RNA bridge/clamp topology explains productive binding",
        fontsize=10,
        fontweight="bold",
        y=0.985,
    )
    return save_figure_all_formats(fig, outdir / "fig23_rna_bridge_clamp_topology_compare", dpi)


def make_fig24_local_pp_competition(
    systems: dict[str, AggregatedSystem], outdir: Path, dpi: int
) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 6.0))
    fig.subplots_adjust(hspace=0.56, wspace=0.46, top=0.83, bottom=0.15)

    plot_summary_metric_bar(
        axes[0, 0],
        systems,
        "rna_topology_summary",
        "inside_dense_only_mean_local_pp_attraction_strength_kjmol",
        "Local PP attraction near RNA",
        r"Mean $-E_{\mathrm{PP}}$ (kJ/mol)",
    )
    add_panel_label_in_margin(fig, axes[0, 0], PANEL_LABELS[0], fontsize=13.0)

    plot_summary_metric_bar(
        axes[0, 1],
        systems,
        "rna_topology_summary",
        "inside_dense_only_mean_local_pp_to_prna_attraction_ratio",
        "PP attraction relative to PRNA",
        r"$-E_{\mathrm{PP}} / -E_{\mathrm{PRNA}}$",
    )
    add_panel_label_in_margin(fig, axes[0, 1], PANEL_LABELS[1], fontsize=13.0)

    plot_summary_metric_bar(
        axes[1, 0],
        systems,
        "rna_topology_summary",
        "inside_dense_only_mean_local_pp_to_bridge_attraction_ratio",
        "PP attraction relative to bridge",
        r"$-E_{\mathrm{PP}} / -E_{\mathrm{bridge}}$",
    )
    add_panel_label_in_margin(fig, axes[1, 0], PANEL_LABELS[2], fontsize=13.0)

    plot_topology_binned_metric(
        axes[1, 1],
        systems,
        "MeanLocalPPtoPRNAAttractionRatio",
        "Competition at matched valency",
        r"$-E_{\mathrm{PP}} / -E_{\mathrm{PRNA}}$",
    )
    add_panel_label_in_margin(fig, axes[1, 1], PANEL_LABELS[3], fontsize=13.0)

    add_system_figure_legend(fig, y=0.92)
    fig.text(
        0.5,
        0.045,
        "Higher local protein-protein attraction around bound RNA would support a self-satisfied protein network that competes with productive RNA bridging.",
        ha="center",
        va="center",
        fontsize=7.1,
        color="#555555",
    )
    fig.suptitle(
        "Figure 24. Local protein-protein competition around RNA",
        fontsize=10,
        fontweight="bold",
        y=0.985,
    )
    return save_figure_all_formats(fig, outdir / "fig24_local_pp_competition_compare", dpi)


def make_fig3_interface(
    systems: dict[str, AggregatedSystem],
    outdir: Path,
    dpi: int,
    residue_min: int,
    residue_max: int,
    smooth_window: int,
) -> list[str]:
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.8), sharex=True, gridspec_kw={"height_ratios": [2.0, 1.0]})
    fig.subplots_adjust(hspace=0.26)
    plot_interface_profile(axes[0], axes[1], systems, residue_min, residue_max, smooth_window)
    add_panel_label(axes[0], "A")
    add_panel_label(axes[1], "B")
    fig.suptitle("Figure 3. Protein residue-wise RNA interaction profile", fontsize=10, fontweight="bold", y=0.995)
    return save_figure_all_formats(fig, outdir / "fig3_interface_compare", dpi)


def make_fig4_dynamics(systems: dict[str, AggregatedSystem], outdir: Path, dpi: int) -> list[str]:
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))
    ax_ct, ax_tau = axes
    plotted = False
    for system in SYSTEM_ORDER:
        table = systems[system].llps_autocorr
        if table is None:
            continue
        x = np.asarray(table.columns["lag_ns"], dtype=float)
        y = np.asarray(table.columns["ct"], dtype=float)
        if x.size == 0 or y.size == 0:
            continue
        ax_ct.plot(x, y, color=system_color(system), linewidth=1.6, label=SYSTEM_DISPLAY[system])
        plotted = True
    if plotted:
        ax_ct.set_title("Contact autocorrelation")
        ax_ct.set_xlabel("Lag (ns)")
        ax_ct.set_ylabel("C(t)")
        ax_ct.grid(True, axis="y")
        ax_ct.legend(loc="best")
    else:
        draw_missing_axis(ax_ct, "Contact autocorrelation\nNot available")

    labels, values, colors = compute_tau_bar_values(systems)
    xpos = np.arange(len(labels))
    ax_tau.bar(xpos, values, color=colors, width=0.6)
    ax_tau.set_xticks(xpos)
    ax_tau.set_xticklabels(labels, rotation=0)
    ax_tau.set_ylabel(r"$\tau$ (ns)")
    ax_tau.set_title("Relaxation time")
    ax_tau.grid(True, axis="y")
    for idx, value in enumerate(values):
        ax_tau.text(idx, value + max(values + [1.0]) * 0.03, f"{value:.2f}", ha="center", va="bottom", fontsize=8)

    for idx, axis in enumerate(axes):
        add_panel_label(axis, PANEL_LABELS[idx])
    fig.suptitle("Figure 4. LLPS dynamics", fontsize=10, fontweight="bold", y=1.04)
    return save_figure_all_formats(fig, outdir / "fig4_dynamics_compare", dpi)


def make_fig5_liquidity(systems: dict[str, AggregatedSystem], outdir: Path, dpi: int) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6))
    fig.subplots_adjust(wspace=0.28, hspace=0.38)
    metric_specs = [
        ("tau_ns", "Contact relaxation time", r"$\tau$ (ns)"),
        ("cv_num_condensates", "Steady-state CV: condensate count", "CV"),
        ("cv_giant_cluster_prochains", "Steady-state CV: protein chains in giant cluster", "CV"),
        ("cv_rg_rna_dense", "Steady-state CV: dense-phase RNA Rg", "CV"),
    ]
    for idx, (metric_key, title, ylabel) in enumerate(metric_specs):
        plot_liquidity_bar(axes.ravel()[idx], systems, metric_key, title, ylabel)
        add_panel_label(axes.ravel()[idx], PANEL_LABELS[idx])
    fig.suptitle("Figure 5. Metrics consistent with condensate liquidity", fontsize=10, fontweight="bold", y=0.995)
    return save_figure_all_formats(fig, outdir / "fig5_liquidity_compare", dpi)


def make_fig6_contact_dynamics(systems: dict[str, AggregatedSystem], outdir: Path, dpi: int) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.8))
    fig.subplots_adjust(wspace=0.32, hspace=0.48, top=0.84)
    panel_axes = axes.ravel()
    plot_lifetime_survival(
        panel_axes[0],
        systems,
        "sticker_contact_lifetimes",
        "Sticker-contact lifetime",
        show_legend=False,
    )
    plot_lifetime_survival(
        panel_axes[1],
        systems,
        "protein_rna_contact_lifetimes",
        "Protein-RNA contact lifetime",
        show_legend=False,
    )
    plot_grouped_summary_bars(
        panel_axes[2],
        systems,
        "neighbor_rewiring",
        [
            ("protein_neighbor_persistence", "Protein"),
            ("rna_neighbor_persistence", "RNA"),
            ("global_neighbor_persistence", "Global"),
        ],
        "Neighbor persistence",
        "Jaccard persistence",
        show_legend=False,
    )
    plot_grouped_summary_bars(
        panel_axes[3],
        systems,
        "interface_turnover",
        [
            ("new_contacts_per_ns", "New"),
            ("lost_contacts_per_ns", "Lost"),
            ("net_change_per_ns", "Net"),
        ],
        "Protein-RNA interface turnover",
        "Contacts/ns",
        show_legend=False,
    )
    for idx, axis in enumerate(panel_axes):
        add_panel_label(axis, PANEL_LABELS[idx])
    fig.suptitle("Figure 6. Contact lifetimes and interface turnover", fontsize=10, fontweight="bold", y=0.995)
    add_system_figure_legend(fig, y=0.94)
    return save_figure_all_formats(fig, outdir / "fig6_llps_contact_dynamics_compare", dpi)


def make_fig7_prna_interface_dynamics(systems: dict[str, AggregatedSystem], outdir: Path, dpi: int) -> list[str]:
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2))
    fig.subplots_adjust(wspace=0.32, top=0.74, bottom=0.24)
    plot_lifetime_survival(
        axes[0],
        systems,
        "protein_rna_contact_lifetimes",
        "Protein-RNA contact lifetime",
        log_x=True,
        show_legend=False,
    )
    plot_grouped_summary_bars(
        axes[1],
        systems,
        "interface_turnover",
        [
            ("new_contacts_per_ns", "New"),
            ("lost_contacts_per_ns", "Lost"),
            ("net_change_per_ns", "Net"),
        ],
        "Protein-RNA interface turnover",
        "Contacts/ns",
        show_legend=False,
    )
    for idx, axis in enumerate(axes):
        add_panel_label(axis, PANEL_LABELS[idx])
    fig.suptitle("Figure 7. Protein-RNA interface stability and turnover", fontsize=10, fontweight="bold", y=0.995)
    add_system_figure_legend(fig, y=0.90)
    return save_figure_all_formats(fig, outdir / "fig7_prna_interface_dynamics_compare", dpi)


def make_figs4_turnover_timeseries(systems: dict[str, AggregatedSystem], outdir: Path, dpi: int) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6))
    fig.subplots_adjust(wspace=0.28, hspace=0.40)
    panel_axes = axes.ravel()
    plot_dynamics_line(panel_axes[0], systems, "neighbor_rewiring", "protein_persistence", "Protein neighbor persistence", "Jaccard")
    plot_dynamics_line(panel_axes[1], systems, "neighbor_rewiring", "rna_persistence", "RNA neighbor persistence", "Jaccard")
    plot_dynamics_line(panel_axes[2], systems, "interface_turnover", "new_per_ns", "New protein-RNA contacts", "Contacts/ns")
    plot_dynamics_line(panel_axes[3], systems, "interface_turnover", "lost_per_ns", "Lost protein-RNA contacts", "Contacts/ns")
    for idx, axis in enumerate(panel_axes):
        add_panel_label(axis, PANEL_LABELS[idx])
    fig.suptitle("Supplementary Figure S4. Time-resolved LLPS contact dynamics", fontsize=10, fontweight="bold", y=0.995)
    return save_figure_all_formats(fig, outdir / "figS4_llps_turnover_timeseries_compare", dpi)


def make_figs1_degree(systems: dict[str, AggregatedSystem], outdir: Path, dpi: int) -> list[str]:
    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    plotted = False
    for system in SYSTEM_ORDER:
        table = systems[system].llps_degree
        if table is None:
            continue
        k = np.asarray(table.columns["k"], dtype=float)
        mean = np.asarray(table.columns["pk_mean"], dtype=float)
        std = np.asarray(table.columns["pk_std"], dtype=float)
        ax.errorbar(
            k,
            mean,
            yerr=std,
            color=system_color(system),
            linewidth=1.2,
            marker="o",
            markersize=3.5,
            capsize=2.5,
            label=SYSTEM_DISPLAY[system],
        )
        plotted = True
    if plotted:
        ax.set_title("Supplementary Figure S1. Degree distribution")
        ax.set_xlabel("k")
        ax.set_ylabel("P(k)")
        ax.grid(True, axis="y")
        ax.legend(loc="best")
    else:
        draw_missing_axis(ax, "Degree distribution\nNot available")
    add_panel_label(ax, "A")
    return save_figure_all_formats(fig, outdir / "figS1_degree_compare", dpi)


def make_figs2_network(systems: dict[str, AggregatedSystem], outdir: Path, dpi: int, topn: int) -> list[str]:
    fig, axes = plt.subplots(3, 2, figsize=(7.2, 8.2))
    metric_builders = [
        ("Node strength", top_network_strength, None),
        ("Top GC edges", top_edges, "gc_edges"),
        ("Top TE edges", top_edges, "te_edges"),
    ]
    for row_idx, (title, builder, selector) in enumerate(metric_builders):
        for col_idx, system in enumerate(SYSTEM_ORDER):
            bundle = systems[system].network_bundle
            if selector is None:
                entries = builder(bundle, topn)
            else:
                entries = builder(bundle, selector, topn)
            plot_ranked_horizontal(
                axes[row_idx, col_idx],
                entries,
                f"{title}: {system}",
                system_color(system),
            )
            if col_idx == 0:
                axes[row_idx, col_idx].set_ylabel("")
    for idx, axis in enumerate(axes.ravel()):
        add_panel_label(axis, PANEL_LABELS[idx])
    fig.suptitle("Supplementary Figure S2. Network / allostery comparison", fontsize=10, fontweight="bold", y=0.995)
    return save_figure_all_formats(fig, outdir / "figS2_network_allostery_compare", dpi)


def make_fig25_dynamic_anchor_network(
    systems: dict[str, AggregatedSystem], outdir: Path, dpi: int, topn: int
) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.4))
    fig.subplots_adjust(hspace=0.46, wspace=0.50, top=0.82, bottom=0.13, left=0.25, right=0.98)
    panels = [
        ("phos", "gc_anchor", "GC anchors: phosphorylated"),
        ("unphos", "gc_anchor", "GC anchors: WT"),
        ("phos", "te_source", "TE source anchors: phosphorylated"),
        ("unphos", "te_source", "TE source anchors: WT"),
    ]
    for idx, (system, metric, title) in enumerate(panels):
        axis = axes.ravel()[idx]
        entries = network_anchor_entries(systems[system], metric, topn)
        plot_ranked_horizontal(axis, entries, title, system_color(system))
        axis.set_xlabel("Anchor score")
        add_panel_label_in_margin(fig, axis, PANEL_LABELS[idx], fontsize=13.0)
    add_system_figure_legend(fig, y=0.91)
    fig.text(
        0.5,
        0.045,
        "GC anchors are incident correlation-strength hubs; TE source anchors are nodes with positive net outgoing transfer entropy.",
        ha="center",
        va="center",
        fontsize=7.2,
        color="#555555",
    )
    fig.suptitle(
        "Figure 25. Dynamic network anchors in protein-RNA condensates",
        fontsize=10,
        fontweight="bold",
        y=0.985,
    )
    return save_figure_all_formats(fig, outdir / "fig25_dynamic_anchor_network_compare", dpi)


def plot_anchor_signature_axis(
    ax: plt.Axes,
    systems: dict[str, AggregatedSystem],
    metric: str,
    molecule: str,
    title: str,
) -> None:
    points_by_system = collect_network_anchor_points(systems, metric, molecule, topn=30)
    plotted = False
    all_indices: list[float] = []
    jitter = {"phos": -0.65, "unphos": 0.65}
    for system in SYSTEM_ORDER:
        points = points_by_system[system]
        if not points:
            continue
        plotted = True
        color = system_color(system)
        x = np.asarray([float(point["index"]) + jitter[system] for point in points], dtype=float)
        y = np.asarray([float(point["relative_score"]) for point in points], dtype=float)
        ranks = np.asarray([float(point["rank"]) for point in points], dtype=float)
        sizes = 78.0 + 250.0 * (1.0 / np.sqrt(np.maximum(ranks, 1.0)))
        all_indices.extend(float(point["index"]) for point in points)
        ax.vlines(x, 0.0, y, color=color, linewidth=0.9, alpha=0.38)
        ax.scatter(
            x,
            y,
            s=sizes,
            color=color,
            edgecolor="white",
            linewidth=0.7,
            alpha=0.92,
            label=SYSTEM_DISPLAY[system],
            zorder=3,
        )
    if not plotted:
        draw_missing_axis(ax, f"{title}\nNot available")
        return
    if molecule == "protein":
        ax.set_xlim(0, 420)
        ax.set_xlabel("Protein residue index")
    else:
        margin = 12.0
        xmin = min(all_indices) - margin if all_indices else 0.0
        xmax = max(all_indices) + margin if all_indices else 1.0
        ax.set_xlim(xmin, xmax)
        ax.set_xlabel("RNA bead / nucleotide index")
    ax.set_ylim(0, 112)
    ax.set_title(title)
    ax.set_ylabel("Relative anchor strength (%)")
    ax.axhline(0.0, color="#222222", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def make_fig26_dynamic_anchor_signature(
    systems: dict[str, AggregatedSystem], outdir: Path, dpi: int
) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.0))
    fig.subplots_adjust(hspace=0.50, wspace=0.35, top=0.80, bottom=0.15, left=0.10, right=0.98)
    panels = [
        ("gc_anchor", "protein", "GC protein-anchor landscape"),
        ("gc_anchor", "rna", "GC RNA-anchor landscape"),
        ("te_source", "protein", "TE source protein anchors"),
        ("te_source", "rna", "TE source RNA anchors"),
    ]
    for idx, (metric, molecule, title) in enumerate(panels):
        axis = axes.ravel()[idx]
        plot_anchor_signature_axis(axis, systems, metric, molecule, title)
        add_panel_label_in_margin(fig, axis, PANEL_LABELS[idx], fontsize=13.0)
    add_system_figure_legend(fig, y=0.90)
    fig.text(
        0.5,
        0.055,
        "Anchor scores are normalized within each system and metric. Bubble size follows rank among top network anchors.",
        ha="center",
        va="center",
        fontsize=7.1,
        color="#555555",
    )
    fig.suptitle(
        "Figure 26. Molecular signature of dynamic network anchors",
        fontsize=10,
        fontweight="bold",
        y=0.985,
    )
    return save_figure_all_formats(fig, outdir / "fig26_dynamic_anchor_signature_compare", dpi)


def make_figs3_tda(systems: dict[str, AggregatedSystem], outdir: Path, dpi: int) -> list[str]:
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharex=True, sharey=True)
    valid_points = [system.tda_points for system in systems.values() if system.tda_points is not None and system.tda_points.size]
    if valid_points:
        stacked = np.concatenate(valid_points, axis=0)
        extent = [float(np.min(stacked[:, 0])), float(np.max(stacked[:, 0])), float(np.min(stacked[:, 1])), float(np.max(stacked[:, 1]))]
    else:
        extent = [0.0, 1.0, 0.0, 1.0]

    artists = []
    for axis, system_key in zip(axes, SYSTEM_ORDER):
        points = systems[system_key].tda_points
        if points is None or points.size == 0:
            draw_missing_axis(axis, f"{SYSTEM_DISPLAY[system_key]}\nNo TDA point cloud")
            continue
        artist = axis.hexbin(
            points[:, 0],
            points[:, 1],
            gridsize=60,
            mincnt=1,
            cmap="viridis",
            extent=extent,
        )
        artists.append(artist)
        axis.set_title(SYSTEM_DISPLAY[system_key])
        axis.set_xlabel("X")
        axis.set_ylabel("Y")
    if artists:
        fig.colorbar(artists[0], ax=axes, fraction=0.046, pad=0.04, label="Counts")
    for idx, axis in enumerate(axes):
        add_panel_label(axis, PANEL_LABELS[idx])
    fig.suptitle("Supplementary Figure S3. Pooled TDA point cloud (XY projection)", fontsize=10, fontweight="bold", y=1.02)
    return save_figure_all_formats(fig, outdir / "figS3_tda_compare", dpi)


def export_system_outputs(system: AggregatedSystem, outdir: Path, options: RuntimeOptions) -> None:
    prefix = system.system
    if system.llps_timeseries is not None:
        save_csv_table(system.llps_timeseries, outdir / f"{prefix}_avg_llps_timeseries.csv")
    if system.llps_autocorr is not None:
        save_csv_table(system.llps_autocorr, outdir / f"{prefix}_avg_llps_autocorr.csv")
    if system.dbscan_timeseries is not None:
        save_csv_table(system.dbscan_timeseries, outdir / f"{prefix}_avg_llps_dbscan_timeseries.csv")
    if system.polymer_scaling is not None:
        save_csv_table(system.polymer_scaling, outdir / f"{prefix}_avg_polymer_scaling.csv")
    if system.llps_degree is not None:
        save_csv_table(system.llps_degree, outdir / f"{prefix}_avg_llps_degree_dist.csv")
    save_summary_text(system, outdir / f"{prefix}_avg_llps_summary.txt")
    if options.export_density_dx and system.density_protein is not None and system.density_protein.has_full_data:
        save_density_dx(system.density_protein, outdir / f"{prefix}_avg_density_protein.dx")
    if options.export_density_dx and system.density_rna is not None and system.density_rna.has_full_data:
        save_density_dx(system.density_rna, outdir / f"{prefix}_avg_density_rna.dx")
    if system.density_protein is not None and system.density_protein.radial_profile_data is not None:
        save_csv_table(system.density_protein.radial_profile_data, outdir / f"{prefix}_avg_radial_density_protein.csv")
    if system.density_rna is not None and system.density_rna.radial_profile_data is not None:
        save_csv_table(system.density_rna.radial_profile_data, outdir / f"{prefix}_avg_radial_density_rna.csv")
    if system.density_shape_metrics:
        save_density_shape_metrics(system, outdir / f"{prefix}_avg_density_shape_metrics.csv")
    if system.nucleic_interactions:
        save_nucleic_interactions(system.nucleic_interactions, outdir / f"{prefix}_avg_nucleic_acid_interactions.csv")
    if system.sidechain_protein:
        save_sidechain_contacts(system.sidechain_protein, outdir / f"{prefix}_avg_protein_sidechain_contacts.csv")
    if system.sidechain_rna:
        save_sidechain_contacts(system.sidechain_rna, outdir / f"{prefix}_avg_rna_sidechain_contacts.csv")
    if options.include_network and system.network_bundle is not None:
        save_network_bundle(system.network_bundle, prefix, outdir)
    if options.include_tda and system.tda_points is not None:
        save_tda_points(system.tda_points, outdir / f"{prefix}_avg_tda_global_pointcloud.csv")
    if system.sticker_contact_lifetimes is not None:
        save_dynamics_bundle(
            system.sticker_contact_lifetimes,
            outdir / f"{prefix}_avg_llps_sticker_contact_lifetimes.csv",
        )
    if system.protein_rna_contact_lifetimes is not None:
        save_dynamics_bundle(
            system.protein_rna_contact_lifetimes,
            outdir / f"{prefix}_avg_llps_protein_rna_contact_lifetimes.csv",
        )
    if system.neighbor_rewiring is not None:
        save_dynamics_bundle(system.neighbor_rewiring, outdir / f"{prefix}_avg_llps_neighbor_rewiring.csv")
    if system.interface_turnover is not None:
        save_dynamics_bundle(system.interface_turnover, outdir / f"{prefix}_avg_llps_interface_turnover.csv")
    if system.dense_msd_protein is not None:
        save_dynamics_bundle(system.dense_msd_protein, outdir / f"{prefix}_avg_llps_dense_msd_protein.csv")
    if system.dense_msd_rna is not None:
        save_dynamics_bundle(system.dense_msd_rna, outdir / f"{prefix}_avg_llps_dense_msd_rna.csv")
    if system.rna_compaction_summary:
        save_key_value_summary(
            system.rna_compaction_summary,
            outdir / f"{prefix}_avg_llps_rna_compaction_binding_summary.txt",
        )
    if system.rna_compaction_binned is not None:
        save_csv_table(
            system.rna_compaction_binned,
            outdir / f"{prefix}_avg_llps_rna_compaction_binding_binned.csv",
        )
    if system.rna_energy_compaction_summary:
        save_key_value_summary(
            system.rna_energy_compaction_summary,
            outdir / f"{prefix}_avg_llps_rna_energy_compaction_summary.txt",
        )
    if system.rna_energy_compaction_binned is not None:
        save_csv_table(
            system.rna_energy_compaction_binned,
            outdir / f"{prefix}_avg_llps_rna_energy_compaction_binned.csv",
        )
    if system.rna_information_summary:
        save_key_value_summary(
            system.rna_information_summary,
            outdir / f"{prefix}_avg_llps_rna_information_summary.txt",
        )
    if system.rna_compaction_state_histogram is not None:
        save_csv_table(
            system.rna_compaction_state_histogram,
            outdir / f"{prefix}_avg_llps_rna_compaction_state_histogram.csv",
        )
    if system.rna_topology_summary:
        save_key_value_summary(
            system.rna_topology_summary,
            outdir / f"{prefix}_avg_llps_rna_topology_summary.txt",
        )
    if system.rna_topology_binned is not None:
        save_csv_table(
            system.rna_topology_binned,
            outdir / f"{prefix}_avg_llps_rna_topology_binding_binned.csv",
        )
    if system.rna_topology_trace is not None:
        save_csv_table(
            system.rna_topology_trace,
            outdir / f"{prefix}_pooled_llps_rna_topology_trace.csv",
        )
    if system.rna_chain_compaction_summary is not None:
        save_csv_table(
            system.rna_chain_compaction_summary,
            outdir / f"{prefix}_avg_llps_rna_chain_compaction_summary.csv",
        )
    if system.rna_intrachain_contacts_summary is not None:
        save_csv_table(
            system.rna_intrachain_contacts_summary,
            outdir / f"{prefix}_avg_llps_rna_intrachain_contacts_summary.csv",
        )
    if system.rna_compaction_trace is not None:
        save_csv_table(
            system.rna_compaction_trace,
            outdir / f"{prefix}_pooled_llps_rna_compaction_trace.csv",
        )
    if system.rna_binding_energy_trace is not None:
        save_csv_table(
            system.rna_binding_energy_trace,
            outdir / f"{prefix}_pooled_llps_rna_binding_energy_trace.csv",
        )


def write_comparison_summary(
    systems: dict[str, AggregatedSystem],
    outdir: Path,
    figure_notes: Sequence[str],
    panel_status: dict[str, str],
    max_time_ns: float | None,
) -> None:
    path = outdir / "comparison_summary.txt"
    common_steady = max(system.steady_start_ns for system in systems.values())
    with path.open("w", encoding="utf-8") as handle:
        handle.write("Comparison summary: phos vs unphos\n")
        handle.write(f"steady_start_ns_common: {common_steady:.6f}\n")
        if max_time_ns is not None:
            handle.write(f"max_time_ns_cap: {max_time_ns:.6f}\n")
        handle.write("\n[Replicates]\n")
        for system_key in SYSTEM_ORDER:
            system = systems[system_key]
            handle.write(f"{system_key} ({SYSTEM_DISPLAY[system_key]}):\n")
            handle.write(f"  total_weighted_frames: {system.total_weight}\n")
            handle.write(f"  steady_start_ns: {system.steady_start_ns:.6f}\n")
            if system.llps_timeseries is not None:
                times = system.llps_timeseries.columns["time_ns"]
                handle.write(
                    f"  llps_time_range_ns: {float(np.min(times)):.6f} -> {float(np.max(times)):.6f}\n"
                )
            if system.dbscan_timeseries is not None:
                times_ps = system.dbscan_timeseries.columns["Time_ps"]
                handle.write(
                    f"  dbscan_time_range_ns: {float(np.min(times_ps)) / 1000.0:.6f} -> {float(np.max(times_ps)) / 1000.0:.6f}\n"
                )
            for rep in system.replicates:
                handle.write(f"  - {rep.replicate_label}: {rep.source_dir} (frame_weight={rep.frame_weight})\n")
                if rep.warnings:
                    for warning in rep.warnings:
                        handle.write(f"      warning: {warning}\n")
            handle.write("  summary_metrics:\n")
            for key in ("percolated", "percolation_time_ns", "tau_ns", "mean_gcc_rg_nm", "mean_gcc_fraction"):
                if key in system.llps_summary:
                    handle.write(f"    {key}: {system.llps_summary[key]}\n")
            if system.warnings:
                handle.write("  system_warnings:\n")
                for warning in sorted(set(system.warnings)):
                    handle.write(f"    - {warning}\n")
            if system.liquidity_metrics:
                handle.write("  liquidity_metrics:\n")
                for key in LIQUIDITY_REPORT_KEYS:
                    if key in system.liquidity_metrics:
                        handle.write(f"    {key}: {system.liquidity_metrics[key]}\n")
            if system.rna_compaction_summary:
                handle.write("  rna_compaction_binding_summary:\n")
                for key in (
                    "samples_all",
                    "samples_inside_dense",
                    "all_frames_spearman_BoundProteinChains_vs_Rg_nm",
                    "all_frames_spearman_BoundProteinChains_vs_EndToEnd_nm",
                    "all_frames_spearman_BoundProteinChains_vs_RelativeShapeAnisotropy",
                    "all_frames_spearman_BoundProteinChains_vs_LongRangeContactFraction",
                    "inside_dense_only_spearman_BoundProteinChains_vs_Rg_nm",
                    "inside_dense_only_spearman_BoundProteinChains_vs_LongRangeContactFraction",
                ):
                    if key in system.rna_compaction_summary:
                        handle.write(f"    {key}: {system.rna_compaction_summary[key]}\n")
            if system.rna_energy_compaction_summary:
                handle.write("  rna_energy_compaction_summary:\n")
                for key in (
                    "samples_all",
                    "samples_inside_dense",
                    "samples_contacting",
                    "all_frames_spearman_AttractionStrength_kJmol_vs_Rg_nm",
                    "all_frames_spearman_AttractionStrength_kJmol_vs_LongRangeContactFraction",
                    "inside_dense_only_spearman_AttractionStrength_kJmol_vs_Rg_nm",
                    "inside_dense_only_spearman_AttractionStrength_kJmol_vs_LongRangeContactFraction",
                    "contacting_only_spearman_AttractionStrength_kJmol_vs_Rg_nm",
                    "contacting_only_spearman_AttractionStrength_kJmol_vs_LongRangeContactFraction",
                ):
                    if key in system.rna_energy_compaction_summary:
                        handle.write(f"    {key}: {system.rna_energy_compaction_summary[key]}\n")
            if system.rna_information_summary:
                handle.write("  rna_information_summary:\n")
                for key in (
                    "inside_dense_only_H_compaction_bits",
                    "inside_dense_only_H_partner_bits",
                    "inside_dense_only_H_partner_given_binding_bits",
                    "inside_dense_only_MI_binding_compaction_bits",
                    "inside_dense_only_MI_binding_compaction_norm",
                    "inside_dense_only_H_compaction_given_binding_bits",
                    "inside_dense_only_transition_entropy_rate_bits_per_ns",
                ):
                    if key in system.rna_information_summary:
                        handle.write(f"    {key}: {system.rna_information_summary[key]}\n")
            if system.rna_topology_summary:
                handle.write("  rna_topology_summary:\n")
                for key in (
                    "inside_dense_only_P_any_bridge_given_bound",
                    "inside_dense_only_P_any_clamp_given_bound",
                    "inside_dense_only_P_longrange_given_bound",
                    "inside_dense_only_P_longrange_given_bridge",
                    "inside_dense_only_P_longrange_given_clamp",
                    "inside_dense_only_bridge_productivity_index",
                    "inside_dense_only_clamp_productivity_index",
                    "inside_dense_only_mean_local_pp_attraction_strength_kjmol",
                    "inside_dense_only_mean_local_pp_to_prna_attraction_ratio",
                    "inside_dense_only_mean_local_pp_to_bridge_attraction_ratio",
                ):
                    if key in system.rna_topology_summary:
                        handle.write(f"    {key}: {system.rna_topology_summary[key]}\n")
            handle.write("\n")

        handle.write("[Panel status]\n")
        for key, status in panel_status.items():
            handle.write(f"{key}: {status}\n")

        if figure_notes:
            handle.write("\n[Figure export notes]\n")
            for note in figure_notes:
                handle.write(f"- {note}\n")


def determine_panel_status(systems: dict[str, AggregatedSystem], options: RuntimeOptions) -> dict[str, str]:
    status: dict[str, str] = {}
    status["fig1_llps_timeseries_compare"] = "ready" if any(system.dbscan_timeseries is not None for system in systems.values()) else "placeholder"
    has_density = any(system.density_protein is not None or system.density_rna is not None for system in systems.values())
    status["fig2_density_polymer_compare"] = "ready" if has_density else "placeholder"
    status["fig3_interface_compare"] = (
        "nucleic"
        if all(system.interface_mode == "nucleic" for system in systems.values())
        else ("sidechain_fallback" if any(system.interface_mode == "sidechain" for system in systems.values()) else "placeholder")
    )
    status["fig4_dynamics_compare"] = "ready" if any(system.llps_autocorr is not None for system in systems.values()) else "placeholder"
    status["fig5_liquidity_compare"] = "ready" if any(system.liquidity_metrics for system in systems.values()) else "placeholder"
    status["fig6_llps_contact_dynamics_compare"] = (
        "ready"
        if any(
            system.sticker_contact_lifetimes is not None
            or system.protein_rna_contact_lifetimes is not None
            or system.neighbor_rewiring is not None
            or system.interface_turnover is not None
            for system in systems.values()
        )
        else "placeholder"
    )
    status["fig7_prna_interface_dynamics_compare"] = (
        "ready"
        if any(system.protein_rna_contact_lifetimes is not None or system.interface_turnover is not None for system in systems.values())
        else "placeholder"
    )
    status["fig8_radial_density_protein_compare"] = (
        "ready" if any(radial_density_table(system, "protein") is not None for system in systems.values()) else "placeholder"
    )
    status["fig9_radial_density_rna_compare"] = (
        "ready" if any(radial_density_table(system, "rna") is not None for system in systems.values()) else "placeholder"
    )
    status["fig10_density_shape_asymmetry_compare"] = (
        "ready" if any(system.density_shape_metrics for system in systems.values()) else "placeholder"
    )
    status["fig11_dense_msd_diffusion_compare"] = (
        "ready"
        if any(system.dense_msd_protein is not None or system.dense_msd_rna is not None for system in systems.values())
        else "placeholder"
    )
    status["fig12_rg_rna_dense_histogram_compare"] = (
        "ready" if any(system.rg_rna_dense_values is not None for system in systems.values()) else "placeholder"
    )
    status["fig13_rna_binding_compaction_binned_compare"] = (
        "ready" if any(system.rna_compaction_binned is not None for system in systems.values()) else "placeholder"
    )
    status["fig14_rna_chain_compaction_compare"] = (
        "ready" if any(system.rna_chain_compaction_summary is not None for system in systems.values()) else "placeholder"
    )
    status["fig15_rna_compaction_binding_trace_compare"] = (
        "ready" if any(system.rna_compaction_trace is not None for system in systems.values()) else "placeholder"
    )
    status["fig16_rna_energy_compaction_compare"] = (
        "ready"
        if any(
            system.rna_energy_compaction_binned is not None or system.rna_binding_energy_trace is not None
            for system in systems.values()
        )
        else "placeholder"
    )
    status["fig17_mechanistic_summary_compare"] = (
        "ready"
        if any(
            system.rna_binding_energy_trace is not None
            or system.dense_msd_protein is not None
            or system.dense_msd_rna is not None
            or system.rg_rna_dense_values is not None
            for system in systems.values()
        )
        else "placeholder"
    )
    status["fig18_compaction_efficiency_compare"] = (
        "ready"
        if any(
            system.rna_binding_energy_trace is not None
            or system.rna_compaction_trace is not None
            or system.rna_compaction_binned is not None
            for system in systems.values()
        )
        else "placeholder"
    )
    status["fig19_dynamic_binding_efficiency_compare"] = (
        "ready"
        if any(
            system.rna_energy_compaction_summary
            or system.rna_compaction_summary
            or system.interface_turnover is not None
            or system.liquidity_metrics
            for system in systems.values()
        )
        else "placeholder"
    )
    status["fig20_main_mechanism_compare"] = (
        "ready"
        if any(
            system.rna_compaction_binned is not None
            or system.liquidity_metrics
            for system in systems.values()
        )
        else "placeholder"
    )
    status["fig21_rna_information_entropy_compare"] = (
        "ready"
        if any(system.rna_information_summary for system in systems.values())
        else "placeholder"
    )
    status["fig22_rna_information_landscape_compare"] = (
        "ready"
        if any(
            system.rna_compaction_state_histogram is not None or system.rna_information_summary
            for system in systems.values()
        )
        else "placeholder"
    )
    status["fig23_rna_bridge_clamp_topology_compare"] = (
        "ready" if any(system.rna_topology_summary for system in systems.values()) else "placeholder"
    )
    status["fig24_local_pp_competition_compare"] = (
        "ready"
        if any(system.rna_topology_summary or system.rna_topology_binned is not None for system in systems.values())
        else "placeholder"
    )
    status["fig25_dynamic_anchor_network_compare"] = (
        "ready"
        if any(
            (
                system.network_anchor_bundle is not None
                and system.network_anchor_bundle.available()
            )
            or (
                system.network_bundle is not None
                and (system.network_bundle.gc_edges or system.network_bundle.te_edges)
            )
            for system in systems.values()
        )
        else ("skipped_by_flag" if not options.include_network_anchors else "placeholder")
    )
    status["fig26_dynamic_anchor_signature_compare"] = status["fig25_dynamic_anchor_network_compare"]
    status["figS1_degree_compare"] = "ready" if any(system.llps_degree is not None for system in systems.values()) else "placeholder"
    status["figS2_network_allostery_compare"] = (
        "ready"
        if any(system.network_bundle is not None for system in systems.values())
        else ("skipped_by_flag" if not options.include_network else "placeholder")
    )
    status["figS3_tda_compare"] = (
        "ready"
        if any(system.tda_points is not None for system in systems.values())
        else ("skipped_by_flag" if not options.include_tda else "placeholder")
    )
    status["figS4_llps_turnover_timeseries_compare"] = (
        "ready"
        if any(system.neighbor_rewiring is not None or system.interface_turnover is not None for system in systems.values())
        else "placeholder"
    )
    return status


def _read_llps_phase_partition_summary(path: Path) -> dict[str, float]:
    metrics: dict[str, float] = {}
    if not path.exists():
        return metrics
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("Metric "):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            name = parts[0]
            try:
                metrics[f"{name}_mean"] = float(parts[1])
                metrics[f"{name}_sem"] = float(parts[2])
                if len(parts) > 3:
                    metrics[f"{name}_blocks"] = float(parts[3])
            except ValueError:
                continue
    return metrics


def _read_llps_key_value_file(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.exists():
        return values
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            if len(parts) == 2:
                values[parts[0]] = parts[1]
    return values


def _as_float(value: object, default: float = math.nan) -> float:
    try:
        if value is None:
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _sem(values: Sequence[float]) -> float:
    clean = [v for v in values if math.isfinite(v)]
    if len(clean) < 2:
        return 0.0
    return float(np.std(clean, ddof=1) / math.sqrt(len(clean)))


def _read_llps_paper_replicate(row: dict[str, str], manifest_dir: Path) -> dict[str, object]:
    out_raw = row.get("output_dir", "").strip()
    output_dir = Path(out_raw)
    if not output_dir.is_absolute():
        output_dir = (manifest_dir / output_dir).resolve()
    result: dict[str, object] = dict(row)
    result["output_dir"] = str(output_dir)
    missing: list[str] = []

    phase_path = output_dir / "llps_phase_partition_summary.dat"
    summary_path = output_dir / "llps_paper_summary.json"
    convergence_path = output_dir / "llps_convergence_report.dat"
    composition_path = output_dir / "llps_composition_timeseries.csv"
    for required in (phase_path, summary_path, convergence_path, composition_path):
        if not required.exists():
            missing.append(required.name)

    result.update(_read_llps_phase_partition_summary(phase_path))
    if summary_path.exists():
        try:
            with summary_path.open("r", encoding="utf-8", errors="replace") as fh:
                for key, value in json.load(fh).items():
                    result[key] = value
        except (OSError, json.JSONDecodeError) as exc:
            missing.append(f"{summary_path.name}:parse_error:{exc}")
    convergence = _read_llps_key_value_file(convergence_path)
    for key, value in convergence.items():
        result[f"convergence_{key}"] = value

    if composition_path.exists():
        with composition_path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
            rows = list(csv.DictReader(fh))
        result["composition_frames"] = len(rows)
        if rows:
            last = rows[-1]
            result["last_dense_fraction"] = _as_float(last.get("DenseFraction"))
            result["last_rna_enrichment"] = _as_float(last.get("RNAEnrichment"))
    else:
        result["composition_frames"] = 0

    status = "ok"
    if missing:
        status = "missing_output"
    elif str(result.get("convergence_stable_dense_phase", "yes")).lower() == "no":
        status = "no_stable_dense_phase"
    elif (
        str(result.get("convergence_insufficient_blocks", "no")).lower() == "yes"
        or str(result.get("convergence_equilibration_drift_warning", "no")).lower() == "yes"
        or str(result.get("convergence_short_trajectory_warning", "no")).lower() == "yes"
    ):
        status = "convergence_warning"
    result["quality_status"] = status
    result["missing_files"] = ";".join(missing)
    return result


def run_llps_paper_ensemble(manifest_path: Path, outdir: Path, dpi: int) -> None:
    if not manifest_path.exists():
        raise SystemExit(f"Missing LLPS ensemble manifest: {manifest_path}")
    required = {
        "condition_id",
        "replica_id",
        "output_dir",
        "temperature",
        "salt",
        "protein_conc",
        "rna_conc",
        "rna_length",
        "label",
    }
    with manifest_path.open("r", encoding="utf-8-sig", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit("Manifest is empty or has no header.")
        missing_columns = required.difference(reader.fieldnames)
        if missing_columns:
            raise SystemExit(
                "Manifest missing required column(s): "
                + ", ".join(sorted(missing_columns))
            )
        replicates = [
            _read_llps_paper_replicate(row, manifest_path.parent)
            for row in reader
            if any((value or "").strip() for value in row.values())
        ]

    comparison_columns = [
        "condition_id",
        "replica_id",
        "label",
        "temperature",
        "salt",
        "protein_conc",
        "rna_conc",
        "rna_length",
        "output_dir",
        "quality_status",
        "missing_files",
        "frames_used",
        "geometry",
        "protein_dense_concentration_mean",
        "protein_dense_concentration_sem",
        "rna_dense_concentration_mean",
        "rna_dense_concentration_sem",
        "protein_partition_mean",
        "protein_partition_sem",
        "rna_partition_mean",
        "rna_partition_sem",
        "dense_protein_rna_ratio_mean",
        "rna_enrichment_mean",
        "dense_fraction_mean",
        "dense_fraction_sem",
        "dense_fraction_drift_slope_per_ns",
        "convergence_insufficient_blocks",
        "convergence_equilibration_drift_warning",
        "convergence_short_trajectory_warning",
        "convergence_stable_dense_phase",
        "composition_frames",
    ]
    comparison_path = outdir / "llps_condition_comparison.csv"
    with comparison_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=comparison_columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(replicates)

    by_condition: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in replicates:
        by_condition[str(row.get("condition_id", ""))].append(row)

    metrics = [
        "protein_dense_concentration_mean",
        "rna_dense_concentration_mean",
        "protein_partition_mean",
        "rna_partition_mean",
        "dense_protein_rna_ratio_mean",
        "rna_enrichment_mean",
        "dense_fraction_mean",
    ]
    phase_columns = [
        "condition_id",
        "label",
        "temperature",
        "salt",
        "protein_conc",
        "rna_conc",
        "rna_length",
        "replicates",
        "ok_replicates",
    ]
    for metric in metrics:
        phase_columns.extend([f"{metric}_mean", f"{metric}_sem"])

    phase_path = outdir / "llps_phase_diagram_summary.csv"
    phase_rows: list[dict[str, object]] = []
    for condition_id, rows in sorted(by_condition.items()):
        first = rows[0]
        out_row: dict[str, object] = {
            "condition_id": condition_id,
            "label": first.get("label", condition_id),
            "temperature": first.get("temperature", ""),
            "salt": first.get("salt", ""),
            "protein_conc": first.get("protein_conc", ""),
            "rna_conc": first.get("rna_conc", ""),
            "rna_length": first.get("rna_length", ""),
            "replicates": len(rows),
            "ok_replicates": sum(1 for r in rows if r.get("quality_status") == "ok"),
        }
        for metric in metrics:
            vals = [_as_float(r.get(metric)) for r in rows]
            vals = [v for v in vals if math.isfinite(v)]
            out_row[f"{metric}_mean"] = float(np.mean(vals)) if vals else math.nan
            out_row[f"{metric}_sem"] = _sem(vals)
        phase_rows.append(out_row)
    with phase_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=phase_columns)
        writer.writeheader()
        writer.writerows(phase_rows)

    report_path = outdir / "llps_replicate_quality_report.md"
    with report_path.open("w", encoding="utf-8") as fh:
        fh.write("# LLPS Replicate Quality Report\n\n")
        fh.write(f"- Manifest: `{manifest_path}`\n")
        fh.write(f"- Replicates: {len(replicates)}\n")
        fh.write(f"- Conditions: {len(by_condition)}\n\n")
        fh.write("| Condition | Replica | Status | Missing files | Drift warning | Blocks |\n")
        fh.write("| --- | --- | --- | --- | --- | --- |\n")
        for row in replicates:
            fh.write(
                f"| {row.get('condition_id', '')} | {row.get('replica_id', '')} | "
                f"{row.get('quality_status', '')} | {row.get('missing_files', '')} | "
                f"{row.get('convergence_equilibration_drift_warning', '')} | "
                f"{row.get('convergence_blocks', '')} |\n"
            )

    if phase_rows:
        labels = [str(row.get("label") or row.get("condition_id")) for row in phase_rows]
        x = np.arange(len(phase_rows), dtype=float)
        width = 0.38
        fig, axes = plt.subplots(1, 2, figsize=(max(7.0, len(labels) * 1.4), 3.6))
        axes[0].bar(
            x - width / 2,
            [float(row.get("protein_dense_concentration_mean", math.nan)) for row in phase_rows],
            width,
            yerr=[float(row.get("protein_dense_concentration_sem", 0.0)) for row in phase_rows],
            label="Protein",
            color="#2f6f9f",
        )
        axes[0].bar(
            x + width / 2,
            [float(row.get("rna_dense_concentration_mean", math.nan)) for row in phase_rows],
            width,
            yerr=[float(row.get("rna_dense_concentration_sem", 0.0)) for row in phase_rows],
            label="RNA",
            color="#b64f6f",
        )
        axes[0].set_ylabel("Dense concentration (nodes/nm$^3$)")
        axes[0].set_xticks(x)
        axes[0].set_xticklabels(labels, rotation=30, ha="right")
        axes[0].legend(frameon=False)
        axes[1].bar(
            x - width / 2,
            [float(row.get("protein_partition_mean", math.nan)) for row in phase_rows],
            width,
            yerr=[float(row.get("protein_partition_sem", 0.0)) for row in phase_rows],
            label="Protein",
            color="#2f6f9f",
        )
        axes[1].bar(
            x + width / 2,
            [float(row.get("rna_partition_mean", math.nan)) for row in phase_rows],
            width,
            yerr=[float(row.get("rna_partition_sem", 0.0)) for row in phase_rows],
            label="RNA",
            color="#b64f6f",
        )
        axes[1].set_ylabel("Partition coefficient")
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(labels, rotation=30, ha="right")
        axes[1].legend(frameon=False)
        fig.tight_layout()
        fig.savefig(outdir / "llps_condition_density_partition.png", dpi=dpi)
        plt.close(fig)

    print(f"Wrote LLPS paper ensemble outputs to: {outdir}")


def main() -> int:
    args = parse_args()
    options = build_runtime_options(args)
    if args.protein_residue_min > args.protein_residue_max:
        raise SystemExit("--protein-residue-min must be <= --protein-residue-max.")
    if args.protein_profile_smooth_window < 1:
        raise SystemExit("--protein-profile-smooth-window must be >= 1.")
    root_dir = Path(args.root_dir).resolve()
    outdir = Path(args.output_dir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    set_publication_style(args.dpi)

    if args.ensemble_manifest:
        run_llps_paper_ensemble(Path(args.ensemble_manifest).resolve(), outdir, args.dpi)
        return 0

    if args.contact_dynamics_from_averaged_dir:
        avg_dir = Path(args.contact_dynamics_from_averaged_dir).resolve()
        systems = load_averaged_interface_systems(avg_dir)
        figure_notes = []
        figure_notes.extend(make_fig6_contact_dynamics(systems, outdir, args.dpi))
        figure_notes.extend(make_fig7_prna_interface_dynamics(systems, outdir, args.dpi))
        for note in figure_notes:
            print(note)
        print(f"Wrote contact-dynamics figures from averaged dynamics CSV files to: {outdir}")
        return 0

    if args.radial_density_from_averaged_dir:
        avg_dir = Path(args.radial_density_from_averaged_dir).resolve()
        systems = load_averaged_radial_density_systems(avg_dir)
        figure_notes = []
        figure_notes.extend(
            make_radial_density_figure(
                systems,
                outdir,
                args.dpi,
                "protein",
                8,
                args.radial_density_smooth_window,
            )
        )
        figure_notes.extend(
            make_radial_density_figure(
                systems,
                outdir,
                args.dpi,
                "rna",
                9,
                args.radial_density_smooth_window,
            )
        )
        for note in figure_notes:
            print(note)
        print(f"Wrote radial density figures from averaged CSV files to: {outdir}")
        return 0

    grouped_dirs = resolve_replicate_dirs(root_dir, args.group)
    systems: dict[str, AggregatedSystem] = {}
    steady_start_arg = args.steady_start_ns
    inferred_steady_start = parse_ns_from_subdir_name(options.steady_state_subdir)
    if steady_start_arg.lower() == "auto" and inferred_steady_start is not None:
        steady_start_arg = f"{inferred_steady_start:.8f}"
    for system in SYSTEM_ORDER:
        replicates = [load_replicate(system, conf_id, directory, options) for conf_id, directory in grouped_dirs[system]]
        systems[system] = aggregate_system(system, replicates, steady_start_arg, options)
        export_system_outputs(systems[system], outdir, options)
        del replicates
        gc.collect()

    figure_notes: list[str] = []
    if not options.include_network:
        figure_notes.append("Network/allostery loading disabled by default to reduce RAM. Use --with-network to enable Supplementary Figure S2.")
    if not options.include_network_anchors:
        figure_notes.append(
            "Low-memory GC/TE anchor loading disabled by default. Use --with-network-anchors to enable Figure 25 without loading full edge maps."
        )
    elif options.include_network_anchors and not options.include_network:
        figure_notes.append(
            "Using low-memory network-anchor mode: GC/TE files are streamed into node scores for Figure 25; full Supplementary Figure S2 remains disabled."
        )
    if not options.include_tda:
        figure_notes.append("TDA loading disabled by default to reduce RAM. Use --with-tda to enable Supplementary Figure S3.")
    if not options.include_density:
        figure_notes.append("Density DX loading was skipped (--skip-density), so Figure 2 is exported as a density placeholder.")
    elif not options.include_radial_density:
        figure_notes.append("Radial density profile calculation was disabled by --no-radial-density.")
    elif not options.include_shape_tensor:
        figure_notes.append("Density shape tensor/asphericity calculation was disabled by --no-shape-tensor.")
    elif not options.export_density_dx:
        figure_notes.append("Averaged 3D density .dx export was skipped to reduce RAM. Use --export-density-dx if you need averaged DX files; radial density CSVs are still exported when DX inputs are available.")
    if options.steady_state_subdir:
        figure_notes.append(
            f"Steady-state analyses prefer files from each replicate subdirectory '{options.steady_state_subdir}', with fallback to the replicate root."
        )
    if args.steady_start_ns.lower() == "auto" and inferred_steady_start is not None:
        figure_notes.append(
            f"Steady-state time-window statistics use t >= {inferred_steady_start:.1f} ns inferred from --steady-state-subdir."
        )
    if options.include_density and options.density_subdir:
        figure_notes.append(
            f"Density DX files were read preferentially from each replicate subdirectory '{options.density_subdir}', with fallback to the replicate root."
        )
    if options.max_time_ns is not None:
        figure_notes.append(
            f"Simulation time-series were truncated at {options.max_time_ns:.1f} ns; after shorter replicates end, later points use the remaining trajectories only."
        )
    figure_notes.extend(make_fig1_timeseries(systems, outdir, args.dpi))
    figure_notes.extend(make_fig2_density_polymer(systems, outdir, args.dpi))
    figure_notes.extend(
        make_fig3_interface(
            systems,
            outdir,
            args.dpi,
            args.protein_residue_min,
            args.protein_residue_max,
            args.protein_profile_smooth_window,
        )
    )
    figure_notes.extend(make_fig4_dynamics(systems, outdir, args.dpi))
    figure_notes.extend(make_fig5_liquidity(systems, outdir, args.dpi))
    figure_notes.extend(make_fig6_contact_dynamics(systems, outdir, args.dpi))
    figure_notes.extend(make_fig7_prna_interface_dynamics(systems, outdir, args.dpi))
    figure_notes.extend(
        make_radial_density_figure(
            systems,
            outdir,
            args.dpi,
            "protein",
            8,
            options.radial_density_smooth_window,
        )
    )
    figure_notes.extend(
        make_radial_density_figure(
            systems,
            outdir,
            args.dpi,
            "rna",
            9,
            options.radial_density_smooth_window,
        )
    )
    figure_notes.extend(make_fig10_shape_asymmetry(systems, outdir, args.dpi))
    figure_notes.extend(make_fig11_dense_msd(systems, outdir, args.dpi))
    figure_notes.extend(make_fig12_rg_rna_dense_histogram(systems, outdir, args.dpi))
    figure_notes.extend(make_fig13_rna_binding_compaction_binned(systems, outdir, args.dpi))
    figure_notes.extend(make_fig14_rna_chain_compaction(systems, outdir, args.dpi))
    figure_notes.extend(make_fig15_rna_compaction_trace(systems, outdir, args.dpi))
    figure_notes.extend(make_fig16_rna_energy_compaction(systems, outdir, args.dpi))
    figure_notes.extend(make_fig17_mechanistic_summary(systems, outdir, args.dpi))
    figure_notes.extend(make_fig18_compaction_efficiency(systems, outdir, args.dpi))
    figure_notes.extend(make_fig19_dynamic_binding_efficiency(systems, outdir, args.dpi))
    figure_notes.extend(make_fig20_main_mechanism(systems, outdir, args.dpi))
    figure_notes.extend(make_fig21_rna_information_entropy(systems, outdir, args.dpi))
    figure_notes.extend(make_fig22_rna_information_landscape(systems, outdir, args.dpi))
    figure_notes.extend(make_fig23_rna_bridge_clamp_topology(systems, outdir, args.dpi))
    figure_notes.extend(make_fig24_local_pp_competition(systems, outdir, args.dpi))
    figure_notes.extend(make_fig25_dynamic_anchor_network(systems, outdir, args.dpi, args.topn))
    figure_notes.extend(make_fig26_dynamic_anchor_signature(systems, outdir, args.dpi))
    figure_notes.extend(make_figs1_degree(systems, outdir, args.dpi))
    figure_notes.extend(make_figs2_network(systems, outdir, args.dpi, args.topn))
    figure_notes.extend(make_figs3_tda(systems, outdir, args.dpi))
    figure_notes.extend(make_figs4_turnover_timeseries(systems, outdir, args.dpi))

    save_liquidity_metrics_compare(systems, outdir / "liquidity_metrics_compare.csv")
    save_density_shape_compare(systems, outdir / "density_shape_metrics_compare.csv")
    save_diffusion_metrics_compare(systems, outdir / "diffusion_metrics_compare.csv")
    save_rg_rna_dense_histogram_compare(systems, outdir / "rg_rna_dense_histogram_compare.csv")
    save_rna_compaction_binned_compare(systems, outdir / "rna_compaction_binned_compare.csv")
    save_rna_energy_compaction_binned_compare(systems, outdir / "rna_energy_compaction_binned_compare.csv")
    save_rna_information_summary_compare(systems, outdir / "rna_information_summary_compare.csv")
    save_rna_topology_summary_compare(systems, outdir / "rna_topology_summary_compare.csv")
    save_rna_chain_compaction_compare(systems, outdir / "rna_chain_compaction_compare.csv")
    save_network_anchor_compare(systems, outdir / "network_anchor_compare.csv", args.topn)
    panel_status = determine_panel_status(systems, options)
    write_comparison_summary(systems, outdir, figure_notes, panel_status, options.max_time_ns)
    print(f"Wrote comparison outputs to: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
